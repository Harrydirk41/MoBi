"""Fit-vision: the model SEES its own predicted-vs-observed curves.

Covers the whole side channel: the Pillow renderer, the series builder, the
ToolResult->Observation->tool_result image-block plumbing (base64 stays OUT of
the streamed/persisted transcript), and the run broadcaster that lets several
pages watch one run at once.
"""

import base64
import queue
import unittest

from pkpd_agent.engines import fit_plot, osp_score


_SERIES = [
    {"study": "Kroboth 1988", "route": "PO", "dose": "1 mg",
     "observed": [[0.5, 0.02], [1, 0.05], [2, 0.04], [4, 0.02], [8, 0.008]],
     "simulated": [[0, 0.0], [0.5, 0.03], [1, 0.052], [2, 0.038], [4, 0.019], [8, 0.007]]},
    {"study": "Smith 1990", "route": "IV", "dose": "2 mg",
     "observed": [[0.1, 0.2], [1, 0.08], [4, 0.015]],
     "simulated": [[0, 0.25], [0.1, 0.19], [1, 0.075], [4, 0.014]]},
]


class TestRenderer(unittest.TestCase):
    def test_render_produces_a_valid_png(self):
        png = fit_plot.render_fit_png(_SERIES, "test")
        self.assertIsInstance(png, bytes)
        self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))     # PNG magic
        self.assertGreater(len(png), 500)

    def test_render_b64_descriptor(self):
        d = fit_plot.render_fit_b64(_SERIES, "t")
        self.assertEqual(d["media_type"], "image/png")
        self.assertTrue(png := base64.b64decode(d["data"]))
        self.assertTrue(png.startswith(b"\x89PNG"))

    def test_nothing_plottable_returns_none(self):
        self.assertIsNone(fit_plot.render_fit_png([]))
        self.assertIsNone(fit_plot.render_fit_b64([]))
        # all concentrations non-positive -> nothing to place on a log axis
        self.assertIsNone(fit_plot.render_fit_png(
            [{"study": "x", "observed": [[1, -1], [2, 0]], "simulated": []}]))

    def test_observed_only_still_renders(self):
        png = fit_plot.render_fit_png(
            [{"study": "obs only", "observed": [[1, 0.1], [2, 0.05]], "simulated": []}])
        self.assertTrue(png and png.startswith(b"\x89PNG"))


class TestOverlaySeries(unittest.TestCase):
    def test_pairs_observed_with_predicted_by_dataset(self):
        observed = [{"dataset": "A", "study": "St A", "route": "PO", "dose": "1 mg",
                     "time_h": [1, 2, 4], "conc_mg_L": [0.1, 0.05, 0.02]}]
        predicted = [{"dataset": "A", "time_h": [0, 1, 2, 4],
                      "pred_conc_mg_L": [0, 0.11, 0.048, 0.021]}]
        ser = osp_score.overlay_series(observed, predicted)
        self.assertEqual(len(ser), 1)
        self.assertEqual(ser[0]["study"], "St A")
        self.assertEqual(len(ser[0]["observed"]), 3)
        self.assertEqual(len(ser[0]["simulated"]), 4)

    def test_unmatched_dataset_has_empty_simulated(self):
        observed = [{"dataset": "A", "time_h": [1], "conc_mg_L": [0.1]}]
        ser = osp_score.overlay_series(observed, [])          # no predictions
        self.assertEqual(ser[0]["simulated"], [])
        self.assertEqual(len(ser[0]["observed"]), 1)

    def test_downsample_caps_dense_curves(self):
        observed = [{"dataset": "A", "time_h": [0], "conc_mg_L": [1.0]}]
        predicted = [{"dataset": "A", "time_h": list(range(500)),
                      "pred_conc_mg_L": [1.0] * 500}]
        ser = osp_score.overlay_series(observed, predicted, cap=60)
        self.assertLessEqual(len(ser[0]["simulated"]), 61)     # cap (+ terminal point)


class TestImageSideChannel(unittest.TestCase):
    """The base64 image must ride a side channel: out of ``content`` (which is
    streamed and persisted), into a tool_result image block at message-build."""

    def test_toolresult_images_not_in_content(self):
        from pkpd_agent.tools.registry import ToolResult
        r = ToolResult.success("ok", gmfe=1.5)
        r.images = [{"media_type": "image/png", "data": "QQ=="}]
        self.assertNotIn("images", r.to_content())            # side channel, not content
        self.assertEqual(r.to_content()["gmfe"], 1.5)

    def test_observe_attaches_image_block(self):
        from pkpd_agent.llm import LLMPolicy
        from pkpd_agent.state import ModelingSession, Observation, Decision
        from pkpd_agent.tools.registry import ToolRegistry
        from pkpd_agent.config import AgentConfig
        pol = LLMPolicy(AgentConfig(), ToolRegistry(), "sys")
        sess = ModelingSession(goal="g")
        pol.observe(sess)                                     # seeds the goal turn
        sess.record(Decision(text="", calls=[]))
        sess.record(Observation(call_id="c1", tool="osp_optimize", ok=True,
                                content={"message": "GMFE 1.5", "gmfe": 1.5},
                                images=[{"media_type": "image/png", "data": "QQ=="}]))
        pol.observe(sess)
        content = pol._messages[-1]["content"]               # last user turn = tool_result(s)
        tr = content[0]
        self.assertEqual(tr["type"], "tool_result")
        blocks = tr["content"]
        self.assertIsInstance(blocks, list)                  # text + image, not a bare string
        self.assertEqual(blocks[0]["type"], "text")
        self.assertEqual(blocks[1]["type"], "image")
        self.assertEqual(blocks[1]["source"]["data"], "QQ==")

    def test_observe_plain_text_when_no_image(self):
        from pkpd_agent.llm import LLMPolicy
        from pkpd_agent.state import ModelingSession, Observation, Decision
        from pkpd_agent.tools.registry import ToolRegistry
        from pkpd_agent.config import AgentConfig
        pol = LLMPolicy(AgentConfig(), ToolRegistry(), "sys")
        sess = ModelingSession(goal="g")
        pol.observe(sess)
        sess.record(Decision(text="", calls=[]))
        sess.record(Observation(call_id="c1", tool="osp_inspect", ok=True,
                                content={"message": "hi"}))
        pol.observe(sess)
        self.assertIsInstance(pol._messages[-1]["content"][0]["content"], str)


class TestBroadcaster(unittest.TestCase):
    """The run queue fans out to many subscribers and replays a backlog, so two
    pages (or two machines) can watch one run without stealing each other's
    events."""

    def test_fanout_and_backlog(self):
        from copilot.server import _LoggingQueue
        q = _LoggingQueue()
        q.put({"type": "meta", "n": 0})                      # before anyone subscribes
        backlog_a, sub_a = q.subscribe()
        q.put({"type": "decision", "n": 1})
        backlog_b, sub_b = q.subscribe()                     # a late joiner
        q.put({"type": "end", "n": 2})
        # A saw event 0 in backlog, and 1 + 2 live
        self.assertEqual([e["n"] for e in backlog_a], [0])
        self.assertEqual(sub_a.get_nowait()["n"], 1)
        self.assertEqual(sub_a.get_nowait()["n"], 2)
        # B joined after event 1, so its backlog has 0 and 1, and it streams 2 live
        self.assertEqual([e["n"] for e in backlog_b], [0, 1])
        self.assertEqual(sub_b.get_nowait()["n"], 2)
        # nothing duplicated: each live queue has exactly one item left then empties
        with self.assertRaises(queue.Empty):
            sub_a.get_nowait()

    def test_unsubscribe_stops_delivery(self):
        from copilot.server import _LoggingQueue
        q = _LoggingQueue()
        _, sub = q.subscribe()
        q.unsubscribe(sub)
        q.put({"type": "end"})
        with self.assertRaises(queue.Empty):
            sub.get_nowait()
        self.assertEqual(q.log[-1]["type"], "end")           # still recorded for replay


if __name__ == "__main__":
    unittest.main()
