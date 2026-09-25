"""LLMPolicy wiring for Anthropic's NATIVE, server-side web tools: when web=True
the request carries web_search/web_fetch, a pause_turn is resumed, and the
model's server-side searches/fetches are recorded on the session for
transparency. The anthropic client is faked — no network, no key."""

import types
import unittest

from pkpd_agent.config import AgentConfig
from pkpd_agent.llm import LLMPolicy
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.state import ModelingSession


def _blk(**kw):
    return types.SimpleNamespace(**kw)


class _FakeMessages:
    def __init__(self, responses, seen):
        self._responses = responses
        self._seen = seen

    def create(self, **kwargs):
        self._seen.append(kwargs)
        return self._responses.pop(0)


class _FakeClient:
    def __init__(self, responses, seen):
        self.messages = _FakeMessages(responses, seen)


def _policy(responses, web=True, model="claude-sonnet-5"):
    cfg = AgentConfig(mock=False)
    cfg.model = model
    p = LLMPolicy(cfg, ToolRegistry(), "sys", web=web)
    p._seen = []
    p._client = _FakeClient(responses, p._seen)
    return p


class TestNativeWebTools(unittest.TestCase):
    def test_web_tools_injected_when_enabled(self):
        resp = _blk(stop_reason="end_turn", content=[_blk(type="text", text="done")])
        p = _policy([resp], web=True)
        p.decide(ModelingSession(goal="g"))
        tools = p._seen[0]["tools"]
        types_ = {t["type"] for t in tools}
        self.assertIn("web_search_20260209", types_)
        self.assertIn("web_fetch_20260209", types_)

    def test_no_web_tools_when_disabled(self):
        resp = _blk(stop_reason="end_turn", content=[_blk(type="text", text="x")])
        p = _policy([resp], web=False)
        p.decide(ModelingSession(goal="g"))
        self.assertEqual(p._seen[0]["tools"], [])

    def test_pause_turn_is_resumed(self):
        r1 = _blk(stop_reason="pause_turn",
                  content=[_blk(type="server_tool_use", name="web_search",
                                input={"query": "alfentanil CYP3A4"})])
        r2 = _blk(stop_reason="end_turn", content=[_blk(type="text", text="found it")])
        p = _policy([r1, r2], web=True)
        step = p.decide(ModelingSession(goal="g"))
        self.assertEqual(len(p._seen), 2)                  # resumed once
        self.assertEqual(step.__class__.__name__, "FinishStep")

    def test_records_searches_and_result_urls(self):
        r1 = _blk(stop_reason="pause_turn", content=[
            _blk(type="server_tool_use", name="web_search", input={"query": "alfentanil PBPK"}),
            _blk(type="web_search_tool_result",
                 content=[_blk(url="https://pubmed/1/"), _blk(url="https://osp/model")]),
        ])
        r2 = _blk(stop_reason="end_turn", content=[
            _blk(type="server_tool_use", name="web_fetch", input={"url": "https://osp/model"}),
            _blk(type="text", text="reproduced the published model"),
        ])
        p = _policy([r1, r2], web=True)
        s = ModelingSession(goal="g")
        p.decide(s)
        log = s.get("web_lookups")
        kinds = [w["kind"] for w in log]
        self.assertIn("search", kinds)
        self.assertIn("results", kinds)
        self.assertIn("fetch", kinds)
        q = next(w["query"] for w in log if w["kind"] == "search")
        self.assertEqual(q, "alfentanil PBPK")
        urls = next(w["urls"] for w in log if w["kind"] == "results")
        self.assertIn("https://osp/model", urls)

    def test_no_recording_when_web_off(self):
        r = _blk(stop_reason="end_turn", content=[_blk(type="text", text="x")])
        p = _policy([r], web=False)
        s = ModelingSession(goal="g")
        p.decide(s)
        self.assertIsNone(s.get("web_lookups"))


if __name__ == "__main__":
    unittest.main()
