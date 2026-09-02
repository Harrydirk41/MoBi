"""Decomposable pipeline: each layer is a Provider with a mode; run() fills the spec and records
per-layer modes. Covers given / data / llm(mocked) providers. Pure."""

import unittest

from pkpd_agent.engines import pipeline as P

PROV = {
    "IL6SecFLS_MaxbyIL1b": {"value_from_reference": 69.0, "from_literature": True, "reference": "F"},
    "kcl_IL6": {"value_from_reference": 5.5, "from_literature": True, "reference": "Y"},
}
LEVELS = {"IL6": 289.0, "IL1b": 1.2, "TNFa": 0.9}


def _ctx(call=None):
    return {"prov": PROV, "levels": LEVELS, "truth": {"IL1b"}, "call": call}


class TestProviders(unittest.TestCase):
    def test_given_and_data_modes_recorded(self):
        providers = {
            "frame": P.given({"objective": "o"}),
            "target": P.given("IL6"),
            "topology": P.data(lambda c: [{"cytokine": "IL1b", "direction": "up"}]),
            "form": P.given({"proliferation_order": "zeroth", "combination": "product", "cap": None}),
        }
        spec = P.run(providers, _ctx())
        self.assertEqual(spec["modes"], {"frame": P.GIVEN, "target": P.GIVEN,
                                         "topology": P.DATA, "form": P.GIVEN})
        self.assertEqual(spec["target"], "IL6")
        self.assertEqual([e["src"] for e in spec["edges"]], ["IL1b"])
        self.assertTrue(spec["edges"][0]["verify"]["in_model"])

    def test_llm_provider_uses_call_boundary(self):
        seen = {}
        def call(system, user):
            seen["system"] = system
            return '{"regulators":[{"cytokine":"IL1b","direction":"up","basis":"x"}]}'
        providers = {
            "frame": P.given({}),
            "target": P.given("IL6"),
            "topology": P.from_llm("SYS", lambda c: "u",
                                   lambda s: __import__("json").loads(s)["regulators"]),
            "form": P.given({"proliferation_order": "zeroth", "combination": "product", "cap": None}),
        }
        spec = P.run(providers, _ctx(call=call))
        self.assertEqual(spec["modes"]["topology"], P.LLM)
        self.assertEqual(seen["system"], "SYS")             # went through the call boundary
        self.assertEqual([e["src"] for e in spec["edges"]], ["IL1b"])

    def test_missing_provider_errors(self):
        with self.assertRaises(ValueError):
            P.run({"frame": P.given({})}, _ctx())


if __name__ == "__main__":
    unittest.main()
