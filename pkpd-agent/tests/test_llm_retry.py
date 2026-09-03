"""The LLM call boundary retries transient overload/5xx errors with backoff (so a long multi-call
agent pass survives a single 529) but propagates non-transient errors immediately. Pure - the
anthropic dependency is stubbed."""

import unittest

from pkpd_agent.engines import llm_tasks as LT


class _Transient(Exception):
    status_code = 529


class TestRetrying(unittest.TestCase):
    def setUp(self):
        self._orig = LT._transient_types
        LT._transient_types = lambda: (_Transient,)     # avoid importing anthropic

    def tearDown(self):
        LT._transient_types = self._orig

    def test_retries_then_succeeds(self):
        n = {"i": 0}

        def flaky():
            n["i"] += 1
            if n["i"] < 3:
                raise _Transient("overloaded")
            return "ok"

        self.assertEqual(LT._retrying(flaky, attempts=6, base=1e-4), "ok")
        self.assertEqual(n["i"], 3)

    def test_non_transient_propagates_immediately(self):
        n = {"i": 0}

        def bad():
            n["i"] += 1
            raise ValueError("bad request")

        with self.assertRaises(ValueError):
            LT._retrying(bad, attempts=6, base=1e-4)
        self.assertEqual(n["i"], 1)                      # no retries on a non-transient error

    def test_exhausts_and_reraises(self):
        def always():
            raise _Transient("still down")

        with self.assertRaises(_Transient):
            LT._retrying(always, attempts=3, base=1e-4)


if __name__ == "__main__":
    unittest.main()
