"""web_tools: open-web lookup for the PBPK agent. Network is injected (stubbed),
so these never touch the wire. Assert the tools return hits, fetch+strip HTML,
and record every lookup on the session for transparency (non-blind provenance)."""

import json
import unittest

from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools import web_tools as W
from pkpd_agent.state import ModelingSession


def _fake_get(pages):
    def get(url, *a, **k):
        for frag, body in pages.items():
            if frag in url:
                return body
        raise RuntimeError("no stub for " + url)
    return get


class TestWebTools(unittest.TestCase):
    def _reg(self, **ctx):
        reg = ToolRegistry()
        W.register_web_tools(reg, None, ctx)
        return reg

    def test_strip_html(self):
        self.assertEqual(W.strip_html("<p>a <b>b</b><script>x=1</script> c</p>"), "a b c")

    def test_search_merges_pubmed_and_web_and_records(self):
        reg = self._reg(
            pubmed=lambda q, k: [{"title": "Alfentanil CYP3A4", "source": "J", "year": "2018",
                                  "pmid": "1", "url": "https://pubmed/1/", "snippet": "cleared by CYP3A4"}],
            searcher=lambda q, k: [{"title": "OSP model", "url": "https://github/x", "snippet": ""}])
        s = ModelingSession(goal="g")
        r = reg.get("osp_web_search").handler({"query": "alfentanil clearance"}, s)
        self.assertTrue(r.ok)
        self.assertEqual(len(r.data["literature"]), 1)
        self.assertEqual(len(r.data["web"]), 1)
        log = s.get("web_lookups")
        self.assertEqual(log[0]["kind"], "search")
        self.assertEqual(log[0]["query"], "alfentanil clearance")

    def test_search_scope_pubmed_only(self):
        reg = self._reg(pubmed=lambda q, k: [{"title": "t", "url": "u"}],
                        searcher=lambda q, k: [{"title": "w", "url": "u2"}])
        s = ModelingSession(goal="g")
        r = reg.get("osp_web_search").handler({"query": "x", "scope": "pubmed"}, s)
        self.assertEqual(len(r.data["literature"]), 1)
        self.assertEqual(len(r.data["web"]), 0)

    def test_fetch_strips_and_records(self):
        reg = self._reg(http_get=_fake_get({"example.org": "<html><body>Hello <b>world</b></body></html>"}))
        s = ModelingSession(goal="g")
        r = reg.get("osp_web_fetch").handler({"url": "https://example.org/p"}, s)
        self.assertTrue(r.ok)
        self.assertEqual(r.data["text"], "Hello world")      # HTML stripped
        self.assertFalse(r.data["truncated"])
        self.assertEqual(s.get("web_lookups")[0]["kind"], "fetch")
        self.assertTrue(s.get("web_lookups")[0]["ok"])

    def test_fetch_caps_long_body(self):
        reg = self._reg(http_get=_fake_get({"big": "x " * 2000}))   # ~4000 chars
        r = reg.get("osp_web_fetch").handler({"url": "https://big/p", "max_chars": 500},
                                             ModelingSession(goal="g"))
        self.assertEqual(len(r.data["text"]), 500)           # capped at the 500 floor
        self.assertTrue(r.data["truncated"])

    def test_fetch_rejects_non_http(self):
        reg = self._reg()
        r = reg.get("osp_web_fetch").handler({"url": "file:///etc/passwd"}, ModelingSession(goal="g"))
        self.assertFalse(r.ok)

    def test_empty_query_errors(self):
        reg = self._reg()
        self.assertFalse(reg.get("osp_web_search").handler({"query": "  "}, ModelingSession(goal="g")).ok)

    def test_pubmed_parse_from_stubbed_eutils(self):
        esearch = json.dumps({"esearchresult": {"idlist": ["111", "222"]}})
        esummary = json.dumps({"result": {"111": {"title": "Paper A", "source": "JPK",
                                                   "pubdate": "2019 Jan"}, "222": {"title": "Paper B",
                                                   "fulljournalname": "CPT", "pubdate": "2020"}}})
        efetch = "Paper A abstract: cleared by CYP3A4.\n\n\nPaper B abstract: renal.\n"
        get = _fake_get({"esearch": esearch, "esummary": esummary, "efetch": efetch})
        hits = W.pubmed_search(get, "alfentanil", k=2)
        self.assertEqual([h["title"] for h in hits], ["Paper A", "Paper B"])
        self.assertEqual(hits[0]["year"], "2019")
        self.assertIn("CYP3A4", hits[0]["snippet"])
        self.assertEqual(hits[0]["url"], "https://pubmed.ncbi.nlm.nih.gov/111/")


if __name__ == "__main__":
    unittest.main()
