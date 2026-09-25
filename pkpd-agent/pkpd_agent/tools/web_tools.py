"""Open-web lookup tools for the PBPK build agent.

A real modeler does not work sealed off from the literature: they look up a
compound's DMPK (which enzyme clears it, measured physchem, transporters) and,
when a PBPK model has already been published, they reproduce and build on it.
These tools give the agent that reach:

  * ``osp_web_search`` - a query -> ranked hits. PubMed (NCBI E-utilities, keyless)
    is the structured literature backend; an optional general-web backend adds
    non-PubMed hits (docs, model repositories). Returns title / source / year /
    url / snippet.
  * ``osp_web_fetch``  - fetch one URL and return its readable text (HTML stripped,
    capped), so the agent can read a specific paper, a DOI page, or a model file.

HONESTY NOTE. Unlike the context-lake tools, this is NOT confined to the handed
materials and NOT leave-one-out: with the open web the agent can reach the target
compound's own published model. That is a deliberate choice - the task is "model
like a real modeler, who can look up prior work" - so a run that used the web is
NOT a blind de-novo grade. The tools record every query/fetch on the session
(``web_lookups``) so the report can state plainly what was consulted; they do not
hide or exclude anything.

The network functions are injectable (``ctx['http_get']`` / ``ctx['searcher']``)
so tests run without touching the network.
"""
from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request

from .registry import Tool, ToolRegistry, ToolResult

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_UA = {"User-Agent": "pkpd-agent/1.0 (PBPK modeling; mailto:noreply@example.org)"}


def _default_get(url: str, timeout: float = 20.0, cap: int = 2_000_000) -> str:
    """GET a URL and return decoded text. Honors HTTPS_PROXY via urllib. Never
    used in tests (they inject a stub)."""
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:      # noqa: S310
        data = r.read(cap)
        enc = r.headers.get_content_charset() or "utf-8"
    return data.decode(enc, "replace")


def strip_html(html: str) -> str:
    html = re.sub(r"(?is)<(script|style|head|noscript)[^>]*>.*?</\1>", " ", html)
    html = re.sub(r"(?s)<[^>]+>", " ", html)
    html = (html.replace("&nbsp;", " ").replace("&amp;", "&")
            .replace("&lt;", "<").replace("&gt;", ">").replace("&#39;", "'")
            .replace("&quot;", '"'))
    return re.sub(r"\s+", " ", html).strip()


def pubmed_search(get, query: str, k: int = 5) -> list[dict]:
    """esearch (relevance) -> esummary; best-effort abstract via efetch. Returns
    up to k {title, source, year, pmid, url, snippet}. Never raises."""
    try:
        q = urllib.parse.urlencode({"db": "pubmed", "term": query, "retmax": k,
                                    "retmode": "json", "sort": "relevance"})
        js = json.loads(get(f"{_EUTILS}/esearch.fcgi?{q}"))
        ids = js.get("esearchresult", {}).get("idlist", []) or []
    except Exception:                                            # noqa: BLE001
        return []
    if not ids:
        return []
    hits: list[dict] = []
    try:
        q2 = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids),
                                     "retmode": "json"})
        summ = json.loads(get(f"{_EUTILS}/esummary.fcgi?{q2}")).get("result", {})
    except Exception:                                            # noqa: BLE001
        summ = {}
    # one efetch for all abstracts, split on the blank-line record separator
    abstracts: dict[str, str] = {}
    try:
        q3 = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids),
                                     "rettype": "abstract", "retmode": "text"})
        blob = get(f"{_EUTILS}/efetch.fcgi?{q3}")
        recs = [r.strip() for r in re.split(r"\n\s*\n\s*\n", blob) if r.strip()]
        for pid, rec in zip(ids, recs):     # order matches the id list
            abstracts[pid] = re.sub(r"\s+", " ", rec)[:700]
    except Exception:                                            # noqa: BLE001
        pass
    for pid in ids:
        rec = summ.get(pid, {}) if isinstance(summ, dict) else {}
        hits.append({
            "title": rec.get("title") or f"PMID {pid}",
            "source": rec.get("fulljournalname") or rec.get("source") or "PubMed",
            "year": (rec.get("pubdate") or "")[:4],
            "pmid": pid,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pid}/",
            "snippet": abstracts.get(pid, ""),
        })
    return hits


def web_search(get, query: str, k: int = 5) -> list[dict]:
    """Best-effort general-web search via DuckDuckGo's keyless lite HTML endpoint.
    Returns {title, url, snippet}. Empty on any failure (the caller still has
    PubMed)."""
    try:
        q = urllib.parse.urlencode({"q": query})
        html = get(f"https://lite.duckduckgo.com/lite/?{q}")
    except Exception:                                            # noqa: BLE001
        return []
    out: list[dict] = []
    for m in re.finditer(r'<a[^>]+href="(http[^"]+)"[^>]*class="result-link"[^>]*>(.*?)</a>',
                         html, re.S):
        url, title = m.group(1), strip_html(m.group(2))
        if title and url:
            out.append({"title": title, "url": url, "snippet": ""})
        if len(out) >= k:
            break
    if not out:            # fallback: any external result anchors
        for m in re.finditer(r'<a[^>]+href="(https?://[^"]+)"[^>]*>(.*?)</a>', html, re.S):
            url, title = m.group(1), strip_html(m.group(2))
            if title and "duckduckgo" not in url and len(title) > 8:
                out.append({"title": title, "url": url, "snippet": ""})
            if len(out) >= k:
                break
    return out


def register_web_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {input?, http_get?, pubmed?, searcher?}. The network functions default
    to real HTTP; tests inject stubs. Every lookup is recorded on the session under
    'web_lookups' for the report (transparency, not exclusion)."""
    get = ctx.get("http_get") or _default_get
    pm = ctx.get("pubmed") or (lambda q, k: pubmed_search(get, q, k))
    ws = ctx.get("searcher") or (lambda q, k: web_search(get, q, k))

    def _record(session, kind, detail):
        log = session.get("web_lookups") or []
        log.append({"kind": kind, **detail})
        session.put("web_lookups", log)

    def do_search(args: dict, session) -> ToolResult:
        query = (args.get("query") or "").strip()
        if not query:
            return ToolResult.error("Provide a 'query' to search for.")
        k = max(1, min(int(args.get("k") or 5), 10))
        scope = (args.get("scope") or "all").lower()
        lit = pm(query, k) if scope in ("all", "pubmed", "literature") else []
        gen = ws(query, k) if scope in ("all", "web") else []
        _record(session, "search", {"query": query, "scope": scope,
                                    "n_pubmed": len(lit), "n_web": len(gen)})
        if not lit and not gen:
            return ToolResult.success(
                "No results (or the search backend was unreachable). Try a "
                "different query, or osp_web_fetch a URL you already know.",
                query=query, literature=[], web=[])
        return ToolResult.success(
            f"{len(lit)} literature + {len(gen)} web hit(s) for '{query}'. Use "
            "osp_web_fetch on a url to read it. NOTE: web lookups are recorded and "
            "make this run non-blind - cite what you use in your final summary.",
            query=query, literature=lit, web=gen)

    def do_fetch(args: dict, session) -> ToolResult:
        url = (args.get("url") or "").strip()
        if not re.match(r"^https?://", url):
            return ToolResult.error("Provide an http(s) 'url' to fetch.")
        cap = max(500, min(int(args.get("max_chars") or 6000), 20000))
        try:
            text = strip_html(get(url))
        except Exception as e:                                   # noqa: BLE001
            _record(session, "fetch", {"url": url, "ok": False})
            return ToolResult.error(f"Could not fetch {url}: {type(e).__name__}")
        _record(session, "fetch", {"url": url, "ok": True, "chars": len(text)})
        return ToolResult.success(
            f"fetched {url} ({len(text)} chars). Extract the DMPK facts / model "
            "structure you need and, if you use them, cite this source.",
            url=url, text=text[:cap], truncated=len(text) > cap)

    registry.register(Tool(
        name="osp_web_search",
        description=(
            "SEARCH the open web / literature for THIS compound's DMPK and prior "
            "modeling: which enzyme(s) clear it, measured physicochemistry, "
            "transporters, and whether a PBPK model has already been published. "
            "Returns ranked hits (title, source, year, url, snippet). 'scope' = "
            "'all' (default) | 'pubmed' | 'web'. This is real open-web access - it "
            "is NOT leave-one-out, so a run that uses it is not a blind de-novo "
            "grade; cite what you rely on."),
        input_schema={"type": "object", "properties": {
            "query": {"type": "string", "description": "search terms, e.g. "
                      "'alfentanil CYP3A4 clearance PBPK'"},
            "scope": {"type": "string", "enum": ["all", "pubmed", "web"]},
            "k": {"type": "integer", "description": "max hits per backend (1-10)"}},
            "required": ["query"]},
        handler=do_search, phase="observe"))

    registry.register(Tool(
        name="osp_web_fetch",
        description=(
            "FETCH one URL and return its readable text (HTML stripped, capped). "
            "Use it to read a paper, a DOI page, or a published model file you "
            "found via osp_web_search or already know. Cite anything you use."),
        input_schema={"type": "object", "properties": {
            "url": {"type": "string"},
            "max_chars": {"type": "integer", "description": "cap on returned text (500-20000)"}},
            "required": ["url"]},
        handler=do_fetch, phase="observe"))
