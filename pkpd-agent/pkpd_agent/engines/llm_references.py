"""Route A for topology drafting: assemble the READING MATERIAL and cache it.

The drafter (``llm_topology.draft_topology``) grounds its edges in whatever text it is
handed. This module builds that text from a paper's own bibliography so the evidence is
fixed, inspectable, and committable with the benchmark - not re-fetched differently on
every run (which is Route B's job, via a web-enabled call).

  * ``parse_references``  - split a References section into individual numbered citations,
                            tolerating the ragged multi-line numbering papers use.
  * ``format_references`` - render citations back to one-per-line, the compact drafter feed.
  * ``fetch_abstracts``   - enrich each citation with its abstract via a PLUGGABLE fetcher
                            (tests stub it; real runs use ``pubmed_fetch``).
  * ``pubmed_fetch``      - a concrete fetcher: NCBI E-utilities title search -> abstract,
                            stdlib-only over the environment's HTTPS proxy, never raises.

General: a numbered-citation text parser plus a title->abstract lookup. Nothing here names
a disease or a model; point it at any paper's reference list.
"""

from __future__ import annotations

import re
import time
import urllib.parse
import urllib.request

# a citation marker: a number then '.' or ')' at the start of a line (the paper prints the
# number alone on its line, or leading the text - both match).
_MARK = re.compile(r"(?m)^[ \t]*(\d{1,3})[.)][ \t]*")


def parse_references(text: str) -> list[dict]:
    """Split a References section into ``[{"n": int, "text": str}]``. Slices the text
    between successive numbered markers and collapses each citation's internal line breaks
    into single spaces. Ignores a leading 'References' heading and any out-of-order markers
    (keeps only a monotonically increasing sequence, so a stray '2020)' year is not a new
    entry)."""
    marks = [(m.start(), int(m.group(1)), m.end()) for m in _MARK.finditer(text)]
    out: list[dict] = []
    expected = 1
    for i, (start, n, end) in enumerate(marks):
        if n != expected:                       # not the next citation in sequence - skip
            continue
        stop = len(text)
        for start2, n2, _ in marks[i + 1:]:     # end at the NEXT expected marker
            if n2 == expected + 1:
                stop = start2
                break
        body = re.sub(r"\s+", " ", text[end:stop]).strip()
        if body:
            out.append({"n": n, "text": body})
            expected += 1
    return out


def extract_references_section(text: str) -> str:
    """Slice the References section out of a full-paper text dump: from the last line that is
    just 'References'/'Bibliography' to the first end-matter heading (Acknowledgements / Author
    contributions / Competing interests) or end of text. Returns '' if no heading is found, so
    the caller can fall back to treating the whole input as the reference list."""
    starts = [m.start() for m in
              re.finditer(r"(?im)^[ \t]*(references|bibliography)[ \t]*$", text)]
    if not starts:
        return ""
    body = text[starts[-1]:]
    end = re.search(r"(?im)^[ \t]*(acknowledge?ments?|author contributions?|"
                    r"competing interests?|declaration of|conflicts? of interest)\b", body)
    return body[: end.start()] if end else body


def format_references(citations: list[dict]) -> str:
    """One citation per line, ``n. text`` - the compact form fed to the drafter."""
    return "\n".join(f"{c['n']}. {c['text']}" +
                     (f"\n   ABSTRACT: {c['abstract']}" if c.get("abstract") else "")
                     for c in citations)


def _title_of(citation_text: str) -> str:
    """Best-effort article title from a citation string: strip the leading author list and the
    trailing journal+volume+year, keep the first sentence between. Good enough to search PubMed
    (a miss just means no abstract for that entry); not a rigorous citation parser."""
    t = citation_text.strip()
    m = re.search(r"\bet\s+al\.?,?\s+", t)                     # authors end at 'et al'
    if m:
        rest = t[m.end():]
    else:                                                      # else consume the author run
        am = re.match(r"^(?:(?:[A-Z][A-Za-z'’-]+|van|de|von|del|der)[, ]+"
                      r"(?:[A-Z]\.[-\s]?)+[, ]*(?:&\s*)?)+", t)
        rest = t[am.end():] if am else t
    ym = list(re.finditer(r"\(\d{4}\)", rest))                 # drop trailing journal+year
    if ym:
        rest = rest[: ym[-1].start()]
    return rest.split(". ")[0].strip(" .,")


_LIGATURES = {"ﬁ": "fi", "ﬂ": "fl", "ﬀ": "ff", "ﬃ": "ffi", "ﬄ": "ffl",
              "–": "-", "—": "-", "’": "'", "‘": "'", "“": '"', "”": '"'}


def _clean(s: str) -> str:
    for k, v in _LIGATURES.items():
        s = s.replace(k, v)
    return s


def _sig_words(s: str) -> set:
    return {w for w in re.findall(r"[a-z0-9]+", _clean(s).lower()) if len(w) > 3}


def _same_paper(want: set, cand: set, min_overlap: float) -> float:
    """Is candidate title ``cand`` the same paper as query title ``want`` (both sig-word sets)?
    Returns a match score (0 = reject). A long, specific query (>=5 words) almost entirely
    contained in the candidate is accepted even if the extracted query title was truncated; a
    short/generic query must match tightly (Jaccard) so an unrelated recent paper repeating a
    few common words does not false-match."""
    if not want or not cand:
        return 0.0
    inter = len(want & cand)
    cont = inter / len(want)
    jac = inter / len(want | cand)
    if jac >= min_overlap or (len(want) >= 5 and cont >= 0.85):
        return jac
    return 0.0


_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get_json(path: str, params: dict, timeout: float):
    import json
    q = urllib.parse.urlencode(params)
    with urllib.request.urlopen(f"{_EUTILS}/{path}?{q}", timeout=timeout) as r:
        return json.load(r)


def pubmed_fetch(title: str, timeout: float = 20.0, min_overlap: float = 0.6) -> str:
    """Look a title up in PubMed and return its abstract text, or '' on any miss/error - a
    fetcher must never break the assembly of the rest of the material. Stdlib urllib over the
    environment's HTTPS proxy.

    Correctness over coverage: PubMed's phrase index is quirky (hyphens like 'IL-6', British
    spelling, ligatures), so a strict phrase query misses real papers, while a loose query
    silently returns an unrelated recent one - a wrong match corrupts the reading material far
    worse than an honest miss. So we search by relevance for a few candidates, then VERIFY:
    accept a candidate only if it contains at least ``min_overlap`` of the query title's
    significant words. Unverifiable -> '' (miss), never a guess."""
    want = _sig_words(title)
    if len(want) < 3:                            # too little to identify a paper -> skip
        return ""
    try:
        res = _get_json("esearch.fcgi", {"db": "pubmed", "term": _clean(title),
                                         "retmax": "5", "sort": "relevance",
                                         "retmode": "json"}, timeout)
        ids = res.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return ""
        summ = _get_json("esummary.fcgi", {"db": "pubmed", "id": ",".join(ids),
                                           "retmode": "json"}, timeout).get("result", {})
        best, best_score = "", 0.0
        for pid in ids:
            cand = _sig_words(summ.get(pid, {}).get("title", ""))
            score = _same_paper(want, cand, min_overlap)
            if score > best_score:
                best, best_score = pid, score
        if not best:                             # no candidate is confidently the same paper
            return ""
        q = urllib.parse.urlencode({"db": "pubmed", "id": best,
                                    "rettype": "abstract", "retmode": "text"})
        with urllib.request.urlopen(f"{_EUTILS}/efetch.fcgi?{q}", timeout=timeout) as r:
            return re.sub(r"\s+", " ", r.read().decode("utf-8", "replace")).strip()
    except Exception:
        return ""


def fetch_abstracts(citations: list[dict], fetch=pubmed_fetch,
                    sleep: float = 0.34) -> list[dict]:
    """Fill each citation's ``abstract`` via ``fetch(title) -> str`` (default PubMed).
    Rate-limited between calls for polite public APIs. Returns the SAME citation dicts,
    mutated in place, so the caller can cache and inspect exactly what was retrieved."""
    for i, c in enumerate(citations):
        c["abstract"] = fetch(_title_of(c["text"])) or ""
        if sleep and i + 1 < len(citations):
            time.sleep(sleep)
    return citations
