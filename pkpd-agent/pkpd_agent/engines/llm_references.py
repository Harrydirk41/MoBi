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


def pubmed_fetch(title: str, timeout: float = 20.0) -> str:
    """Look a title up in PubMed (E-utilities esearch) and return its abstract text (efetch).
    Stdlib urllib over the environment's HTTPS proxy. Returns '' on any miss or error - a
    fetcher must never break the assembly of the rest of the material."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    try:
        q = urllib.parse.urlencode({"db": "pubmed", "term": title + "[Title]",
                                    "retmax": "1", "retmode": "json"})
        with urllib.request.urlopen(f"{base}/esearch.fcgi?{q}", timeout=timeout) as r:
            import json
            ids = json.load(r).get("esearchresult", {}).get("idlist", [])
        if not ids:                             # retry without the [Title] field restriction
            q = urllib.parse.urlencode({"db": "pubmed", "term": title,
                                        "retmax": "1", "retmode": "json"})
            with urllib.request.urlopen(f"{base}/esearch.fcgi?{q}", timeout=timeout) as r:
                import json
                ids = json.load(r).get("esearchresult", {}).get("idlist", [])
        if not ids:
            return ""
        q = urllib.parse.urlencode({"db": "pubmed", "id": ids[0],
                                    "rettype": "abstract", "retmode": "text"})
        with urllib.request.urlopen(f"{base}/efetch.fcgi?{q}", timeout=timeout) as r:
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
