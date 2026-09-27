"""One source of truth for walling the OSP answer source off from the native web
tools, used by BOTH the decision policy (``llm.py``) and the literature-reading
Route-B call (``llm_tasks.py``).

Why a module and not one list: Anthropic's ``blocked_domains`` matches web_search
and web_fetch DIFFERENTLY (verified against the server-tools domain-filtering docs):

  * web_search  - subpaths ARE honored: ``github.com/Open-Systems-Pharmacology``
                  blocks that org's repos while the rest of GitHub stays searchable.
  * web_fetch   - matches on the DOMAIN ONLY; "an entry that includes a path never
                  matches a web fetch URL." So a path-scoped entry is silently
                  IGNORED by web_fetch, leaving the OSP repos fetchable — a real
                  leak of the benchmark answers.

Therefore the fetch blocklist must name the BARE HOSTS (github.com,
raw.githubusercontent.com). That over-blocks non-OSP GitHub for direct fetches,
which is the accepted cost of the inviolable "never see the answer source"
guarantee: the agent can still DISCOVER non-OSP GitHub via web_search and read the
snippets search returns.
"""

from typing import Any

# The benchmark ANSWERS live in the OSP model library: its website, its GitHub org
# (Open-Systems-Pharmacology/OSP-PBPK-Model-Library and the <Compound>-Model repos),
# and that org's raw file host.
_OSP_SITE = "open-systems-pharmacology.org"
_OSP_GITHUB_ORG = "Open-Systems-Pharmacology"
_GITHUB_HOSTS = ("github.com", "raw.githubusercontent.com")


def search_blocklist() -> list[str]:
    """Path-scoped: block ONLY the OSP org path on GitHub (keeping the rest of GitHub
    searchable for genuine literature) plus the OSP website (subdomains auto-included).
    NOTE: a personal MIRROR of a model outside this org is not caught by an org-path
    block — an inherent limit of source-blocking, not a bug here."""
    return [f"{h}/{_OSP_GITHUB_ORG}" for h in _GITHUB_HOSTS] + [_OSP_SITE]


def fetch_blocklist() -> list[str]:
    """Bare hosts: web_fetch ignores any entry with a path, so path-scoping the OSP
    org would be a no-op and the answer repos would stay fetchable. Block the whole
    GitHub hosts for direct fetch, plus the OSP site."""
    return list(_GITHUB_HOSTS) + [_OSP_SITE]


def web_tools(model: str | None, want_fetch: bool = True,
              max_uses: int | None = None) -> list[dict[str, Any]]:
    """Build Anthropic's native server-side web tools with the answer source walled
    off. Newer Claude-family models get the dynamic-filtering variants + web_fetch;
    older ones get basic web_search only (web_fetch is not offered there, so there is
    no unblockable fetch path). ``max_uses`` caps searches/fetches when given."""
    m = (model or "").lower()
    new = any(k in m for k in ("opus-5", "opus-4", "sonnet-5", "sonnet-4-6"))

    def _cap(t: dict) -> dict:
        if max_uses is not None:
            t["max_uses"] = max_uses
        return t

    if new:
        tools = [_cap({"type": "web_search_20260209", "name": "web_search",
                       "blocked_domains": search_blocklist()})]
        if want_fetch:
            tools.append(_cap({"type": "web_fetch_20260209", "name": "web_fetch",
                               "blocked_domains": fetch_blocklist()}))
        return tools
    # older models: basic search only, path-scoped block applies correctly
    return [_cap({"type": "web_search_20250305", "name": "web_search",
                  "blocked_domains": search_blocklist()})]
