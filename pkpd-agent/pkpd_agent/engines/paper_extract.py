"""Y3: pull the machine-readable inputs out of the paper - and be honest about the rest.

Two kinds of per-model input come from the paper, not the model file:
  * clinical reference data (trial ACR/DAS numbers)  - these live in the TEXT/tables
  * the global-sensitivity top-params (the GSA list)  - these live in a FIGURE (Fig 9)

The text/table part is reliably extractable (``pdf_text`` / ``find_param_tokens``). The
figure part is NOT: the GSA parameter labels are vector-drawn, so PDF text extraction gets
almost none of them (empirically 2 of ~20). That input needs a VISION step - render the
figure page to an image and read it (a vision-capable model, or a person) - which is
exactly how a project's spec.json:gsa_top was produced. ``render_page`` does
the rendering; the reading is the vision step, deliberately not faked here with brittle OCR.
"""

from __future__ import annotations

import re


def pdf_text(path: str) -> str:
    import pymupdf
    doc = pymupdf.open(path)
    return "\n".join(f"\n===== PAGE {i + 1} =====\n{p.get_text()}"
                     for i, p in enumerate(doc))


def pdf_page_text(path: str, page: int) -> str:
    import pymupdf
    return pymupdf.open(path)[page - 1].get_text()


def render_page(path: str, page: int, out_png: str, zoom: float = 3.0) -> str:
    """Render one page to PNG for the vision step (reading a figure like the GSA chart)."""
    import pymupdf
    doc = pymupdf.open(path)
    pix = doc[page - 1].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
    pix.save(out_png)
    return out_png


def find_param_tokens(path: str) -> list[str]:
    """All parameter-name-like tokens that are SELECTABLE text in the PDF. Useful for
    tables/appendices; will miss labels that are drawn as vector graphics (e.g. figure
    axes), which is why the GSA list is not reliably auto-extractable."""
    toks = set(re.findall(r"[A-Za-z]+_[A-Za-z0-9_]+", pdf_text(path)))
    return sorted(t for t in toks if len(t) > 3)


_TRIPLE = r"[^\n]*?(\d{1,3}(?:\.\d)?)\s+(\d{1,3}(?:\.\d)?)\s+(\d{1,3}(?:\.\d)?)"


def _acr_row(drugs=None):
    """Compile the 'drug ... n n n' row regex, anchored on the given drug names (or any
    capitalized drug-like token when none are given). Model-agnostic: no drug is baked in."""
    drug_alt = "|".join(re.escape(d) for d in drugs) if drugs else r"[A-Z][A-Za-z]{2,}"
    return re.compile(rf"({drug_alt})\b{_TRIPLE}")


def find_clinical_acr(path: str, drugs=None) -> list[dict]:
    """Best-effort scan for 'drug ... ACR20 ACR50 ACR70' triples in the text. Pass the
    project's drug names (e.g. from tasks.json) to anchor the rows; with none given it
    matches any capitalized drug-like token. Flattened PDF tables are noisy, so treat
    this as a first pass to verify by eye, not ground truth."""
    row = _acr_row(drugs)
    out = []
    for m in row.finditer(pdf_text(path)):
        a20, a50, a70 = float(m.group(2)), float(m.group(3)), float(m.group(4))
        if a20 >= a50 >= a70 and a20 <= 100:          # ACR monotonicity sanity filter
            out.append({"drug": m.group(1), "ACR20": a20, "ACR50": a50, "ACR70": a70})
    return out
