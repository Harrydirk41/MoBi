"""Render an observed-vs-simulated fit overlay as a PNG the MODEL can see.

This is the "fit-vision" channel: after a fit, the agent is handed a semilog
plot of its own predicted curves against the observed data, so it can reason
about SHAPE the scalar GMFE hides - a missing distribution phase, a wrong
terminal slope, an over/undershot Cmax, a shifted tmax. PBPK curves span orders
of magnitude, so the y-axis is log10 (the standard concentration-time view).

Deliberately dependency-light: Pillow + stdlib only (no matplotlib, no numpy),
because the container has Pillow but not matplotlib, and a fit overlay is a few
polylines and markers - it does not need a plotting stack. Returns raw PNG bytes
(or ``None`` when there is nothing positive to plot); the caller base64-encodes
it into a tool_result image block.
"""

from __future__ import annotations

import base64
import io
import math

_PALETTE = [
    (37, 99, 235), (220, 38, 38), (5, 150, 105), (217, 119, 6),
    (124, 58, 237), (13, 148, 136), (190, 24, 93), (100, 116, 139),
]
_BG = (255, 255, 255)
_AXIS = (55, 65, 81)
_GRID = (226, 232, 240)
_INK = (17, 24, 39)


def _font(size: int = 12):
    """A TrueType font if one is on the box, else Pillow's bitmap default (which
    works fully headless and ignores the size)."""
    from PIL import ImageFont
    for path in ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
                 "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"):
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    try:
        return ImageFont.load_default()
    except Exception:                                    # noqa: BLE001
        return None


def _finite_xy(pts):
    """Keep (t, c) pairs with a finite time and a strictly positive conc (log-y)."""
    out = []
    for p in pts or []:
        try:
            t, c = float(p[0]), float(p[1])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(t) and math.isfinite(c) and c > 0:
            out.append((t, c))
    return out


def _nice_pow10_ticks(lo: float, hi: float):
    """Decade tick values spanning [lo, hi] in log10 concentration."""
    a, b = math.floor(lo), math.ceil(hi)
    return [10.0 ** k for k in range(int(a), int(b) + 1)]


def _fmt_conc(v: float) -> str:
    if v >= 1000 or (0 < v < 0.01):
        return f"{v:.0e}".replace("e-0", "e-").replace("e+0", "e")
    if v >= 1:
        return f"{v:g}"
    return f"{v:g}"


def render_fit_png(series: list[dict], title: str = "", *,
                   width: int = 760, height: int = 520,
                   max_series: int = 8) -> "bytes | None":
    """Draw observed (markers) vs simulated (lines) on a shared semilog-y axis.

    ``series`` entries: ``{study, route, dose, observed:[[t,c],...],
    simulated:[[t,c],...]}``. Returns PNG bytes, or ``None`` if nothing is
    plottable (e.g. all concentrations non-positive)."""
    from PIL import Image, ImageDraw

    ser = [s for s in (series or []) if s][:max_series]
    # collect finite/positive points to set the data ranges
    prepared = []
    tmin = tmax = None
    lymin = lymax = None
    for i, s in enumerate(ser):
        obs = _finite_xy(s.get("observed"))
        sim = _finite_xy(s.get("simulated"))
        if not obs and not sim:
            continue
        for t, c in obs + sim:
            tmin = t if tmin is None else min(tmin, t)
            tmax = t if tmax is None else max(tmax, t)
            ly = math.log10(c)
            lymin = ly if lymin is None else min(lymin, ly)
            lymax = ly if lymax is None else max(lymax, ly)
        prepared.append((i, s, obs, sim))
    if not prepared or tmin is None or lymin is None:
        return None
    if tmax <= tmin:
        tmax = tmin + 1.0
    if lymax - lymin < 0.5:                               # pad a near-flat range
        lymin -= 0.5
        lymax += 0.5

    ml, mr, mt, mb = 66, 16, 30 if title else 14, 40
    legend_h = 16 * min(len(prepared), max_series) + 8
    plot_w = width - ml - mr
    plot_h = height - mt - mb - legend_h
    img = Image.new("RGB", (width, height), _BG)
    d = ImageDraw.Draw(img)
    f = _font(12)
    fb = _font(13)

    def sx(t):
        return ml + (t - tmin) / (tmax - tmin) * plot_w

    def sy(c):
        ly = math.log10(c)
        return mt + (lymax - ly) / (lymax - lymin) * plot_h

    if title:
        d.text((ml, 8), title[:90], fill=_INK, font=fb)

    # y grid + labels (decades)
    for v in _nice_pow10_ticks(lymin, lymax):
        ly = math.log10(v)
        if ly < lymin - 1e-9 or ly > lymax + 1e-9:
            continue
        y = mt + (lymax - ly) / (lymax - lymin) * plot_h
        d.line([(ml, y), (ml + plot_w, y)], fill=_GRID, width=1)
        lbl = _fmt_conc(v)
        d.text((ml - 6 - d.textlength(lbl, font=f), y - 6), lbl, fill=_AXIS, font=f)

    # x ticks (nice-ish: 5 intervals)
    for k in range(6):
        t = tmin + (tmax - tmin) * k / 5.0
        x = sx(t)
        d.line([(x, mt), (x, mt + plot_h)], fill=_GRID, width=1)
        lbl = f"{t:g}"
        d.text((x - d.textlength(lbl, font=f) / 2, mt + plot_h + 5), lbl, fill=_AXIS, font=f)

    # axis box + titles
    d.rectangle([ml, mt, ml + plot_w, mt + plot_h], outline=_AXIS, width=1)
    d.text((ml + plot_w / 2 - 18, mt + plot_h + 22), "time (h)", fill=_AXIS, font=f)

    # series: observed markers + simulated polyline, same color
    for slot, (idx, s, obs, sim) in enumerate(prepared):
        col = _PALETTE[idx % len(_PALETTE)]
        if len(sim) >= 2:
            d.line([(sx(t), sy(c)) for t, c in sim], fill=col, width=2, joint="curve")
        for t, c in obs:
            x, y = sx(t), sy(c)
            d.ellipse([x - 3, y - 3, x + 3, y + 3], fill=col, outline=_BG)
        # legend row
        ly0 = mt + plot_h + mb + slot * 16
        d.line([(ml, ly0 + 7), (ml + 22, ly0 + 7)], fill=col, width=2)
        d.ellipse([ml + 9, ly0 + 4, ml + 15, ly0 + 10], fill=col, outline=_BG)
        tag = str(s.get("study") or f"study {idx+1}")
        extra = " · ".join(x for x in (str(s.get("route") or ""), str(s.get("dose") or "")) if x)
        lbl = tag + (f"  ({extra})" if extra else "")
        d.text((ml + 28, ly0), lbl[:80], fill=_INK, font=f)

    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


def render_fit_b64(series: list[dict], title: str = "", **kw) -> "dict | None":
    """render_fit_png + base64, packaged as an image descriptor ready for a
    tool_result image block: ``{media_type, data}``. ``None`` when unplottable."""
    png = render_fit_png(series, title, **kw)
    if png is None:
        return None
    return {"media_type": "image/png", "data": base64.b64encode(png).decode("ascii")}
