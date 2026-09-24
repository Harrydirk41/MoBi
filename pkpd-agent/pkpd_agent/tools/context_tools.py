"""Context-data-lake tools: the agent builds its OWN givens from raw materials.

The realistic modeler workflow (and the honest analogue of a literature-retrieval
agent) is: you are handed a report and a set of data files, and you extract the
physicochemical / in-vitro inputs yourself, recording each with its source. That
is what these tools support -- instead of receiving a pre-digested
``literature_physicochemical`` list, the agent:

  * ``osp_read_report``  - reads the handed report (background prose, the
    physicochemical table, the enzyme identity, the clinical study descriptions);
  * ``osp_list_studies`` / ``osp_read_study`` - inspects the raw data files;
  * ``osp_record_givens`` - declares the parameters it extracted, each with a
    value, unit, source citation, and a provenance tag (``given`` = read from the
    report, ``judged`` = inferred/assumed). The recorded list is written into the
    task's ``literature_physicochemical`` so the rest of the pipeline (measured
    ranges, the optimizer's fix-vs-fit split) uses the agent's OWN givens.

Extraction is confined to the handed materials -- there is no open-web search --
so the leave-one-out wall holds: the agent can never pull the target compound's
own published model off the internet.
"""
from __future__ import annotations

import csv
import glob
import json
import os

from .registry import Tool, ToolRegistry, ToolResult

_PROVENANCE = ("given", "judged")


def _manifest(data_dir: str) -> list[dict]:
    man = os.path.join(data_dir, "_manifest.json")
    if os.path.isfile(man):
        try:
            return json.load(open(man, encoding="utf-8"))
        except (OSError, ValueError):
            pass
    return [{"file": os.path.basename(f)}
            for f in sorted(glob.glob(os.path.join(data_dir, "*.csv")))]


def _read_points(path: str, cap: int = 500) -> list[list]:
    pts = []
    try:
        with open(path, encoding="utf-8") as fh:
            for i, row in enumerate(csv.reader(fh)):
                if i == 0 or len(row) < 2:
                    continue
                try:
                    pts.append([float(row[0]), float(row[1])])
                except ValueError:
                    continue
    except OSError:
        return []
    return pts[:cap]


def _match_study(entries: list[dict], study: str) -> dict | None:
    s = (study or "").lower().strip()
    for e in entries:
        if s and (s == (e.get("study") or "").lower()
                  or s in e.get("file", "").lower()):
            return e
    return None


def register_context_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {report_path, data_dir, input}. Any of report_path/data_dir may be
    None (then the corresponding tool reports that nothing was handed)."""
    report_path = ctx.get("report_path")
    data_dir = ctx.get("data_dir")
    # keep the SAME object the caller passed (even an empty {}), so record_givens
    # writes back into the shared task input that the loop tools also read.
    inp: dict = ctx["input"] if ctx.get("input") is not None else {}

    # -- read the handed report ----------------------------------------- #
    def read_report(args: dict, session) -> ToolResult:
        if not report_path or not os.path.isfile(report_path):
            return ToolResult.error(
                "No report was handed for this task (givens may be pre-digested; "
                "call osp_inspect instead).")
        md = open(report_path, encoding="utf-8").read()
        section = (args.get("section") or "").strip().lower()
        if section:
            # return only paragraphs/tables under the first heading that matches
            lines, out, keep = md.splitlines(), [], False
            for ln in lines:
                h = ln.lstrip().startswith("#")
                if h:
                    keep = section in ln.lower()
                if keep:
                    out.append(ln)
            body = "\n".join(out).strip() or md
        else:
            body = md
        return ToolResult.success(
            "the report you were handed - extract the physicochemical/in-vitro "
            "givens (and the metabolizing enzyme identity) yourself, then call "
            "osp_record_givens. The chosen methods and fitted values are redacted.",
            report_markdown=body[:20000],
            truncated=len(body) > 20000)

    # -- inspect the raw data files ------------------------------------- #
    def list_studies(args: dict, session) -> ToolResult:
        if not data_dir or not os.path.isdir(data_dir):
            return ToolResult.error("No data files were handed for this task.")
        entries = _manifest(data_dir)
        out = []
        for e in entries:
            p = os.path.join(data_dir, e.get("file", ""))
            pts = _read_points(p, cap=10000)
            if not pts:
                continue
            ts = [t for t, _ in pts]
            out.append({"study": e.get("study") or e["file"],
                        "route": e.get("route", ""), "dose": e.get("dose", ""),
                        "n_points": len(pts),
                        "t_range": [min(ts), max(ts)] if ts else None,
                        "file": e["file"]})
        return ToolResult.success(
            f"{len(out)} clinical study curve(s) were handed as raw "
            "time,concentration files. Use osp_read_study to see the points; the "
            "route and dose define each simulated administration.",
            studies=out)

    def read_study(args: dict, session) -> ToolResult:
        if not data_dir or not os.path.isdir(data_dir):
            return ToolResult.error("No data files were handed for this task.")
        entries = _manifest(data_dir)
        e = _match_study(entries, args.get("study", ""))
        if not e:
            return ToolResult.error(
                f"no study matching {args.get('study')!r}; call osp_list_studies.")
        pts = _read_points(os.path.join(data_dir, e["file"]))
        return ToolResult.success(
            f"{e.get('study')} - {e.get('route','?')} {e.get('dose','?')} "
            f"({len(pts)} points, time[h], concentration)",
            study=e.get("study"), route=e.get("route"), dose=e.get("dose"),
            points=pts)

    # -- the agent's own parameter table (provenance-tagged) ------------ #
    def record_givens(args: dict, session) -> ToolResult:
        raw = args.get("givens")
        if not isinstance(raw, list) or not raw:
            return ToolResult.error(
                "pass givens=[{parameter, value, unit, source, provenance}] - the "
                "physicochemical/in-vitro values you extracted from the report.")
        cleaned, problems = [], []
        for g in raw:
            if not isinstance(g, dict) or not g.get("parameter"):
                problems.append(f"skipped malformed entry: {g!r}")
                continue
            prov = (g.get("provenance") or "given").lower()
            if prov not in _PROVENANCE:
                prov = "given"
            entry = {"parameter": g.get("parameter"),
                     "value": g.get("value"), "unit": g.get("unit", ""),
                     "source": g.get("source", ""), "provenance": prov}
            if entry["value"] is None and not entry["source"]:
                problems.append(f"'{entry['parameter']}' has no value and no source")
            cleaned.append(entry)
        # write the agent's givens into the task so measured-range / fix-vs-fit
        # logic downstream uses the SAME list a pre-digested task would carry.
        gd = inp.setdefault("given_data", {})
        gd["literature_physicochemical"] = cleaned
        session.put("extracted_givens", cleaned)
        n_given = sum(1 for e in cleaned if e["provenance"] == "given")
        return ToolResult.success(
            f"recorded {len(cleaned)} givens ({n_given} read from the report, "
            f"{len(cleaned) - n_given} judged). These now feed the model's "
            "fix-vs-fit split. Revise by calling osp_record_givens again.",
            recorded=cleaned, warnings=problems)

    registry.register(Tool(
        name="osp_read_report",
        description=(
            "READ the report you were handed for THIS compound (background, the "
            "physicochemical/in-vitro data table, the metabolizing enzyme named in "
            "the prose, and the clinical study descriptions). The modeling "
            "strategy, chosen methods and fitted values are redacted - you extract "
            "the inputs yourself. Optional 'section' returns just one heading "
            "(e.g. 'in vitro', 'data', 'background')."),
        input_schema={"type": "object", "properties": {
            "section": {"type": "string",
                        "description": "optional heading substring to return only that section"}}},
        handler=read_report, phase="observe"))

    registry.register(Tool(
        name="osp_list_studies",
        description=(
            "LIST the raw clinical data files you were handed - each is a "
            "digitized time,concentration curve with a route and dose. Use these "
            "to know which administrations the model must reproduce."),
        input_schema={"type": "object", "properties": {}},
        handler=list_studies, phase="observe"))

    registry.register(Tool(
        name="osp_read_study",
        description=(
            "READ one study's raw time,concentration points (by study name from "
            "osp_list_studies)."),
        input_schema={"type": "object", "properties": {
            "study": {"type": "string"}}, "required": ["study"]},
        handler=read_study, phase="observe"))

    registry.register(Tool(
        name="osp_record_givens",
        description=(
            "RECORD the physicochemical/in-vitro givens you extracted from the "
            "report, as a provenance-tagged table. Each entry: {parameter, value, "
            "unit, source (citation), provenance ('given'=read from the report, "
            "'judged'=inferred/assumed}}. This becomes the model's known-input "
            "list, so record measured values (lipophilicity, fu, pKa, solubility, "
            "in-vitro Km/Vmax) here BEFORE optimizing - the optimizer fixes what "
            "you mark measured and fits only the rest."),
        input_schema={"type": "object", "properties": {
            "givens": {"type": "array", "items": {"type": "object", "properties": {
                "parameter": {"type": "string"}, "value": {"type": ["number", "string", "null"]},
                "unit": {"type": "string"}, "source": {"type": "string"},
                "provenance": {"type": "string", "enum": list(_PROVENANCE)}},
                "required": ["parameter"]}}},
            "required": ["givens"]},
        handler=record_givens, phase="observe"))
