"""Rapport HTML autonome du mode individu (aucune ressource externe)."""

from __future__ import annotations

import datetime as dt
import html
from pathlib import Path

from jinja2 import Environment, PackageLoader, select_autoescape

from .. import __version__
from ..classify import THREE, HaplotypeEffect, hgvs_c, parse_diff
from ..indiv import IndivResult, SampleResult

# géométrie des pistes protéiques (unités du viewBox SVG)
WIDTH, LEFT, RIGHT = 1000, 70, 20
TRUNCATING = {"stop_gained", "frameshift", "start_lost"}
KIND_LABEL = {
    "missense": "faux-sens",
    "inframe_indel": "indel en phase",
    "stop_gained": "stop gagné",
    "frameshift": "décalage du cadre",
    "start_lost": "start perdu",
    "stop_lost": "stop perdu",
    "masked": "masqué en cis",
}
STATUS = {
    "good": ("✓", "Aucune copie altérée"),
    "warning": ("!", "Une copie altérée, une copie normale"),
    "critical": ("✕", "Toutes les copies altérées"),
    "serious": ("?", "Conclusion dépendante de la phase"),
}


def _x(pos: float, length: int) -> float:
    span = WIDTH - LEFT - RIGHT
    return LEFT + (min(max(pos, 1), length) - 1) / max(length - 1, 1) * span


def _status(altered: int, total: int) -> str:
    return "good" if altered == 0 else "critical" if altered == total else "warning"


def _copy_view(copy: HaplotypeEffect, index: int, length: int, ref_cds: str) -> dict:
    marks = []
    for c in copy.changes:
        marks.append({"pos": c.position, "kind": c.kind, "label": c.hgvs.removeprefix("p."),
                      "tip": f"{c.hgvs} — {KIND_LABEL.get(c.kind, c.kind)}"})
    for d in copy.masked:
        pd = parse_diff(d)
        name = hgvs_c(d, ref_cds)
        marks.append({"pos": pd.codon, "kind": "masked", "label": name,
                      "tip": f"{name} (codon {pd.codon}) — sans effet propre : en aval de la "
                             f"troncature sur la même copie"})
    marks.sort(key=lambda m: m["pos"])
    # étiquettes : on alterne la hauteur quand deux marques sont proches
    last_x, level = -1e9, 0
    for m in marks:
        m["x"] = round(_x(m["pos"], length), 1)
        level = (level + 1) % 2 if m["x"] - last_x < 110 else 0
        m["level"], last_x = level, m["x"]
    end = copy.truncated_at if copy.truncated_at is not None else length
    return {
        "index": index,
        "label": copy.label,
        "altered": copy.altered,
        "truncated_at": copy.truncated_at,
        "compensated": copy.compensated_frameshift,
        "end_x": round(_x(max(end, 1), length), 1) if end else LEFT,
        "marks": marks,
        "changes": [{"hgvs": c.hgvs, "kind": KIND_LABEL.get(c.kind, c.kind), "pos": c.position}
                    for c in copy.changes],
        "masked": [hgvs_c(d, ref_cds) for d in copy.masked],
    }


def _sample_view(name: str, res: SampleResult, length: int, ref_cds: str) -> dict:
    scenarios = []
    for s in res.scenarios:
        v = s.verdict
        hyp = [f"{a} {r} {b}" for a, b, r in s.relations] + list(s.scenario.assumptions)
        scenarios.append({
            "title": f"Scénario {s.scenario.index + 1}" if len(res.scenarios) > 1 else None,
            "hypotheses": hyp,
            "status": _status(v.altered_copies, len(v.copies)),
            "summary": v.summary[:1].upper() + v.summary[1:],
            "copies": [_copy_view(c, i, length, ref_cds) for i, c in enumerate(v.copies, 1)],
        })
    level = "serious" if not res.conclusive else scenarios[0]["status"]
    return {
        "name": name,
        "anchor": "s-" + "".join(ch if ch.isalnum() else "-" for ch in name),
        "status": level,
        "notes": res.audit.notes,
        "pairs": [{"a": p.a.label, "b": p.b.label, "relation": p.relation, "reason": p.reason}
                  for p in res.audit.pairs],
        "units": res.audit.units,
        "scenarios": scenarios,
    }


def _ticks(length: int) -> list[dict]:
    step = next(s for s in (10, 25, 50, 100, 200, 250, 500, 1000, 2000) if length / s <= 10)
    return [{"pos": p, "x": round(_x(p, length), 1)} for p in [1, *range(step, length + 1, step)]]


def render_indiv(result: IndivResult, out: Path, source: str = "") -> Path:
    from ..classify import translate

    ref_protein = translate(result.ref_cds)
    length = len(ref_protein.rstrip("*"))
    env = Environment(loader=PackageLoader("gemprot.report", "templates"),
                      autoescape=select_autoescape(["html", "j2"]))
    page = env.get_template("indiv.html.j2").render(
        transcript=result.transcript,
        length=length,
        source=source,
        version=__version__,
        date=dt.date.today().isoformat(),
        samples=[_sample_view(n, r, length, result.ref_cds) for n, r in result.samples.items()],
        ticks=_ticks(length),
        geom={"width": WIDTH, "left": LEFT, "right": RIGHT, "line_end": WIDTH - RIGHT},
        status=STATUS,
        truncating=TRUNCATING,
        three=THREE,
        escape=html.escape,
    )
    out.write_text(page, encoding="utf-8")
    return out
