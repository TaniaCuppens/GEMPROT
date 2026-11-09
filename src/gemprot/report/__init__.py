"""Self-contained HTML report for the sample mode (no external resources)."""

from __future__ import annotations

import datetime as dt
import html
from pathlib import Path

from jinja2 import Environment, PackageLoader, select_autoescape

from .. import __version__
from ..annotate import Annotations, CopyVariant, VariantAnnotation, copy_variants
from ..classify import THREE, HaplotypeEffect, hgvs_c, parse_diff
from ..indiv import IndivResult, SampleResult
from ..phasing import SOURCE_LABELS
from ..scope import Outside

# protein track geometry (SVG viewBox units)
WIDTH, LEFT, RIGHT = 1000, 70, 20
TRUNCATING = {"stop_gained", "frameshift", "start_lost"}
KIND_LABEL = {
    "missense": "missense",
    "inframe_indel": "in-frame indel",
    "stop_gained": "stop gained",
    "frameshift": "frameshift",
    "start_lost": "start lost",
    "stop_lost": "stop lost",
    "masked": "masked in cis",
}
VARIANT_STATUS = {
    "applied": "isolated effect = effect on this copy",
    "masked": "no effect of its own: downstream of a truncation on the same copy",
    "combined": "combined with another variant on the same copy",
    "no protein effect": "synonymous",
    "outside CDS": "outside the coding sequence: GEMPROT does not evaluate it",
}
STATUS = {
    "good": ("✓", "No copy with a changed sequence"),
    "warning": ("!", "One copy with a changed sequence, one unchanged copy"),
    "critical": ("✕", "All copies have a changed sequence"),
    "serious": ("?", "Conclusion depends on the phase"),
}


def _x(pos: float, length: int) -> float:
    span = WIDTH - LEFT - RIGHT
    return LEFT + (min(max(pos, 1), length) - 1) / max(length - 1, 1) * span


def _status(altered: int, total: int) -> str:
    return "good" if altered == 0 else "critical" if altered == total else "warning"


def _num(v: float | None, fmt: str = ".2f") -> str:
    return "—" if v is None else format(v, fmt)


def _evidence(a: VariantAnnotation) -> str:
    """Short database summary: "ClinVar Pathogenic ★★ · AlphaMissense 0.87"."""
    parts = []
    if a.clinvar:
        parts.append(f"ClinVar {a.clinvar.significance} {'★' * a.clinvar.stars}".rstrip())
    if a.scores and a.scores.alphamissense is not None:
        parts.append(f"AlphaMissense {_num(a.scores.alphamissense)}")
    return " · ".join(parts)


def _variant_row(cv: CopyVariant, a: VariantAnnotation) -> dict:
    iso, sc, cl = a.isolated, a.scores or None, a.clinvar
    return {
        "label": a.label,
        "hgvsc": iso.hgvsc if iso else None,
        "isolated": (iso.hgvsp if iso and iso.hgvsp else (iso.consequence.replace("_", " ") if iso else "—")),
        "status": cv.status,
        "status_tip": VARIANT_STATUS[cv.status],
        "effect_on_copy": cv.effect_on_copy,
        "misleading": cv.misleading,
        "clinvar": cl and {"sig": cl.significance, "stars": cl.stars, "id": cl.variation_id,
                           "pathogenic": cl.pathogenic, "conditions": ", ".join(cl.conditions[:3])},
        "am": sc and sc.alphamissense is not None and {"score": _num(sc.alphamissense), "cls": sc.alphamissense_class},
        "revel": _num(sc.revel) if sc else "—",
        "cadd": _num(sc.cadd_phred, ".1f") if sc else "—",
    }


def _warnings(index: int, cvs: list[CopyVariant], ann: Annotations) -> list[str]:
    out = []
    for cv in cvs:
        a = ann.variants[cv.key]
        ev = _evidence(a)
        who = f"{a.label} ({ev})" if ev else a.label
        if cv.status == "masked":
            out.append(f"Copy {index}: {cv.effect_on_copy} masks {who}. "
                       f"Its annotation describes the variant alone, not this copy.")
        elif cv.status == "combined":
            what = f"into {cv.effect_on_copy}" if cv.effect_on_copy else "with another variant"
            out.append(f"Copy {index}: {who} is combined {what}. "
                       f"No database evaluates this combined change.")
        elif cv.status == "outside CDS" and a.clinvar and a.clinvar.pathogenic:
            out.append(f"Copy {index}: {who} is outside the coding sequence (splicing?). GEMPROT does not evaluate it.")
    return out


def _copy_view(copy: HaplotypeEffect, index: int, length: int, ref_cds: str,
               alleles: dict[int, int] | None = None, ann: Annotations | None = None) -> dict:
    cvs = copy_variants(copy, alleles, ann) if ann is not None and alleles is not None else []
    marks = []
    for c in copy.changes:
        ev = [_evidence(ann.variants[cv.key]) for cv in cvs
              if (cv.status == "applied" and ann.variants[cv.key].isolated.hgvsp
                  and ann.variants[cv.key].isolated.hgvsp.split("ext")[0] == c.hgvs.split("ext")[0])
              or cv.effect_on_copy == c.hgvs and cv.status == "combined"]
        ev = [e for e in ev if e]
        tip = f"{c.hgvs} — {KIND_LABEL.get(c.kind, c.kind)}" + (f" · {' ; '.join(ev)}" if ev else "")
        marks.append({"pos": c.position, "kind": c.kind, "label": c.hgvs.removeprefix("p."), "tip": tip})
    for d in copy.masked:
        pd = parse_diff(d)
        name = hgvs_c(d, ref_cds)
        marks.append({"pos": pd.codon, "kind": "masked", "label": name,
                      "tip": f"{name} (codon {pd.codon}) — no effect of its own: downstream of the "
                             f"truncation on the same copy"})
    marks.sort(key=lambda m: m["pos"])
    # labels: alternate the height when two marks are near
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
        "variants": [_variant_row(cv, ann.variants[cv.key]) for cv in cvs],
        "warnings": _warnings(index, cvs, ann) if ann is not None else [],
    }


def _sample_view(name: str, res: SampleResult, length: int, ref_cds: str,
                 ann: Annotations | None = None, outside: list[Outside] = ()) -> dict:
    scenarios = []
    for s in res.scenarios:
        v = s.verdict
        hyp = [f"{a} {r} {b}" for a, b, r in s.relations] + list(s.scenario.assumptions)
        scenarios.append({
            "title": f"Scenario {s.scenario.index + 1}" if len(res.scenarios) > 1 else None,
            "hypotheses": hyp,
            "status": _status(v.altered_copies, len(v.copies)),
            "summary": v.summary[:1].upper() + v.summary[1:],
            "copies": [_copy_view(c, i, length, ref_cds, al, ann)
                       for i, (c, al) in enumerate(zip(v.copies, s.scenario.copies), 1)],
        })
        scenarios[-1]["warnings"] = [w for c in scenarios[-1]["copies"] for w in c["warnings"]]
    level = "serious" if not res.conclusive else scenarios[0]["status"]
    sources, evidence, phase_notes = [], {}, []
    if res.phase:
        sources = res.phase.sources
        evidence = {(e.a, e.b): e for e in res.phase.pairs}
        phase_notes = res.phase.notes + [
            f"{e.label_a} / {e.label_b}: " + ", ".join(f"{SOURCE_LABELS[k]} {r}" for k, r in e.by_source.items())
            + f" — kept: {e.final} ({SOURCE_LABELS[e.retained]})"
            for e in res.phase.disagreements
        ]

    def pair(p) -> dict:
        e = evidence.get((min(p.a.record, p.b.record), max(p.a.record, p.b.record)))
        return {"a": p.a.label, "b": p.b.label, "relation": p.relation, "reason": p.reason,
                "by_source": [(e.by_source.get(k, "—") if e else "—") for k in sources],
                "status": e.status if e else ""}

    return {
        "name": name,
        "anchor": "s-" + "".join(ch if ch.isalnum() else "-" for ch in name),
        "status": level,
        "notes": res.audit.notes,
        "pairs": [pair(p) for p in res.audit.pairs],
        "sources": [SOURCE_LABELS[k] for k in sources],
        "phase_notes": phase_notes,
        "units": res.audit.units,
        "outside": [{"label": site.record.label(1), "consequence": site.consequence, "splice": site.splice}
                    for site in outside if site.record.genotypes[name].carried],
        "scenarios": scenarios,
        "n_warnings": max(len(sc["warnings"]) for sc in scenarios),
    }


def _fit(text: str, n: int) -> str:
    """Cut a label to n characters (empty if fewer than 5 characters fit)."""
    if n < 5:
        return ""
    return text if len(text) <= n else text[: n - 1].rstrip() + "…"


def _domains(ann: Annotations | None, length: int) -> list[dict]:
    if ann is None or ann.protein is None:
        return []
    out = []
    for i, d in enumerate(ann.protein.domains):
        x1, x2 = _x(d.start, length), _x(min(d.end, length), length)
        out.append({"x": round(x1, 1), "w": round(max(x2 - x1, 2), 1), "name": d.name, "acc": d.accession,
                    "tip": f"{d.name} ({d.accession}): {d.start}–{d.end}",
                    "label": _fit(d.name, int((x2 - x1 - 8) / 6.5)), "alt": i % 2})
    return out


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
        samples=[_sample_view(n, r, length, result.ref_cds, result.annotations, result.outside)
                 for n, r in result.samples.items()],
        annotated=result.annotations is not None,
        protein=result.annotations.protein if result.annotations else None,
        ann_notes=result.annotations.notes if result.annotations else [],
        domains=_domains(result.annotations, length),
        ticks=_ticks(length),
        geom={"width": WIDTH, "left": LEFT, "right": RIGHT, "line_end": WIDTH - RIGHT},
        status=STATUS,
        truncating=TRUNCATING,
        three=THREE,
        escape=html.escape,
    )
    out.write_text(page, encoding="utf-8")
    return out
