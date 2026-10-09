"""Mode individu : audit de la phase → scénarios → haplo → effet de chaque copie."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

from .annotate import Annotations, annotate
from .classify import SampleVerdict, classify_pseudo_samples, translate
from .config import Config
from .haplotypes import TranscriptHaplotypes, parse_haplo_json, run_haplo
from .phase import SampleAudit, Scenario, Vcf, audit_sample, read_vcf, scenarios, write_scenario_vcf


@dataclass
class ScenarioResult:
    scenario: Scenario
    verdict: SampleVerdict
    relations: list[tuple[str, str, str]]  # (variant a, variant b, cis/trans) pour les paires de phase inconnue


@dataclass
class SampleResult:
    audit: SampleAudit
    scenarios: list[ScenarioResult]

    @property
    def conclusive(self) -> bool:
        """Toutes les hypothèses de phase mènent au même nombre de copies altérées."""
        return len({s.verdict.altered_copies for s in self.scenarios}) == 1


@dataclass
class IndivResult:
    transcript: str
    ref_cds: str
    samples: dict[str, SampleResult]
    annotations: Annotations | None = None
    vcf: Vcf | None = field(default=None, repr=False)


def prepare(vcf_path: Path, samples: Sequence[str] = ()) -> tuple:
    vcf = read_vcf(vcf_path)
    names = list(samples) or vcf.samples
    audits = {s: audit_sample(vcf, s) for s in names}
    plan = {s: scenarios(a) for s, a in audits.items()}
    return vcf, audits, plan


def assemble(
    tx: TranscriptHaplotypes,
    audits: dict[str, SampleAudit],
    plan: dict[str, list[Scenario]],
) -> dict[str, SampleResult]:
    verdicts = classify_pseudo_samples(tx)
    results = {}
    for sample, audit in audits.items():
        unknown = [p for p in audit.pairs if p.relation == "inconnue"]
        results[sample] = SampleResult(
            audit,
            [
                ScenarioResult(
                    sc,
                    verdicts[sample][sc.index],
                    [(p.a.label, p.b.label, sc.relation(p.a, p.b)) for p in unknown],
                )
                for sc in plan[sample]
            ],
        )
    return results


def run(
    vcf_path: Path,
    transcript: str,
    cfg: Config,
    workdir: Path,
    samples: Sequence[str] = (),
    with_annotations: bool = True,
    online: bool = True,
) -> IndivResult:
    vcf, audits, plan = prepare(vcf_path, samples)
    pseudo = write_scenario_vcf(vcf, plan, workdir / "scenarios.vcf")
    txs = {t.transcript_id: t for t in parse_haplo_json(run_haplo(pseudo, workdir / "haplo.json", cfg))}
    if transcript not in txs:
        known = ", ".join(sorted(txs)) or "aucun"
        raise KeyError(f"{transcript} absent de la sortie haplo (transcrits : {known})")
    tx = txs[transcript]
    ann = None
    if with_annotations:
        length = len(translate(tx.ref_cds).rstrip("*"))
        ann = annotate(vcf, transcript, cfg, workdir, length, online)
    return IndivResult(transcript, tx.ref_cds, assemble(tx, audits, plan), ann, vcf)
