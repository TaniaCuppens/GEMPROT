"""Annotations des variants (ClinVar, AlphaMissense, REVEL, CADD) et domaines de la protéine.

Les bases de données évaluent un variant **isolé**. GEMPROT dit si cet effet isolé est bien celui
que produit la copie qui le porte : un variant peut être masqué (en aval d'une troncature en cis)
ou combiné à un autre (même codon, décalage compensé). Dans ces deux cas, son score ne décrit pas
la protéine réellement produite, et le rapport le signale.

Sources :
- effet isolé de chaque allèle : `vep` (cache RefSeq, transcrit demandé) ;
- ClinVar et dbNSFP : `tabix` sur les fichiers partagés (score du transcrit Ensembl MANE) ;
- domaines : UniProt (via l'identifiant RefSeq NP_) puis InterPro (Pfam), mis en cache (JSON).
"""

from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .classify import HaplotypeEffect, ProteinChange
from .config import Config
from .phase import Vcf
from .tools import ToolError, run_tool

VariantKey = tuple[int, int]  # (indice du variant dans le VCF, allèle)

TRUNCATING = ("stop_gained", "frameshift", "start_lost")

# Seuils AlphaMissense (Cheng et al., Science 2023)
AM_CLASSES = ((0.564, "probablement pathogène"), (0.34, "ambigu"), (0.0, "probablement bénin"))
CLINVAR_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_single_submitter": 1,
}


# --- modèles ------------------------------------------------------------------------------


@dataclass(frozen=True)
class IsolatedEffect:
    """Effet du variant seul sur le transcrit (VEP)."""

    hgvsc: str | None
    hgvsp: str | None  # notation GEMPROT : p.Trp7*, p.Pro137Hisfs*8 ; None si non codant
    consequence: str
    protein_start: int | None


@dataclass(frozen=True)
class ClinVarEntry:
    variation_id: str
    significance: str
    review: str
    conditions: tuple[str, ...]

    @property
    def stars(self) -> int:
        return CLINVAR_STARS.get(self.review, 0)

    @property
    def pathogenic(self) -> bool:
        s = self.significance.lower()
        return "pathogenic" in s and "conflicting" not in s


@dataclass(frozen=True)
class Scores:
    alphamissense: float | None = None
    revel: float | None = None
    cadd_phred: float | None = None

    @property
    def alphamissense_class(self) -> str | None:
        if self.alphamissense is None:
            return None
        return next(label for cut, label in AM_CLASSES if self.alphamissense >= cut)


@dataclass
class VariantAnnotation:
    label: str
    isolated: IsolatedEffect | None = None
    clinvar: ClinVarEntry | None = None
    scores: Scores | None = None


@dataclass(frozen=True)
class Domain:
    accession: str
    name: str
    start: int
    end: int


@dataclass
class ProteinInfo:
    protein_id: str  # RefSeq NP_
    uniprot: str | None = None
    length: int | None = None
    domains: list[Domain] = field(default_factory=list)


@dataclass
class Annotations:
    variants: dict[VariantKey, VariantAnnotation]
    protein: ProteinInfo | None = None
    notes: list[str] = field(default_factory=list)


# --- effet isolé : VEP ---------------------------------------------------------------------


def sites_id(key: VariantKey) -> str:
    return f"v{key[0]}_{key[1]}"


def write_sites_vcf(vcf: Vcf, out: Path) -> Path:
    """Un allèle par ligne, sans échantillons, ID = clé GEMPROT (pour relire la sortie VEP)."""
    with open(out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, rec in enumerate(vcf.records):
            for a, alt in enumerate(rec.alts, 1):
                if alt in ("*", "."):
                    continue
                fh.write(f"{rec.chrom}\t{rec.pos}\t{sites_id((i, a))}\t{rec.ref}\t{alt}\t.\t.\t.\n")
    return out


def run_vep(sites: Path, out_json: Path, cfg: Config) -> Path:
    if cfg.fasta is None:
        raise ToolError("vep a besoin de la FASTA (champ `fasta` de la config)")
    cmd = [
        "vep", "--cache", "--offline",
        *(["--refseq"] if cfg.vep.cache_type == "refseq" else []),
        "--dir_cache", cfg.vep.cache_dir,
        "--cache_version", str(cfg.vep.cache_version),
        "--assembly", cfg.assembly,
        "--fasta", cfg.fasta,
        "--hgvs", "--protein", "--mane", "--json", "--no_stats", "--force_overwrite",
        "-i", sites, "-o", out_json,
    ]
    binds = (*cfg.binds, cfg.vep.cache_dir.resolve(), cfg.fasta.parent, sites.parent, out_json.parent)
    run_tool(cmd, cfg.vep.container, binds)
    return out_json


def short_hgvsp(hgvsp: str | None) -> str | None:
    """`NP_000302.1:p.Trp7Ter` → `p.Trp7*` (notation du reste de GEMPROT)."""
    if not hgvsp:
        return None
    p = hgvsp.split(":", 1)[-1].replace("%3D", "=")
    return re.sub(r"(?<!p\.)Ter", "*", p)  # p.Ter254Argext*16 : le stop d'origine garde Ter


@dataclass
class VepResult:
    effects: dict[str, IsolatedEffect]  # par ID du VCF de sites
    protein_id: str | None = None
    mane_ensembl: str | None = None  # ENST MANE correspondant (pour dbNSFP)


def parse_vep(path: Path, transcript: str) -> VepResult:
    res = VepResult({})
    with open(path) as fh:
        for line in fh:
            v = json.loads(line)
            for t in v.get("transcript_consequences", []):
                if t.get("transcript_id") != transcript:
                    continue
                res.effects[v["id"]] = IsolatedEffect(
                    hgvsc=t["hgvsc"].split(":", 1)[-1] if t.get("hgvsc") else None,
                    hgvsp=short_hgvsp(t.get("hgvsp")),
                    consequence=t["consequence_terms"][0],
                    protein_start=t.get("protein_start"),
                )
                res.protein_id = res.protein_id or t.get("protein_id")
                res.mane_ensembl = res.mane_ensembl or t.get("mane_select")
    return res


# --- ClinVar et dbNSFP : tabix ----------------------------------------------------------------


def _bare(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def tabix(path: Path, regions: Iterable[str], cfg: Config, header: bool = False) -> list[str]:
    cmd = ["tabix", *(["-H"] if header else []), path, *([] if header else regions)]
    out = run_tool(cmd, cfg.vep.container, (*cfg.binds, path.parent)).stdout
    return [line for line in out.splitlines() if line]


def _info(field_str: str) -> dict[str, str]:
    return dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in field_str.split(";"))


def parse_clinvar(lines: Iterable[str]) -> dict[tuple[str, int, str, str], ClinVarEntry]:
    out = {}
    for line in lines:
        if line.startswith("#"):
            continue
        chrom, pos, vid, ref, alt, _, _, info = line.split("\t")[:8]
        f = _info(info)
        if "CLNSIG" not in f:
            continue
        out[(_bare(chrom), int(pos), ref, alt)] = ClinVarEntry(
            variation_id=vid,
            significance=f["CLNSIG"].replace("_", " "),
            review=f.get("CLNREVSTAT", ""),
            conditions=tuple(c.replace("_", " ") for c in f.get("CLNDN", "").split("|")
                             if c and c not in ("not_provided", "not_specified")),
        )
    return out


def _float(value: str) -> float | None:
    return None if value in ("", ".") else float(value)


def parse_dbnsfp(header: str, lines: Iterable[str], mane_ensembl: str | None) -> dict[tuple, Scores]:
    """Score du transcrit MANE : colonne `Ensembl_transcriptid`, sinon colonne `MANE` = Select."""
    cols = {name: i for i, name in enumerate(header.lstrip("#").split("\t"))}
    enst = (mane_ensembl or "").split(".")[0]
    out = {}
    for line in lines:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        txs = f[cols["Ensembl_transcriptid"]].split(";")
        mane = f[cols["MANE"]].split(";")
        idx = txs.index(enst) if enst in txs else (mane.index("Select") if "Select" in mane else None)
        if idx is None:
            continue

        def pick(col: str) -> float | None:
            if col not in cols:  # colonne absente (ex. fixtures publiques sans REVEL/CADD)
                return None
            vals = f[cols[col]].split(";")
            return _float(vals[idx] if len(vals) > idx else vals[0])

        key = (_bare(f[cols["chr"]]), int(f[cols["pos(1-based)"]]), f[cols["ref"]], f[cols["alt"]])
        out[key] = Scores(pick("AlphaMissense_score"), pick("REVEL_score"), pick("CADD_phred"))
    return out


# --- domaines : UniProt + InterPro (cache JSON par protéine) ------------------------------------------


def _get_json(url: str, timeout: int = 30):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r) if r.status == 200 else None


def _uniprot_for(protein_id: str) -> tuple[str, int] | None:
    q = urllib.parse.quote(f"xref:refseq-{protein_id} AND reviewed:true")
    d = _get_json(f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields=accession,length&format=json")
    hits = (d or {}).get("results", [])
    return (hits[0]["primaryAccession"], hits[0]["sequence"]["length"]) if hits else None


def _pfam_domains(uniprot: str) -> list[Domain]:
    d = _get_json(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{uniprot}")
    out = []
    for r in (d or {}).get("results", []):
        m = r["metadata"]
        for p in r["proteins"]:
            for loc in p["entry_protein_locations"]:
                for frag in loc["fragments"]:
                    out.append(Domain(m["accession"], m["name"], int(frag["start"]), int(frag["end"])))
    return sorted(out, key=lambda x: x.start)


def protein_info(protein_id: str, cache_dir: Path | None, online: bool = True) -> ProteinInfo | None:
    """UniProt (revu) correspondant au NP_ RefSeq, et ses domaines Pfam ; None si indisponible.

    Cache : un fichier JSON par protéine (pas de SQLite : le home est souvent sur NFS).
    """
    cached = cache_dir / f"{protein_id}.json" if cache_dir else None
    if cached and cached.exists():
        return protein_info_from_json(json.loads(cached.read_text()))
    if not online:
        return None
    hit = _uniprot_for(protein_id)
    info = ProteinInfo(protein_id)
    if hit:
        info.uniprot, info.length = hit
        info.domains = _pfam_domains(info.uniprot)
    if cached:
        cached.parent.mkdir(parents=True, exist_ok=True)
        tmp = cached.with_suffix(".tmp")
        tmp.write_text(json.dumps(protein_info_to_json(info), indent=1))
        tmp.replace(cached)
    return info


def protein_info_to_json(info: ProteinInfo) -> dict:
    return {"protein_id": info.protein_id, "uniprot": info.uniprot, "length": info.length,
            "domains": [d.__dict__ for d in info.domains]}


def protein_info_from_json(d: dict) -> ProteinInfo:
    return ProteinInfo(d["protein_id"], d.get("uniprot"), d.get("length"), [Domain(**x) for x in d.get("domains", [])])


# --- assemblage ---------------------------------------------------------------------------


def combine(
    vcf: Vcf,
    vep: VepResult,
    clinvar: dict[tuple, ClinVarEntry],
    scores: dict[tuple, Scores],
    protein: ProteinInfo | None,
    protein_length: int | None = None,
) -> Annotations:
    """Partie pure : réunit les sources pour chaque allèle du VCF."""
    ann = Annotations({}, protein)
    for i, rec in enumerate(vcf.records):
        for a, alt in enumerate(rec.alts, 1):
            key = (_bare(rec.chrom), rec.pos, rec.ref, alt)
            ann.variants[(i, a)] = VariantAnnotation(
                label=rec.label(a),
                isolated=vep.effects.get(sites_id((i, a))),
                clinvar=clinvar.get(key),
                scores=scores.get(key),
            )
    if protein and protein.uniprot and protein_length and protein.length != protein_length:
        ann.notes.append(
            f"UniProt {protein.uniprot} ({protein.length} aa) ne correspond pas à la protéine du transcrit "
            f"({protein_length} aa) : positions des domaines approximatives"
        )
    elif protein is None or not protein.uniprot:
        ann.notes.append("domaines protéiques indisponibles (hors ligne ou protéine absente d'UniProt)")
    return ann


def annotate(
    vcf: Vcf,
    transcript: str,
    cfg: Config,
    workdir: Path,
    protein_length: int | None = None,
    online: bool = True,
) -> Annotations:
    vep = parse_vep(run_vep(write_sites_vcf(vcf, workdir / "sites.vcf"), workdir / "sites.vep.json", cfg),
                    transcript)
    regions = sorted({f"{_bare(r.chrom)}:{r.pos}-{r.pos + len(r.ref) - 1}" for r in vcf.records})
    clinvar, scores, notes = {}, {}, []
    if (path := cfg.annotations.get("clinvar")) and path.exists():
        clinvar = parse_clinvar(tabix(path, regions, cfg))
    else:
        notes.append("ClinVar non configuré")
    if (path := cfg.annotations.get("dbnsfp")) and path.exists():
        header = tabix(path, [], cfg, header=True)[-1]
        scores = parse_dbnsfp(header, tabix(path, regions, cfg), vep.mane_ensembl)
    else:
        notes.append("dbNSFP non configuré (AlphaMissense, REVEL, CADD)")
    protein = None
    if vep.protein_id:
        try:
            protein = protein_info(vep.protein_id, cfg.annotations.get("protein_cache"), online)
        except OSError:
            protein = None
    ann = combine(vcf, vep, clinvar, scores, protein, protein_length)
    ann.notes[:0] = notes
    return ann


# --- lecture par copie --------------------------------------------------------------------


@dataclass(frozen=True)
class CopyVariant:
    key: VariantKey
    status: str  # appliqué | masqué | combiné | sans effet protéique | hors CDS
    effect_on_copy: str | None  # changement de la copie qui contient ce variant (si combiné)

    @property
    def misleading(self) -> bool:
        """L'effet isolé (celui des bases de données) n'est pas celui de cette copie."""
        return self.status in ("masqué", "combiné")


def _same_change(isolated: str, ours: str) -> bool:
    if isolated == ours:
        return True
    # stop perdu : VEP connaît la longueur de l'extension, haplo non (CDS sans 3'UTR)
    return "ext" in isolated and isolated.split("ext")[0] == ours.split("ext")[0]


def _span(change: ProteinChange) -> tuple[int, int]:
    """Positions couvertes : p.Asp167_Glu168delins… → (167, 168)."""
    m = re.match(r"p\.[A-Z][a-z]{2}\d+_[A-Z][a-z]{2}(\d+)", change.hgvs)
    return change.position, int(m.group(1)) if m else change.position


def copy_variants(copy: HaplotypeEffect, alleles: dict[int, int], ann: Annotations) -> list[CopyVariant]:
    """Pour chaque allèle porté par la copie : son effet isolé est-il celui de la copie ?"""
    out = []
    for rec, allele in sorted(alleles.items()):
        if not allele:
            continue
        key = (rec, allele)
        iso = ann.variants.get(key, VariantAnnotation("")).isolated
        if iso is None or iso.protein_start is None:
            out.append(CopyVariant(key, "hors CDS", None))
        elif iso.hgvsp and any(_same_change(iso.hgvsp, c.hgvs) for c in copy.changes):
            out.append(CopyVariant(key, "appliqué", None))
        elif copy.truncated_at is not None and iso.protein_start > copy.truncated_at + 1:
            cut = min((c for c in copy.changes if c.kind in TRUNCATING), key=lambda c: c.position)
            out.append(CopyVariant(key, "masqué", cut.hgvs))
        elif cover := [c for c in copy.changes if _span(c)[0] <= iso.protein_start <= _span(c)[1]]:
            out.append(CopyVariant(key, "combiné", cover[0].hgvs))
        elif not iso.hgvsp or iso.hgvsp.endswith("="):
            out.append(CopyVariant(key, "sans effet protéique", None))
        else:
            out.append(CopyVariant(key, "combiné", None))  # effet isolé absent de la copie
    return out
