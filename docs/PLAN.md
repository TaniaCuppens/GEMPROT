# GEMPROT 2 — plan (v2: positioning, phasing module, single container)

## Context

GEMPROT (Cuppens et al., *Bioinformatics* 35(14):2492, 2019): from a phased VCF, show the two
proteins of a sample, because variants interpreted one at a time can mislead
(variants in cis that combine or cancel each other). The Perl code no longer runs. We keep the idea.

**Review of existing tools (Oct. 2026)** — each tool covers one part, none covers the full chain:
Haplosaurus (`haplo`, Ensembl) translates each haplotype but gives raw JSON. BCFtools/csq gives
phase-aware consequences, one variant at a time. SmartPhase and gnomAD
variant co-occurrence (Guo 2023) give cis/trans for one pair. WhatsHap/HiPhase/SHAPEIT5 do the phasing.
ProHap is for proteomics. **Position of GEMPROT 2**: a cis/trans interpretation layer:
where does the phase come from and how much do we trust it → what are the two proteins → are both
copies affected? All of this in a readable report. We do not write the translation again (haplo).

## Done

- Repository: `docs/` (paper, `legacy_review.md`, `PLAN.md`),
  `legacy/GEMPROT/`, `config/gemprot.yaml`, `tests/data/gjb2_r75q_35delG.grch38.vcf`.
- Shared resources connected (no copy): VEP RefSeq 115 cache (links in `cache/vep`),
  ClinVar 2026-04, dbNSFP 5.3.1a, gnomAD v4.1 sites, SpliceAI, FASTA GRCh38,
  ensembl-vep 115.2 container.
- `haplo` validated on GJB2: CIS = truncated + REF. TRANS = truncated + 75R>Q.

## New: phasing module (`phase.py`)

1. **Audit of the input phase**, from v1: for each pair of variants in the gene → `cis` / `trans`
   / `unknown`, from `GT` *and* `PS`. Two `0|1` in different `PS` blocks = **not phased
   with each other** (haplo treats them as phased, which is wrong → we correct it). Unphased genotypes, `0|.` and haploid genotypes are handled.
2. **Phase sources**, used when they are available:
   - **trio / pedigree** (parent VCF + `.ped`): Mendelian transmission, almost certain.
   - **reads** (BAM/CRAM, short or long reads): WhatsHap, which combines reads + pedigree natively.
   - **population**: SHAPEIT5 with `--use-PS` (keeps the read phase and connects the blocks) for
     common variants. **gnomAD co-occurrence** (cis/trans probability for each pair) for rare variants.
3. **Combination of the evidence** (answer to "can we combine them"): for each pair, keep each opinion
   with its source and its confidence, then:
   - hierarchy: trio > reads > population (statistical).
   - agreement between sources → high confidence. Conflict → reported (switch error, swapped sample…).
   - the trio and the population also **connect `PS` blocks** that the reads did not join.
   - the trio checks the switch errors of the read phasing.
4. **Report**: for each pair, the source(s) and the confidence. If the phase stays unknown, show
   the **two scenarios** ("if cis → …, if trans → …") and do not choose one.

Data to get: gnomAD v4 variant co-occurrence table (not in our gnomAD v4.1 data,
which has only the `joint` sites). Genetic map + GRCh38 reference panel for SHAPEIT5 (option, later).

## New: distribution — one image, tools included, data outside

- **One `gemprot` image** that contains all the tools: Python + `gemprot` package, ensembl-vep 115
  (`vep`, `haplo`), bcftools/htslib, WhatsHap, HiPhase, SHAPEIT5. Built from a
  `Dockerfile` (micromamba + bioconda, versions pinned in `environment.lock`).
- The **data do not go in the image** (VEP cache 24 GB, dbNSFP, gnomAD: too large, and already
  shared here). `gemprot setup`: either point to existing data (lab case, through the
  config), or download the minimum (VEP RefSeq cache, ClinVar) into a selected folder.
- **Build**: no `--fakeroot` on this server (no `/etc/subuid` entry) → **GitHub Actions** builds the
  image and publishes it on GHCR. Here: `singularity pull docker://ghcr.io/…/gemprot:<version>`.
  Docker users use it directly. In the lab, we can put it with the shared containers.
- The Python code calls the tools **directly** when it runs in the image, and through
  `singularity exec` otherwise (development here) — one function `run_tool()` in `tools.py`.

## Architecture

```
src/gemprot/
  cli.py         # typer: indiv, pop, phase, setup
  config.py      # YAML (data paths, optional container)
  tools.py       # run_tool(): direct or singularity exec; subprocess with a list, never through a shell
  extract.py     # bcftools: MANE transcript region, samples, chr/no-chr normalization
  phase.py       # GT/PS audit, WhatsHap/SHAPEIT5 runs, gnomAD co-occurrence, evidence fusion
  haplotypes.py  # haplo → JSON → ProteinHaplotype (split by reliable phase blocks)
  engine.py      # minimal in-house translation (pyfaidx), cross-check of haplo
  classify.py    # haplo notation → HGVS; combined effects (same codon, compensated frameshift,
                 #   variant masked by a stop/frameshift in cis); verdict "0, 1 or 2 affected copies"
  properties.py  # physico-chemical properties (spec: Fonction_final.pm:459)
  annotate.py    # ClinVar, dbNSFP (AlphaMissense, REVEL, CADD), InterPro domains (JSON cache)
  structure.py   # AlphaFold DB, 3D distances between variants in cis
  popstats.py    # haplotypes × groups, Fisher/χ², permutations, FDR
  report/templates/  # Jinja2 → self-contained HTML + TSV + FASTA
Dockerfile, environment.yml, .github/workflows/image.yml
```

Convention: chromosomes **without `chr`** internally (VEP cache, ClinVar, test VCF). Conversion at
input/output when the user FASTA or VCF uses `chr`.

## Order of work

1. Skeleton: `pyproject.toml`, `environment.yml`, `config.py`, `tools.py`, `extract.py`,
   `haplotypes.py`, `classify.py`. First pytest test = GJB2 cis/trans.
2. `phase.py` **audit** part (GT/PS, cis/trans scenarios) + tests: different PS blocks, unphased, `0|.`.
3. PRNP and HFE cases + synthetic cases (same codon, compensated frameshift, start lost, − strand, chrX, multi-allelic).
4. Indiv mode → HTML report (two tracks, phase status, verdict).
5. `Dockerfile` + GitHub Actions → GHCR image. Check `singularity pull` + `gemprot indiv` in the image.
   Repository: branch **`v2`** of `TaniaCuppens/GEMPROT` (decided on 2026-10-09). The shared VEP container does not have
   `bcftools` (only `tabix`): `--region` works only after the gemprot image is built.
   **Written**: `Dockerfile` (micromamba), `environment.image.yml` (pinned versions: vep 115.2, bcftools 1.21,
   WhatsHap 2.8, HiPhase 1.7.0, SHAPEIT5 5.1.1), `config/gemprot.image.yaml` (data under `/data`),
   `.github/workflows/image.yml` (build, tool check, pytest in the image, publication to
   `ghcr.io/taniacuppens/gemprot`). Same environment installed locally: tests + `gemprot indiv -r`
   (bcftools) on GJB2 OK without a container. **Done** (2026-10-14): Actions OK (tools + tests in the image),
   `ghcr.io/taniacuppens/gemprot:v2` public (visibility inherited from the repository), `singularity pull` here (1.1 GB),
   then `gemprot indiv -r … --html` in the image with the mounted data: same results.
6. Annotations (ClinVar, dbNSFP, InterPro) — **done** (`annotate.py`): isolated effect of each allele from `vep`
   (HGVSc/HGVSp, CDS position), ClinVar and dbNSFP through `tabix` (scores of the Ensembl MANE transcript),
   UniProt + Pfam domains through the API (JSON cache for each protein, no SQLite: home is on NFS).
   GEMPROT adds: for each variant of each copy, "applied / masked / combined / outside CDS".
   Alert when the annotation of a variant does not describe the protein that the copy makes.
7. `phase.py` **phasing** part: WhatsHap (reads ± trio), then gnomAD co-occurrence, then SHAPEIT5. Evidence fusion.
   **7a done** (`phasing.py`): separate sources — trio (Mendelian transmission calculated from the
   parent genotypes + `.ped`), reads (WhatsHap without `--ped`), input VCF — then fusion by
   union-find with parity in the order trio > reads > VCF. Conflicts for each pair are reported, indirect
   contradictions are noted, blocks connected by another source (trio that joins two read blocks → "inferred").
   The result is written again as GT/PS, then the usual audit runs. CLI `--bam` (repeatable), `--ped`. Report: one opinion for each
   source and each pair. Test case: NAT2, CEPH 1463 trio from 1000G 30x (NA12878: c.341C + c.481T in trans of
   c.803. 32 scenarios without phase → 1. Trio and reads agree, and agree with the 1000G panel).
   To do: 7b gnomAD co-occurrence (rare pairs), 7c SHAPEIT5 (common variants, `--use-PS`).
   Done (2026-11-06): scenarios and pairs only include the variants that change the protein (`scope.py`,
   after the phase fusion). The report lists the other variants separately and flags splicing variants (not modeled).
   NAT2 without phase: NA12891 256 → 4 scenarios, NA12878 32 → 2. `--all-variants` gives the old behavior.
8. Pop mode + stats. Graphic output: **mirrored lolliplot** group A (top) / group B (bottom) in SVG
   in the pop HTML report (idea from `lolliplot_update`,
   trackViewer → written again without R): MANE transcript, height = number of carrier **haplotypes**,
   variants masked in cis shown differently (hollow lollipop), shapes/colors of the indiv report (no
   hash colors), InterPro domains (step 6), plus a track "samples with 0 / 1 / 2 copies
   with a changed sequence" for each group.
9. AlphaFold. 10. README + real GRCh38 examples (phased 1000G 30x, regional extract).

## Verification

- `pytest`: GJB2 (cis → 1 normal copy. trans → 2 affected copies). PRNP 178N/129M vs 178N/129V.
  Same codon → 1 amino acid. Two `0|1` variants in different PS blocks → "unknown phase" + 2 scenarios.
- Cross-check `classify` / `engine.py` / `bcftools csq` on the simple cases.
- Phasing: on a 1000G trio (e.g. CEPH family NA12878/NA12891/NA12892, GRCh38 30x), read phasing
  alone vs trio: agreement rate. The fusion must report the conflicts.
- Image: `singularity pull`, then `singularity exec gemprot.sif gemprot indiv …` on the GJB2 case with
  nothing installed on the user side, data mounted with `-B`.
- HTML report opened offline: no network request when it shows.
