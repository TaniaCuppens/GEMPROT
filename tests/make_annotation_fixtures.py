"""Freezes the annotation sources of the test VCFs in tests/data/annot/ (needs the container and the network).

    python tests/make_annotation_fixtures.py [case …]   # all cases by default

For each case: VEP output of the sites, ClinVar lines, dbNSFP lines (useful columns only),
UniProt protein + Pfam domains.
"""

import json
import sys
import tempfile
from pathlib import Path

from gemprot import annotate as an
from gemprot import config
from gemprot.phase import read_vcf

DATA = Path(__file__).parent / "data"
OUT = DATA / "annot"
CASES = {
    "gjb2_phase_cases": "NM_004004.6",
    "gjb2_r75q_35delG": "NM_004004.6",
    "prnp_d178n_m129v": "NM_000311.5",
    "prnp_synthetic": "NM_000311.5",
    "hfe_c282y_h63d": "NM_000410.4",
    "nat2_trio": "NM_000015.3",
}
# REVEL and CADD are excluded: non-commercial licenses, no redistribution in a public repository.
# AlphaMissense is under CC BY 4.0.
DBNSFP_COLS = ["chr", "pos(1-based)", "ref", "alt", "Ensembl_transcriptid", "MANE", "AlphaMissense_score"]


def main() -> None:
    cfg = config.load()
    OUT.mkdir(exist_ok=True)
    header = an.tabix(cfg.annotations["dbnsfp"], [], cfg, header=True)[-1]
    cols = header.lstrip("#").split("\t")
    keep = [cols.index(c) for c in DBNSFP_COLS]
    for name in sys.argv[1:] or CASES:
        tx = CASES[name]
        vcf = read_vcf(DATA / f"{name}.grch38.vcf")
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as tmp:
            sites = an.write_sites_vcf(vcf, Path(tmp) / "sites.vcf")
            out = an.run_vep(sites, Path(tmp) / "sites.vep.json", cfg)
            (OUT / f"{name}.vep.json").write_text(out.read_text())
        vep = an.parse_vep(OUT / f"{name}.vep.json", tx)
        regions = sorted({f"{an._bare(r.chrom)}:{r.pos}-{r.pos + len(r.ref) - 1}" for r in vcf.records})
        clin = an.tabix(cfg.annotations["clinvar"], regions, cfg)
        (OUT / f"{name}.clinvar.vcf").write_text("".join(line + "\n" for line in clin))
        rows = ["#" + "\t".join(DBNSFP_COLS)]
        for line in an.tabix(cfg.annotations["dbnsfp"], regions, cfg):
            f = line.split("\t")
            rows.append("\t".join(f[i] for i in keep))
        (OUT / f"{name}.dbnsfp.tsv").write_text("\n".join(rows) + "\n")
        info = an.protein_info(vep.protein_id, None)
        (OUT / f"{name}.protein.json").write_text(json.dumps(an.protein_info_to_json(info), indent=1) + "\n")
        print(name, len(vep.effects), "effets,", len(clin), "ClinVar,", len(rows) - 1, "dbNSFP,", info.uniprot)


if __name__ == "__main__":
    main()
