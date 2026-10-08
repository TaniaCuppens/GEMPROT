"""Régénère les VCF de test PRNP et HFE (GRCh38).

Coordonnées réelles prises dans ClinVar 2026-04 ; les individus sont fictifs. Les cas synthétiques
utilisent la séquence de référence (PRNP : un seul exon codant, brin +, c.1 = 20:4699221).

    python tests/make_test_vcfs.py           # écrit dans tests/data/
"""

from pathlib import Path

from gemprot import config
from gemprot.reference import Fasta

DATA = Path(__file__).parent / "data"
PRNP_C1 = 4699221
HEADER = [
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
]


def write(path, chrom, length, records, samples):
    """records : (pos, id, ref, alt, {échantillon: GT}) ; un échantillon absent est 0|0."""
    with open(path, "w") as fh:
        fh.write("\n".join(HEADER + [f"##contig=<ID={chrom},length={length}>"]) + "\n")
        fh.write("\t".join("#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split() + samples) + "\n")
        for pos, vid, ref, alt, gts in sorted(records):
            row = [chrom, str(pos), vid, ref, alt, ".", "PASS", ".", "GT:PS"]
            fh.write("\t".join(row + [f"{gts.get(s, '0|0')}:1" for s in samples]) + "\n")


def main():
    fa = Fasta(config.load().fasta)

    def g(c):  # PRNP c. → position génomique
        return PRNP_C1 + c - 1

    def b(c):
        return fa.fetch("20", g(c), g(c))

    # Monari 1994 : D178N avec 129M en cis → insomnie fatale familiale ; avec 129V → Creutzfeldt-Jakob
    assert b(385) == "A" and b(532) == "G"
    write(DATA / "prnp_d178n_m129v.grch38.vcf", "20", 64444167, [
        (g(385), "rs1799990", "A", "G", {"IFF": "0|1", "MCJ": "1|0"}),
        (g(532), "rs74315403", "G", "A", {"IFF": "1|0", "MCJ": "1|0"}),
    ], ["IFF", "MCJ"])

    # Uguen 2017 : HFE C282Y / H63D
    write(DATA / "hfe_c282y_h63d.grch38.vcf", "6", 170805979, [
        (26090951, "rs1799945", "C", "G", {"C282Y_H63D_TRANS": "0|1", "C282Y_H63D_CIS": "1|0"}),
        (26092913, "rs1800562", "G", "A", {"C282Y_HOMO": "1|1", "C282Y_H63D_TRANS": "1|0", "C282Y_H63D_CIS": "1|0"}),
    ], ["C282Y_HOMO", "C282Y_H63D_TRANS", "C282Y_H63D_CIS"])

    assert b(1) == "A" and b(20) == "G" and b(760) == "T" and b(500) == "A" and b(503) == "A"
    write(DATA / "prnp_synthetic.grch38.vcf", "20", 64444167, [
        (g(1), "start", "A", "G", {"START_PERDU": "1|0"}),
        (g(20), "stop20", "G", "A", {"STOP_GAGNE": "1|0"}),
        (g(384), "del129", b(384) + "ATG", b(384), {"DEL_EN_PHASE": "1|0"}),
        (g(385), "multi", "A", "G,C", {"MULTI": "1|2"}),
        (g(399), "fs400", b(399) + b(400), b(399), {"FS_STOP": "1|0"}),
        (g(409), "ins409", b(409), b(409) + "A", {"FS_STOP": "1|0"}),  # en aval du stop : masqué
        (g(499), "fs500", b(499) + "A", b(499), {"FS_COMPENSE": "1|0"}),
        (g(503), "ins503", "A", "AC", {"FS_COMPENSE": "1|0"}),
        (g(532), "c532", "G", "A", {"MEME_CODON_CIS": "1|0", "MEME_CODON_TRANS": "1|0"}),
        (g(534), "c534", "C", "A", {"MEME_CODON_CIS": "1|0", "MEME_CODON_TRANS": "0|1"}),
        (g(760), "stoplost", "T", "C", {"STOP_PERDU": "1|0"}),
    ], ["START_PERDU", "STOP_GAGNE", "DEL_EN_PHASE", "MULTI", "FS_STOP", "FS_COMPENSE",
        "MEME_CODON_CIS", "MEME_CODON_TRANS", "STOP_PERDU"])


if __name__ == "__main__":
    main()
