"""Fige le cas de phasage NAT2 (trio CEPH 1463 de 1000 Genomes, 30x GRCh38) dans tests/data/.

    python tests/make_phasing_fixtures.py   # conteneur + réseau requis

NAT2 : profil d'acétylation déterminé par la combinaison de SNP faux-sens sur chaque copie.
NA12878 porte c.341T>C (p.Ile114Thr), c.481C>T (p.Leu161=) et c.803 (p.Arg268Lys) hétérozygotes.

- nat2_trio.grch38.vcf : génotypes du panel phasé 1000G, phase retirée (`/`), INFO retirée ;
- phasing/nat2_trio.panel.vcf : les mêmes, avec la phase du panel (SHAPEIT + trio) ;
- phasing/nat2_trio.whatshap.vcf : WhatsHap sur les reads seuls (CRAM 1000G, région extraite) ;
- nat2_trio.ped.
Les BAM régionaux ne sont pas gardés (régénérés par ce script).
"""

import re
import tempfile
from pathlib import Path

from gemprot import config
from gemprot.phasing import run_whatshap
from gemprot.tools import run_tool

DATA = Path(__file__).parent / "data"
REGION = "chr8:18399000-18402000"
PANEL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
         "20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr8.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
CRAMS = {
    "NA12878": "https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram",
    "NA12891": "https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989341/NA12891.final.cram",
    "NA12892": "https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989342/NA12892.final.cram",
}
PED = "CEPH1463\tNA12878\tNA12891\tNA12892\t2\t0\nCEPH1463\tNA12891\t0\t0\t1\t0\nCEPH1463\tNA12892\t0\t0\t2\t0\n"


def main() -> None:
    cfg = config.load()
    (DATA / "phasing").mkdir(exist_ok=True)
    (DATA / "nat2_trio.ped").write_text(PED)
    with tempfile.TemporaryDirectory(dir=Path.cwd()) as tmp:
        tmp = Path(tmp)

        def tool(*cmd):
            return run_tool(list(cmd), cfg.vep.container, (*cfg.binds, tmp, cfg.fasta.parent)).stdout

        panel = tool("bcftools", "view", "-r", REGION, "-s", ",".join(CRAMS), "-a", "-c", "1", PANEL)
        lines = [line for line in panel.splitlines()
                 if not line.startswith("##INFO") and not line.startswith("##bcftools")]
        rows = []
        for line in lines:
            if line.startswith("#"):
                rows.append(line)
                continue
            f = line.split("\t")
            f[7] = "."
            rows.append("\t".join(f))
        (DATA / "phasing" / "nat2_trio.panel.vcf").write_text("\n".join(rows) + "\n")
        unphased = [r if r.startswith("#") else re.sub(r"(\d)\|(\d)", r"\1/\2", r) for r in rows]
        vcf = DATA / "nat2_trio.grch38.vcf"
        vcf.write_text("\n".join(unphased) + "\n")

        bams = []
        for sample, url in CRAMS.items():
            bam = tmp / f"{sample}.bam"
            tool("samtools", "view", "-T", cfg.fasta, "-b", "-o", bam, "-X", url, url + ".crai", REGION)
            tool("samtools", "index", bam)
            bams.append(bam)
        tmp_vcf = tmp / "in.vcf"
        tmp_vcf.write_text(vcf.read_text())
        out = run_whatshap(tmp_vcf, bams, list(CRAMS), tmp / "whatshap.vcf", cfg)
        kept = [line for line in out.read_text().splitlines() if not line.startswith("##commandline")]
        (DATA / "phasing" / "nat2_trio.whatshap.vcf").write_text("\n".join(kept) + "\n")
    print("ok :", vcf, "et tests/data/phasing/")


if __name__ == "__main__":
    main()
