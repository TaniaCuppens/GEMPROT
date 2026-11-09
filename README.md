# GEMPROT 2

GEMPROT 2 interprets variants in cis and in trans. From a VCF (phased or unphased), it gives the two
protein copies of each sample, the phase status of the variants, and a self-contained HTML report.

Version 2 is in development (branch `v2`). Original version: Cuppens *et al.*,
*Bioinformatics* 35(14):2492, 2019 (branch `master`).

## Image

```bash
singularity pull gemprot.sif docker://ghcr.io/taniacuppens/gemprot:v2
singularity exec -B /path/to/data:/data gemprot.sif \
    gemprot indiv patient.vcf.gz -t NM_004004.6 --html report.html
```

The image does not contain the data. For the expected layout under `/data`, see
`config/gemprot.image.yaml`. You can also give your own config (`--config`, `GEMPROT_CONFIG`).

## Development

```bash
conda env create -f environment.yml
cp config/gemprot.example.yaml config/gemprot.yaml   # set the paths to your data
pytest
```
