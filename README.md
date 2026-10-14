# GEMPROT 2

Interprétation cis/trans des variants : à partir d'un VCF (phasé ou non), les deux copies de la
protéine de chaque individu, le statut de phase des variants et un rapport HTML autonome.

Version 2 en cours de développement (branche `v2`). Version d'origine : Cuppens *et al.*,
*Bioinformatics* 35(14):2492, 2019 (branche `master`).

## Image

```bash
singularity pull gemprot.sif docker://ghcr.io/taniacuppens/gemprot:v2
singularity exec -B /chemin/vers/données:/data gemprot.sif \
    gemprot indiv patient.vcf.gz -t NM_004004.6 --html rapport.html
```

Les données ne sont pas dans l'image. Disposition attendue sous `/data` : voir
`config/gemprot.image.yaml`, ou fournir sa propre config (`--config`, `GEMPROT_CONFIG`).

## Développement

```bash
conda env create -f environment.yml
cp config/gemprot.example.yaml config/gemprot.yaml   # adapter les chemins
pytest
```
