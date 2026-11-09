# gemprot image: Python package + vep/haplo, bcftools/htslib, WhatsHap, HiPhase, SHAPEIT5.
# The data (VEP cache, FASTA, ClinVar, dbNSFP) stay outside the image. See config/gemprot.image.yaml.
FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.source="https://github.com/TaniaCuppens/GEMPROT" \
      org.opencontainers.image.description="GEMPROT 2: the two proteins of a sample from a phased VCF" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.image.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /opt/gemprot/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /opt/gemprot/src
COPY --chown=$MAMBA_USER:$MAMBA_USER config/gemprot.image.yaml /opt/gemprot/gemprot.yaml
RUN pip install --no-deps --no-cache-dir /opt/gemprot

# singularity exec does not use the micromamba entrypoint. Set PATH explicitly.
ENV PATH=/opt/conda/bin:$PATH \
    GEMPROT_CONFIG=/opt/gemprot/gemprot.yaml
