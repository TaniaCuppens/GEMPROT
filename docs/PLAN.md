# GEMPROT 2 — plan (v2 : positionnement, module phasage, conteneur unique)

## Contexte

GEMPROT (Cuppens et al., *Bioinformatics* 35(14):2492, 2019) : à partir d'un VCF phasé, montrer les
deux protéines d'un individu, parce que des variants interprétés un par un peuvent tromper
(variants en cis qui se combinent ou se neutralisent). Le code Perl ne tourne plus ; on garde l'idée.

**Revue de l'existant (oct. 2026)** — chaque outil couvre un morceau, aucun la chaîne complète :
Haplosaurus (`haplo`, Ensembl) traduit par haplotype mais rend un JSON brut ; BCFtools/csq donne des
conséquences tenant compte de la phase, variant par variant ; SmartPhase et gnomAD
variant co-occurrence (Guo 2023) disent cis/trans pour une paire ; WhatsHap/HiPhase/SHAPEIT5 phasent ;
ProHap vise la protéomique. **Positionnement de GEMPROT 2** : couche d'interprétation cis/trans —
d'où vient la phase et à quel point on y croit → quelles sont les deux protéines → les deux copies
sont-elles atteintes ? — dans un rapport lisible. On ne réécrit pas la traduction (haplo).

## Déjà fait

- Dépôt : `docs/` (papier, `legacy_review.md`, `PLAN.md`),
  `legacy/GEMPROT/`, `config/gemprot.yaml`, `tests/data/gjb2_r75q_35delG.grch38.vcf`.
- Ressources partagées branchées (aucune copie) : cache VEP RefSeq 115 (liens dans `cache/vep`),
  ClinVar 2026-04, dbNSFP 5.3.1a, gnomAD v4.1 sites, SpliceAI, FASTA GRCh38,
  conteneur ensembl-vep 115.2.
- `haplo` validé sur GJB2 : CIS = tronquée + REF ; TRANS = tronquée + 75R>Q.

## Nouveau : module phasage (`phase.py`)

1. **Audit de la phase d'entrée**, dès la v1 : pour chaque paire de variants du gène → `cis` / `trans`
   / `inconnu`, en lisant `GT` *et* `PS`. Deux `0|1` dans des blocs `PS` différents = **non phasés
   entre eux** (haplo les traite à tort comme phasés → on corrige). Génotypes non phasés, `0|.`, haploïdes gérés.
2. **Sources de phase**, utilisées selon ce qui est fourni :
   - **trio / pedigree** (VCF des parents + `.ped`) : transmission mendélienne, quasi certaine ;
   - **reads** (BAM/CRAM, short ou long reads) : WhatsHap, qui combine nativement reads + pedigree ;
   - **population** : SHAPEIT5 avec `--use-PS` (garde la phase des reads et relie les blocs) pour
     les variants communs ; **gnomAD co-occurrence** (probabilité cis/trans par paire) pour les rares.
3. **Combinaison des preuves** (réponse à « peut-on les combiner ») : par paire, on garde chaque avis
   avec sa source et sa confiance, puis
   - hiérarchie : trio > reads > population (statistique) ;
   - accord entre sources → confiance haute ; désaccord → signalé (erreur de switch, échantillon inversé…) ;
   - le trio et la population servent aussi à **relier des blocs `PS`** que les reads n'ont pas joints ;
   - le trio contrôle les erreurs de switch du phasage par reads.
4. **Rapport** : pour chaque paire, source(s) et confiance ; si la phase reste inconnue, on affiche
   les **deux scénarios** (« si cis → …, si trans → … ») au lieu d'en choisir un.

Données à récupérer : table gnomAD v4 variant co-occurrence (absente de nos données gnomAD v4.1,
qui ne contient que les sites `joint`) ; carte génétique + panel de référence GRCh38 pour SHAPEIT5 (option, plus tard).

## Nouveau : distribution — une image, outils inclus, données à part

- **Une image `gemprot`** contenant tous les outils : Python + paquet `gemprot`, ensembl-vep 115
  (`vep`, `haplo`), bcftools/htslib, WhatsHap, HiPhase, SHAPEIT5. Construite depuis un
  `Dockerfile` (micromamba + bioconda, versions figées dans `environment.lock`).
- Les **données ne vont pas dans l'image** (cache VEP 24 Go, dbNSFP, gnomAD : trop gros, et déjà
  partagées ici). `gemprot setup` : soit pointer vers des données existantes (cas du labo, via la
  config), soit télécharger le minimum (cache VEP RefSeq, ClinVar) dans un dossier choisi.
- **Construction** : pas de `--fakeroot` sur ce serveur (aucune entrée `/etc/subuid`) → l'image est
  construite par **GitHub Actions** et publiée sur GHCR ; ici : `singularity pull docker://ghcr.io/…/gemprot:<version>`.
  Les utilisateurs Docker l'utilisent directement. Au labo on pourra la déposer avec les conteneurs partagés.
- Le code Python appelle les outils **directement** quand il tourne dans l'image, et via
  `singularity exec` sinon (développement ici) — une seule fonction `run_tool()` dans `tools.py`.

## Architecture

```
src/gemprot/
  cli.py         # typer : indiv, pop, phase, setup
  config.py      # YAML (chemins des données, conteneur éventuel)
  tools.py       # run_tool() : direct ou singularity exec ; subprocess en liste, jamais via shell
  extract.py     # bcftools : région du transcrit MANE, échantillons, normalisation chr/sans chr
  phase.py       # audit GT/PS, lancement WhatsHap/SHAPEIT5, gnomAD co-occurrence, fusion des preuves
  haplotypes.py  # haplo → JSON → ProteinHaplotype (découpé par blocs de phase fiables)
  engine.py      # traduction maison minimale (pyfaidx), contrôle croisé de haplo
  classify.py    # notation haplo → HGVS ; effets combinés (même codon, frameshift compensé,
                 #   variant masqué par un stop/frameshift en cis) ; verdict « 0, 1 ou 2 copies atteintes »
  properties.py  # propriétés physico-chimiques (spec : Fonction_final.pm:459)
  annotate.py    # ClinVar, dbNSFP (AlphaMissense, REVEL, CADD), domaines InterPro (cache JSON)
  structure.py   # AlphaFold DB, distances 3D entre variants en cis
  popstats.py    # haplotypes × groupes, Fisher/χ², permutations, FDR
  report/templates/  # Jinja2 → HTML autonome + TSV + FASTA
Dockerfile, environment.yml, .github/workflows/image.yml
```

Convention : chromosomes **sans `chr`** en interne (cache VEP, ClinVar, VCF test) ; conversion à
l'entrée/sortie quand la FASTA ou le VCF utilisateur est en `chr`.

## Ordre de réalisation

1. Squelette : `pyproject.toml`, `environment.yml`, `config.py`, `tools.py`, `extract.py`,
   `haplotypes.py`, `classify.py` ; premier test pytest = GJB2 cis/trans.
2. `phase.py` partie **audit** (GT/PS, scénarios cis/trans) + tests : blocs PS différents, non phasé, `0|.`.
3. Cas PRNP et HFE + cas synthétiques (même codon, frameshift compensé, start perdu, brin −, chrX, multi-allélique).
4. Mode indiv → rapport HTML (deux pistes, statut de phase, verdict).
5. `Dockerfile` + GitHub Actions → image GHCR ; vérifier `singularity pull` + `gemprot indiv` dans l'image.
   Dépôt : branche **`v2`** de `TaniaCuppens/GEMPROT` (décidé le 2026-10-09). Le conteneur VEP partagé n'a pas
   `bcftools` (seulement `tabix`) : `--region` ne marche qu'une fois l'image gemprot construite.
   **Écrit** : `Dockerfile` (micromamba), `environment.image.yml` (versions figées : vep 115.2, bcftools 1.21,
   WhatsHap 2.8, HiPhase 1.7.0, SHAPEIT5 5.1.1), `config/gemprot.image.yaml` (données sous `/data`),
   `.github/workflows/image.yml` (construction, contrôle des outils, pytest dans l'image, publication
   `ghcr.io/taniacuppens/gemprot`). Même environnement installé localement : tests + `gemprot indiv -r`
   (bcftools) sur GJB2 OK sans conteneur. **Fait** (2026-10-14) : Actions OK (outils + tests dans l'image),
   `ghcr.io/taniacuppens/gemprot:v2` public (visibilité héritée du dépôt), `singularity pull` ici (1,1 Go)
   puis `gemprot indiv -r … --html` dans l'image avec les données montées : mêmes résultats.
6. Annotations (ClinVar, dbNSFP, InterPro) — **fait** (`annotate.py`) : effet isolé de chaque allèle par `vep`
   (HGVSc/HGVSp, position CDS), ClinVar et dbNSFP par `tabix` (scores du transcrit Ensembl MANE),
   UniProt + domaines Pfam par API (cache JSON par protéine, pas de SQLite : home sur NFS).
   Apport GEMPROT : pour chaque variant de chaque copie, « appliqué / masqué / combiné / hors CDS » ;
   alerte quand l'annotation d'un variant ne décrit pas la protéine produite.
7. `phase.py` partie **phasage** : WhatsHap (reads ± trio), puis gnomAD co-occurrence, puis SHAPEIT5 ; fusion des preuves.
8. Mode pop + stats. Sortie graphique : **lolliplot miroir** groupe A (haut) / groupe B (bas) en SVG
   dans le rapport HTML pop (idée reprise de `lolliplot_update`,
   trackViewer → réécrit sans R) : transcrit MANE, hauteur = nombre de **haplotypes** porteurs,
   variants masqués en cis distingués (lollipop creux), formes/couleurs du rapport indiv (pas de
   couleurs par hachage), domaines InterPro (étape 6), plus une piste « individus à 0 / 1 / 2 copies
   altérées » par groupe.
9. AlphaFold ; 10. README + exemples GRCh38 réels (1000G 30x phasé, extrait régional).

## Vérification

- `pytest` : GJB2 (cis → 1 copie normale ; trans → 2 copies atteintes) ; PRNP 178N/129M vs 178N/129V ;
  même codon → 1 acide aminé ; deux variants `0|1` dans des PS différents → « phase inconnue » + 2 scénarios.
- Contrôle croisé `classify` / `engine.py` / `bcftools csq` sur les cas simples.
- Phasage : sur un trio 1000G (ex. famille CEPH NA12878/NA12891/NA12892, GRCh38 30x), la phase par
  reads seule vs trio : taux d'accord ; la fusion doit signaler les désaccords.
- Image : `singularity pull` puis `singularity exec gemprot.sif gemprot indiv …` sur le cas GJB2 sans
  rien d'installé côté utilisateur, données montées par `-B`.
- Rapport HTML ouvert hors ligne : aucune requête réseau à l'affichage.
