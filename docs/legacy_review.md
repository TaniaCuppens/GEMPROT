# Revue du code GEMPROT v1 (Perl, 2016-2018)

Code de référence : `legacy/GEMPROT/` (clone de https://github.com/TaniaCuppens/GEMPROT, commit 9b09fbc).
Sert de **spécification** pour GEMPROT 2 : ce qu'il faut garder (l'idée, les sorties) et ce qu'il faut éviter.

## Pourquoi v1 ne tourne plus

- `pfam.xfam.org` est fermé (Pfam intégré à InterPro). `Fonction_final.pm:1169` interroge Pfam à
  chaque lancement, même si le gène est en cache ; la page d'erreur 403 commence par `<!--`, passe le
  test `m/^\</` et fait planter `XML::LibXML`. Les boucles `while` (l.1225, 1257) n'ont ni limite ni délai.
- `lysine.univ-brest.fr` (CSS/JS/header) et `med-laennec.univ-brest.fr` (application web) sont hors ligne.
- CCDS, ClinVar 2018 (ancien format INFO) codés en dur dans la config.

## Bugs qui faussent les résultats

| Où | Problème | Effet |
|---|---|---|
| `Indiv.pm:320` | `$chr = "chr";` au lieu de `"chr" . $chr` | En mode indiv avec un VCF `chr1`, aucun variant extrait, ClinVar vide |
| `Fonction_final.pm:821,834,996,1010` | `$assembly !~ [..]` compare à `ARRAY(0x…)` | `chr` toujours ajouté dans le lien MuPIT |
| `Indiv.pm:141` | idem, dans un `defined(...)` | `--assembly` accepte n'importe quoi |
| `Fonction_final.pm:1613,1824` | `my $i_inf` masque la variable extérieure | Insertions avec REF > 1 nt placées au mauvais endroit |
| `Fonction_final.pm:1694,1874,1921` | test sur `$i`, affectation sur `$i_inf` | Phase enregistrée à une autre position que celle testée |
| `Fonction_final.pm:1960` | `$REF_nt[$i]` en double (hap2 seulement) | Délétions hap2 affichées `AA>A` |
| `Indiv.pm:842` | imprime une référence de hash | `hap_<gene>.txt` contient `HASH(0x…)` |
| `Fonction_final.pm:1520` | `0\|.` → manquant sur les deux haplotypes | Faux « missing » |
| `Fonction_final.pm:1797` | `$info[1]` indéfini pour un génotype haploïde | chrX/Y masculins faux |
| `haplotype_translation_final.pl:110` | bloc `END` crée `tmp/OK` même après échec | L'app web voit un plantage comme un succès |

## Robustesse / performance

- Commandes shell construites par concaténation (injection possible, chemins avec espaces cassés).
- Mode pop en O(échantillons × variants) : `grep` sur tout le VCF puis relecture complète par échantillon.
- `zgrep` sur tout ClinVar par gène (→ tabix).
- Choix interactif du CCDS avec `alarm(10)` : inutilisable en batch.

## À garder comme spécification

- Statuts : phasé (bleu hap1 / rouge hap2), non phasé (gris), manquant — `Fonction_final.pm:509-782`
- Types : synonyme, faux-sens, stop, élongation, start-loss, insertion, délétion, frameshift (+ nouveau stop)
- Propriétés physico-chimiques des acides aminés — `Fonction_final.pm:459`
- Rapport indiv (tableau échantillon × hap1/hap2) et résumé des haplotypes (n chromosomes, fréquence, individus) — `Indiv.pm:765-877`
- Mode pop par groupes — `Pop.pm`
- Sorties FASTA par haplotype (`--fasta`)
