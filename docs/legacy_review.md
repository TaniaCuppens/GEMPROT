# Review of the GEMPROT v1 code (Perl, 2016-2018)

Reference code: `legacy/GEMPROT/` (clone of https://github.com/TaniaCuppens/GEMPROT, commit 9b09fbc).
This code is the **specification** for GEMPROT 2: what to keep (the idea, the outputs) and what to avoid.

## Why v1 no longer runs

- `pfam.xfam.org` is closed (Pfam is now part of InterPro). `Fonction_final.pm:1169` queries Pfam at
  each run, even when the gene is in the cache. The 403 error page starts with `<!--`, passes the
  test `m/^\</` and makes `XML::LibXML` crash. The `while` loops (l.1225, 1257) have no limit and no delay.
- `lysine.univ-brest.fr` (CSS/JS/header) and `med-laennec.univ-brest.fr` (web application) are offline.
- CCDS and ClinVar 2018 (old INFO format) are hard-coded in the config.

## Bugs that make the results wrong

| Where | Problem | Effect |
|---|---|---|
| `Indiv.pm:320` | `$chr = "chr";` instead of `"chr" . $chr` | In indiv mode with a `chr1` VCF, no variant is extracted and ClinVar is empty |
| `Fonction_final.pm:821,834,996,1010` | `$assembly !~ [..]` compares with `ARRAY(0x…)` | `chr` is always added to the MuPIT link |
| `Indiv.pm:141` | same problem, in a `defined(...)` | `--assembly` accepts any value |
| `Fonction_final.pm:1613,1824` | `my $i_inf` hides the outer variable | Insertions with REF > 1 nt go to the wrong position |
| `Fonction_final.pm:1694,1874,1921` | test on `$i`, assignment to `$i_inf` | The phase is stored at a position different from the tested one |
| `Fonction_final.pm:1960` | `$REF_nt[$i]` twice (hap2 only) | hap2 deletions show as `AA>A` |
| `Indiv.pm:842` | prints a hash reference | `hap_<gene>.txt` contains `HASH(0x…)` |
| `Fonction_final.pm:1520` | `0\|.` → missing on the two haplotypes | False "missing" |
| `Fonction_final.pm:1797` | `$info[1]` is undefined for a haploid genotype | Male chrX/Y results are wrong |
| `haplotype_translation_final.pl:110` | the `END` block creates `tmp/OK` even after a failure | The web app sees a crash as a success |

## Robustness and performance

- Shell commands are built by concatenation (injection is possible, paths with spaces break).
- Pop mode is O(samples × variants): `grep` on the full VCF, then a full read again for each sample.
- `zgrep` on all of ClinVar for each gene (→ tabix).
- Interactive CCDS selection with `alarm(10)`: not usable in batch.

## Keep as specification

- Statuses: phased (blue hap1 / red hap2), unphased (grey), missing — `Fonction_final.pm:509-782`
- Types: synonymous, missense, stop, elongation, start-loss, insertion, deletion, frameshift (+ new stop)
- Physico-chemical properties of the amino acids — `Fonction_final.pm:459`
- Indiv report (table of sample × hap1/hap2) and haplotype summary (n chromosomes, frequency, samples) — `Indiv.pm:765-877`
- Pop mode by groups — `Pop.pm`
- FASTA output for each haplotype (`--fasta`)
