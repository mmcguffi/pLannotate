# Annotation differences for manual review

This report compares the refactored implementation and its current database
bundle with the vanilla pLannotate 1.2.5 controls. It was generated on
2026-08-15 with the toolchain recorded in `regression-context.json`, and
replaces the 2026-06-21 report, which predates the GenBank qualifier work.

Regenerate it with `python tools/annotation_controls.py compare`, which writes
the same comparison to `artifacts/annotation-controls/`.

## Summary

- 31 mode/plasmid combinations were compared.
- All 31 differ, because every annotation now carries qualifiers 1.2.5 never
  emitted. That is the intended change, not a regression.
- No plasmid sequence, topology, feature location, or feature type changed, and
  no annotation was added. The only feature-set change is the removal of eight
  duplicate `11.0` `misc_feature` annotations across four detailed-mode outputs.
- Only three qualifiers changed value rather than appearing: `note` on every
  feature, and `identity` and `match_length` on a handful.

| Plasmid | Regular | Detailed | Linear | Detailed + linear |
| --- | --- | --- | --- | --- |
| RF0G-IodoY | qualifiers (17) | qualifiers (22) | qualifiers (17) | — |
| pACYC184 | qualifiers (10) | qualifiers (14) | qualifiers (11) | — |
| pBTK562 | qualifiers (11) | **15 → 13** | qualifiers (12) | — |
| pCA-mTmG | qualifiers (29) | **38 → 36** | qualifiers (29) | — |
| pCMVR8.74 | database data (27) | database data (40) | database data (29) | — |
| pPAGFP-C | qualifiers (14) | **23 → 21** | qualifiers (14) | — |
| pSC101 | qualifiers (10) | qualifiers (10) | qualifiers (10) | — |
| pTN7-pa1-GFP-Kan | qualifiers (22) | **30 → 28** | qualifiers (23) | — |
| pUC19 | qualifiers (11) | qualifiers (14) | qualifiers (11) | — |
| pXampl3 | Rfam data (19) | Rfam data (21) | Rfam data (18) | Rfam data (20) |

Numbers in parentheses are annotation counts. Bold count changes require the
most attention.

## Feature-set changes

Only detailed mode changes the feature set, and only by removing duplicates:

- `pBTK562`: removes two duplicate FPbase `11.0` `misc_feature` entries at
  3671–4385. The retained `GFPmut3` CDS covers the same interval.
- `pCA-mTmG`: removes two duplicate FPbase `11.0` `misc_feature` entries at
  5092–5803. The retained fluorescent-protein CDS covers the same interval.
- `pPAGFP-C`: removes two duplicate FPbase `11.0` `misc_feature` entries at
  618–1329. The retained `mPA-GFP` CDS covers 612–1329.
- `pTN7-pa1-GFP-Kan`: removes two duplicate FPbase `11.0` `misc_feature`
  entries at 4041–4755. The retained `GFPmut3` CDS covers the same interval.

Vanilla 1.2.5 creates these duplicates because circular sequences are searched
twice, FPbase identifier `11.0` does not join to metadata key `11`, and the
resulting untyped hits evade its detailed-mode overlap grouping. SQLite resolves
the identifier and type, so the refactored overlap filter removes the duplicates.

## New qualifiers

583 annotations were emitted across the 31 controls. Each count below is how
many of them gained that qualifier; a qualifier is absent where it does not
apply, never blank.

| Qualifier | Features | Source |
| --- | ---: | --- |
| `annotator`, `subject_start`, `subject_end` | 583 | every hit |
| `btop` | 557 | BLAST and DIAMOND alignment traces |
| `domain`, `host_range` | 155 | curated marker and origin tables |
| `reference` | 136 | curated marker and origin tables |
| `selection_marker`, `selection_agent` | 97 | curated marker table |
| `copy_number_class`, `copy_number_note` | 58 | curated origin table |
| `copy_number` | 36 | curated origin table |
| `structure` | 26 | Infernal consensus structure, WUSS notation |

`btop` and `structure` partition the 583 exactly: a covariance-model hit has no
alignment trace, and no other source has a consensus structure.

## Changed qualifier values

- `note`, all 583: 1.2.5 wrote the literal string `pLannotate` into `/note` and
  had nowhere to put the feature description. The description now occupies
  `/note`, and the tool name moved to `/annotator`.
- `identity`, 15 annotations: every one is an Rfam hit. Infernal reports no
  identity, so 1.2.5 recorded a hardcoded `100.0`; the value is now the
  alignment's average posterior probability. `tRNA` reads 97.0, `RNAI` 93.0,
  and the two `AAC AAD leader` hits 91.0 and 94.0. A covariance model scores
  structure rather than base identity, so these are confidences, not
  base counts, and none of them indicates a worse hit than 1.2.5 reported.
- `match_length`, 7 annotations:
  - `5S rRNA` in all four `pXampl3` modes, 97.5% → 96.7%. The Rfam model length
    is read from the model's `CLEN` rather than from the aligned span, giving
    120 nt where 1.2.5 recorded 119.
  - `vpu (fragment)` in the three `pCMVR8.74` modes, 39.5% → 39.0%, from the
    Swiss-Prot change below.

## Changed CSV data

- Rfam `sseqid`, every control containing an ncRNA: 1.2.5 reported the SQLite
  row number (`1`, `3`, `5`) where the refactor reports the stable Rfam
  accession (`RF00106`, `RF02912`, `RF00005`). The `Description` column already
  carried the accession in both.
- `pCMVR8.74`, all three modes: the `vpu` fragment remains at 7860–7956, but
  Swiss-Prot changes from accession `P69699` to `P05919`, identity 87.5% to
  90.9%, and source length 243 to 246 nt. Its description is updated.
- Descriptions throughout are updated by the newer FPbase, Swiss-Prot, and Rfam
  releases; identifier casing changed for `GFPmut3` and `mPA-GFP`, the `neo`
  description drops the old `nptII` alias, and `MARCKSL1` drops the old `MRP`
  alias. None of these changes an annotation's location, type, or score.

## Database manifest

The comparison reports one context warning: the installed database bundle does
not match `current-database-manifest.json`. The difference is packaging only —
every search index is byte-identical, and the manifest differs in build-date
strings and in how the description SQLite databases are split per source.
