# Annotation differences for manual review

This historical report compares the refactored implementation and its current
database bundle with the full vanilla pLannotate 1.2.5 mode matrix. It was regenerated
on 2026-08-23 for #83 after enabling the nested-feature policy in detailed mode. The
current product has one annotation behavior and its active controls use only the
equivalent legacy detailed outputs; the broader matrix below is retained as audit
history. The toolchain and database differences from the frozen control are recorded
below.

This archived matrix is the review artifact produced for #83. The current
`python tools/annotation_controls.py compare` command intentionally covers only the
single supported annotation behavior and therefore does not reproduce this matrix.

## Summary

- 31 mode/plasmid combinations were compared.
- All 31 differ because every retained annotation now carries qualifiers 1.2.5
  never emitted. That is an intended change, not a regression.
- No plasmid sequence, topology, retained feature location, or retained feature
  type changed, and no annotation was added. Detailed mode removes 16 calls:
  eight duplicate FPbase records and eight contained fragments rejected by the
  nested-feature policy.
- Only three qualifiers changed value rather than appearing: `note` on every
  feature, and `identity` and `match_length` on a handful.

| Plasmid | Regular | Detailed | Linear | Detailed + linear |
| --- | --- | --- | --- | --- |
| RF0G-IodoY | qualifiers (17) | qualifiers (22) | qualifiers (17) | — |
| pACYC184 | qualifiers (10) | qualifiers (14) | qualifiers (11) | — |
| pBTK562 | qualifiers (11) | **15 → 13** | qualifiers (12) | — |
| pCA-mTmG | qualifiers (29) | **38 → 33** | qualifiers (29) | — |
| pCMVR8.74 | database data (27) | **40 → 38** | database data (29) | — |
| pPAGFP-C | qualifiers (14) | **23 → 18** | qualifiers (14) | — |
| pSC101 | qualifiers (10) | qualifiers (10) | qualifiers (10) | — |
| pTN7-pa1-GFP-Kan | qualifiers (22) | **30 → 28** | qualifiers (23) | — |
| pUC19 | qualifiers (11) | qualifiers (14) | qualifiers (11) | — |
| pXampl3 | Rfam data (19) | Rfam data (21) | Rfam data (18) | Rfam data (20) |

Numbers in parentheses are annotation counts. Bold count changes require the
most attention.

## Feature-set changes

Only detailed mode changes the feature set.

Eight removals are duplicate FPbase records:

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

The other eight removals are short fragments contained by larger annotations:

- `pCA-mTmG`: `VP1_SV40` CDS (57/1086 nt), a 17/200 nt CMV promoter,
  and a 36/978 nt CMV IE94 promoter.
- `pCMVR8.74`: a 36/978 nt CMV IE94 promoter and a 19/537 nt mCMV promoter.
- `pPAGFP-C`: a 17/200 nt CMV promoter, a 36/978 nt CMV IE94 promoter,
  and a 17/286 nt CMV enhancer.

These are exact but very low-coverage fragments. Exact non-CDS elements at
least 30% complete are retained, as are near-complete children, boundary
crossings, structured ncRNAs, compound/same-kind children, and all calls not
actually contained by another hit. This is why the T5 promoter in the pTN7
detailed control survives while the tiny CMV-family fragments above do not.

## New qualifiers

575 annotations were emitted across the 31 controls. Each count below is how
many of them gained that qualifier; a qualifier is absent where it does not
apply, never blank.

| Qualifier | Features | Source |
| --- | ---: | --- |
| `annotator`, `subject_start`, `subject_end` | 575 | every hit |
| `btop` | 549 | BLAST and DIAMOND alignment traces |
| `domain`, `host_range` | 155 | curated marker and origin tables |
| `reference` | 136 | curated marker and origin tables |
| `selection_marker`, `selection_agent` | 97 | curated marker table |
| `copy_number_class`, `copy_number_note` | 58 | curated origin table |
| `copy_number` | 36 | curated origin table |
| `structure` | 26 | Infernal consensus structure, WUSS notation |

`btop` and `structure` partition the 575 exactly: a covariance-model hit has no
alignment trace, and no other source has a consensus structure.

For translated DIAMOND hits, `subject_start` and `subject_end` are emitted in
nucleotide-equivalent units so they use the same unit as the database feature
length. This is a qualifier representation change only; query locations do not
move.

## Changed qualifier values

- `note`, all 575 retained annotations: 1.2.5 wrote the literal string
  `pLannotate` into `/note` and
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

The comparison reports two context warnings. The installed database bundle does
not match `current-database-manifest.json`; that difference is packaging only —
every search index is byte-identical, and the manifest differs in build-date
strings and in how the description SQLite databases are split per source. The
local BLAST and DIAMOND patch versions also differ from the frozen control
(`2.16.0` vs `2.17.0`, and `2.2.3` vs `2.1.24`); Infernal remains `1.1.5`.
