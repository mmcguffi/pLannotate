# Annotation differences for manual review

This report compares the refactored implementation and its current database
bundle with the vanilla pLannotate 1.2.5 controls. It was generated on
2026-06-21 with the toolchain recorded in `regression-context.json`.

## Summary

- 31 mode/plasmid combinations were compared.
- 13 combinations are unchanged.
- 18 combinations have a CSV or GenBank difference.
- No plasmid sequence, topology, feature location, or feature type changed,
  except that four detailed-mode outputs removed two duplicate `11.0`
  `misc_feature` annotations each.
- Most remaining changes are expected consequences of newer FPbase,
  Swiss-Prot, and Rfam metadata.

| Plasmid | Regular | Detailed | Linear | Detailed + linear |
| --- | --- | --- | --- | --- |
| RF0G-IodoY | unchanged (17) | unchanged (22) | unchanged (17) | — |
| pACYC184 | unchanged (10) | unchanged (14) | unchanged (11) | — |
| pBTK562 | metadata (11) | **15 → 13** | metadata (12) | — |
| pCA-mTmG | unchanged (29) | **38 → 36** | metadata (29) | — |
| pCMVR8.74 | database data (27) | database data (40) | database data (29) | — |
| pPAGFP-C | metadata (14) | **23 → 21** | metadata (14) | — |
| pSC101 | unchanged (10) | unchanged (10) | unchanged (10) | — |
| pTN7-pa1-GFP-Kan | metadata (22) | **30 → 28** | metadata (23) | — |
| pUC19 | unchanged (11) | unchanged (14) | unchanged (11) | — |
| pXampl3 | Rfam data (19) | Rfam data (21) | Rfam data (18) | Rfam data (20) |

Numbers in parentheses are annotation counts. Bold count changes require the
most attention.

## Feature-set changes

Only detailed mode changes the feature set:

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

## Same features with changed database data

- `pCMVR8.74`, all three modes:
  - The `vpu` fragment remains at 7860–7956, but Swiss-Prot changes from
    accession `P69699` to `P05919`, identity 87.5% to 90.9%, source length 243
    to 246 nt, and match length 39.506% to 39.024%. Its description is updated.
  - Both `gpt` fragments retain their locations and scores; descriptions are
    updated by the newer Swiss-Prot release.
- `pXampl3`, regular, detailed, linear, and detailed+linear modes:
  - The `5S rRNA` remains at 1435–1551. Rfam model length changes from 119 to
    120 nt, changing percent match from 97.479% to 96.667%. This also updates
    the GenBank feature qualifier, but not its location or type.

## CSV metadata-only changes

- `pBTK562` regular/detailed/linear and `pTN7-pa1-GFP-Kan`
  regular/detailed/linear: `GFPmut3` keeps the same feature, location, and match;
  its FPbase identifier casing and description are updated.
- `pPAGFP-C` regular/detailed/linear: `mPA-GFP` keeps the same feature,
  location, and match; identifier casing and the FPbase description are updated.
- `pTN7-pa1-GFP-Kan` regular/detailed: the `neo` description drops the old
  `nptII` alias; its annotation is otherwise unchanged.
- `pCA-mTmG` linear: the `MARCKSL1` description drops the old `MRP` alias; its
  annotation is otherwise unchanged.

## Fully unchanged controls

All requested modes are unchanged for `RF0G-IodoY`, `pACYC184`, `pSC101`, and
`pUC19`. Regular `pCA-mTmG` is also unchanged.
