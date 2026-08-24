# Nested-feature curation policy

This policy defines which nested annotations may coexist. Its central rule is that
biological validity and display redundancy are different questions: a child can be a
real functional element even when showing both the parent and child is visually noisy.

The policy is deliberately conservative about deletion. It automates suppression only
for weak contained fragments and preserves whole-feature nesting unless a curated
decision demonstrates a broken record or false match.

## Outcomes

Every nested parent/child pair receives exactly one outcome:

1. `keep` — both annotations describe independently useful biological functions.
2. `suppress_child` — the child alignment is a false positive, or an explicit display
   policy says that an already-curated composite parent is the only useful label.
3. `replace_parent` — the parent record is a composite or mislabeled feature and should
   be renamed, retyped, split, or replaced in the source database. Valid children stay.
4. `trim_parent` — primary/source evidence provides exact canonical bounds excluding
   contaminating flanking sequence. This is a database correction, not a runtime guess.
5. `replace_child` — the overlap is real, but hit ranking selected the wrong homolog or
   strain accession. Rerank the child while preserving the interval.
6. `relabel_child` — the overlap supports generic source provenance but not the specific
   protein accession or a functional full-length CDS claim.
7. `review` — evidence is insufficient. Review fails open: keep the call until a curator
   records a decision, so annotation does not silently lose real biology.

`replace_parent` and `trim_parent` change the database bundle. `suppress_child` is the
only outcome implemented as a runtime filter.

## Rules, allowlists, and denylists

This is not a pure whitelist or blacklist. The default behavior is an ordered rule
engine: it keeps biologically credible nesting and suppresses weak, low-coverage
contained fragments. The packaged CSVs add narrowly scoped exceptions:

- a `keep` or other non-suppression row in `nested_feature_overrides.csv` acts like a
  relationship whitelist: that exact parent/child pair survives even if a later
  general rule would suppress it;
- a `suppress_child` row in `nested_feature_overrides.csv` acts like a relationship
  blacklist, but only while that exact child is a fragment strictly contained by that
  exact parent; and
- a row in `feature_suppressions.csv` is a global, source-pinned blacklist entry. It
  removes that exact `(db, sseqid)` from annotation results in every context; and
- a row in `composite_reference_regions.csv` identifies bases of a larger reference
  that come from a known component. It suppresses only fragments whose aligned
  subject window has no more than three bases outside those curated regions; and
- a row in `fragment_suppression_regions.csv` suppresses only a reviewed low-identity
  fragment from one source interval and expected source-record length. Whole and
  higher-identity matches to that accession remain detectable.

The current pair table is mostly whitelist-like because the audit showed that many
nested relationships are legitimate. Actions such as `replace_parent`,
`replace_child`, and `relabel_child` preserve the call at runtime and record future
database or presentation work; they do not currently rewrite an annotation.

## Automated rules

Rules run in this order.

### 1. Suppress component-only fragments of composite references

A real sequence match can still support the wrong label. When a larger reference is
a known composite, a fragment aligning almost entirely to one of its embedded
components is evidence for the component, not for the complete composite record.

The curated `composite_reference_regions.csv` table records the component interval in
the reference's one-based, inclusive, nucleotide-equivalent subject coordinates. A
fragment is suppressed when no more than three aligned subject bases lie outside the
union of those regions. Whole hits are never affected, and fragments containing four
or more bases of sequence specific to the larger record fail open. This rule does not
require the component itself to win annotation ranking, so its result is stable when
other overlapping hits change.

This is deliberately source-pinned rather than a global coverage threshold. Promoter,
packaging-signal, LTR, origin, and other fragments derived from CDSs remain untouched
unless their exact source record has a reviewed component interval.

### 2. Suppress reviewed low-specificity source-region fragments

Some translated searches repeatedly assign a large protein's name to a short,
noncoding interval even when there is no single containing parent feature. Do not
globally blacklist that protein: a full match or a strong fragment elsewhere may be
real. Instead, `fragment_suppression_regions.csv` pins the exact source accession,
expected subject length, aligned subject interval, and maximum identity of the known
artifact. At least 90% of the aligned subject interval must overlap its curated bounds.
This tolerates modest endpoint variation without capturing a hit substantially shifted
to another part of the source record. Every other geometry or stronger match fails
open.

This is intentionally a manual escape hatch, not a learned identity cutoff. Add a row
only after reproducing a persistent false label and recording why the alignment lacks
biological support.

### 3. Curated accession-pair decisions win among relationship rules

A decision is keyed by the complete tuple
`(parent_db, parent_sseqid, child_db, child_sseqid)`. Names are never keys: names are
not unique and can drift independently of sequence accessions.

A child is suppressed only while its query interval is contained by the matching
parent interval. The same child accession remains detectable elsewhere in the plasmid.
The source-level composite-region rule runs first because a fragment carrying no
sequence evidence beyond a known embedded component cannot recover the larger label
merely from a contextual pair exception.

### 4. Keep structured ncRNAs

Keep nested `ncRNA` calls by default, including calls from custom Infernal sources.
Covariance models detect structured RNA families rather than ordinary nucleotide
motifs. Pair-specific evidence can still suppress a demonstrated exception.

### 5. Keep whole components of gene-level compound cassettes

A whole promoter, CDS, terminator, intron, or poly(A) signal inside a parent typed as
`gene` is a functional component of a compound cassette. Keep both levels. This covers
the `HIS3MX6`, `bleMX6`, `hphMX6`, `kanMX`, `natMX6`, and `patMX4` families without
hardcoding every parent/child combination.

### 6. Preserve exact short functional elements

Short promoters, operators, recombination sites, repeats, tRNAs, and other non-CDS
elements can remain functional even when an installed reference includes extra flank.
Keep a non-CDS fragment when identity is at least 98% and it covers at least 30% of
the child reference. This coverage floor prevents the generic fragment rule from
discarding exact FRT, CRISPR-repeat, and tRNA sequence while still
suppressing tiny motifs such as a 14 bp match to a 324 bp terminator reference.

### 7. Preserve near-complete fragments

A record classified as a fragment can still represent nearly all of the child. Keep a
child covering at least 80% of its reference when either:

- it is a same-source DNA match at 95% or greater identity; or
- it has at least 70% identity and E-value at most `1e-10`.

This protects genuine clipped features such as signal sequences, introns, viral locus
segments, and source-derived components. E-value corroborates coverage and identity;
it never rescues a short low-coverage match on its own.

### 8. Preserve high-confidence CDS-derived fragments

A partial CDS match is not automatically false. Promoters, packaging signals, LTRs,
introns, UTRs, and origins can be cut directly from coding loci. Automatically retain
a different-kind CDS fragment when it spans at least 90 nt (30 amino acids) at 95% or
greater amino-acid identity. Label it explicitly as a fragment or source provenance;
do not imply that the construct encodes a complete functional protein.

The threshold is deliberately conservative. Shorter or more divergent relationships
need an accession-pair decision. Curated decisions also distinguish a genuine locus
relationship reported under the wrong protein accession from a false translated frame.

### 9. Keep strong boundary extensions

A child is not truly contained when its aligned fragment reaches within 3 bp of a
parent boundary and the corresponding unaligned end of the child reference points
through that boundary. Keep such a call when the aligned interval is at least 90 nt,
identity is at least 95%, the reference tail is at least 10% of the child and 30 nt,
and E-value is at most `1e-10`. Strand determines whether the subject prefix or suffix
continues through each query edge. This is an edge-clipped feature candidate, not a
guess based merely on the child reference being longer than the parent.

E-value supports this decision but cannot establish biological meaning by itself. A
wrong-frame homolog can be statistically decisive; conversely, a short exact functional
element can have a less impressive E-value. Accession, frame, boundaries, and parent
semantics remain primary evidence.

### 10. Suppress unsupported contained fragment noise

Suppress a child when all of the following are true:

- the child is wholly contained by a retained parent of a different annotation kind;
- the parent is a whole, high-confidence call;
- the child is a fragment under pLannotate's existing fragment definition;
- it does not pass the exact-element, near-complete, boundary-extension, or
  high-confidence CDS-derived rules; and
- the pair has no curated `keep` exception.

This is the main general rule. In the current audit, 352 of 478 nested calls are
fragments; 310 cover at most 50% of the child reference. For example, the 14 bp
`tbb-2 terminator` hit inside the 45 bp T5 record covers only 4.3% of its 324 bp
reference and should disappear without a feature-specific blacklist.

Do not replace this rule with a raw minimum length. Short complete elements such as
operators and promoters are real; subject coverage distinguishes them from short
fragments of much larger records.

### 11. Keep whole children by default

A match covering the complete child reference is positive sequence evidence. Keep it
by default even when its function in the construct is uncertain or its label is
redundant with the parent. Those concerns require `relabel_child`, `replace_child`,
parent correction, or a display-collapse policy—not deletion of the sequence match.

Different feature types are precisely what nesting is intended to recover. Valid
examples include:

- regulatory RNAs inside replication origins;
- operators or upstream activating sequences inside promoter/operator regions;
- structured RNAs inside viral LTR and packaging regions;
- promoter, CDS, and terminator parts inside a named expression cassette; and
- introns inside genes or engineered expression regions.

Same-type whole children can represent a fusion component or a redundant homolog.
Preserve the relationship and resolve redundancy through ranking or presentation.
A curated pair may still override this rule when independent evidence demonstrates a
broken source record or a specific false match.

### 12. Bounds changes require independent coordinates

Never infer new parent bounds merely by subtracting a child interval. A child may
overlap a functional parent rather than sit in contaminating flank.

Use `trim_parent` only when a primary paper, vendor map, or authoritative source record
gives exact canonical parent coordinates and shows that the current database record
contains extra sequence. Rebuild the search sequence and metadata together. Runtime
coordinate clipping would make alignment statistics and subject coverage dishonest.

### 13. Composite parents are corrected, not treated as false sequence matches

If the parent intentionally contains independently functional children but its name or
type implies an atomic feature, use `replace_parent`:

- rename it as a composite region;
- retype it as `regulatory_region`, `gene`, or another appropriate container; or
- split it into canonical child records and remove the composite record.

Whether the UI later collapses a composite and its children is a presentation policy,
not an annotation-accuracy filter.

## Whole-child adjudication test

For each non-fragment child, record evidence for the following questions:

1. **Sequence:** Is the child near-full-length and high-identity?
2. **Function:** Can the child perform its named function in this parent context?
3. **Source:** Does the installed parent description explicitly mention the child?
4. **External evidence:** Does a primary paper, vendor map, or authoritative database
   place the child in the construct?
5. **Parent semantics:** Is the parent an umbrella/composite region, or is it intended
   to be one atomic feature?
6. **Bounds:** Are exact alternative parent coordinates documented independently?

Apply the following decision table:

| Evidence | Outcome |
|---|---|
| Weak/fragment child inside a whole parent | `suppress_child` by the general rule |
| Whole child with supported independent function | `keep` |
| Whole child supported by sequence/source, but parent is misleadingly atomic | `replace_parent`; keep child |
| Child lies outside independently documented canonical parent bounds | `trim_parent`; keep child |
| Whole child with uncertain function but a real sequence match | `keep`, optionally `relabel_child` |
| Whole child demonstrated to be a false match or broken source record | tuple-specific correction |
| Evidence unresolved | `keep`; optionally queue non-blocking review |

## Anchor decisions

### ColE1-family origins containing RNAI: `keep`

Keep the RNAI calls within the installed `ori`, `p15A_ori`, `ColA_ori`,
`CloDF13_ori`, and `RSF_ori` records. RNAI is an antisense regulator encoded within
ColE1-like replicons; RNAII supplies the primer precursor and RNAI inhibits primer
formation. The nested ncRNA conveys a function that the generic origin label does not.

Evidence:

- Lin-Chao and Bremer, *J. Bacteriol.* 1987, PMID 2434459:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC211922/>
- Eguchi and Tomizawa, *J. Mol. Biol.* 1990:
  <https://doi.org/10.1016/0022-2836(90)90230-J>

### T5 promoter containing lac operator: `replace_parent`; keep lacO

The exact 17/17 bp lac-operator call is not an alignment false positive. The installed
SnapGene record describes itself as a T5 promoter "with embedded lac operator," and
published work describes the engineered sequence as a hybrid T5-lac promoter containing
a T5-like promoter and lac operator.

Preferred correction: rename/redefine the parent as an engineered `T5-lac promoter` or
`T5 promoter/operator`. Keep the lac-operator child because it explains LacI/IPTG
regulation. If product policy requires one visible label, record that separately as a
display-collapse decision rather than claiming the sequence match is false.

The inverse relationship is not equivalent. On pUC19, the complete 17 bp lac operator
causes a 19/45 bp `T5_promoter` fragment spanning subject positions 21–39. Positions
21–37 are exactly the embedded lac operator and the remaining two bases do not provide
meaningful T5-promoter evidence. The composite-region rule therefore suppresses that
fragment while the raw-policy opt-out retains it for audit. A fragment covering the
distinctive T5 portion of the reference remains detectable.

Evidence:

- Ivanov et al., *Microbiologica* 1990, PMID 2352484:
  <https://pubmed.ncbi.nlm.nih.gov/2352484/>

The unrelated 14 bp `tbb-2 terminator` fragment in the same record is
`suppress_child` under the general fragment rule.

### Lac regulatory region producing a penA fragment: source-region suppression

Raw annotation of pUC19 and pCA-mTmG repeatedly reports the same
`swissprot:Q02940` penA fragment over lac regulatory DNA. The alignment covers only
96/939 nucleotide-equivalent subject bases (10.2%), at 68.8% amino-acid identity and
E-value `2.51e-4`. Translating its 32 aligned query codons produces two internal stop
codons, so this interval cannot encode a continuous penA CDS fragment.

Suppress only Q02940 fragments aligning to subject bases 580–675 at no more than 70%
identity, with at least 90% of the aligned subject interval overlapping those bounds,
and only while the installed subject length remains 939. A whole penA call, a
higher-identity fragment, a substantially shifted subject region, or a changed source
record survives. This belongs in `fragment_suppression_regions.csv`, not the global
feature blacklist and not a lac-parent pair override: the hit spans several regulatory
features and is not strictly contained by any one of them.

### Origin-associated P03845/P03846/P03851 calls: global suppression

These three Swiss-Prot records are exact source-pinned suppressions, not a blanket
rule against proteins in origins. Local DIAMOND reproduction against the installed
bundle shows:

- `ori` bases 249–371 align to P03845 but translate as
  `AISRVLPGWTQDDSYRIRRSGRAERGVRAHSPAWSERPTPN`, with no start methionine;
- `CloDF13_ori` bases 264–401 align at 100% identity to P03846 but translate as
  `AFYRAFPGWTQVNSYRIRRSSRAERGVLAYSPAWSERPTPSRDTSV`, again with no start
  methionine; and
- the reverse-strand `ori` interval 333→85 aligns at 100% identity to P03851 but
  translates from `HEPPVQPD...`, with no start methionine.

The underlying sequence homology is real, but these stop-bounded frames do not encode
the asserted ORFs in the origin records. In contrast, P03849 in `oriV` is a whole
feature call and remains annotated.

### ISS: legacy compatibility suppression pending source adjudication

The exact `snapgene:ISS` suppression preserves the behavior of the former Python
blacklist. It is deliberately source-pinned and is not evidence for a general rule
against insertion-sequence annotations. Its current rationale is compatibility and
non-informative overlap behavior, not a completed biological source review; keep it in
the record-correction queue rather than using it as a precedent for new suppressions.

### CMV intron containing T7 promoter: `keep`

The child is an exact 19/19 bp internal match, so the sequence relationship is real
even though the installed CMV-intron description does not mention T7. Keep it under
the whole-child rule. Further source review may improve the representation:

- if its documented intron intentionally contains a functional T7 promoter, `keep` and
  correct the parent description;
- if the T7 sequence is cloning/vector flank outside canonical intron coordinates,
  `trim_parent` and rebuild the record; or
- if the motif is not expected to function in this context, keep the match but label
  that uncertainty separately from sequence detection.

Do not globally suppress `T7_promoter`; it must remain detectable elsewhere.

## Runtime curation table

Global source records known to be false positives live in
`plannotate/data/data/feature_suppressions.csv`, keyed by exact `(db, sseqid)`. This
replaces the former accession-only Python blacklist and prevents a suppression in one
database from leaking to an unrelated record with the same identifier in another.

Reviewed low-specificity fragment intervals live in
`plannotate/data/data/fragment_suppression_regions.csv`. They are unary source-record
rules like composite regions, but they describe a reproduced false alignment rather
than a known embedded component.

| Column | Meaning |
|---|---|
| `db`, `sseqid`, `name` | exact source record and its documentary display name |
| `subject_length` | expected nucleotide-equivalent source length; drift fails open |
| `region_start`, `region_end` | one-based inclusive false-alignment subject bounds |
| `max_identity` | highest fragment identity covered by the reviewed decision |
| `rationale` | why this alignment is an unsupported feature label |
| `source` | stable decision identifier |

Nested parent/child exceptions live in
`plannotate/data/data/nested_feature_overrides.csv`; their keys and actions differ from
global record suppression. Semicolon-separated parent or child accessions compactly
represent one decision that applies to several exact keys.

Component intervals within composite references live in
`plannotate/data/data/composite_reference_regions.csv`. Unlike pair overrides, these
are unary source-record facts: suppression does not depend on a second annotation
surviving overlap resolution. Coordinates are one-based and inclusive in the
nucleotide-equivalent subject coordinate system emitted as `sstart`/`send`.

| Column | Meaning |
|---|---|
| `db`, `sseqid` | exact source record containing the composite interval |
| `region_start`, `region_end` | one-based inclusive component bounds in that record |
| `component_db`, `component_sseqid` | exact source record that explains the interval |
| `rationale` | why the interval does not independently support the composite label |
| `source` | stable decision identifier |

For same-source SnapGene records, `curation_pins.py check` also extracts both sequences
and verifies that the declared interval still equals the complete component sequence.
DIAMOND reports protein subject coordinates to the policy in nucleotide-equivalent
units: amino-acid bounds `a..b` become `(a-1)*3+1 .. b*3`. The pin checker currently
fails closed on non-SnapGene composite rows until equivalent sequence validation is
implemented; do not add a protein-region row without extending that validator.

Composite-region and fragment-region suppression are part of the nested policy.
The raw-policy opt-out retains those fragments for auditing.

The nested override table has this schema:

| Column | Meaning |
|---|---|
| `parent_db` | source database of the containing annotation |
| `parent_sseqid` | containing source accession |
| `child_db` | source database of the nested annotation |
| `child_sseqid` | nested source accession |
| `status` | audit bin: `good`, `bad`, or `review` |
| `action` | `keep`, `suppress_child`, or a source-correction action |
| `rationale` | concise biological/technical reason |
| `source` | stable decision identifier |

The production filter uses only `suppress_child`; correction actions
remain visible guidance and do not silently rewrite database accessions. Pair
overrides apply only when the exact pair occurs in strict containment geometry. The
runtime considers fragment children only, so a curated suppression of a whole child
fails open. If several whole, higher-scoring parents contain a child, suppression
occurs only when every applicable decision says `suppress_child`. Any keep or review
result wins.

Use `apply_nested_policy=False` in the Python API or
`--keep-nested-fragments` on `plannotate batch` to retain raw contained fragment calls
for curation and auditing.

## Maintainer workflow

The audit outputs and runtime policy tables have different ownership. Do not edit a
generated decisions report as though it were policy.

| File | Ownership | Runtime input |
|---|---|---|
| `docs/nested-feature-audit.csv` | generated snapshot of raw nested calls | no |
| `docs/nested-feature-audit.md` | generated human-readable audit snapshot | no |
| `docs/nested-feature-decisions.csv` | generated by applying the current policy to the audit | no |
| `plannotate/data/data/nested_feature_overrides.csv` | manually curated exact pair decisions | yes |
| `plannotate/data/data/composite_reference_regions.csv` | manually curated component intervals inside composite references | yes |
| `plannotate/data/data/fragment_suppression_regions.csv` | manually curated low-specificity source fragment intervals | yes |
| `plannotate/data/data/feature_suppressions.csv` | manually curated global source-record suppressions | yes |
| `figures/nested-feature-audit-linear/` | generated local inspection viewer; gitignored | no |

Regenerate the complete audit whenever the packaged database bundle, search-tool
versions or options, hit normalization, overlap resolution, fragment detection, or
source-feature metadata change. The audit deliberately disables the nested policy so
it records raw candidates instead of feeding the previous decisions back into itself:

```bash
plannotate setupdb
python tools/nested_feature_audit.py \
  --csv docs/nested-feature-audit.csv \
  --markdown docs/nested-feature-audit.md
```

Review the raw calls, then edit only the packaged curation tables as needed. Use exact
database/accession keys, add a concise rationale and stable decision source, and
prefer a `keep` or `review` decision when the biological evidence is uncertain.
Validate that every curated accession still exists in the installed bundle:

```bash
python tools/curation_pins.py check
```

Regenerate the derived decisions report after changing the audit, executable rules,
or any curation table:

```bash
python tools/nested_feature_rules.py
```

For visual inspection, regenerate the local linear viewer. This command recomputes
the decisions from the raw audit, writes the same
`docs/nested-feature-decisions.csv`, and builds the gitignored viewer:

```bash
python tools/nested_feature_figures.py --genes-d3 ~/projects/genes-d3
python -m http.server 8000
```

Then open
`http://127.0.0.1:8000/figures/nested-feature-audit-linear/`. Finish by
running the unit and integration suites and inspecting annotation-control diffs. CI
checks the committed files and accession pins, but does not rerun the exhaustive
database audit; keeping the generated snapshot current is an explicit maintainer
step.

As a follow-up when the source bundle is rebuilt, the audit tooling should enumerate
records whose complete sequence occurs inside another record and emit their derived
intervals as review candidates. Candidate generation must remain separate from the
manually curated region table: substring identity finds possible composite records but
does not establish that the larger label is uninformative on the shared interval.

## Audit completeness warnings

The exhaustive report separately lists source records for which annotation did not
recover an exact full-length self-hit. Those records are still useful as parent
sequences for discovering nested matches, but the missing self-hit is a database/search
provenance warning, not evidence that every child is false. Do not automate deletion
from this warning alone. Fix or retire the source record, rerun the audit, and use an
accession-pair decision only when the child itself has been adjudicated.

## Rollout plan and decision

1. Implement contained-fragment suppression and measure its effect on annotation
   controls. This addresses most accumulated noise without biological pair decisions.
2. Add the version-pinned pair table and seed only adjudicated `keep`/`suppress_child`
   rows.
3. Review whole-child calls only to improve labels, ranking, and parent records; do not
   make review a prerequisite for retaining them.
4. Correct database records marked `replace_parent` or `trim_parent`, rebuild the
   bundle, and rerun the nested-feature audit.
5. Make nested behavior the default while retaining a temporary legacy-flat escape
   hatch for one release. Remove `--detailed` only after control diffs and the curated
   known-suppression and record-correction queues are clean.

PR #84 intentionally supersedes step 5 at maintainer direction: nested behavior is the
only supported behavior, without a one-release legacy-flat window. The remaining
`replace_parent`, `replace_child`, and `relabel_child` decisions are recorded source
curation follow-ups; because those actions fail open at runtime, they do not block the
mode consolidation.
