# Collection scope and preservation

This organized archive extends the existing private repository
`goodcarp/forced-route-2026-09`. Its original source bytes are preserved
in place or at the explicit provenance locations below. The new snapshot brings the later core-line passes and the relevant
clock-line source campaign together without altering the active Desktop
research checkout.

## What was collected

| Source | Current archive location | Treatment |
|---|---|---|
| Existing repository at `a37d92e6053a5c9de03cc078d65634cf01807f19` | Original paths and explicit `provenance/` copies | All 601 tracked files preserved byte for byte. The old root README and three documents revised by later collaborator commits have named original-version copies. |
| Core-line research | `core-line/outputs/euler-ns-transfer/` and `core-line/work/` | Reader-facing pass packages plus relevant analytical/code/result working files, with source-relative paths and hashes. |
| Clock-line release-response work | `clock-line/source-campaign/external/alpoge-buckmaster-2026-09-08/` | Original notes, statement audits, scripts and distinct variants. |
| Log-clock development | `clock-line/source-campaign/deepest-think/DTC-2026-09-06/` | `s3close/`, its named `lower/`, `write/`, `gaps/` inputs, the far/near kernel and referee work, and the explicitly referenced `sharp/`, `bfg/`, `verify-astra/`, `toys/` dependency subtrees; stopped and refuted variants remain historical sources. |
| Publication scope reviews | `docs/` | New synthesis and corrections; these do not silently revise historical proofs. |
| Explicitly referenced external Euler build | `clock-line/replication-snapshot-2026-09-08/` | Seven bounded log, launcher and dependency/input files captured at 14:13:51Z; no compiled cache or claimed completed replay. |
| Large numerical fields | Local paths listed in the manifest; separate GitHub release assets | Byte-preserved, independently hashed, excluded from Git blobs. |

The collection deliberately excludes unrelated account/session stores,
other research campaigns, runtime caches, compiled build products and redundant
external PDF/HTML fetches. The source URLs and pins remain in the curated
source records, and earlier external paper copies already preserved in
`strained-core/` are left intact. An omission is not a finding that the
omitted external source is irrelevant to the mathematics.

## Machine-readable record

- [Preservation baseline](../manifests/preserved-base.json): original
  commit, file hashes and explicit original-version preservation locations.
- [Collected source inventory](../manifests/collected-sources.json): each
  copied file's original source-relative path, byte count, SHA-256 and
  Git-versus-release storage classification; also records the snapshot
  time, source roots and intentional fetch omissions.
- [Integrity and local-link report](../manifests/archive-check.json):
  verification of those claims, with historical link warnings kept
  separate from failures in the new reader guides.
- [Assembly comparison](../manifests/clock-line-assembly-comparison.json): the
  original curated, collected and live S3 assembly were checked to have
  identical bytes at the recorded comparison. Later source additions are
  recorded separately in the collected-source snapshot.
- [Large data guide](data-assets.md): release retrieval and exact local
  placement.
- [Pass12 integration record](../manifests/pass12-integration.json): the
  intervening main revision, naming migration, original-version copies,
  intentional top-status changes and frozen new source collection.
- [Pass12 bounded checks](../manifests/pass12-checks.json): historical and
  current package manifests plus four programs rerun in a temporary copy.

The collector hashes each source before copying and verifies both the
copied bytes and the source again afterward. It does not claim a single
atomic filesystem transaction across a live research session. The current
numbered package is added only after the research owner supplies its final
snapshot flag. Files and source hashes, rather than an inferred story
about session authorship, define the collected evidence.

Internal working names, referee labels and historical “proved” headings are
preserved as provenance. They do not imply external peer review or assign
priority. The [current scope review](cross-review-clock-line.md) and
[evidence-status guide](evidence-status.md) distinguish the levels of
support and the open mathematical gaps.

The [dated external-build addendum](external-replication-update.md) records
the separately captured log's status and limitations. The current
[expanded cross-review](cross-review-clock-line-expanded.md) is pinned to
the collected live-draft hash, so later edits do not inherit its review.

## Pass12 integration without rewriting the first snapshot

The original release points to `3dc241bc8f22abbd10fe6bb25d81566f5d3f0d5e`.
Before adding pass12, the local clone was fast-forwarded through the
collaborators' commits to `722cadb3a749874ce0235f54a20c367f7649b9e0`.
Their `codex`→`core-line` and `claude`→`clock-line` naming, front-page
rewrites and new clock sources were retained. The original release and
its two asset URLs were not changed.

Those upstream commits changed the event record, publication page and
historical strained-core README. Their original baseline bytes now live
under `provenance/original-base-a37d92e/`, with the exact locations in
the preservation manifest. The original working contribution-map file
edited during the naming change is similarly retained under
`provenance/first-snapshot-3dc241b/`; its current version is identified
as upstream Git content. These copies were recovered from the pinned
Git objects and checked against their earlier hashes, without reverting
the collaborators' current documents.

Pass2–11 canonical package bytes remain unchanged. Pass12 was copied
only after its owner-frozen flag; the broad live campaign was not
recollected. Its separately reviewed pmax revision is in pass12's own
provenance. The [dated scope update](pass12-scope-update.md) records how
the new results relate to the intervening clock addendum's claims.
