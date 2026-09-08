# External Euler replication: dated build snapshot

Read-only inspection and file capture: **2026-09-08 14:13:51 UTC**.
The historical [replication README](../replication/README.md) is unchanged.
Its explicitly referenced scratchpad build was located, and the
[captured build log](../clock-line/replication-snapshot-2026-09-08/build.log)
has progressed to **[8740/9776]**, beyond the earlier [4085/4141] entry.
Lake's changing denominator is a dynamically discovered task count, not a
fixed completion percentage.

The checked-out external repository HEAD is
`d0124689230b58b4f86e7b90ac59de06404b3b6b`. Its tracked source tree was clean;
only the build log, empty `nohup.out`, and local launcher were untracked.
The bounded snapshot preserves the launcher, toolchain, Lake dependency
manifest, Lake configuration, comparator configuration and the PrintAxioms
**input script** alongside the log. The inventory gives their byte counts
and hashes. The full external checkout and compiled dependency cache are
not duplicated here.

The log records an earlier unsupported `-j` option and exit 1 at 06:27:18Z,
followed by the documented restart at 06:34:47Z. At capture there is **no
exit marker for that restart, no successful final build message, and no
PrintAxioms output**. `Challenge.olean` existed while `Solution.olean` did
not. No comparator or nanoda replay output was found in the targeted build
directory. The process-list query was unavailable in the sandbox; the
advancing log supports continuing activity but is not an observed process
completion status.

Accordingly this is an improved provenance record of an **incomplete
independent build**, not a completed Lean theorem/axiom verification or a
comparator replay. No build, axiom command or imported Lean code was run by
the archival collector, and no claim is made about later live progress.

Log snapshot SHA-256:
`eea285232e5fe8b13c281897026ba60384b1930cc7ff1c91d8dad7e9d5d7e217`.
All seven captured files are covered by the
[collected-source manifest](../manifests/collected-sources.json).
