# Complete numerical endpoint arrays

The NPZ files contain final complex angular coefficient arrays and their
finite-cylinder MAC grid. They preserve phase information for independent
endpoint comparisons. They are numerical checkpoints, not continuous
whole-space solutions or formal certificates.

Large arrays are retained in the local disk archive. The unified GitHub
repository may distribute them as release assets rather than ordinary git
objects; consult its asset inventory for exact download locations. Verify
downloads with the separate `SHA256SUMS.txt` in this directory.

The text/code archive manifest intentionally excludes NPZ payloads. Its
inclusion of this directory's checksum manifest binds the expected payloads
without requiring a large-file checkout for the short verification runner.
