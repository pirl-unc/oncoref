# Bounded assignment regression fixture

This fixture was produced by ORFdata 0.42.2 at commit
`985a93e59c34b1c59dee76ac52b044af4a1c32c8`, using its existing source fixtures and
native importer. The compressed 218 KiB SQLite file exercises the real public
query API and source-specific summary index without downloading the full 24.7 GB
release. `receipt.json` pins both the compressed and expanded bytes.

The original source records are bounded factual slices from the public Pfizer
canine deposit. The producer fixture adds **synthetic** summary quantities,
acquisition ordering, shared memberships and a dangling endpoint. These additions
exercise ambiguity and missingness; they are not biological observations. The
summary representative is `A0A5F4CET6`; its native group ID is `1`. Temporary file
URLs inside the fixture record its original local inputs and are never fetched.
See `PRODUCER-NOTICE.md` for primary-source identity and licenses. This fixture
is not registered as a published dataset and makes no whole-source completeness
claim. Tests install a separate explicit fixture registry with its own byte pins.

Regenerate with a checkout of the pinned producer and a Python 3.11+ environment
containing ORFdata's development dependencies:

```bash
python tests/data/pfizer-assignments/build.py /path/to/orfdata
```

The script imports producer test construction helpers only; no production reader
imports test code, and no network or vendor-database access occurs. New temporary
source URLs change binary hashes, so regeneration updates the receipt as well.
