# Isolated shard caching

Set `ONCOREF_SHARD_CACHE_DIR` before loading a sharded dataset to place its
concatenated Parquet cache and signature in a separate directory:

```sh
ONCOREF_SHARD_CACHE_DIR=/tmp/oncoref-replay-cache python scripts/regenerate_cta_walkthrough.py
```

The default remains `~/.cache/oncoref/shard_cache`. An unset or empty override
uses that default. Tilde expansion is supported. The override neither reads,
migrates nor prunes the default cache. It is separate from downloaded source
matrices, reference assets and the in-memory dataframe cache; set it before the
first dataset load in a process.

Cache signatures identify the source directory and each shard's name, byte size
and modification time. These are cache invalidation hints, not scientific
content receipts. An unreadable cache warns and rebuilds from the authoritative
CSV shards. A failed cache write warns but still returns the loaded dataset.
