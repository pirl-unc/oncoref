# Licensed under the Apache License, Version 2.0.
import bz2
import gzip
import io
import lzma
import tarfile
import zipfile

import pandas as pd
import pytest

from oncoref._table_io import read_csv
from oncoref.expression_builders import read_source_expression_matrix
from oncoref.load_dataset import _read_csv_for_cache


def table_bytes(text, compression):
    data = text.encode("utf-8")
    if compression == "gzip":
        return gzip.compress(data)
    if compression == "bz2":
        return bz2.compress(data)
    if compression == "xz":
        return lzma.compress(data)
    if compression == "zstd":
        return pytest.importorskip("zstandard").ZstdCompressor().compress(data)
    if compression == "zip":
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr("table.tsv", data)
        return buffer.getvalue()
    if compression and compression.startswith("tar"):
        buffer = io.BytesIO()
        mode = "w:" + compression.removeprefix("tar").lstrip(".")
        with tarfile.open(fileobj=buffer, mode=mode) as archive:
            member = tarfile.TarInfo("table.tsv")
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
        return buffer.getvalue()
    assert compression is None
    return data


@pytest.mark.parametrize(
    "compression", [None, "gzip", "bz2", "xz", "zstd", "zip", "tar", "tar.gz", "tar.bz2", "tar.xz"]
)
@pytest.mark.parametrize("filename", ["override", "override.tsv", "override.tsv.gz"])
@pytest.mark.parametrize("chunked", [False, True])
def test_table_reads_use_content_and_preserve_parsing_options(
    tmp_path, compression, filename, chunked
):
    path = tmp_path / filename
    path.write_bytes(table_bytes("gene\tsample\ng1\t0\ng2\tNA\ng3\t\n", compression))
    kwargs = {"sep": "\t", "dtype": str, "keep_default_na": False}
    if chunked:
        with read_csv(path, chunksize=1, **kwargs) as chunks:
            frame = pd.concat(chunks, ignore_index=True)
    else:
        frame = read_csv(path, **kwargs)
    assert frame.to_dict("list") == {"gene": ["g1", "g2", "g3"], "sample": ["0", "NA", ""]}


@pytest.mark.parametrize("compression,filename", [("gzip", "override"), (None, "override.gz")])
def test_shared_dataset_readers_use_content(tmp_path, compression, filename):
    path = tmp_path / filename
    path.write_bytes(table_bytes("gene,sample\ng1,0\ng2,2\n", compression))

    source = read_source_expression_matrix(path, sep=",")
    assert source.to_dict("list") == {"gene": ["g1", "g2"], "sample": ["0", "2"]}
    assert source.attrs["row_id_col"] == "gene"
    cached = _read_csv_for_cache(path, "fixture")
    assert cached.to_dict("list") == {"gene": ["g1", "g2"], "sample": [0, 2]}


def test_content_detection_preserves_dataset_specific_na_policy(tmp_path):
    path = tmp_path / "override"
    path.write_bytes(gzip.compress(b"alias,official_symbol\nNA,XK\nNaN,SCN11A\n,EMPTY\n"))
    frame = _read_csv_for_cache(path, "ncbi-symbol-synonyms")
    assert frame["alias"].iloc[:2].tolist() == ["NA", "NaN"]
    assert pd.isna(frame["alias"].iloc[2])


def test_corrupt_gzip_raises_without_falling_back_to_plain_text(tmp_path):
    path = tmp_path / "override"
    path.write_bytes(gzip.compress(b"gene\tsample\ng1\t0\n")[:-8])
    with pytest.raises(EOFError):
        read_csv(path, sep="\t")
