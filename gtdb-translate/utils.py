"""Shared utility functions for gtdb_translate."""

import gzip
import json
from pathlib import Path
from typing import Any, Dict, Union

import msgpack
import zstandard as zstd

# ---------------------------------------------------------------------------
# GTDB lineage parsing
# ---------------------------------------------------------------------------

RANK_PREFIXES = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def parse_gtdb_lineage(lineage: str) -> dict:
    """Parse a GTDB semicolon-delimited lineage string into a dict of ranks.

    Example input:  ``"d__Bacteria;p__Bacillota;c__Bacilli;..."``
    Returns:        ``{"domain": "Bacteria", "phylum": "Bacillota", ...}``
    """
    result = {}
    for field in lineage.split(";"):
        prefix = field[0]
        name = field[3:]  # skip "X__"
        rank = RANK_PREFIXES.get(prefix, prefix)
        result[rank] = name
    return result


# ---------------------------------------------------------------------------
# JSON persistence (lightweight / human-readable, for single dicts)
# ---------------------------------------------------------------------------

def save_json(d: dict, path: Union[str, Path]) -> None:
    """Save a dictionary to a JSON file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(d, f, indent=2)


def load_json(path: Union[str, Path]) -> dict:
    """Load a dictionary from a JSON file."""
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Bundle persistence  (msgpack + zstandard)
# ---------------------------------------------------------------------------

def save_bundle(data: Dict[str, Any], path: Union[str, Path]) -> None:
    """Serialize *data* to a msgpack + zstandard compressed file.

    Parameters
    ----------
    data : dict
        Arbitrary dict whose values must be msgpack-serialisable
        (dicts, lists, strings, ints, floats, None).
    path : str or Path
        Destination file (conventionally ``*.msgpack.zst``).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = msgpack.packb(data, use_bin_type=True)
    compressor = zstd.ZstdCompressor(level=9)
    compressed = compressor.compress(raw)
    with open(path, "wb") as f:
        f.write(compressed)


def load_bundle(path: Union[str, Path]) -> Dict[str, Any]:
    """Load a msgpack + zstandard compressed file.

    Returns
    -------
    dict
        The deserialised bundle contents.
    """
    path = Path(path)
    with open(path, "rb") as f:
        compressed = f.read()
    decompressor = zstd.ZstdDecompressor()
    raw = decompressor.decompress(compressed)
    return msgpack.unpackb(raw, raw=False)


def load_legacy_gzip_json(path: Union[str, Path]) -> list:
    """Load the legacy ``translation_dicts_rXXX.json.gz`` format.

    Returns the raw list ``[ncbi_name_to_gtdb_mode,
    ncbi_id_to_ncbi_scientific, gtdb_name_to_lineage]``.
    """
    with gzip.open(path, "rt") as f:
        return json.load(f)
