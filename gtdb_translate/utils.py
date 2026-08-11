"""Shared utility functions for gtdb_translate."""

import gzip
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple, Union

import msgpack
import zstandard as zstd


RANK_PREFIXES = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

RANK_TO_PREFIX = {v: k for k, v in RANK_PREFIXES.items()}

RANK_ORDER = ("domain", "phylum", "class", "order", "family", "genus", "species")

_PREFIX_RE = re.compile(r"^([dpcofgs])__")


def parse_gtdb_lineage(lineage: str) -> dict:
    """Parse a GTDB semicolon-delimited lineage string into a dict of ranks.

    Example input:  ``"d__Bacteria;p__Bacillota;c__Bacilli;..."``
    Returns:        ``{"domain": "Bacteria", "phylum": "Bacillota", ...}``
    """
    result = {}
    for field in lineage.split(";"):
        prefix = field[0]
        name = field[3:]
        rank = RANK_PREFIXES.get(prefix, prefix)
        result[rank] = name
    return result


def split_rank_prefix(name: str) -> Tuple[Optional[str], str]:
    """Split a possibly-prefixed taxon name into ``(prefix_letter, bare_name)``.

    ``"s__Escherichia coli"`` -> ``("s", "Escherichia coli")``
    ``"Escherichia coli"``    -> ``(None, "Escherichia coli")``

    Only the standard single-letter GTDB rank prefixes (``d p c o f g s``)
    are recognised; anything else (including ``sk__``, which SILVA uses
    for superkingdom) is treated as unprefixed.
    """
    if not isinstance(name, str):
        return None, name
    stripped = name.strip()
    m = _PREFIX_RE.match(stripped)
    if m:
        return m.group(1), stripped[3:]
    return None, stripped


def strip_rank_prefix(name: str) -> str:
    """Return *name* with any leading rank prefix (``d__`` … ``s__``) removed."""
    _, bare = split_rank_prefix(name)
    return bare


def detect_column_format(
    values: Iterable[str],
    sample_size: int = 100,
) -> Dict[str, Any]:
    """Infer separator / prefix / lineage conventions from sample *values*.

    This is a lightweight heuristic used to pick sensible defaults (e.g.
    for the CLI) — it is not a guarantee of correctness, and any setting
    it infers can always be overridden explicitly.

    Parameters
    ----------
    values : iterable of str
        Raw column values to inspect (e.g. ``df[col].dropna()``).
    sample_size : int
        Maximum number of values to sample.

    Returns
    -------
    dict
        * ``sep`` (str) — best-guess separator between multiple entries
          in a single cell (one of ``";"``, ``"|"``, ``","``).
        * ``has_prefix`` (bool) — whether tokens carry a rank prefix
          like ``s__``.
        * ``is_full_lineage`` (bool) — whether cells look like multi-rank
          lineages rather than single taxon names.
        * ``from_taxids`` (bool) — whether values are plain integers
          (NCBI tax IDs) rather than names.
    """
    sample = [
        str(v).strip() for v in values if isinstance(v, str) and v.strip()
    ][:sample_size]

    if not sample:
        return {
            "sep": "|",
            "has_prefix": False,
            "is_full_lineage": False,
            "from_taxids": False,
        }

    candidates = (";", "|", ",")
    best_sep, best_score = "|", -1.0
    scores = {}
    for cand in candidates:
        token_counts = [len(v.split(cand)) for v in sample]
        multi_frac = sum(1 for t in token_counts if t > 1) / len(sample)
        avg_tokens = sum(token_counts) / len(sample)
        scores[cand] = (multi_frac, avg_tokens)
        score = multi_frac * avg_tokens
        if score > best_score:
            best_score, best_sep = score, cand

    tokens = [t.strip() for v in sample for t in v.split(best_sep) if t.strip()]
    if not tokens:
        tokens = sample

    from_taxids = sum(t.isdigit() for t in tokens) / len(tokens) > 0.8

    prefix_hits = sum(1 for t in tokens if split_rank_prefix(t)[0] is not None)
    has_prefix = (prefix_hits / len(tokens)) > 0.5

    avg_tokens_best = scores[best_sep][1]
    is_full_lineage = avg_tokens_best >= 3

    return {
        "sep": best_sep,
        "has_prefix": has_prefix,
        "is_full_lineage": is_full_lineage,
        "from_taxids": from_taxids,
    }


def resolve_votes(counter: Mapping[str, int]) -> Tuple[str, int, float]:
    """Pick the winning target from a vote counter, with support statistics.

    Selection is a plain argmax over vote counts -- no minimum-count or
    minimum-purity threshold is applied, because vote count tracks how
    many genomes a taxon has rather than how good the evidence is, and
    thresholding would preferentially delete rare taxa.  Callers that
    care can filter on the returned statistics instead.

    Ties are broken by name so that two builds from the same input
    produce byte-identical bundles.

    Parameters
    ----------
    counter : mapping of str to int
        Candidate target name -> number of votes.

    Returns
    -------
    tuple of (str, int, float)
        ``(target, total_votes, purity)`` where *purity* is the winner's
        share of all votes as a fraction in ``[0, 1]`` -- 1.0 when the
        vote was unanimous.  Rounded to four decimal places so that
        bundles stay compact and byte-identical across rebuilds.

        A unanimous single vote gives ``(target, 1, 1.0)`` -- which is
        why both numbers are reported and not just the purity.
    """
    target = max(counter.items(), key=lambda kv: (kv[1], kv[0]))[0]
    total = sum(counter.values())
    purity = round(counter[target] / total, 4) if total else 0.0
    return target, total, purity


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
