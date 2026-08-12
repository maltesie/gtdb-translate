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

#: Matches any short prefix-shaped token ending in a double underscore --
#: ``sk__`` (SILVA superkingdom), ``k__`` (Greengenes/MetaPhlAn kingdom),
#: ``D_0__`` (SILVA 132 via QIIME 2), ``t__`` (MetaPhlAn strain).  Used to
#: strip prefixes the rank-letter regex above does not recognise, so that
#: unfamiliar schemes still reach the dictionaries as bare names.
_GENERIC_PREFIX_RE = re.compile(r"^[A-Za-z0-9_]{1,4}__")

#: Tokens that carry no taxonomic information in any source.  These are
#: pipeline placeholders rather than taxa, so they are never looked up.
#: Names that merely *look* vague (``uncultured bacterium`` and friends)
#: are deliberately absent: they are often real NCBI taxa with genomes
#: behind them, and the purity threshold separates the usable ones from
#: the arbitrary ones far better than a name list can.
PLACEHOLDER_TOKENS = frozenset(
    {"", "na", "n/a", "nan", "none", "null", "unassigned", "ambiguous_taxa"}
)


def is_placeholder(token: str) -> bool:
    """Return ``True`` if *token* is a pipeline placeholder, not a taxon."""
    return str(token).strip().lower() in PLACEHOLDER_TOKENS


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
    carry a rank.  Other prefix-shaped tokens -- ``sk__`` (SILVA
    superkingdom), ``k__`` (Greengenes/MetaPhlAn), ``D_0__`` (SILVA 132
    via QIIME 2), ``t__`` (MetaPhlAn strain) -- are still stripped, but
    return ``None`` for the rank:

    ``"sk__Bacteria"``       -> ``(None, "Bacteria")``
    ``"D_3__Bacillales"``    -> ``(None, "Bacillales")``
    """
    if not isinstance(name, str):
        return None, name
    stripped = name.strip()
    m = _PREFIX_RE.match(stripped)
    if m:
        return m.group(1), stripped[3:].strip()
    m = _GENERIC_PREFIX_RE.match(stripped)
    if m:
        # Prefix-shaped but not a known rank letter: strip it so the bare
        # name is usable, and leave the rank unknown so callers fall back
        # to positional assignment.
        return None, stripped[m.end():].strip()
    return None, stripped


def strip_rank_prefix(name: str) -> str:
    """Return *name* with any leading rank prefix (``d__`` … ``s__``) removed."""
    _, bare = split_rank_prefix(name)
    return bare


def _cell_is_lineage(
    value: str,
    sep: str,
    resolve,
    lineage_dict: Mapping[str, str],
) -> bool:
    """Return ``True`` if *value* looks like a lineage when split on *sep*.

    A lineage is recognised by its taxa nesting: each resolvable token
    should be an ancestor of the next.  Tokens that resolve to nothing --
    placeholders, ``uncultured``, ranks the dictionary does not cover --
    are skipped rather than counted against the cell, so ``NA``-padded and
    noisy lineages still pass.

    One failed link is tolerated, because vote-derived translations are
    occasionally inconsistent along a chain, but at least one satisfied
    link is required.  Without that floor a two-token cell would pass with
    no evidence at all, and shallow lineages are common in amplicon data.
    """
    tokens = [t for t in str(value).split(sep) if t.strip()]
    if len(tokens) < 2:
        return False

    resolved = []
    for token in tokens:
        _, bare = split_rank_prefix(token)
        if not bare or bare.lower() in PLACEHOLDER_TOKENS:
            continue
        target = resolve(bare)
        if target:
            resolved.append(target)
    if len(resolved) < 2:
        return False

    satisfied = failed = 0
    for parent, child in zip(resolved, resolved[1:]):
        if parent in lineage_dict.get(child, "").split(";"):
            satisfied += 1
        else:
            failed += 1
    return satisfied >= 1 and failed <= 1


def detect_column_format(
    values: Iterable[str],
    sample_size: int = 100,
    resolve=None,
    lineage_dict: Optional[Mapping[str, str]] = None,
    multi_sep: Optional[str] = None,
    candidates: Tuple[str, ...] = (";", "|", ","),
    min_lineage_fraction: float = 0.5,
) -> Dict[str, Any]:
    """Infer separator / prefix / lineage conventions from sample *values*.

    When *resolve* and *lineage_dict* are supplied, a column counts as
    holding lineages only if more than *min_lineage_fraction* of its cells
    pass :func:`_cell_is_lineage` for some candidate separator.  This is
    far more reliable than counting tokens, which cannot tell a lineage
    from a cell holding several taxon names.  Without them the function
    falls back to a token-count heuristic so it stays usable standalone.

    Parameters
    ----------
    values : iterable of str
        Raw column values to inspect.
    sample_size : int
        Maximum number of values to sample.
    resolve : callable, optional
        Maps a bare taxon name to a prefixed GTDB taxon, or ``None``.
    lineage_dict : mapping, optional
        Prefixed GTDB taxon -> full lineage, used for the ancestry test.
    multi_sep : str, optional
        Separator already reserved for multiple entries per cell; excluded
        from the candidates so the two cannot collide.
    candidates : tuple of str
        Separators to consider.
    min_lineage_fraction : float
        Share of cells that must nest for the column to count as lineages.

    Returns
    -------
    dict
        * ``sep`` (str) -- best-guess separator between ranks.
        * ``has_prefix`` (bool) -- whether tokens carry a rank prefix.
        * ``is_full_lineage`` (bool) -- whether cells are lineages.
        * ``from_taxids`` (bool) -- whether values are plain integers.
        * ``lineage_fraction`` (float) -- share of cells that nested, or
          ``None`` when no resolver was given.  Reported so that "nothing
          nested" is distinguishable from "no separator found".
    """
    sample = [
        str(v).strip() for v in values if isinstance(v, str) and v.strip()
    ][:sample_size]

    usable = tuple(c for c in candidates if c != multi_sep)

    if not sample or not usable:
        return {
            "sep": ";",
            "has_prefix": False,
            "is_full_lineage": False,
            "from_taxids": False,
            "lineage_fraction": None,
        }

    best_sep = usable[0]
    lineage_fraction = None
    is_full_lineage = False

    if resolve is not None and lineage_dict is not None:
        scores = {}
        for cand in usable:
            hits = sum(
                _cell_is_lineage(v, cand, resolve, lineage_dict)
                for v in sample
            )
            scores[cand] = hits / len(sample)
        best_sep = max(scores, key=lambda c: (scores[c], -usable.index(c)))
        lineage_fraction = scores[best_sep]
        is_full_lineage = lineage_fraction > min_lineage_fraction
        if not is_full_lineage:
            best_sep = ";" if ";" in usable else usable[0]
    else:
        best_score = -1.0
        token_counts = {}
        for cand in usable:
            counts = [len(v.split(cand)) for v in sample]
            multi_frac = sum(1 for t in counts if t > 1) / len(sample)
            avg_tokens = sum(counts) / len(sample)
            token_counts[cand] = avg_tokens
            score = multi_frac * avg_tokens
            if score > best_score:
                best_score, best_sep = score, cand
        is_full_lineage = token_counts[best_sep] >= 3

    tokens = [
        t.strip() for v in sample for t in v.split(best_sep) if t.strip()
    ] or sample

    from_taxids = sum(t.isdigit() for t in tokens) / len(tokens) > 0.8
    prefix_hits = sum(
        1 for t in tokens if _PREFIX_RE.match(t) or _GENERIC_PREFIX_RE.match(t)
    )
    has_prefix = (prefix_hits / len(tokens)) > 0.5

    return {
        "sep": best_sep,
        "has_prefix": has_prefix,
        "is_full_lineage": is_full_lineage,
        "from_taxids": from_taxids,
        "lineage_fraction": lineage_fraction,
    }


def resolve_votes(counter: Mapping[str, int]) -> Tuple[str, int, float]:
    """Pick the winning target from a vote counter, with support statistics.

    Selection is a plain argmax over vote counts -- no minimum-count or
    minimum-purity threshold is applied at build time, because vote count
    tracks how many genomes a taxon has rather than how good the evidence
    is, and thresholding would preferentially delete rare taxa.  Filtering
    happens at lookup time instead, where the caller can choose it.

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
