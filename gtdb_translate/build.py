"""Build translation bundles from raw GTDB metadata + NCBI ``names.dmp``.

This module is the only place that reads the raw GTDB metadata TSVs.  It
makes a *single* streaming pass over each file and collects every vote
needed by all downstream translators:

* NCBI name -> GTDB taxon (flat, see :func:`_finalise_ncbi`)
* SILVA token -> GTDB taxon (rank-keyed, see :func:`_finalise_silva`)
* GTDB taxon -> its full lineage

Keeping construction here means :mod:`gtdb_translate.ncbi` and
:mod:`gtdb_translate.silva` are pure consumers of a finished bundle and
never need pandas or the multi-gigabyte metadata files.

Why the two name dictionaries are keyed differently
---------------------------------------------------
NCBI names are unique across ranks -- a genus and a family never share a
bare name -- so a flat ``name -> taxon`` dict is unambiguous.  SILVA
invents labels outside that system (``Incertae Sedis`` appears at several
ranks, and genus tokens such as ``Lachnospiraceae NK4A136 group`` embed a
family name), so its dict is keyed by ``(rank, token)``.  Rank-keying
also means SILVA needs no string normalisation to avoid rank promotion:
a genus token is only ever looked up against genus-rank votes.
"""

from __future__ import annotations

import csv
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple, Union

from .utils import RANK_ORDER, resolve_votes

logger = logging.getLogger(__name__)

# The GTDB metadata TSVs carry very long free-text fields.
csv.field_size_limit(1 << 31)

#: Metadata column holding the GTDB lineage.
GTDB_COLUMN = "gtdb_taxonomy"
#: Metadata column holding the NCBI lineage.
NCBI_COLUMN = "ncbi_taxonomy"
#: Metadata column holding the NCBI organism name.
NCBI_ORGANISM_COLUMN = "ncbi_organism_name"

#: Metadata columns holding SILVA classifications of recovered rRNA genes.
#: Both have the same layout, so votes from them are pooled by default.
SILVA_COLUMNS: Tuple[str, ...] = (
    "ssu_silva_taxonomy",
    "lsu_silva_23s_taxonomy",
)

#: Values GTDB uses for "no data" in the SILVA columns.
_MISSING = {"", "none", "n/a", "na", "null"}

#: SILVA domains that correspond to something in GTDB.  Anything else is
#: a spurious BLAST hit to a eukaryotic rRNA sequence -- these rows carry
#: plant, fungal or metazoan lineages attached to ordinary bacterial
#: genomes and must not contribute votes.
_PROKARYOTIC_DOMAINS = {"Bacteria", "Archaea"}

#: SILVA tokens that carry no taxonomic information.  Skipped when
#: building votes and rejected at lookup time (see
#: :data:`gtdb_translate.silva.REJECT_TOKENS`, which is a superset).
_UNINFORMATIVE_TOKENS = {
    "incertae sedis",
    "unknown family",
    "uncultured",
    "unidentified",
    "metagenome",
    "ambiguous_taxa",
}

#: SILVA ranks that are voted on.  Index 6 of a SILVA path is the
#: reference sequence's *organism name* (often strain-level, e.g.
#: ``"Klebsiella pneumoniae subsp. pneumoniae DSM 30104"``), not a rank,
#: so it is excluded here and routed to the NCBI dictionary at lookup
#: time instead.
SILVA_VOTED_RANKS: Tuple[str, ...] = RANK_ORDER[:6]

_SPP_RE = re.compile(r"\s+spp?\.?$")


# --------------------------------------------------------------------------
# vote collection
# --------------------------------------------------------------------------


class _Votes:
    """Mutable accumulator for a single streaming pass over the metadata."""

    def __init__(self) -> None:
        # NCBI name -> GTDB taxon -> count.  Includes organism-name votes.
        self.ncbi: Dict[str, Dict[str, int]] = defaultdict(
            lambda: defaultdict(int)
        )
        # Same keys, but only votes from the rank-aligned lineage zip.
        self.rank_aligned: Dict[str, Dict[str, int]] = defaultdict(
            lambda: defaultdict(int)
        )
        # rank -> SILVA token -> GTDB taxon -> count.
        self.silva: Dict[str, Dict[str, Dict[str, int]]] = {
            rank: defaultdict(lambda: defaultdict(int))
            for rank in SILVA_VOTED_RANKS
        }
        # Prefixed GTDB taxon -> its lineage, in first-seen order.
        self.lineages: Dict[str, str] = {}
        self.n_rows = 0
        self.n_silva_rows = 0
        self.n_silva_skipped_domain = 0


def _iter_rows(path: Union[str, Path]) -> Iterator[dict]:
    """Stream a metadata TSV row by row.

    Uses :mod:`csv` rather than pandas: the files have >100 columns and
    hundreds of thousands of rows, only five columns are needed, and
    streaming keeps peak memory flat.
    """
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def _accumulate_ncbi(votes: _Votes, row: dict, gtdb_lineage: List[str]) -> None:
    """Collect NCBI name votes from one metadata row."""
    organism = row.get(NCBI_ORGANISM_COLUMN) or ""
    stripped = organism.strip()
    if stripped and " " in stripped and not _SPP_RE.search(stripped):
        votes.ncbi[organism][gtdb_lineage[-1]] += 1

    ncbi_lineage = str(row.get(NCBI_COLUMN) or "").split(";")
    for ii, (ncbi_tax, gtdb_tax) in enumerate(zip(ncbi_lineage, gtdb_lineage)):
        # ``len(...) == 3`` means a bare prefix such as "g__" (rank absent).
        if len(ncbi_tax) == 3 or len(gtdb_tax) == 3:
            continue
        bare = ncbi_tax[3:]
        # A single-word "species" is a genus name that NCBI placed at the
        # species position; voting on it would map a genus to a species.
        if ii == 6 and " " not in bare.strip():
            continue
        votes.ncbi[bare][gtdb_tax] += 1
        votes.rank_aligned[bare][gtdb_tax] += 1


def _accumulate_silva(
    votes: _Votes,
    row: dict,
    gtdb_lineage: List[str],
    silva_columns: Sequence[str],
) -> None:
    """Collect SILVA token votes from one metadata row.

    Votes are cast per ``(rank, token)`` pair against the GTDB taxon at
    the *same* rank.  SILVA paths in the metadata are a fixed seven
    fields (domain..genus plus an organism name), so positions 0-5 align
    directly with GTDB's first six ranks.
    """
    for column in silva_columns:
        value = str(row.get(column) or "").strip()
        if value.lower() in _MISSING:
            continue

        tokens = [t.strip() for t in value.split(";")]
        if not tokens or tokens[0] not in _PROKARYOTIC_DOMAINS:
            # Eukaryotic hit: the SSU/LSU BLAST matched a plant, fungal or
            # metazoan reference sequence.  Nothing here is usable.
            votes.n_silva_skipped_domain += 1
            continue

        votes.n_silva_rows += 1
        for idx, rank in enumerate(SILVA_VOTED_RANKS):
            if idx >= len(tokens) or idx >= len(gtdb_lineage):
                break
            token = tokens[idx]
            gtdb_tax = gtdb_lineage[idx]
            if not token or token.lower() in _UNINFORMATIVE_TOKENS:
                continue
            if len(gtdb_tax) <= 3:
                continue
            votes.silva[rank][token][gtdb_tax] += 1


def collect_votes(
    metadata_paths: Sequence[Union[str, Path]],
    silva_columns: Sequence[str] = SILVA_COLUMNS,
) -> _Votes:
    """Make one streaming pass over each metadata TSV and collect all votes."""
    votes = _Votes()
    for path in metadata_paths:
        logger.info("Parsing metadata: %s", path)
        for row in _iter_rows(path):
            gtdb_lineage = str(row.get(GTDB_COLUMN) or "").split(";")
            if len(gtdb_lineage) < 2:
                continue
            votes.n_rows += 1

            # Record every GTDB taxon's lineage.  The keys of this dict
            # double as the set of all GTDB taxa, used later to add
            # identity mappings without a second pass over the file.
            for ii, gtdb_tax in enumerate(gtdb_lineage):
                votes.lineages[gtdb_tax] = ";".join(gtdb_lineage[: ii + 1])

            _accumulate_ncbi(votes, row, gtdb_lineage)
            if silva_columns:
                _accumulate_silva(votes, row, gtdb_lineage, silva_columns)

    logger.info(
        "Read %d genomes (%d usable SILVA classifications, "
        "%d eukaryotic hits discarded)",
        votes.n_rows,
        votes.n_silva_rows,
        votes.n_silva_skipped_domain,
    )
    return votes


# --------------------------------------------------------------------------
# finalisation
# --------------------------------------------------------------------------


def _finalise_silva(
    votes: _Votes,
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, list]]]:
    """Resolve SILVA votes into rank-keyed translation and support dicts.

    Support values are ``[total_votes, purity]`` with purity a fraction
    in ``[0, 1]``.
    """
    translations: Dict[str, Dict[str, str]] = {}
    support: Dict[str, Dict[str, list]] = {}

    for rank in SILVA_VOTED_RANKS:
        rank_votes = votes.silva.get(rank) or {}
        if not rank_votes:
            continue
        rank_map: Dict[str, str] = {}
        rank_support: Dict[str, list] = {}
        for token, counter in rank_votes.items():
            target, total, purity = resolve_votes(counter)
            rank_map[token] = target
            rank_support[token] = [total, purity]
        translations[rank] = rank_map
        support[rank] = rank_support
        n_pure = sum(1 for s in rank_support.values() if s[1] == 1.0)
        logger.info(
            "  SILVA %-7s %6d tokens (%d unanimous)",
            rank + ":",
            len(rank_map),
            n_pure,
        )

    return translations, support


def _finalise_ncbi(
    votes: _Votes,
    names_dmp_path: Optional[Union[str, Path]],
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, list]]:
    """Resolve NCBI votes and expand them across ``names.dmp`` synonyms.

    Returns ``(ncbi_name_to_gtdb, ncbi_id_to_scientific, support)``.

    Support obeys one invariant: ``support[k]`` describes the mapping
    stored at ``ncbi_name_to_gtdb[k]``.  Synonyms therefore inherit their
    representative's statistics, because they are stored with the
    representative's target.  Entries with no votes anywhere behind them
    -- GTDB taxa mapping to themselves, and representatives that were
    never voted on -- are absent rather than given a fabricated score.
    """
    ncbi_name_to_gtdb: Dict[str, str] = {}
    support: Dict[str, list] = {}
    for name, counter in votes.ncbi.items():
        target, total, purity = resolve_votes(counter)
        ncbi_name_to_gtdb[name] = target
        support[name] = [total, purity]

    rank_aligned: Dict[str, str] = {
        name: resolve_votes(counter)[0]
        for name, counter in votes.rank_aligned.items()
    }

    # Prefer the rank-aligned answer for single-word capitalised names, so
    # a species-level organism-name vote cannot override a genus mapping.
    for name in list(ncbi_name_to_gtdb):
        if " " not in name and name[0:1].isupper() and name in rank_aligned:
            ncbi_name_to_gtdb[name] = rank_aligned[name]

    # Every GTDB taxon maps to itself.  ``votes.lineages`` is keyed by
    # every prefixed GTDB taxon in first-seen order, so this replaces the
    # second pass over the metadata files that earlier versions made.
    for gtdb_tax in votes.lineages:
        bare = gtdb_tax[3:]
        if bare and bare not in ncbi_name_to_gtdb:
            ncbi_name_to_gtdb[bare] = gtdb_tax

    ncbi_id_to_scientific: Dict[str, str] = {}
    if names_dmp_path is None:
        logger.info(
            "No names.dmp given: skipping synonym expansion "
            "(tax-ID translation will be unavailable)"
        )
        return ncbi_name_to_gtdb, ncbi_id_to_scientific, support

    name_to_id, id_to_scientific = _parse_names_dmp(names_dmp_path)

    protected: Set[str] = set(rank_aligned)

    def _representative(name: str) -> str:
        """Map a name onto the scientific name of its tax ID."""
        if name in name_to_id:
            return id_to_scientific.get(name_to_id[name], name)
        return name

    # Collapse onto representative names, letting rank-aligned mappings win.
    rep_to_gtdb: Dict[str, str] = {}
    rep_support: Dict[str, list] = {}
    for name in protected:
        if name in ncbi_name_to_gtdb:
            rep = _representative(name)
            rep_to_gtdb[rep] = ncbi_name_to_gtdb[name]
            if name in support:
                rep_support[rep] = support[name]
    for name, gtdb_name in ncbi_name_to_gtdb.items():
        if name in protected:
            continue
        rep = _representative(name)
        if rep not in rep_to_gtdb:
            rep_to_gtdb[rep] = gtdb_name
            if name in support:
                rep_support[rep] = support[name]

    # Expand across every names.dmp name.  A synonym is stored with its
    # *representative's* target, so it must also carry the
    # representative's support -- keeping its own tally here would
    # describe a mapping that is not the one stored.  The invariant this
    # maintains is that ``support[k]`` always describes
    # ``ncbi_name_to_gtdb[k]``.
    expanded: Dict[str, str] = {}
    for name, taxid in name_to_id.items():
        rep = id_to_scientific.get(taxid, name)
        target = rep_to_gtdb.get(rep, "none")
        expanded[name] = target
        if name == rep:
            continue
        # ``"none"`` entries are treated as misses at lookup time and can
        # never surface, so they get no support regardless.
        stats = None if target == "none" else rep_support.get(rep)
        if stats is None:
            support.pop(name, None)
        else:
            # Shared, not copied: synonyms of one representative point at
            # the same list, which is what keeps this affordable.
            support[name] = stats
    expanded.update(rep_to_gtdb)

    # Representative names take the support of the vote that reached them.
    support.update(rep_support)

    return (
        expanded,
        {str(k): v for k, v in id_to_scientific.items()},
        support,
    )


def _parse_names_dmp(
    path: Union[str, Path],
) -> Tuple[Dict[str, int], Dict[int, str]]:
    """Parse NCBI ``names.dmp`` into name->taxid and taxid->scientific name.

    Names shared by more than one tax ID are dropped: they cannot be
    resolved without additional context and would otherwise attach an
    arbitrary lineage to an ambiguous string.
    """
    logger.info("Parsing names.dmp: %s", path)
    name_to_id: Dict[str, int] = {}
    id_to_scientific: Dict[int, str] = {}
    ambiguous: Set[str] = set()

    with open(path) as fh:
        for line in fh:
            parts = [x.strip() for x in line.split("|")]
            if len(parts) < 4:
                continue
            name_class = parts[3]
            if name_class not in (
                "scientific name",
                "synonym",
                "equivalent name",
            ):
                continue
            name = parts[1]
            if not name:
                continue
            taxid = int(parts[0])
            if name in name_to_id and name_to_id[name] != taxid:
                ambiguous.add(name)
            name_to_id[name] = taxid
            if name_class == "scientific name":
                id_to_scientific[taxid] = name

    if ambiguous:
        logger.info(
            "Dropped %d ambiguous names.dmp name(s) shared by multiple taxids",
            len(ambiguous),
        )
        for name in ambiguous:
            del name_to_id[name]

    return name_to_id, id_to_scientific


# --------------------------------------------------------------------------
# public entry point
# --------------------------------------------------------------------------


def build_bundle(
    metadata_paths: Sequence[Union[str, Path]],
    names_dmp_path: Optional[Union[str, Path]] = None,
    changelog_path: Optional[Union[str, Path]] = None,
    version: str = "r226",
    silva_columns: Sequence[str] = SILVA_COLUMNS,
) -> dict:
    """Build a complete translation bundle from raw source files.

    Parameters
    ----------
    metadata_paths : sequence of str or Path
        GTDB metadata TSVs, e.g. ``["bac120_metadata_r226.tsv",
        "ar53_metadata_r226.tsv"]``.  Each is read exactly once.
    names_dmp_path : str or Path, optional
        NCBI ``names.dmp``.  Without it, synonym expansion and tax-ID
        translation are unavailable but everything else still builds.
    changelog_path : str or Path, optional
        ``gtdb-taxid-changelog.csv``.  If given, a forward translator is
        built and included.
    version : str
        GTDB release label stored in the bundle.
    silva_columns : sequence of str
        Metadata columns to pool SILVA votes from.  Pass a single-element
        sequence to use only SSU, or an empty sequence to skip SILVA.

    Returns
    -------
    dict
        The bundle, ready for :func:`gtdb_translate.utils.save_bundle`.
    """
    votes = collect_votes(metadata_paths, silva_columns=silva_columns)

    ncbi_name_to_gtdb, ncbi_id_to_scientific, ncbi_support = _finalise_ncbi(
        votes, names_dmp_path
    )
    silva_name_to_gtdb, silva_support = _finalise_silva(votes)

    bundle = {
        "version": version,
        "ncbi_name_to_gtdb": ncbi_name_to_gtdb,
        "ncbi_name_to_gtdb_support": ncbi_support,
        "ncbi_id_to_scientific": ncbi_id_to_scientific,
        "gtdb_name_to_lineage": votes.lineages,
        "silva_name_to_gtdb": silva_name_to_gtdb,
        "silva_name_to_gtdb_support": silva_support,
        "forward_trans_dict": {},
        "forward_rank_dicts": {},
    }

    if changelog_path is not None:
        from .forward import ForwardTranslator

        logger.info("Building forward translator from: %s", changelog_path)
        forward = ForwardTranslator().build(changelog_path)
        bundle["forward_trans_dict"] = forward.translation_dict
        bundle["forward_rank_dicts"] = forward.rank_translation_dicts

    return bundle
