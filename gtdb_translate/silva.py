"""Translate SILVA taxonomy to GTDB.

SILVA is largely NCBI-derived, so most of the work is already done by
:class:`~gtdb_translate.ncbi.NCBITranslator`.  What SILVA adds is a layer
of labels that have no NCBI equivalent -- ``SAR92 clade``, ``SUP05
cluster``, ``NS4 marine group``, ``Lachnospiraceae NK4A136 group``,
``Escherichia-Shigella`` -- and those are resolved from a dictionary
built by voting SILVA classifications of GTDB genomes against those
genomes' GTDB lineages (see :mod:`gtdb_translate.build`).

:class:`SILVATranslator` *composes* an :class:`NCBITranslator` rather
than inheriting from it, for three reasons:

* Ordering.  SILVA-specific resolution has to run *before* the NCBI
  genus fallback, not after it.  Inheritance would put it after.
* Rank awareness.  SILVA lookups are rank-scoped; NCBI lookups are not.
  The two ``translate`` methods are deliberately not substitutable.
* Memory.  Both views share one loaded bundle instead of two copies.

Usage
-----

.. code-block:: python

    from gtdb_translate import SILVATranslator

    t = SILVATranslator.default()
    t.translate(["Bacteria;Bacillota;Bacilli;Lactobacillales;"
                 "Lachnospiraceae;Lachnospiraceae NK4A136 group"],
                full_lineage=True)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

from .ncbi import NCBITranslator
from .utils import RANK_ORDER, RANK_PREFIXES, RANK_TO_PREFIX, split_rank_prefix

logger = logging.getLogger(__name__)

NO_TRANSLATION = "no_translation"

#: Tokens that must never resolve to anything.
#:
#: ``Chloroplast`` and ``Mitochondria`` matter most: SILVA places them
#: inside Cyanobacteriota and Rickettsiales respectively, so without this
#: guard an organellar read would translate to a plausible-looking
#: bacterial lineage -- a wrong answer is worse than no answer.
#:
#: ``Incertae Sedis`` is rejected because it names a placement problem
#: rather than a taxon, and appears at several ranks with no shared
#: meaning.
REJECT_TOKENS = frozenset(
    {
        "",
        "ambiguous_taxa",
        "chloroplast",
        "incertae sedis",
        "metagenome",
        "mitochondria",
        "uncultured",
        "uncultured bacterium",
        "uncultured archaeon",
        "uncultured organism",
        "unidentified",
        "unknown family",
        "wrong",
    }
)

#: Organellar clades.  SILVA nests these inside real bacterial lineages,
#: so they are rejected at whole-lineage level (see
#: :meth:`SILVATranslator.translate_lineage`) rather than per token.
_ORGANELLAR = frozenset({"chloroplast", "mitochondria"})

#: Domains SILVA uses that have no GTDB counterpart.
_NON_PROKARYOTIC_DOMAINS = frozenset({"eukaryota", "unclassified"})

_REJECT_PREFIXES = ("uncultured ", "unknown ", "unclassified ")


def _is_rejected(token: str) -> bool:
    """Return ``True`` if *token* carries no usable taxonomic information."""
    low = token.strip().lower()
    if low in REJECT_TOKENS:
        return True
    return low.startswith(_REJECT_PREFIXES)


def _split_composite(token: str) -> List[str]:
    """Split a SILVA composite genus into its components.

    ``"Escherichia-Shigella"`` -> ``["Escherichia", "Shigella"]``.  Only
    applied when every component looks like a capitalised genus name, so
    hyphenated single genera (``"CL500-3"``, ``"P3OB-42"``,
    ``"Tychonema CCAP 1459-11B"``) are left alone.
    """
    if "-" not in token or " " in token:
        return []
    parts = [p for p in token.split("-") if p]
    if len(parts) < 2:
        return []
    if all(p[:1].isupper() and p.isalpha() for p in parts):
        return parts
    return []


class SILVATranslator:
    """Translate SILVA taxonomy names and lineages to GTDB.

    Parameters
    ----------
    ncbi : NCBITranslator
        A loaded NCBI translator.  Supplies the species-level and
        fallback lookups, and shares the bundle's lineage dictionary.
    silva_name_to_gtdb : dict, optional
        ``{rank: {token: prefixed GTDB taxon}}``.  Taken from the bundle.
    silva_support : dict, optional
        ``{rank: {token: [votes, purity]}}`` with purity in ``[0, 1]``.

    Attributes
    ----------
    ncbi : NCBITranslator
        The wrapped NCBI translator, available for direct use.
    """

    def __init__(
        self,
        ncbi: NCBITranslator,
        silva_name_to_gtdb: Optional[Dict[str, Dict[str, str]]] = None,
        silva_support: Optional[Dict[str, Dict[str, list]]] = None,
    ) -> None:
        self.ncbi = ncbi
        self.silva_name_to_gtdb = silva_name_to_gtdb or {}
        self.silva_support = silva_support or {}

    # -- construction ------------------------------------------------------

    @classmethod
    def load(cls, path: Union[str, Path]) -> "SILVATranslator":
        """Load from a ``.msgpack.zst`` bundle."""
        ncbi = NCBITranslator.load(path)
        return cls.from_ncbi(ncbi)

    @classmethod
    def default(
        cls,
        version: Optional[str] = None,
        cache_dir: Optional[Union[str, Path]] = None,
        force_download: bool = False,
    ) -> "SILVATranslator":
        """Load the pre-built bundle, downloading from GitHub if needed."""
        ncbi = NCBITranslator.default(
            version=version, cache_dir=cache_dir, force_download=force_download
        )
        return cls.from_ncbi(ncbi)

    @classmethod
    def from_ncbi(cls, ncbi: NCBITranslator) -> "SILVATranslator":
        """Wrap an already-loaded :class:`NCBITranslator`.

        Use this to get both translators from a single loaded bundle
        instead of loading the (large) name dictionaries twice.
        """
        return cls(
            ncbi,
            silva_name_to_gtdb=ncbi.silva_name_to_gtdb,
            silva_support=ncbi.silva_support,
        )

    @property
    def version(self) -> str:
        """GTDB release this translator targets."""
        return self.ncbi.version

    # -- single-token lookup ----------------------------------------------

    def lookup_token(
        self,
        token: str,
        rank: Optional[str] = None,
        genus_fallback: bool = False,
    ) -> Optional[Tuple[str, Optional[list]]]:
        """Resolve one SILVA token to a prefixed GTDB taxon.

        Resolution order, stopping at the first hit:

        1. Reject list -- organellar and placeholder tokens fail here
           rather than matching something plausible further down.
        2. The SILVA dictionary at *rank*.  Direct evidence from GTDB
           genomes, so it outranks every string heuristic.
        3. The NCBI dictionary, exact match only.
        4. Components of a composite genus (``Escherichia-Shigella``),
           tried against both dictionaries.
        5. The NCBI genus fallback, only if *genus_fallback* is set.

        Because the SILVA dictionary is rank-keyed, no ``group``/suffix
        stripping is needed: ``Lachnospiraceae NK4A136 group`` is looked
        up against genus-rank votes and cannot be promoted to the family
        ``Lachnospiraceae``.

        Parameters
        ----------
        token : str
            A single SILVA taxon name, with or without a rank prefix.
        rank : str, optional
            Which rank the token sits at.  Without it, all ranks in the
            SILVA dictionary are tried, deepest first.
        genus_fallback : bool
            Enable the NCBI binomial/genus fallback as a last resort.

        Returns
        -------
        tuple of (str, list or None), or None
            ``(prefixed GTDB taxon, [votes, purity])`` with purity a
            fraction in ``[0, 1]``.  Support is ``None`` when the hit did
            not come from a vote-derived entry.
        """
        prefix, bare = split_rank_prefix(str(token).strip())
        if prefix is not None and rank is None:
            rank = RANK_PREFIXES.get(prefix)

        if _is_rejected(bare):
            return None

        hit = self._lookup_silva(bare, rank)
        if hit is not None:
            return hit

        hit = self._lookup_ncbi(bare, genus_fallback=False)
        if hit is not None:
            return hit

        for part in _split_composite(bare):
            hit = self._lookup_silva(part, rank)
            if hit is not None:
                return hit
            hit = self._lookup_ncbi(part, genus_fallback=False)
            if hit is not None:
                return hit

        if genus_fallback:
            return self._lookup_ncbi(bare, genus_fallback=True)

        return None

    def _lookup_silva(
        self, token: str, rank: Optional[str]
    ) -> Optional[Tuple[str, Optional[list]]]:
        """Look *token* up in the SILVA dictionary, at *rank* if known."""
        ranks: Sequence[str]
        if rank is not None:
            ranks = (rank,)
        else:
            # Deepest first: a token that exists at several ranks is most
            # informative at the most specific one.
            ranks = tuple(reversed(RANK_ORDER))

        for candidate_rank in ranks:
            rank_map = self.silva_name_to_gtdb.get(candidate_rank)
            if not rank_map:
                continue
            target = rank_map.get(token)
            if target:
                support = self.silva_support.get(candidate_rank, {}).get(token)
                return target, support
        return None

    def _lookup_ncbi(
        self, token: str, genus_fallback: bool
    ) -> Optional[Tuple[str, Optional[list]]]:
        """Look *token* up in the wrapped NCBI translator."""
        result = self.ncbi.lookup_name(token, genus_fallback=genus_fallback)
        if result is None:
            return None
        target, matched_key = result
        return target, self.ncbi.support_for(matched_key)

    # -- lineage translation ----------------------------------------------

    def parse_lineage(self, lineage: str, sep: str = ";") -> List[Tuple[str, str]]:
        """Parse a SILVA lineage into ``(rank, token)`` pairs.

        Handles both bare SILVA paths (``"Bacteria;Bacillota;..."``) and
        prefixed ones as emitted by QIIME/DADA2 (``"d__Bacteria;..."``).
        Where a token carries a prefix its rank is taken from that;
        otherwise the rank is assigned by position, which is safe because
        SILVA paths are rank-ordered from domain down.

        A seventh field is treated as ``species``.  In GTDB metadata that
        field holds the reference sequence's organism name, often with a
        strain suffix; the species-rank lookup routes it to NCBI, whose
        binomial fallback strips such suffixes.
        """
        raw = [t.strip() for t in str(lineage).split(sep)]
        while raw and not raw[-1]:
            raw.pop()

        parsed: List[Tuple[str, str]] = []
        for idx, token in enumerate(raw):
            prefix, bare = split_rank_prefix(token)
            if prefix is not None:
                rank = RANK_PREFIXES.get(prefix)
            elif idx < len(RANK_ORDER):
                rank = RANK_ORDER[idx]
            else:
                rank = None
            if rank is None or not bare:
                continue
            parsed.append((rank, bare))
        return parsed

    def translate_lineage(
        self,
        lineage: str,
        sep: str = ";",
        genus_fallback: bool = False,
    ) -> Tuple[Optional[str], Optional[list]]:
        """Translate one SILVA lineage to a current GTDB lineage.

        Walks the lineage from its lowest rank upward and returns the
        stored GTDB lineage of the first token that resolves, so the
        result is always a real lineage from the current release.

        Returns
        -------
        tuple of (str or None, list or None)
            ``(gtdb_lineage, [votes, purity])`` with purity in ``[0, 1]``.
        """
        parsed = self.parse_lineage(lineage, sep=sep)
        if not parsed:
            return None, None

        # A eukaryotic path is a spurious rRNA match; refuse it outright
        # rather than letting a lower token match a bacterial name.
        if parsed[0][0] == "domain" and parsed[0][1].lower() in _NON_PROKARYOTIC_DOMAINS:
            return None, None

        # An organellar read is rejected as a whole lineage, not just at
        # the offending token.  SILVA nests Chloroplast inside
        # Cyanobacteriota and Mitochondria inside Rickettsiales, so
        # skipping the token alone would let the walk continue upward and
        # return a real -- but badly wrong -- free-living bacterial
        # lineage for what is actually plant or host DNA.
        if any(token.lower() in _ORGANELLAR for _, token in parsed):
            return None, None

        for rank, token in reversed(parsed):
            hit = self.lookup_token(
                token, rank=rank, genus_fallback=genus_fallback
            )
            if hit is None:
                continue
            target, support = hit
            resolved = self.ncbi.gtdb_name_to_lineage.get(target)
            if resolved:
                return resolved, support
        return None, None

    # -- batch API ---------------------------------------------------------

    def translate(
        self,
        entries: Sequence[str],
        sep: str = ";",
        full_lineage: bool = True,
        genus_fallback: bool = False,
        multi_sep: Optional[str] = None,
        with_support: bool = False,
    ):
        """Translate a batch of SILVA entries.

        Parameters
        ----------
        entries : sequence of str
            SILVA lineages, or bare taxon names when *full_lineage* is
            ``False``.
        sep : str
            Separator between ranks within a lineage (default ``";"``).
        full_lineage : bool
            Treat entries as lineages (default) rather than single names.
        genus_fallback : bool
            Enable the NCBI binomial/genus fallback.
        multi_sep : str, optional
            Splits an entry into several independent values, translated
            separately and rejoined.
        with_support : bool
            Also return per-entry ``(purity, votes)`` columns, with
            purity a fraction in ``[0, 1]``.

        Returns
        -------
        list of str, or tuple of three lists
            Translations, or ``(translations, purity, votes)`` when
            *with_support* is set.  Failed lookups give
            ``"no_translation"`` and empty support strings.
        """
        translations: List[str] = []
        purities: List[str] = []
        supports: List[str] = []

        for entry in entries:
            if not isinstance(entry, str) or not entry.strip():
                translations.append(NO_TRANSLATION)
                purities.append("")
                supports.append("")
                continue

            parts = entry.split(multi_sep) if multi_sep else [entry]
            results: List[str] = []
            part_purity: List[str] = []
            part_votes: List[str] = []

            for part in parts:
                part = part.strip()
                if not part:
                    results.append(NO_TRANSLATION)
                    part_purity.append("")
                    part_votes.append("")
                    continue

                if full_lineage:
                    resolved, support = self.translate_lineage(
                        part, sep=sep, genus_fallback=genus_fallback
                    )
                else:
                    hit = self.lookup_token(
                        part, genus_fallback=genus_fallback
                    )
                    if hit is None:
                        resolved, support = None, None
                    else:
                        target, support = hit
                        resolved = target[3:]

                results.append(resolved if resolved else NO_TRANSLATION)
                if support:
                    part_votes.append(str(support[0]))
                    part_purity.append(str(support[1]))
                else:
                    part_votes.append("")
                    part_purity.append("")

            joiner = multi_sep or ""
            if all(r == NO_TRANSLATION for r in results):
                translations.append(NO_TRANSLATION)
                purities.append("")
                supports.append("")
            else:
                translations.append(
                    joiner.join(results) if multi_sep else results[0]
                )
                purities.append(
                    joiner.join(part_purity) if multi_sep else part_purity[0]
                )
                supports.append(
                    joiner.join(part_votes) if multi_sep else part_votes[0]
                )

        if with_support:
            return translations, purities, supports
        return translations

    # -- introspection -----------------------------------------------------

    def __len__(self) -> int:
        return sum(len(d) for d in self.silva_name_to_gtdb.values())

    def __repr__(self) -> str:
        per_rank = ", ".join(
            f"{rank}={len(self.silva_name_to_gtdb[rank])}"
            for rank in RANK_ORDER
            if rank in self.silva_name_to_gtdb
        )
        return (
            f"SILVATranslator(version={self.version!r}, "
            f"{len(self)} SILVA tokens"
            f"{': ' + per_rank if per_rank else ''})"
        )
