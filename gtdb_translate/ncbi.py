"""Translate NCBI taxonomy names and tax IDs to GTDB taxonomy.

This module provides :class:`NCBITranslator`, which wraps the three
translation dictionaries built from GTDB metadata and NCBI ``names.dmp``,
plus an optional :class:`~gtdb_translate.forward.ForwardTranslator` for
resolving renamed GTDB species across releases.

This module does not handle SILVA taxonomy strings; SILVA-classified
input needs to be converted to NCBI-style names/lineages separately
before translation.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

from .forward import ForwardTranslator
from .utils import (
    detect_column_format,
    load_bundle,
    load_legacy_gzip_json,
    save_bundle,
    strip_rank_prefix,
)

logger = logging.getLogger(__name__)


class NCBITranslator:
    """Translate NCBI taxonomy to GTDB.

    The translator holds three dictionaries (built from GTDB metadata TSVs
    and NCBI ``names.dmp``) and an optional :class:`ForwardTranslator`:

    * ``ncbi_name_to_gtdb`` — any NCBI name → best GTDB taxon
      (prefixed, e.g. ``"s__Bacillus subtilis"``).
    * ``ncbi_id_to_scientific`` — NCBI tax-ID (str) → scientific name.
    * ``gtdb_name_to_lineage`` — prefixed GTDB taxon → partial lineage
      string.
    * ``forward`` — optional :class:`ForwardTranslator` targeting the same
      GTDB version.

    Parameters
    ----------
    version : str
        GTDB release version this translator targets (e.g. ``"r226"``).
    """

    def __init__(self, version: str = "r226") -> None:
        self.version = version
        self.ncbi_name_to_gtdb: Dict[str, str] = {}
        self.ncbi_id_to_scientific: Dict[str, str] = {}
        self.gtdb_name_to_lineage: Dict[str, str] = {}
        self.forward: Optional[ForwardTranslator] = None
        # Vote statistics, present only for vote-derived entries.  Synonym
        # expansions and GTDB-identity mappings are not votes and are
        # deliberately absent rather than given a fabricated score.
        self.ncbi_support: Dict[str, list] = {}
        # Carried through so SILVATranslator can share one loaded bundle.
        self.silva_name_to_gtdb: Dict[str, Dict[str, str]] = {}
        self.silva_support: Dict[str, Dict[str, list]] = {}

    def support_for(self, name: Optional[str]) -> Optional[list]:
        """Return ``[votes, purity]`` for *name*, or ``None``.

        *purity* is the winning share of the vote as a fraction in
        ``[0, 1]``; 1.0 means the vote was unanimous.

        The statistics always describe the mapping actually stored for
        *name*.  A ``names.dmp`` synonym is stored with its
        representative scientific name's target, so it reports that
        representative's votes and purity -- the evidence genuinely is
        the representative's.

        ``None`` means no votes lie behind the mapping at all: a GTDB
        taxon mapping to itself, or a name whose representative was never
        voted on.
        """
        if name is None:
            return None
        return self.ncbi_support.get(name)

    @classmethod
    def build(
        cls,
        metadata_paths,
        names_dmp_path=None,
        changelog_path=None,
        version: str = "r226",
        silva_columns=None,
    ) -> "NCBITranslator":
        """Build a translator from raw GTDB + NCBI files.

        Thin wrapper around :func:`gtdb_translate.build.build_bundle`,
        which does the actual single-pass parsing.  Construction lives
        there so that this module stays a pure consumer of a finished
        bundle.
        """
        from .build import SILVA_COLUMNS, build_bundle

        bundle = build_bundle(
            metadata_paths=metadata_paths,
            names_dmp_path=names_dmp_path,
            changelog_path=changelog_path,
            version=version,
            silva_columns=(
                SILVA_COLUMNS if silva_columns is None else silva_columns
            ),
        )
        return cls.from_bundle(bundle)

    def translate(
        self,
        entries: Sequence[str],
        sep: str = "|",
        full_lineage: bool = False,
        genus_fallback: bool = False,
        multi_sep: Optional[str] = None,
        with_support: bool = False,
    ):
        """Translate NCBI names to GTDB.

        Parameters
        ----------
        entries : sequence of str
            Each entry is one or more taxon names joined by *sep*
            (non-full-lineage mode), or one or more full lineage
            strings (full-lineage mode; multiple lineages per entry
            require *multi_sep*).
        sep : str
            Separator between multiple names within a single entry
            when *full_lineage* is ``False``. When *full_lineage* is
            ``True``, this is instead the separator between ranks
            within a single lineage.
        full_lineage : bool
            If ``True``, treat each entry (or, if *multi_sep* is given,
            each *multi_sep*-separated part of an entry) as a complete
            lineage and return the best matching GTDB lineage for it.
        genus_fallback : bool
            If ``True``, fall back to binomial and genus-level
            lookups when an exact match fails.  If ``False`` (default),
            only try exact match, bracket removal, and parenthetical
            removal.
        multi_sep : str, optional
            Only used when *full_lineage* is ``True``. If given, splits
            each entry into one or more independent lineages on
            *multi_sep*, translates them separately, and re-joins the
            results with *multi_sep* (an entry where every lineage
            fails becomes a single ``"no_translation"``). If ``None``
            (default), each entry is treated as a single lineage.

        with_support : bool
            If ``True``, also return per-entry purity and vote-count
            columns.  Purity is a fraction in ``[0, 1]``.  Values are
            empty for entries with no votes behind them (GTDB taxa that
            map to themselves, and names whose representative was never
            voted on).

        Returns
        -------
        list of str, or tuple of three lists
            Translated entries, or ``(translations, purity, votes)``
            when *with_support* is set.  ``"no_translation"`` when no
            mapping was found.
        """
        NO_TRANS = "no_translation"
        translations = [NO_TRANS] * len(entries)
        purities = [""] * len(entries)
        supports = [""] * len(entries)

        def _fmt(stats):
            """Render one entry's support stats as (purity, votes) strings."""
            if not stats:
                return "", ""
            return str(stats[1]), str(stats[0])

        for i, entry in enumerate(entries):
            if not isinstance(entry, str):
                continue

            if full_lineage:
                lineages = entry.split(multi_sep) if multi_sep else [entry]
                results, part_pur, part_vot = [], [], []
                for lineage in lineages:
                    lineage = lineage.strip()
                    if not lineage:
                        results.append(NO_TRANS)
                        part_pur.append("")
                        part_vot.append("")
                        continue
                    result, stats = self._translate_single_lineage(
                        lineage, sep=sep, genus_fallback=genus_fallback
                    )
                    results.append(result if result is not None else NO_TRANS)
                    pur, vot = _fmt(stats)
                    part_pur.append(pur)
                    part_vot.append(vot)
                if all(r == NO_TRANS for r in results):
                    continue
                if multi_sep:
                    translations[i] = multi_sep.join(results)
                    purities[i] = multi_sep.join(part_pur)
                    supports[i] = multi_sep.join(part_vot)
                else:
                    translations[i] = results[0]
                    purities[i] = part_pur[0]
                    supports[i] = part_vot[0]
            else:
                taxa = entry.split(sep)[::-1]
                parts, part_pur, part_vot = [], [], []
                for tax in taxa:
                    hit = self.lookup_name(tax, genus_fallback=genus_fallback)
                    if hit is not None and hit[0] != "none":
                        parts.append(hit[0][3:])
                        pur, vot = _fmt(self.support_for(hit[1]))
                    else:
                        parts.append(NO_TRANS)
                        pur, vot = "", ""
                    part_pur.append(pur)
                    part_vot.append(vot)
                parts = parts[::-1]
                if all(p == NO_TRANS for p in parts):
                    continue
                translations[i] = sep.join(parts)
                purities[i] = sep.join(part_pur[::-1])
                supports[i] = sep.join(part_vot[::-1])

        if with_support:
            return translations, purities, supports
        return translations

    def _translate_single_lineage(
        self,
        lineage: str,
        sep: str,
        genus_fallback: bool,
    ) -> Tuple[Optional[str], Optional[list]]:
        """Translate a single full lineage string.

        Tries the lowest (most specific) rank first and works up toward
        domain, returning the first GTDB match found. Returns ``None``
        if no rank in the lineage resolves to a known GTDB taxon.

        The lineage does not need to start at domain -- a phylum-first
        lineage (e.g. one that's already had its domain/kingdom prefix
        stripped) is tried at every rank it does contain. Non-prokaryotic
        input isn't special-cased; it's expected to simply fail to match
        anything in ``ncbi_name_to_gtdb`` (which is built exclusively
        from prokaryotic GTDB/NCBI data) and fall through to ``None``.
        """
        taxa = lineage.split(sep)[::-1]
        is_domain_anchored = taxa[-1].endswith("Bacteria") or taxa[-1].endswith("Archaea")
        for ii, tax in enumerate(taxa):
            hit = self.lookup_name(tax, genus_fallback=genus_fallback)
            if hit is None:
                continue
            gtdb_name, matched_key = hit
            if gtdb_name in self.gtdb_name_to_lineage:
                best = self.gtdb_name_to_lineage[gtdb_name]
                if is_domain_anchored and best.count(";") >= len(taxa) - ii:
                    best = ";".join(best.split(";")[: len(taxa) - ii])
                return best, self.support_for(matched_key)
        return None, None

    def translate_ids(
        self,
        entries: Sequence[str],
        sep: str = "|",
        full_lineage: bool = False,
        with_support: bool = False,
    ):
        """Translate NCBI tax IDs to GTDB names.

        Parameters
        ----------
        entries : sequence of str
            Each entry is one or more NCBI tax IDs joined by *sep*.
        sep : str
            Separator between multiple IDs within a single entry.
        full_lineage : bool
            If ``True``, return full GTDB lineage strings.
        with_support : bool
            If ``True``, also return purity and vote-count columns.

        Returns
        -------
        list of str, or tuple of three lists
        """
        NO_TRANS = "no_translation"
        translations = [NO_TRANS] * len(entries)
        purities = [""] * len(entries)
        supports = [""] * len(entries)
        for i, taxs in enumerate(entries):
            if not isinstance(taxs, str):
                continue
            parts, part_pur, part_vot = [], [], []
            for taxid in taxs.split(sep):
                sci = self.ncbi_id_to_scientific.get(taxid.strip())
                if sci is None:
                    parts.append(NO_TRANS)
                    part_pur.append("")
                    part_vot.append("")
                    continue
                gtdb_name = self.ncbi_name_to_gtdb.get(sci, "none")
                if gtdb_name == "none":
                    parts.append(NO_TRANS)
                    part_pur.append("")
                    part_vot.append("")
                    continue
                stats = self.support_for(sci)
                part_pur.append(str(stats[1]) if stats else "")
                part_vot.append(str(stats[0]) if stats else "")
                if full_lineage and gtdb_name in self.gtdb_name_to_lineage:
                    gtdb_name = self.gtdb_name_to_lineage[gtdb_name]
                elif not full_lineage:
                    gtdb_name = strip_rank_prefix(gtdb_name)
                parts.append(gtdb_name)
            if all(p == NO_TRANS for p in parts):
                translations[i] = NO_TRANS
            else:
                translations[i] = sep.join(parts)
                purities[i] = sep.join(part_pur)
                supports[i] = sep.join(part_vot)

        if with_support:
            return translations, purities, supports
        return translations

    def lookup_name(
        self, name: str, genus_fallback: bool = False
    ) -> Optional[Tuple[str, str]]:
        """Look up a single NCBI name with progressive fallbacks.

        Tries, in order:

        0. Strip bare ``sp.``/``sp``/``spp.``/``spp`` suffix
           (``Pseudomonas sp.`` -> ``Pseudomonas``)
        1. Exact match
        2. Bracket removal  (``[Clostridium]`` -> ``Clostridium``)
        3. Parenthetical removal  (``Klebsiella pneumoniae (resistant)``
           -> ``Klebsiella pneumoniae``)
        4. Replace the first ``_`` with a space (``Escherichia_coli``
           -> ``Escherichia coli``). Simple and doesn't catch every
           underscore-joined format, but exact match on the untouched
           name is always tried first, so this doesn't affect names
           that already match as given.

        If *genus_fallback* is ``True``, also tries:

        5. Binomial (first two words) -- strips strain IDs like
           ``Acinetobacter baumannii AB03``
           -> ``Acinetobacter baumannii``
           Skipped when the second word is ``sp.``, ``sp``, ``spp.``,
           or ``spp`` (these are genus-level, not true binomials).
        6. Genus only  (``Bombilactobacillus sp.``
           -> ``Bombilactobacillus``)

        At each step, ``"none"`` values in the dictionary are treated as
        misses so that the fallback chain continues.

        Returns
        -------
        tuple of (str, str), or None
            ``(prefixed GTDB taxon, matched dictionary key)``.  The key
            is returned so callers can retrieve the vote statistics
            behind the hit via :meth:`support_for`.
        """
        import re

        def _get(key: str) -> Optional[Tuple[str, str]]:
            """Return (value, key) if the key exists and isn't 'none'."""
            val = self.ncbi_name_to_gtdb.get(key)
            if val is not None and val != "none":
                return val, key
            return None

        name = strip_rank_prefix(name.strip())
        name = re.sub(r"\s+spp?\.?$", "", name.strip())

        result = _get(name)
        if result is not None:
            return result

        if "[" in name:
            cleaned = name.replace("[", "").replace("]", "")
            result = _get(cleaned)
            if result is not None:
                return result
            name = cleaned

        if "(" in name:
            stripped = re.sub(r"\([^)]*\)", "", name).strip()
            stripped = re.sub(r"\s+", " ", stripped)
            result = _get(stripped)
            if result is not None:
                return result
            name = stripped

        if "_" in name:
            name = name.replace("_", " ", 1)
            result = _get(name)
            if result is not None:
                return result

        if not genus_fallback:
            return None

        _SP_WORDS = {"sp.", "sp", "spp.", "spp"}
        words = name.split()
        if len(words) > 2 and words[1] not in _SP_WORDS:
            binomial = f"{words[0]} {words[1]}"
            result = _get(binomial)
            if result is not None:
                return result

        if len(words) >= 1:
            result = _get(words[0])
            if result is not None:
                return result

        return None

    def _lookup_name(
        self, name: str, genus_fallback: bool = False
    ) -> Optional[str]:
        """Return only the GTDB taxon from :meth:`lookup_name`."""
        result = self.lookup_name(name, genus_fallback=genus_fallback)
        return result[0] if result else None

    def to_bundle(self) -> dict:
        """Return this translator's contents as a bundle dict."""
        return {
            "version": self.version,
            "ncbi_name_to_gtdb": self.ncbi_name_to_gtdb,
            "ncbi_name_to_gtdb_support": self.ncbi_support,
            "ncbi_id_to_scientific": self.ncbi_id_to_scientific,
            "gtdb_name_to_lineage": self.gtdb_name_to_lineage,
            "silva_name_to_gtdb": self.silva_name_to_gtdb,
            "silva_name_to_gtdb_support": self.silva_support,
            "forward_trans_dict": (
                self.forward.translation_dict if self.forward else {}
            ),
            "forward_rank_dicts": (
                self.forward.rank_translation_dicts if self.forward else {}
            ),
        }

    @classmethod
    def from_bundle(cls, data: dict) -> "NCBITranslator":
        """Construct a translator from an in-memory bundle dict."""
        obj = cls(version=data.get("version", "unknown"))
        obj.ncbi_name_to_gtdb = data["ncbi_name_to_gtdb"]
        obj.ncbi_support = data.get("ncbi_name_to_gtdb_support", {})
        obj.ncbi_id_to_scientific = data["ncbi_id_to_scientific"]
        obj.gtdb_name_to_lineage = data["gtdb_name_to_lineage"]
        obj.silva_name_to_gtdb = data.get("silva_name_to_gtdb", {})
        obj.silva_support = data.get("silva_name_to_gtdb_support", {})

        fwd_dict = data.get("forward_trans_dict", {})
        fwd_rank_dicts = data.get("forward_rank_dicts", {})
        if fwd_dict or fwd_rank_dicts:
            obj.forward = ForwardTranslator()
            obj.forward._trans_dict = fwd_dict
            obj.forward._rank_trans_dicts = fwd_rank_dicts
        return obj

    def save(self, path: Union[str, Path]) -> None:
        """Save all dictionaries as a single ``.msgpack.zst`` bundle.

        The bundle contains:

        * ``ncbi_name_to_gtdb`` and ``ncbi_name_to_gtdb_support``
        * ``ncbi_id_to_scientific``
        * ``gtdb_name_to_lineage``
        * ``silva_name_to_gtdb`` and ``silva_name_to_gtdb_support``
        * ``forward_trans_dict`` / ``forward_rank_dicts``
        * ``version``
        """
        save_bundle(self.to_bundle(), path)

    @classmethod
    def load(cls, path: Union[str, Path]) -> "NCBITranslator":
        """Load a translator from a ``.msgpack.zst`` bundle.

        Also supports the legacy ``translation_dicts_rXXX.json.gz``
        format, which carries neither support statistics nor a forward
        or SILVA dictionary.
        """
        path = Path(path)

        if path.suffixes[-2:] == [".json", ".gz"] or path.suffix == ".gz":
            dicts = load_legacy_gzip_json(path)
            obj = cls()
            obj.ncbi_name_to_gtdb = dicts[0]
            obj.ncbi_id_to_scientific = {str(k): v for k, v in dicts[1].items()}
            obj.gtdb_name_to_lineage = dicts[2]
            return obj

        return cls.from_bundle(load_bundle(path))

    @classmethod
    def default(
        cls,
        version: Optional[str] = None,
        cache_dir: Optional[Union[str, Path]] = None,
        force_download: bool = False,
    ) -> "NCBITranslator":
        """Load the pre-built translator, downloading from GitHub if needed.

        The bundle is fetched from
        `github.com/maltesie/gtdb-translate/releases <https://github.com/maltesie/gtdb-translate/releases>`_.

        Parameters
        ----------
        version : str, optional
            GTDB release version (e.g. ``"r226"``).  If ``None`` (the
            default), the latest release is resolved automatically.
        cache_dir : str or Path, optional
            Override default cache directory (``~/.cache/gtdb_translate/``).
        force_download : bool
            Re-download even if a cached file exists.

        Returns
        -------
        NCBITranslator
            Ready-to-use translator targeting *version*.
        """
        from .data import ensure_bundle

        bundle_path = ensure_bundle(
            version=version, cache_dir=cache_dir, force=force_download
        )
        return cls.load(bundle_path)

    def detect_column(
        self,
        df,
        sample_rows: int = 100,
        sep: str = "|",
    ) -> Optional[str]:
        """Heuristically detect which column in *df* contains translatable names.

        Scores each string column by how many of its sampled values have
        at least one taxon covered by the translation dictionary.

        Parameters
        ----------
        df : pandas DataFrame
            Input table.
        sample_rows : int
            Number of rows to sample for scoring.
        sep : str
            Separator for multi-valued cells.

        Returns
        -------
        str or None
            Column name with the highest coverage, or ``None`` if
            no column scores above zero.
        """
        import pandas as pd

        best_col = None
        best_score = 0
        sample = df.head(sample_rows)
        for col in sample.columns:
            if not pd.api.types.is_string_dtype(sample[col]):
                continue
            hits = 0
            total = 0
            for val in sample[col].dropna():
                total += 1
                names = str(val).split(sep)
                if any(n in self.ncbi_name_to_gtdb for n in names):
                    hits += 1
            score = hits / total if total > 0 else 0
            if score > best_score:
                best_score = score
                best_col = col
        if best_score > 0:
            logger.info(
                "Auto-detected column '%s' (%.0f%% coverage)",
                best_col,
                best_score * 100,
            )
        return best_col

    def detect_format(
        self,
        df,
        column: str,
        sample_rows: int = 100,
    ) -> Dict[str, object]:
        """Infer separator / prefix / lineage conventions for *column*.

        Thin wrapper around :func:`gtdb_translate.utils.detect_column_format`
        that samples *column* from *df*. See that function for the meaning
        of the returned dict.

        Parameters
        ----------
        df : pandas DataFrame
            Input table.
        column : str
            Name of the column to inspect.
        sample_rows : int
            Number of rows to sample.
        """
        values = df[column].dropna().astype(str).tolist()[:sample_rows]
        detected = detect_column_format(values, sample_size=sample_rows)
        logger.info("Auto-detected format for '%s': %s", column, detected)
        return detected

    def __len__(self) -> int:
        return len(self.ncbi_name_to_gtdb)

    def __repr__(self) -> str:
        fwd = f", forward={len(self.forward)} translations" if self.forward else ""
        return (
            f"NCBITranslator(version={self.version!r}, "
            f"{len(self.ncbi_name_to_gtdb)} NCBI mappings{fwd})"
        )
