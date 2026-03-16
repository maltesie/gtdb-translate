"""Translate NCBI/SILVA taxonomy names and tax IDs to GTDB taxonomy.

This module provides :class:`NCBITranslator`, which wraps the three
translation dictionaries built from GTDB metadata and NCBI ``names.dmp``,
plus an optional :class:`~gtdb_translate.forward.ForwardTranslator` for
resolving renamed GTDB species across releases.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Union

from .forward import ForwardTranslator
from .utils import load_bundle, load_legacy_gzip_json, save_bundle

logger = logging.getLogger(__name__)


class NCBITranslator:
    """Translate NCBI taxonomy to GTDB.

    The translator holds three dictionaries (built from GTDB metadata TSVs
    and NCBI ``names.dmp``) and an optional :class:`ForwardTranslator`:

    * ``ncbi_name_to_gtdb`` — any NCBI / SILVA name → best GTDB taxon
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

    # ------------------------------------------------------------------
    # Building from raw data
    # ------------------------------------------------------------------
    def build(
        self,
        metadata_paths: Sequence[Union[str, Path]],
        names_dmp_path: Union[str, Path],
        changelog_path: Optional[Union[str, Path]] = None,
    ) -> "NCBITranslator":
        """Build all translation dicts from raw GTDB + NCBI files.

        Parameters
        ----------
        metadata_paths : sequence of str/Path
            Paths to GTDB metadata TSVs (e.g.
            ``["bac120_metadata_r226.tsv", "ar53_metadata_r226.tsv"]``).
        names_dmp_path : str or Path
            Path to NCBI ``names.dmp``.
        changelog_path : str or Path, optional
            Path to ``gtdb-taxid-changelog.csv``.  If provided, a
            :class:`ForwardTranslator` is built and included in the bundle.

        Returns
        -------
        self
        """
        import pandas as pd  # deferred — only needed for building

        # -- Step 1: NCBI/SILVA name → GTDB name (vote-based) -------------
        ncbi_votes: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        gtdb_name_to_lineage: Dict[str, str] = {}

        for fpath in metadata_paths:
            logger.info("Parsing metadata: %s", fpath)
            df = pd.read_csv(fpath, sep="\t", low_memory=False)
            for row in df.itertuples(index=False):
                ncbi_organism_name = getattr(row, "ncbi_organism_name", "")
                ncbi_lineage = str(getattr(row, "ncbi_taxonomy", "")).split(";")
                silva_lineage = str(getattr(row, "ssu_silva_taxonomy", "")).split(";")
                silva_23s_lineage = str(getattr(row, "lsu_silva_23s_taxonomy", "")).split(";")
                gtdb_lineage = str(getattr(row, "gtdb_taxonomy", "")).split(";")
                unfiltered_ncbi_last = str(
                    getattr(row, "ncbi_taxonomy_unfiltered", "")
                ).split(";")[-1]

                # organism name → species-level GTDB taxon
                ncbi_votes[ncbi_organism_name][gtdb_lineage[-1]] += 1

                # rank-by-rank alignment
                for ii, (ncbi_tax, gtdb_tax) in enumerate(
                    zip(ncbi_lineage, gtdb_lineage)
                ):
                    gtdb_name_to_lineage[gtdb_tax] = ";".join(
                        gtdb_lineage[: ii + 1]
                    )
                    # fall back to unfiltered species if filtered is empty
                    if (
                        ii == 6
                        and len(ncbi_tax) == 3
                        and len(unfiltered_ncbi_last) > 3
                        and unfiltered_ncbi_last.startswith("s__")
                    ):
                        ncbi_tax = unfiltered_ncbi_last

                    if not (len(ncbi_tax) == 3 or len(gtdb_tax) == 3):
                        ncbi_votes[ncbi_tax[3:]][gtdb_tax] += 1

                    # SILVA 16S
                    if len(silva_lineage) == 7 or (
                        ii + 1 < len(silva_lineage) < 7
                    ):
                        silva_tax = silva_lineage[ii].split("str.")[0].strip()
                        ncbi_votes[silva_tax][gtdb_tax] += 1
                    # SILVA 23S fallback
                    elif len(silva_23s_lineage) == 7 or (
                        ii + 1 < len(silva_23s_lineage) < 7
                    ):
                        silva_tax = (
                            silva_23s_lineage[ii].split("str.")[0].strip()
                        )
                        ncbi_votes[silva_tax][gtdb_tax] += 1

        # majority vote → single best GTDB name per NCBI name
        ncbi_name_to_gtdb: Dict[str, str] = {
            ncbi: max(counter.items(), key=lambda x: x[1])[0]
            for ncbi, counter in ncbi_votes.items()
        }

        # ensure every GTDB name (without prefix) maps to itself
        for fpath in metadata_paths:
            df = pd.read_csv(fpath, sep="\t", low_memory=False)
            for row in df.itertuples(index=False):
                gtdb_lineage = str(getattr(row, "gtdb_taxonomy", "")).split(";")
                for gtdb_tax in gtdb_lineage:
                    bare = gtdb_tax[3:]
                    if bare not in ncbi_name_to_gtdb:
                        ncbi_name_to_gtdb[bare] = gtdb_tax

        # -- Step 2: NCBI names.dmp — synonyms & scientific names ----------
        logger.info("Parsing names.dmp: %s", names_dmp_path)
        ncbi_name_to_id: Dict[str, int] = {}
        ncbi_id_to_scientific: Dict[int, str] = {}
        with open(names_dmp_path) as fh:
            for line in fh:
                parts = [x.strip() for x in line.split("|")]
                if len(parts) < 4:
                    continue
                taxid = int(parts[0])
                name = parts[1].split("(")[0].strip()
                ncbi_name_to_id[name] = taxid
                if parts[3] == "scientific name":
                    ncbi_id_to_scientific[taxid] = name

        # Resolve synonyms through scientific names
        ncbi_rep_to_gtdb = {}
        for ncbi_name, gtdb_name in ncbi_name_to_gtdb.items():
            if ncbi_name in ncbi_name_to_id:
                rep = ncbi_id_to_scientific.get(ncbi_name_to_id[ncbi_name], ncbi_name)
            else:
                rep = ncbi_name
            ncbi_rep_to_gtdb[rep] = gtdb_name

        expanded: Dict[str, str] = {}
        for ncbi_name, taxid in ncbi_name_to_id.items():
            rep = ncbi_id_to_scientific.get(taxid, ncbi_name)
            expanded[ncbi_name] = ncbi_rep_to_gtdb.get(rep, "none")
        expanded.update(ncbi_rep_to_gtdb)

        # Store with string keys for serialisation
        id_to_sci_str = {str(k): v for k, v in ncbi_id_to_scientific.items()}

        self.ncbi_name_to_gtdb = expanded
        self.ncbi_id_to_scientific = id_to_sci_str
        self.gtdb_name_to_lineage = gtdb_name_to_lineage

        # -- Step 3 (optional): Forward translator -------------------------
        if changelog_path is not None:
            logger.info("Building forward translator from: %s", changelog_path)
            self.forward = ForwardTranslator()
            self.forward.build(changelog_path)
        else:
            self.forward = None

        return self

    # ------------------------------------------------------------------
    # Translation — names
    # ------------------------------------------------------------------
    @staticmethod
    def sanitize_lineages(
        lineages: Iterable[str],
        lineage_sep: str = ";",
        check_merge_species: bool = False,
        replace_symbols: Optional[Dict[str, str]] = None,
        check_remove_sk: bool = False,
    ) -> List[str]:
        """Clean up raw lineage strings before translation.

        Parameters
        ----------
        lineages : iterable of str
            Raw lineage strings.
        lineage_sep : str
            Delimiter between ranks (default ``";"``)
        check_merge_species : bool
            If the species name does not start with the genus, prepend genus.
        replace_symbols : dict, optional
            Characters to replace (e.g. ``{"_": " "}``).
        check_remove_sk : bool
            Remove the ``sk__`` (superkingdom) rank if present.
        """
        sanitized = []
        for lineage in lineages:
            taxa = [t.strip() for t in lineage.split(lineage_sep)]
            check_species = taxa[-1].startswith("s__")
            if check_remove_sk and taxa[0].startswith("sk__"):
                taxa = [taxa[0]] + taxa[2:]
            taxa = [
                "__".join(t.split("__")[1:]) if "__" in t else t for t in taxa
            ]
            if replace_symbols is not None:
                for old, new in replace_symbols.items():
                    taxa = [t.replace(old, new) for t in taxa]
            if (
                check_merge_species
                and check_species
                and not taxa[-1].startswith(taxa[-2])
            ):
                taxa[-1] = f"{taxa[-2]} {taxa[-1]}"
            sanitized.append(lineage_sep.join(taxa))
        return sanitized

    def translate(
        self,
        entries: Sequence[str],
        sep: str = "|",
        full_lineage: bool = False,
    ) -> List[str]:
        """Translate NCBI/SILVA names to GTDB.

        Parameters
        ----------
        entries : sequence of str
            Each entry is one or more taxon names joined by *sep*,
            or a full lineage string if *full_lineage* is ``True``.
        sep : str
            Separator between multiple names within a single entry.
        full_lineage : bool
            If ``True``, treat each entry as a complete lineage and
            return the best matching GTDB lineage.

        Returns
        -------
        list of str
            Translated entries.  ``"no_translation"`` when no mapping
            was found.
        """
        translations = ["no_translation"] * len(entries)
        for i, entry in enumerate(entries):
            if not isinstance(entry, str):
                continue
            taxa = entry.split(sep)[::-1]
            if full_lineage and not (
                taxa[-1].endswith("Bacteria") or taxa[-1].endswith("Archaea")
            ):
                continue
            for ii, tax in enumerate(taxa):
                gtdb_name = self._lookup_name(tax)
                if full_lineage:
                    best = "no_gtdb"
                    if gtdb_name in self.gtdb_name_to_lineage:
                        best = self.gtdb_name_to_lineage[gtdb_name]
                    if best != "no_gtdb":
                        if best.count(";") >= len(taxa) - ii:
                            best = ";".join(
                                best.split(";")[: len(taxa) - ii]
                            )
                        translations[i] = best
                        break
                else:
                    translation = (
                        gtdb_name[3:]
                        if (gtdb_name is not None and gtdb_name != "none")
                        else ""
                    )
                    if translations[i] == "no_translation":
                        translations[i] = translation
                    else:
                        translations[i] = translation + sep + translations[i]
        return translations

    def translate_ids(
        self,
        entries: Sequence[str],
        sep: str = "|",
        full_lineage: bool = False,
    ) -> List[str]:
        """Translate NCBI tax IDs to GTDB names.

        Parameters
        ----------
        entries : sequence of str
            Each entry is one or more NCBI tax IDs joined by *sep*.
        sep : str
            Separator between multiple IDs within a single entry.
        full_lineage : bool
            If ``True``, return full GTDB lineage strings.

        Returns
        -------
        list of str
        """
        translations = ["no_translation"] * len(entries)
        for i, taxs in enumerate(entries):
            if not isinstance(taxs, str):
                continue
            parts = []
            for taxid in taxs.split(sep):
                sci = self.ncbi_id_to_scientific.get(taxid)
                if sci is None:
                    continue
                gtdb_name = self.ncbi_name_to_gtdb.get(sci, "none")
                if full_lineage and gtdb_name in self.gtdb_name_to_lineage:
                    gtdb_name = self.gtdb_name_to_lineage[gtdb_name]
                if gtdb_name == "none":
                    gtdb_name = ""
                parts.append(gtdb_name)
            translations[i] = sep.join(parts)
        return translations

    def _lookup_name(self, name: str) -> Optional[str]:
        """Look up a single NCBI name, trying bracket-removal as fallback."""
        if name in self.ncbi_name_to_gtdb:
            return self.ncbi_name_to_gtdb[name]
        if "[" in name:
            cleaned = name.replace("[", "").replace("]", "")
            if cleaned in self.ncbi_name_to_gtdb:
                return self.ncbi_name_to_gtdb[cleaned]
        return None

    # ------------------------------------------------------------------
    # Persistence — bundle format
    # ------------------------------------------------------------------
    def save(self, path: Union[str, Path]) -> None:
        """Save all translation dicts as a single ``.msgpack.zst`` bundle.

        The bundle contains:

        * ``ncbi_name_to_gtdb``
        * ``ncbi_id_to_scientific``
        * ``gtdb_name_to_lineage``
        * ``forward_trans_dict`` (if a forward translator is attached)
        * ``version``
        """
        data = {
            "version": self.version,
            "ncbi_name_to_gtdb": self.ncbi_name_to_gtdb,
            "ncbi_id_to_scientific": self.ncbi_id_to_scientific,
            "gtdb_name_to_lineage": self.gtdb_name_to_lineage,
            "forward_trans_dict": (
                self.forward.translation_dict if self.forward else {}
            ),
        }
        save_bundle(data, path)

    @classmethod
    def load(cls, path: Union[str, Path]) -> "NCBITranslator":
        """Load a translator from a ``.msgpack.zst`` bundle.

        Also supports the legacy ``translation_dicts_rXXX.json.gz``
        format (forward translator will not be available in that case).
        """
        path = Path(path)

        # Legacy format detection
        if path.suffixes[-2:] == [".json", ".gz"] or path.suffix == ".gz":
            dicts = load_legacy_gzip_json(path)
            obj = cls()
            obj.ncbi_name_to_gtdb = dicts[0]
            # Legacy stored int keys; normalise to str
            obj.ncbi_id_to_scientific = {
                str(k): v for k, v in dicts[1].items()
            }
            obj.gtdb_name_to_lineage = dicts[2]
            return obj

        data = load_bundle(path)
        obj = cls(version=data.get("version", "unknown"))
        obj.ncbi_name_to_gtdb = data["ncbi_name_to_gtdb"]
        obj.ncbi_id_to_scientific = data["ncbi_id_to_scientific"]
        obj.gtdb_name_to_lineage = data["gtdb_name_to_lineage"]
        fwd_dict = data.get("forward_trans_dict", {})
        if fwd_dict:
            obj.forward = ForwardTranslator()
            obj.forward._trans_dict = fwd_dict
        return obj

    # ------------------------------------------------------------------
    # Convenience — auto-download
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # Column auto-detection (scaffold)
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return len(self.ncbi_name_to_gtdb)

    def __repr__(self) -> str:
        fwd = f", forward={len(self.forward)} translations" if self.forward else ""
        return (
            f"NCBITranslator(version={self.version!r}, "
            f"{len(self.ncbi_name_to_gtdb)} NCBI mappings{fwd})"
        )
