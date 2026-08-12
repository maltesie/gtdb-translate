"""Forward-translate GTDB names that were renamed across releases.

Uses the changelog from https://github.com/shenwei356/gtdb-taxdump to build a
directed acyclic graph (DAG) of species name transitions, and collects
rank-level votes for higher taxonomy renames (phylum, class, order, family,
genus) from the same genome-level data.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Optional, Set, Union

import networkx as nx

from .taxonomy import GTDBTaxonomy
from .utils import (
    RANK_ORDER,
    RANK_PREFIXES,
    RANK_TO_PREFIX,
    load_json,
    save_json,
    split_rank_prefix,
)

logger = logging.getLogger(__name__)

RANK_POSITIONS = {
    0: "domain",
    1: "phylum",
    2: "class",
    3: "order",
    4: "family",
    5: "genus",
}


class ForwardTranslator:
    """Translate outdated GTDB names to their current equivalents.

    Build with :meth:`build` (parses the gtdb-taxdump changelog into a
    species-level DAG plus per-rank translation dicts), then look up
    names with :meth:`translate` / :meth:`translate_many` (species) or
    :meth:`translate_rank` (higher ranks). Persist with :meth:`save` /
    :meth:`load`.

    Parameters
    ----------
    taxonomy : GTDBTaxonomy, optional
        If provided, only translations whose *target* name exists in the
        taxonomy are kept.
    """

    def __init__(self, taxonomy: Optional[GTDBTaxonomy] = None) -> None:
        self.taxonomy = taxonomy
        self._dag: Optional[nx.DiGraph] = None
        self._trans_dict: Dict[str, str] = {}
        self._rank_trans_dicts: Dict[str, Dict[str, str]] = {}

    def build(self, changelog_path: Union[str, Path]) -> "ForwardTranslator":
        """Parse the gtdb-taxdump changelog and build translation dicts.

        Builds both:

        * A species-level DAG (from genome entries where the species name
          changed between versions).
        * Per-rank translation dicts (from *all* genome lineage changes,
          including those where only higher ranks changed).

        Parameters
        ----------
        changelog_path : str or Path
            Path to ``gtdb-taxid-changelog.csv``.

        Returns
        -------
        self
            For method-chaining convenience.
        """
        dag = nx.DiGraph()
        current_taxid: Optional[str] = None
        chain: list[str] = []
        species_set: Set[str] = set()
        prev_lineage: list[str] = []

        rank_votes: Dict[str, Dict[str, Dict[str, int]]] = {
            name: defaultdict(lambda: defaultdict(int))
            for name in RANK_POSITIONS.values()
        }

        with open(changelog_path) as fh:
            reader = csv.DictReader(fh, delimiter=",")
            for row in reader:
                rank = row.get("rank", "")
                if rank != "no rank":
                    continue

                taxid = row["taxid"]
                change = row.get("change", "")
                lineage = row.get("lineage", "")
                parts = lineage.split(";")

                if taxid != current_taxid:
                    self._commit_chain(dag, chain)
                    chain = []
                    species_set = set()
                    prev_lineage = []
                    current_taxid = taxid

                if len(parts) < 2:
                    continue

                if (
                    change in ("NEW", "CHANGE_LIN_TAX")
                    and len(prev_lineage) >= 7
                    and len(parts) >= 7
                ):
                    for pos, rank_name in RANK_POSITIONS.items():
                        if pos < len(prev_lineage) and pos < len(parts):
                            old = prev_lineage[pos]
                            new = parts[pos]
                            if old != new:
                                rank_votes[rank_name][old][new] += 1

                if change == "DELETE":
                    prev_lineage = []
                else:
                    prev_lineage = parts

                species = parts[-2]

                if species in species_set:
                    continue
                species_set.add(species)

                version = row.get("version", "")
                node = f"{version}|{species}"

                if change in ("NEW", "CHANGE_LIN_TAX"):
                    chain.append(node)
                elif change == "DELETE":
                    chain = []

            self._commit_chain(dag, chain)

        self._dag = dag
        self._trans_dict = self._compute_translations(dag)
        self._rank_trans_dicts = self._compute_rank_translations(rank_votes)

        n_rank = sum(len(d) for d in self._rank_trans_dicts.values())
        logger.info(
            "Built %d species translations + %d rank translations",
            len(self._trans_dict),
            n_rank,
        )
        return self

    @staticmethod
    def _commit_chain(dag: nx.DiGraph, chain: list[str]) -> None:
        """Add a chain of version|species nodes to the DAG."""
        if len(chain) > 1:
            for i in range(len(chain) - 1):
                for n in (chain[i], chain[i + 1]):
                    if not dag.has_node(n):
                        dag.add_node(n, species=n.split("|")[1], count=1)
                    else:
                        dag.nodes[n]["count"] += 1
                if not dag.has_edge(chain[i], chain[i + 1]):
                    dag.add_edge(chain[i], chain[i + 1], count=1)
                else:
                    dag[chain[i]][chain[i + 1]]["count"] += 1
        elif len(chain) == 1:
            n = chain[0]
            if not dag.has_node(n):
                dag.add_node(n, species=n.split("|")[1], count=1)
            else:
                dag.nodes[n]["count"] += 1

    @staticmethod
    def _best_translation(dag: nx.DiGraph, node: str) -> str:
        """Walk the DAG greedily to find the most likely current name."""
        successors = dag[node]
        if len(successors) == 0:
            return node
        total_out = sum(successors[x]["count"] for x in successors)
        if total_out / dag.nodes[node]["count"] < 0.5:
            return node
        best_next = max(successors, key=lambda x: successors[x]["count"])
        return ForwardTranslator._best_translation(dag, best_next)

    def _compute_translations(self, dag: nx.DiGraph) -> Dict[str, str]:
        """Derive a species→species translation dict from the finished DAG."""
        raw: Dict[str, str] = {}
        for node in dag.nodes():
            translation = self._best_translation(dag, node)
            species = node.split("|")[1]
            species_trans = translation.split("|")[1]
            if species_trans == species:
                continue
            if " " not in species or " " not in species_trans:
                continue
            if self.taxonomy is not None and species_trans not in self.taxonomy:
                continue
            raw[species] = species_trans

        resolved: Dict[str, str] = {}
        for src, dst in raw.items():
            seen = {src}
            while dst in raw and dst not in seen:
                seen.add(dst)
                dst = raw[dst]
            if src != dst:
                resolved[src] = dst

        return resolved

    def _compute_rank_translations(
        self,
        rank_votes: Dict[str, Dict[str, Dict[str, int]]],
    ) -> Dict[str, Dict[str, str]]:
        """Compute per-rank translation dicts from collected votes.

        For each rank, takes the majority-vote target for each old name,
        optionally filters against the taxonomy, and resolves transitive
        chains.
        """
        result: Dict[str, Dict[str, str]] = {}
        taxonomy_values: Dict[str, Set[str]] = {}

        if self.taxonomy is not None:
            for rank_name in RANK_POSITIONS.values():
                taxonomy_values[rank_name] = self.taxonomy.unique_values(
                    rank_name
                )

        for rank_name, votes in rank_votes.items():
            if not votes:
                continue

            raw: Dict[str, str] = {}
            for old_name, targets in votes.items():
                best = max(targets.items(), key=lambda x: x[1])
                new_name = best[0]
                if new_name == old_name:
                    continue
                if (
                    self.taxonomy is not None
                    and rank_name in taxonomy_values
                    and new_name not in taxonomy_values[rank_name]
                ):
                    continue
                raw[old_name] = new_name

            resolved: Dict[str, str] = {}
            for src, dst in raw.items():
                seen = {src}
                while dst in raw and dst not in seen:
                    seen.add(dst)
                    dst = raw[dst]
                if src != dst:
                    resolved[src] = dst

            if resolved:
                result[rank_name] = resolved
                logger.info(
                    "  %s: %d translations", rank_name, len(resolved)
                )

        return result

    def translate(self, species: str) -> str:
        """Return the current name for *species*, or *species* itself.

        Accepts either a bare species name (``"Escherichia coli"``) or a
        rank-prefixed one (``"s__Escherichia coli"``) — the lookup is
        always done on the bare name, and if the input was prefixed the
        same prefix is restored on the output.
        """
        prefix, bare = split_rank_prefix(species)
        mapped = self._trans_dict.get(bare)
        if mapped is None:
            return species
        return f"{prefix}__{mapped}" if prefix else mapped

    def translate_many(self, species_names: Iterable[str]) -> list[str]:
        """Translate a list of species names."""
        return [self.translate(s) for s in species_names]

    def translate_rank(self, name: str, rank: str) -> str:
        """Forward-translate a name at a specific taxonomic rank.

        Use this when the input is a bare genus, family, order, etc. that
        may not exist in the current GTDB.  Once you have the current name,
        look up its full lineage from the taxonomy or lineage dict.

        Parameters
        ----------
        name : str
            The (potentially outdated) taxon name.
        rank : str
            One of ``"phylum"``, ``"class"``, ``"order"``, ``"family"``,
            ``"genus"``.

        Returns
        -------
        str
            The current name, or *name* itself if no mapping exists.
        """
        prefix, bare = split_rank_prefix(name)
        rank_dict = self._rank_trans_dicts.get(rank, {})
        mapped = rank_dict.get(bare)
        if mapped is None:
            return name
        return f"{prefix}__{mapped}" if prefix else mapped

    def translate_lineage(
        self,
        lineage: str,
        gtdb_name_to_lineage: Dict[str, str],
        sep: str = ";",
        lineage_fallback: bool = True,
    ) -> Optional[str]:
        """Forward-translate a GTDB lineage to its current form.

        Works bottom-up from the lowest rank: if a taxon exists in the
        current GTDB lineage dict, its stored lineage is returned
        directly (capturing any higher-rank renames). If the lowest
        rank is not found, it is forward-mapped first, then looked up
        again. Falls back to progressively higher ranks.

        Accepts lineages whose ranks are either explicitly prefixed
        (``"d__Bacteria;p__Firmicutes;...;s__Escherichia coli"``) or bare
        (``"Bacteria;Firmicutes;...;Escherichia coli"``). Bare lineages
        are assumed to be given in standard domain-to-species order,
        anchored at domain (so a partial bare lineage is assumed to be
        missing its *lowest* ranks, e.g. a 5-token bare lineage is
        domain through family). A lineage may also mix prefixed and bare
        tokens.

        Every returned lineage is guaranteed to come from
        *gtdb_name_to_lineage* and therefore be in the current GTDB.

        Parameters
        ----------
        lineage : str
            A GTDB lineage, prefixed or bare (see above).
        gtdb_name_to_lineage : dict
            Maps prefixed GTDB names (e.g. ``"s__Escherichia coli"``)
            to their full current lineage string.
        sep : str
            Separator between ranks (default ``";"``).
        lineage_fallback : bool
            When ``True`` (the default), work up through higher ranks if
            the lowest one cannot be resolved.  When ``False``, only the
            lowest rank is tried.

        Returns
        -------
        str or None
            The current lineage from the dict, or ``None`` if nothing
            could be resolved.
        """
        raw_parts = [p.strip() for p in lineage.split(sep) if p.strip()]
        n = len(raw_parts)
        positional_ranks = RANK_ORDER[:n] if n <= len(RANK_ORDER) else [None] * n

        resolved = []
        for idx, part in enumerate(raw_parts):
            letter, bare = split_rank_prefix(part)
            if letter is not None:
                rank = RANK_PREFIXES.get(letter)
            else:
                rank = positional_ranks[idx]
                letter = RANK_TO_PREFIX.get(rank)
            if not bare or rank is None or letter is None:
                continue
            resolved.append((letter, rank, bare))

        for letter, rank, name in reversed(resolved):
            prefixed = f"{letter}__{name}"

            if prefixed in gtdb_name_to_lineage:
                return gtdb_name_to_lineage[prefixed]

            if rank == "species":
                mapped = self.translate(name)
            else:
                mapped = self.translate_rank(name, rank)

            if mapped != name:
                mapped_prefixed = f"{letter}__{mapped}"
                if mapped_prefixed in gtdb_name_to_lineage:
                    return gtdb_name_to_lineage[mapped_prefixed]

            if not lineage_fallback:
                return None

        return None

    def has_translation(self, species: str) -> bool:
        """Return ``True`` if *species* has a known forward mapping."""
        return species in self._trans_dict

    @property
    def translation_dict(self) -> Dict[str, str]:
        """The full old-name → new-name species mapping (read-only copy)."""
        return dict(self._trans_dict)

    @property
    def rank_translation_dicts(self) -> Dict[str, Dict[str, str]]:
        """Per-rank translation dicts (read-only copy)."""
        return {k: dict(v) for k, v in self._rank_trans_dicts.items()}

    def __len__(self) -> int:
        return len(self._trans_dict)

    def __repr__(self) -> str:
        rank_counts = ", ".join(
            f"{r}={len(d)}" for r, d in self._rank_trans_dicts.items()
        )
        rank_info = f", rank translations: {rank_counts}" if rank_counts else ""
        return (
            f"ForwardTranslator({len(self._trans_dict)} species translations"
            f"{rank_info})"
        )

    def save(self, path: Union[str, Path]) -> None:
        """Save translation dictionaries to a JSON file.

        Saves both the species-level dict and the per-rank dicts.
        The DAG is *not* serialised.
        """
        data = {
            "species": self._trans_dict,
            "ranks": self._rank_trans_dicts,
        }
        save_json(data, path)

    @classmethod
    def load(
        cls,
        path: Union[str, Path],
        taxonomy: Optional[GTDBTaxonomy] = None,
    ) -> "ForwardTranslator":
        """Load previously saved translation dictionaries.

        Supports both the new format (with rank dicts) and the old format
        (species-only flat dict).

        Parameters
        ----------
        path : str or Path
            JSON file written by :meth:`save`.
        taxonomy : GTDBTaxonomy, optional
            Optionally attach a taxonomy for downstream queries.
        """
        obj = cls(taxonomy=taxonomy)
        raw = load_json(path)

        if isinstance(raw, dict) and "species" in raw:
            obj._trans_dict = raw["species"]
            obj._rank_trans_dicts = raw.get("ranks", {})
        else:
            obj._trans_dict = raw
            obj._rank_trans_dicts = {}

        return obj
