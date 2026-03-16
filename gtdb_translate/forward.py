"""Forward-translate GTDB species names that were renamed across releases.

Uses the changelog from https://github.com/shenwei356/gtdb-taxdump to build a
directed acyclic graph (DAG) of name transitions, then picks the best path to
the current name.
"""

from __future__ import annotations

import csv
import warnings
from pathlib import Path
from typing import Dict, Iterable, Optional, Set, Union

import networkx as nx

from .taxonomy import GTDBTaxonomy
from .utils import load_json, save_json


class ForwardTranslator:
    """Translate outdated GTDB species names to their current equivalents.

    The translator is built in two steps:

    1. :meth:`build` — parse the gtdb-taxdump changelog CSV and construct an
       internal translation DAG.
    2. :meth:`translate` / :meth:`translate_many` — look up one or more species
       names.

    The resulting translation dictionary can be persisted with :meth:`save` and
    later restored with :meth:`load`, so the (potentially large) changelog does
    not need to be re-parsed every time.

    Parameters
    ----------
    taxonomy : GTDBTaxonomy, optional
        If provided, only translations whose *target* species exists in the
        taxonomy are kept.  This filters out intermediate renames that are
        themselves now outdated.
    """

    def __init__(self, taxonomy: Optional[GTDBTaxonomy] = None) -> None:
        self.taxonomy = taxonomy
        self._dag: Optional[nx.DiGraph] = None
        self._trans_dict: Dict[str, str] = {}

    # ------------------------------------------------------------------
    # Building the DAG
    # ------------------------------------------------------------------
    def build(self, changelog_path: Union[str, Path]) -> "ForwardTranslator":
        """Parse the gtdb-taxdump changelog and build the translation dict.

        Parameters
        ----------
        changelog_path : str or Path
            Path to ``gtdb-taxid-changelog.csv`` (comma-separated, with
            columns: ``taxid, version, change, change-value, rank, lineage``
            and optionally ``name``).

        Returns
        -------
        self
            For method-chaining convenience.
        """
        dag = nx.DiGraph()
        current_taxid: Optional[str] = None
        chain: list[str] = []
        species_set: Set[str] = set()

        with open(changelog_path) as fh:
            reader = csv.DictReader(fh, delimiter=",")
            for row in reader:
                rank = row.get("rank", "")
                if rank != "no rank":
                    continue

                taxid = row["taxid"]
                if taxid != current_taxid:
                    self._commit_chain(dag, chain)
                    chain = []
                    species_set = set()
                    current_taxid = taxid

                lineage = row.get("lineage", "")
                parts = lineage.split(";")
                if len(parts) < 2:
                    continue
                species = parts[-2]  # second-to-last field is species

                if species in species_set:
                    continue
                species_set.add(species)

                version = row.get("version", "")
                node = f"{version}|{species}"
                change = row.get("change", "")

                if change in ("NEW", "CHANGE_LIN_TAX"):
                    chain.append(node)
                elif change == "DELETE":
                    chain = []

            # flush last chain
            self._commit_chain(dag, chain)

        self._dag = dag
        self._trans_dict = self._compute_translations(dag)
        return self

    # ------------------------------------------------------------------
    # Internal DAG helpers
    # ------------------------------------------------------------------
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
            # Only keep translations that look like valid binomial names
            if " " not in species or " " not in species_trans:
                continue
            # If a taxonomy is available, only keep translations that land in it
            if self.taxonomy is not None and species_trans not in self.taxonomy:
                continue
            raw[species] = species_trans

        # Resolve transitive chains  (A→B, B→C  ⇒  A→C)
        resolved: Dict[str, str] = {}
        for src, dst in raw.items():
            seen = {src}
            while dst in raw and dst not in seen:
                seen.add(dst)
                dst = raw[dst]
            resolved[src] = dst

        return resolved

    # ------------------------------------------------------------------
    # Translation interface
    # ------------------------------------------------------------------
    def translate(self, species: str) -> str:
        """Return the current name for *species*, or *species* itself."""
        return self._trans_dict.get(species, species)

    def translate_many(self, species_names: Iterable[str]) -> list[str]:
        """Translate a list of species names."""
        return [self.translate(s) for s in species_names]

    def has_translation(self, species: str) -> bool:
        """Return ``True`` if *species* has a known forward mapping."""
        return species in self._trans_dict

    @property
    def translation_dict(self) -> Dict[str, str]:
        """The full old-name → new-name mapping (read-only copy)."""
        return dict(self._trans_dict)

    def __len__(self) -> int:
        return len(self._trans_dict)

    def __repr__(self) -> str:
        return f"ForwardTranslator({len(self._trans_dict)} translations)"

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------
    def save(self, path: Union[str, Path]) -> None:
        """Save the translation dictionary to a JSON file.

        Only the final dictionary is saved — the DAG is *not* serialised.
        To rebuild from the changelog, call :meth:`build` again.
        """
        save_json(self._trans_dict, path)

    @classmethod
    def load(cls, path: Union[str, Path], taxonomy: Optional[GTDBTaxonomy] = None) -> "ForwardTranslator":
        """Load a previously saved translation dictionary.

        Parameters
        ----------
        path : str or Path
            JSON file written by :meth:`save`.
        taxonomy : GTDBTaxonomy, optional
            Optionally attach a taxonomy for downstream queries.
        """
        obj = cls(taxonomy=taxonomy)
        obj._trans_dict = load_json(path)
        return obj
