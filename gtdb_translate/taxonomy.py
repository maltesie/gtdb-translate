"""Load and query GTDB taxonomy (e.g. bac120_taxonomy_r226.tsv)."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

from .utils import parse_gtdb_lineage

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")


@dataclass
class GTDBTaxonomy:
    """In-memory index of a GTDB taxonomy TSV file.

    The TSV is expected to have two columns (no header):
        genome_id <tab> d__…;p__…;c__…;o__…;f__…;g__…;s__…

    Attributes
    ----------
    species_to_lineage : dict[str, dict]
        Maps each species name to its full parsed lineage dict.
    species_list : list[str]
        Ordered list of unique species (insertion order).
    species_to_index : dict[str, int]
        Maps species name → index in ``species_list``.
    """

    species_to_lineage: Dict[str, dict] = field(default_factory=dict)
    species_list: List[str] = field(default_factory=list)
    species_to_index: Dict[str, int] = field(default_factory=dict)

    @classmethod
    def from_tsv(cls, path: Union[str, Path]) -> "GTDBTaxonomy":
        """Load a GTDB taxonomy TSV (e.g. ``bac120_taxonomy_r226.tsv``).

        Parameters
        ----------
        path : str or Path
            Path to the tab-separated taxonomy file.
        """
        taxonomy = cls()
        with open(path) as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if len(row) < 2:
                    continue
                lineage = parse_gtdb_lineage(row[1])
                species = lineage.get("species", "")
                if not species or species in taxonomy.species_to_lineage:
                    continue
                taxonomy.species_to_lineage[species] = lineage
                taxonomy.species_to_index[species] = len(taxonomy.species_list)
                taxonomy.species_list.append(species)
        return taxonomy

    def __contains__(self, species: str) -> bool:
        return species in self.species_to_lineage

    def __len__(self) -> int:
        return len(self.species_list)

    def get_lineage(self, species: str) -> Optional[dict]:
        """Return the full lineage dict for *species*, or ``None``."""
        return self.species_to_lineage.get(species)

    def get_rank(self, species: str, rank: str) -> Optional[str]:
        """Return a single rank value (e.g. ``"phylum"``) for *species*."""
        lineage = self.species_to_lineage.get(species)
        if lineage is None:
            return None
        return lineage.get(rank)

    def species_at_rank(self, rank: str, value: str) -> List[str]:
        """Return all species that share *value* at the given *rank*."""
        return [
            sp
            for sp, lin in self.species_to_lineage.items()
            if lin.get(rank) == value
        ]

    def unique_values(self, rank: str) -> Set[str]:
        """Return the set of unique values observed at *rank*."""
        return {lin[rank] for lin in self.species_to_lineage.values() if rank in lin}
