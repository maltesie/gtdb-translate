"""Validate a freshly built bundle against the metadata it was built from.

This is a *test on train*: the same genomes that cast the votes are used
to score them, so the absolute numbers are optimistic and are not an
estimate of accuracy on new data.  That is deliberate -- the purpose is
to expose artefacts of the voting scheme, which show up clearly even on
training data:

* A rank whose accuracy is well below 100% has vote conflicts that the
  argmax is resolving badly, or an alignment bug.
* Ambiguous tokens should land in the low-purity bucket.  If accuracy is
  flat across purity buckets, purity is not measuring anything and is
  not worth reporting to users.
* Coverage below the fraction of rows carrying a SILVA classification
  means tokens are being dropped somewhere they shouldn't be.

The SILVA and NCBI paths are scored independently against each genome's
true GTDB lineage, and their agreement with each other is reported as a
third figure.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

from .build import GTDB_COLUMN, NCBI_COLUMN, SILVA_COLUMNS, _iter_rows
from .ncbi import NCBITranslator
from .silva import SILVATranslator
from .utils import RANK_ORDER

logger = logging.getLogger(__name__)

#: Ranks scored by the self-test.  Species is excluded: SILVA's seventh
#: field is a reference organism name rather than a species assignment,
#: so scoring it would measure the NCBI fallback, not the SILVA votes.
SCORED_RANKS: Tuple[str, ...] = RANK_ORDER[:6]

#: Purity buckets used to check that the statistic separates good
#: mappings from bad ones.  Bounds are half-open on the 0-1 purity
#: scale, with the top bucket extending past 1.0 so that exactly
#: unanimous mappings land in it.
_PURITY_BUCKETS: Tuple[Tuple[str, float, float], ...] = (
    ("unanimous (1.0)", 1.0, 1.01),
    ("high (0.8-1.0)", 0.8, 1.0),
    ("mixed (0.5-0.8)", 0.5, 0.8),
    ("low (<0.5)", 0.0, 0.5),
)


@dataclass
class PathResult:
    """Scores for one translation path (SILVA or NCBI)."""

    name: str
    n_input: int = 0
    n_translated: int = 0
    n_with_support: int = 0
    rank_correct: Counter = field(default_factory=Counter)
    rank_total: Counter = field(default_factory=Counter)
    purity_correct: Counter = field(default_factory=Counter)
    purity_total: Counter = field(default_factory=Counter)

    @property
    def coverage(self) -> float:
        return self.n_translated / self.n_input if self.n_input else 0.0

    @property
    def support_coverage(self) -> float:
        """Share of successful translations that report vote statistics."""
        return (
            self.n_with_support / self.n_translated
            if self.n_translated
            else 0.0
        )

    def rank_accuracy(self, rank: str) -> Optional[float]:
        total = self.rank_total[rank]
        return self.rank_correct[rank] / total if total else None


@dataclass
class SelfTestResult:
    """Outcome of :func:`run_self_test`."""

    silva: PathResult
    ncbi: PathResult
    n_rows: int = 0
    n_comparable: int = 0
    agreement: Counter = field(default_factory=Counter)
    agreement_total: Counter = field(default_factory=Counter)

    @property
    def passed(self) -> bool:
        """Whether the run looks free of gross voting artefacts.

        The thresholds are loose on purpose.  This checks for a broken
        build, not for good accuracy -- on training data anything much
        below these values indicates a bug rather than a hard dataset.
        """
        if self.silva.n_input == 0:
            return True
        if self.silva.coverage < 0.90:
            return False
        for rank in SCORED_RANKS:
            acc = self.silva.rank_accuracy(rank)
            if acc is not None and acc < 0.80:
                return False
        return True


def _split_bare(lineage: str) -> List[str]:
    """Split a prefixed GTDB lineage into bare names per rank."""
    return [part[3:] if len(part) > 3 else "" for part in lineage.split(";")]


def _score(
    result: PathResult,
    predicted: Optional[str],
    truth_bare: Sequence[str],
    purity: Optional[float],
) -> Optional[List[str]]:
    """Score one prediction against the truth, updating *result* in place."""
    result.n_input += 1
    if not predicted:
        return None

    result.n_translated += 1
    if purity is not None:
        result.n_with_support += 1
    pred_bare = _split_bare(predicted)

    correct_here = 0
    scored_here = 0
    for idx, rank in enumerate(SCORED_RANKS):
        if idx >= len(pred_bare) or idx >= len(truth_bare):
            break
        if not pred_bare[idx] or not truth_bare[idx]:
            continue
        scored_here += 1
        result.rank_total[rank] += 1
        if pred_bare[idx] == truth_bare[idx]:
            result.rank_correct[rank] += 1
            correct_here += 1

    if purity is not None and scored_here:
        for label, low, high in _PURITY_BUCKETS:
            if low <= purity < high:
                result.purity_total[label] += scored_here
                result.purity_correct[label] += correct_here
                break

    return pred_bare


def run_self_test(
    translator: NCBITranslator,
    metadata_paths: Sequence[Union[str, Path]],
    silva_columns: Sequence[str] = SILVA_COLUMNS,
    limit: Optional[int] = None,
) -> SelfTestResult:
    """Score a built bundle against the metadata it came from.

    Parameters
    ----------
    translator : NCBITranslator
        The freshly built translator.
    metadata_paths : sequence of str or Path
        The same metadata TSVs used to build it.  Read once, streaming.
    silva_columns : sequence of str
        SILVA columns to evaluate.  The first one carrying data for a
        given genome is used, matching the pooling done at build time.
    limit : int, optional
        Stop after this many rows.  Useful for a quick check on a large
        release.

    Returns
    -------
    SelfTestResult
    """
    silva = SILVATranslator.from_ncbi(translator)
    result = SelfTestResult(
        silva=PathResult("SILVA"), ncbi=PathResult("NCBI")
    )

    for path in metadata_paths:
        for row in _iter_rows(path):
            if limit is not None and result.n_rows >= limit:
                break

            truth = str(row.get(GTDB_COLUMN) or "")
            if len(truth.split(";")) < 2:
                continue
            result.n_rows += 1
            truth_bare = _split_bare(truth)

            silva_value = ""
            for column in silva_columns:
                value = str(row.get(column) or "").strip()
                if value and value.lower() != "none":
                    silva_value = value
                    break

            silva_pred = None
            silva_bare = None
            if silva_value:
                predicted, support = silva.translate_lineage(silva_value)
                silva_pred = predicted
                silva_bare = _score(
                    result.silva,
                    predicted,
                    truth_bare,
                    support[1] if support else None,
                )

            ncbi_value = str(row.get(NCBI_COLUMN) or "").strip()
            ncbi_bare = None
            if ncbi_value:
                predicted, support = translator._translate_single_lineage(
                    ncbi_value, sep=";", genus_fallback=False
                )
                ncbi_bare = _score(
                    result.ncbi,
                    predicted,
                    truth_bare,
                    support[1] if support else None,
                )

            if silva_bare and ncbi_bare:
                result.n_comparable += 1
                for idx, rank in enumerate(SCORED_RANKS):
                    if idx >= len(silva_bare) or idx >= len(ncbi_bare):
                        break
                    if not silva_bare[idx] or not ncbi_bare[idx]:
                        continue
                    result.agreement_total[rank] += 1
                    if silva_bare[idx] == ncbi_bare[idx]:
                        result.agreement[rank] += 1

    return result


def format_report(result: SelfTestResult) -> str:
    """Render a :class:`SelfTestResult` as a plain-text report."""
    lines: List[str] = []
    add = lines.append

    add("")
    add("Self-test (scored on the build data -- see module docstring)")
    add("-" * 62)
    add(f"  genomes read:            {result.n_rows:,}")
    for path in (result.silva, result.ncbi):
        add(
            f"  {path.name} input / translated: "
            f"{path.n_input:,} / {path.n_translated:,} "
            f"({path.coverage:.1%} coverage, "
            f"{path.support_coverage:.1%} with support stats)"
        )

    add("")
    add(f"  {'rank':<9}{'SILVA acc':>11}{'NCBI acc':>11}{'agreement':>12}")
    for rank in SCORED_RANKS:
        silva_acc = result.silva.rank_accuracy(rank)
        ncbi_acc = result.ncbi.rank_accuracy(rank)
        agree_total = result.agreement_total[rank]
        agree = (
            result.agreement[rank] / agree_total if agree_total else None
        )
        add(
            f"  {rank:<9}"
            f"{('-' if silva_acc is None else f'{silva_acc:.1%}'):>11}"
            f"{('-' if ncbi_acc is None else f'{ncbi_acc:.1%}'):>11}"
            f"{('-' if agree is None else f'{agree:.1%}'):>12}"
        )

    if result.silva.purity_total:
        add("")
        add("  SILVA accuracy by reported purity:")
        for label, _, _ in _PURITY_BUCKETS:
            total = result.silva.purity_total[label]
            if not total:
                continue
            acc = result.silva.purity_correct[label] / total
            add(f"    {label:<18} {acc:>7.1%}   ({total:,} rank calls)")

    add("")
    add(f"  result: {'PASS' if result.passed else 'CHECK'}")
    if not result.passed:
        add(
            "  Coverage or per-rank accuracy is below the sanity "
            "threshold; inspect the voting step before publishing."
        )
    add("")
    return "\n".join(lines)
