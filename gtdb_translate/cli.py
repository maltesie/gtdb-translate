"""Command-line interface for gtdb-translate.

Subcommands
-----------
build    Build a translation bundle from raw GTDB + NCBI files.
ncbi     Batch-translate NCBI names or tax IDs to GTDB.
silva    Batch-translate SILVA taxonomy to GTDB.
forward  Forward-translate old GTDB names to the current release.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

from .utils import RANK_PREFIXES


def _build(args: argparse.Namespace) -> None:
    from .build import SILVA_COLUMNS
    from .ncbi import NCBITranslator

    silva_columns = SILVA_COLUMNS
    if args.no_silva:
        silva_columns = ()
    elif args.silva_columns:
        silva_columns = tuple(args.silva_columns)

    translator = NCBITranslator.build(
        metadata_paths=args.metadata,
        names_dmp_path=args.names_dmp,
        changelog_path=args.changelog,
        version=args.version,
        silva_columns=silva_columns,
    )

    out = args.output or f"gtdb_translate_{args.version}.msgpack.zst"
    translator.save(out)
    print(f"Bundle saved to {out}")
    print(f"  NCBI mappings:        {len(translator.ncbi_name_to_gtdb):,}")
    print(f"  NCBI ID -> scientific: {len(translator.ncbi_id_to_scientific):,}")
    print(f"  GTDB lineage entries: {len(translator.gtdb_name_to_lineage):,}")
    n_silva = sum(len(d) for d in translator.silva_name_to_gtdb.values())
    print(f"  SILVA tokens:         {n_silva:,}")
    if translator.forward:
        print(f"  Forward translations: {translator.forward}")

    if not args.self_test:
        return

    from .selftest import format_report, run_self_test

    result = run_self_test(
        translator,
        metadata_paths=args.metadata,
        silva_columns=silva_columns,
        limit=args.self_test_limit,
        min_purity=args.min_purity,
    )
    print(format_report(result))
    if not result.passed:
        sys.exit(2)


def _load_translator(args: argparse.Namespace):
    """Load an NCBITranslator from bundle or auto-download."""
    from .ncbi import NCBITranslator

    if args.bundle:
        return NCBITranslator.load(args.bundle)
    return NCBITranslator.default(
        version=args.version, force_download=args.force
    )


def _read_table(path: str):
    """Read a CSV or TSV into a DataFrame."""
    import pandas as pd

    sep_in = "," if path.endswith(".csv") else "\t"
    return pd.read_csv(path, dtype=str, sep=sep_in)


def _ncbi(args: argparse.Namespace) -> None:
    import pandas as pd

    translator = _load_translator(args)
    df = _read_table(args.in_file)

    column_name = args.column_name
    if column_name is None:
        column_name = translator.detect_column(df, sep=args.multi_sep or "|")
        if column_name is None:
            print(
                "ERROR: Could not auto-detect the column to translate. "
                "Use --column_name to specify it.",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Auto-detected column: '{column_name}'")

    detected = translator.detect_format(
        df, column_name, multi_sep=args.multi_sep
    )
    from_taxids = args.from_taxids if args.from_taxids is not None else detected["from_taxids"]
    full_lineage = args.full_lineage if args.full_lineage is not None else detected["is_full_lineage"]
    multi_sep = args.multi_sep
    lineage_sep = args.lineage_sep if args.lineage_sep is not None else (
        detected["sep"] if full_lineage else ";"
    )

    auto_notes = []
    if args.from_taxids is None:
        auto_notes.append(f"from_taxids={from_taxids}")
    if args.full_lineage is None:
        auto_notes.append(_lineage_note(full_lineage, detected))
    if args.lineage_sep is None and full_lineage and not from_taxids:
        auto_notes.append(f"lineage_sep={lineage_sep!r}")
    if auto_notes:
        print(f"Auto-detected format: {', '.join(auto_notes)}")

    column = df[column_name]

    if from_taxids:
        # Always request lineages: the output column carries the full
        # lineage regardless of the input shape.
        translations, purity, support = translator.translate_ids(
            column, sep=multi_sep, full_lineage=True, with_support=True,
            min_purity=args.min_purity,
        )
    elif full_lineage:
        translations, purity, support = translator.translate(
            column, sep=lineage_sep, full_lineage=True,
            genus_fallback=args.genus_fallback, multi_sep=multi_sep,
            with_support=True, min_purity=args.min_purity,
            lineage_fallback=args.lineage_fallback,
        )
    else:
        translations, purity, support = translator.translate(
            column, sep=multi_sep, full_lineage=False,
            genus_fallback=args.genus_fallback, with_support=True,
            min_purity=args.min_purity,
        )

    # tax-ID translation already returns lineages when asked; names do not.
    returns_lineage = True if from_taxids else full_lineage
    _attach_outputs(
        df, args, translator, translations, returns_lineage,
        multi_sep, lineage_sep,
    )
    _attach_support(df, args, purity, support)

    _apply_empty_on_fail(df, args)

    df.to_csv(args.out_file, index=False)
    print(f"Output written to {args.out_file}")

NO_TRANS = "no_translation"

#: Rank prefixes tried, most specific first, when resolving a bare GTDB
#: name back to its full lineage.
_RANK_PREFIXES = ("s__", "g__", "f__", "o__", "c__", "p__", "d__")


def _lineage_note(full_lineage, detected) -> str:
    """Render the full_lineage decision with the evidence behind it.

    Reporting the nesting fraction distinguishes "nothing in this column
    nests" from "no separator was found at all" -- the former is what a
    bundle that does not cover the input taxonomy looks like.
    """
    fraction = detected.get("lineage_fraction")
    if fraction is None:
        return f"full_lineage={full_lineage}"
    return f"full_lineage={full_lineage} ({fraction:.0%} of cells nest)"


def _attach_support(df, args, purity, support) -> None:
    """Add purity and support columns unless the user opted out."""
    if not args.report_purity:
        return
    df[args.out_column_name + "_purity"] = purity
    df[args.out_column_name + "_support"] = support


def _map_parts(values, multi_sep, per_part):
    """Apply *per_part* to each part of every entry, preserving joins.

    An entry that is ``no_translation``, or whose parts all fail, stays
    ``no_translation`` so the failure marker never turns into an empty
    cell halfway through the pipeline.
    """
    out = []
    for entry in values:
        if not isinstance(entry, str) or not entry or entry == NO_TRANS:
            out.append(NO_TRANS)
            continue
        parts = entry.split(multi_sep) if multi_sep else [entry]
        mapped = [per_part(p.strip()) for p in parts]
        if all(m in ("", NO_TRANS) for m in mapped):
            out.append(NO_TRANS)
        elif multi_sep:
            out.append(multi_sep.join(m or NO_TRANS for m in mapped))
        else:
            out.append(mapped[0] or NO_TRANS)
    return out


def _lineage_of_name(translator, name):
    """Resolve a bare GTDB name to its full lineage, deepest rank first."""
    if not name or name == NO_TRANS:
        return NO_TRANS
    for prefix in _RANK_PREFIXES:
        lineage = translator.gtdb_name_to_lineage.get(prefix + name)
        if lineage:
            return lineage
    return NO_TRANS


def _lowest_of_lineage(lineage, lineage_sep):
    """Take the most specific rank of a lineage, without its prefix."""
    if not lineage or lineage == NO_TRANS:
        return NO_TRANS
    from .utils import strip_rank_prefix

    return strip_rank_prefix(lineage.split(lineage_sep)[-1])


def _apply_empty_on_fail(df, args) -> None:
    """Blank the failure marker in every column that can carry it."""
    if not args.empty_on_fail:
        return
    for column in (args.out_column_name, args.out_column_name + "_lowest"):
        if column in df.columns:
            df[column] = df[column].replace(NO_TRANS, "")


def _attach_outputs(
    df, args, translator, direct, is_lineage, multi_sep, lineage_sep
) -> None:
    """Write the lineage column, and the lowest-rank column if requested.

    *direct* is whatever the translator returned: already a lineage when
    the input was lineage-shaped, otherwise bare taxon names.  The output
    column always ends up carrying the full lineage either way, so the
    shape of the output no longer depends on the shape of the input.
    """
    if is_lineage:
        lineages = direct
        lowest = _map_parts(
            direct, multi_sep, lambda p: _lowest_of_lineage(p, lineage_sep)
        )
    else:
        lineages = _map_parts(
            direct, multi_sep, lambda p: _lineage_of_name(translator, p)
        )
        lowest = direct

    df[args.out_column_name] = lineages
    if args.output_lowest_rank:
        df[args.out_column_name + "_lowest"] = lowest


def _silva(args: argparse.Namespace) -> None:
    from .silva import SILVATranslator

    translator = _load_translator(args)
    silva = SILVATranslator.from_ncbi(translator)

    df = _read_table(args.in_file)

    column_name = args.column_name
    if column_name is None:
        column_name = translator.detect_column(df, sep=args.multi_sep or "|")
        if column_name is None:
            print(
                "ERROR: Could not auto-detect the column to translate. "
                "Use --column_name to specify it.",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Auto-detected column: '{column_name}'")

    detected = translator.detect_format(
        df, column_name, multi_sep=args.multi_sep
    )
    full_lineage = (
        args.full_lineage
        if args.full_lineage is not None
        else detected["is_full_lineage"]
    )
    lineage_sep = args.lineage_sep or detected["sep"] or ";"

    auto_notes = []
    if args.full_lineage is None:
        auto_notes.append(_lineage_note(full_lineage, detected))
    if args.lineage_sep is None:
        auto_notes.append(f"lineage_sep={lineage_sep!r}")
    if auto_notes:
        print(f"Auto-detected format: {', '.join(auto_notes)}")

    translations, purity, support = silva.translate(
        df[column_name],
        sep=lineage_sep,
        full_lineage=full_lineage,
        genus_fallback=args.genus_fallback,
        multi_sep=args.multi_sep,
        with_support=True,
        min_purity=args.min_purity,
        lineage_fallback=args.lineage_fallback,
    )
    _attach_outputs(
        df, args, translator, translations, full_lineage,
        args.multi_sep, lineage_sep,
    )
    _attach_support(df, args, purity, support)

    n_ok = sum(1 for t in translations if t != NO_TRANS)
    print(f"SILVA translation: {n_ok}/{len(translations)} entries translated")

    _apply_empty_on_fail(df, args)

    df.to_csv(args.out_file, index=False)
    print(f"Output written to {args.out_file}")


def _forward_resolver(translator):
    """Resolve a bare name the way the forward step would.

    A name already in the current GTDB maps to itself; an outdated one is
    forward-mapped first.  Using this rather than the NCBI dictionary is
    what lets lineage detection work on a column of outdated names.
    """
    fwd = translator.forward
    lineage_dict = translator.gtdb_name_to_lineage

    def resolve(name):
        for prefix in _RANK_PREFIXES:
            if prefix + name in lineage_dict:
                return prefix + name
        if fwd is None:
            return None
        if " " in name:
            mapped = fwd.translate(name)
            if mapped != name and "s__" + mapped in lineage_dict:
                return "s__" + mapped
        for prefix in _RANK_PREFIXES[1:]:
            rank = RANK_PREFIXES[prefix[0]]
            mapped = fwd.translate_rank(name, rank)
            if mapped != name and prefix + mapped in lineage_dict:
                return prefix + mapped
        return None

    return resolve


def _detect_forward_column(translator, df, sample_rows=100):
    """Pick the column most likely to hold GTDB names to forward-translate.

    Scores each string column by the share of sampled values containing at
    least one name the forward step could act on: a name already in the
    current GTDB, or one the forward translator knows how to update.
    Matching against both is what separates this from
    :meth:`NCBITranslator.detect_column` -- an outdated GTDB name is
    absent from the NCBI dictionary, so scoring on that alone would rank
    the most relevant column lowest.

    Values are tokenised on both the multi-value and the rank separator,
    so full-lineage columns score on their individual ranks.
    """
    import pandas as pd

    from .utils import split_rank_prefix

    fwd = translator.forward
    known = set(fwd.translation_dict)
    for rank_map in fwd.rank_translation_dicts.values():
        known.update(rank_map)

    lineage_dict = translator.gtdb_name_to_lineage

    def _recognised(token):
        _, bare = split_rank_prefix(token.strip())
        if not bare:
            return False
        if bare in known:
            return True
        return any(prefix + bare in lineage_dict for prefix in _RANK_PREFIXES)

    best_col, best_score = None, 0.0
    sample = df.head(sample_rows)
    for col in sample.columns:
        if not pd.api.types.is_string_dtype(sample[col]):
            continue
        hits = total = 0
        for value in sample[col].dropna():
            total += 1
            tokens = re.split(r"[|;,]", str(value))
            if any(_recognised(t) for t in tokens):
                hits += 1
        score = hits / total if total else 0.0
        if score > best_score:
            best_col, best_score = col, score

    if best_col is not None:
        print(
            f"Auto-detected column: '{best_col}' "
            f"({best_score:.0%} coverage)"
        )
    return best_col


def _forward(args: argparse.Namespace) -> None:
    import pandas as pd

    from .utils import split_rank_prefix

    translator = _load_translator(args)

    if translator.forward is None:
        print(
            "ERROR: The loaded bundle does not contain a forward "
            "translator. Rebuild the bundle with --changelog.",
            file=sys.stderr,
        )
        sys.exit(1)

    df = _read_table(args.in_file)
    column_name = args.column_name

    if column_name is None:
        column_name = _detect_forward_column(translator, df)
        if column_name is None:
            print(
                "ERROR: Could not auto-detect the column to translate. "
                "Use --column_name to specify it.",
                file=sys.stderr,
            )
            sys.exit(1)
    elif column_name not in df.columns:
        print(
            f"ERROR: Column '{column_name}' not found in {args.in_file}. "
            f"Available: {', '.join(df.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    detected = translator.detect_format(
        df, column_name, multi_sep=args.multi_sep,
        resolve=_forward_resolver(translator),
    )
    full_lineage = args.full_lineage if args.full_lineage is not None else detected["is_full_lineage"]
    lineage_sep = args.lineage_sep if args.lineage_sep is not None else (
        detected["sep"] if full_lineage else ";"
    )

    auto_notes = []
    if args.full_lineage is None:
        auto_notes.append(_lineage_note(full_lineage, detected))
    if args.lineage_sep is None:
        auto_notes.append(f"lineage_sep={lineage_sep!r}")
    if auto_notes:
        print(f"Auto-detected format: {', '.join(auto_notes)}")

    fwd = translator.forward
    lineage_dict = translator.gtdb_name_to_lineage

    RANK_PREFIXES = [
        ("s__", "species"),
        ("g__", "genus"),
        ("f__", "family"),
        ("o__", "order"),
        ("c__", "class"),
        ("p__", "phylum"),
        ("d__", "domain"),
    ]

    def forward_one(name: str) -> str:
        """Forward-translate a single GTDB name.

        1. Check if it's already in the current GTDB (try all rank prefixes)
        2. If not, forward-map it, then check again
        3. Return the current name or 'no_translation'
        """
        if not isinstance(name, str) or not name.strip():
            return NO_TRANS

        _, name = split_rank_prefix(name.strip())

        for prefix, rank in RANK_PREFIXES:
            if (prefix + name) in lineage_dict:
                return name

        if " " in name:
            mapped = fwd.translate(name)
            if mapped != name and ("s__" + mapped) in lineage_dict:
                return mapped
        for prefix, rank in RANK_PREFIXES[1:]:
            mapped = fwd.translate_rank(name, rank)
            if mapped != name and (prefix + mapped) in lineage_dict:
                return mapped

        return NO_TRANS

    multi_sep = args.multi_sep

    def forward_entry(entry: str) -> str:
        """Forward-translate an entry that may contain multiple names."""
        if not isinstance(entry, str) or not entry.strip():
            return NO_TRANS
        # multi_sep=None means one value per cell.
        parts = entry.split(multi_sep) if multi_sep else [entry]
        if full_lineage:
            translated = []
            for part in parts:
                part = part.strip()
                if not part:
                    translated.append(NO_TRANS)
                    continue
                result = fwd.translate_lineage(
                    part, lineage_dict, sep=lineage_sep,
                    lineage_fallback=args.lineage_fallback,
                )
                translated.append(result if result is not None else NO_TRANS)
        else:
            translated = [forward_one(p.strip()) for p in parts]
        if all(t == NO_TRANS for t in translated):
            return NO_TRANS
        return (multi_sep or "").join(translated)

    forwarded = [forward_entry(v) for v in df[column_name]]
    _attach_outputs(
        df, args, translator, forwarded, full_lineage,
        multi_sep, lineage_sep,
    )

    def _bare_entry(entry):
        if not isinstance(entry, str) or not entry.strip():
            return entry
        parts = entry.split(multi_sep) if multi_sep else [entry]
        return (multi_sep or "").join(
            split_rank_prefix(p.strip())[1] for p in parts
        )

    bare_input = df[column_name].map(_bare_entry)
    n_total = len(df)
    n_translated = sum(1 for v in forwarded if v != NO_TRANS)
    n_unchanged = sum(
        1 for v, b in zip(forwarded, bare_input) if v == b
    )
    n_mapped = n_translated - n_unchanged
    print(f"Forward mapping: {n_total} entries")
    print(f"  Already current: {n_unchanged}")
    print(f"  Forward-mapped:  {n_mapped}")
    print(f"  No translation:  {n_total - n_translated}")

    _apply_empty_on_fail(df, args)

    df.to_csv(args.out_file, index=False)
    print(f"Output written to {args.out_file}")


def _add_output_args(parser: argparse.ArgumentParser) -> None:
    """Add the shared output-shape flags to a subparser.

    The main output column always carries the full GTDB lineage.  The
    lowest-rank column is the opt-in extra.
    """
    parser.add_argument(
        "--lineage_fallback",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="When the lowest rank of a lineage does not resolve, work "
             "up through higher ranks and return the first that does "
             "(default: enabled). Pass --no-lineage_fallback to fail the "
             "entry instead, so that a rejected translation disappears "
             "rather than coming back as a shallower rank. A rank counts "
             "as failing whether it was missing or rejected.",
    )
    parser.add_argument(
        "--output_lowest_rank",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Add an <out_column>_lowest column holding just the most "
             "specific translated rank, without its rank prefix "
             "(default: disabled; the main column always holds the full "
             "lineage)",
    )
    parser.add_argument(
        "--empty_on_fail",
        action="store_true",
        help="Output empty string instead of 'no_translation' for failed lookups",
    )


def _add_report_args(parser: argparse.ArgumentParser) -> None:
    """Add the shared purity flags to a subparser."""
    parser.add_argument(
        "--min_purity",
        type=float,
        default=0.0,
        help="Reject a translation when the winning share of the vote "
             "falls below this fraction (default: 0.0, i.e. no filtering; "
             "0.5 is a reasonable starting point). Mappings with no votes "
             "behind them are never filtered.",
    )
    parser.add_argument(
        "--report_purity",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add <out_column>_purity and <out_column>_support columns "
             "giving the winning share of votes (a fraction in 0-1, "
             "1.0 = unanimous) and the total vote count behind each "
             "mapping (default: enabled). Values are empty where no "
             "votes lie behind the mapping.",
    )


def _add_bundle_args(parser: argparse.ArgumentParser) -> None:
    """Add --bundle, --version, --force to a subparser."""
    parser.add_argument(
        "--bundle",
        default=None,
        help="Path to a local .msgpack.zst bundle (downloads latest if omitted)",
    )
    parser.add_argument(
        "--version",
        default=None,
        help="GTDB version to download (default: latest release)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download the bundle even if cached locally",
    )


def main(argv: list[str] | None = None) -> None:
    """Main CLI entry point."""
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s: %(message)s"
    )

    parser = argparse.ArgumentParser(
        prog="gtdb-translate",
        description="Translate taxonomy names across GTDB releases and from NCBI.",
    )
    subparsers = parser.add_subparsers(dest="command")

    p_build = subparsers.add_parser(
        "build",
        help="Build a translation bundle from raw GTDB + NCBI files.",
    )
    p_build.add_argument(
        "--metadata",
        nargs="+",
        required=True,
        help="GTDB metadata TSV(s) (e.g. bac120_metadata_r226.tsv ar53_metadata_r226.tsv)",
    )
    p_build.add_argument(
        "--names_dmp",
        default=None,
        help="Path to NCBI names.dmp (optional; without it, synonym "
             "expansion and tax-ID translation are unavailable)",
    )
    p_build.add_argument(
        "--silva_columns",
        nargs="+",
        default=None,
        help="Metadata columns to pool SILVA votes from (default: "
             "ssu_silva_taxonomy lsu_silva_23s_taxonomy)",
    )
    p_build.add_argument(
        "--no_silva",
        action="store_true",
        help="Skip the SILVA dictionary entirely",
    )
    p_build.add_argument(
        "--self_test",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Score the finished bundle against the metadata it was "
             "built from (default: enabled)",
    )
    p_build.add_argument(
        "--min_purity",
        type=float,
        default=0.0,
        help="Purity threshold the self-test scores at; match it to the "
             "threshold you expect users to run with (default: 0.0)",
    )
    p_build.add_argument(
        "--self_test_limit",
        type=int,
        default=None,
        help="Only score this many genomes in the self-test",
    )
    p_build.add_argument(
        "--changelog",
        default=None,
        help="Path to gtdb-taxid-changelog.csv (optional, for forward translation)",
    )
    p_build.add_argument(
        "--version",
        default="r226",
        help="GTDB version label (default: r226)",
    )
    p_build.add_argument(
        "--output", "-o",
        default=None,
        help="Output bundle path (default: gtdb_translate_<version>.msgpack.zst)",
    )

    p_ncbi = subparsers.add_parser(
        "ncbi",
        help="Batch-translate NCBI names or tax IDs to GTDB.",
    )
    p_ncbi.add_argument("--in_file", required=True, help="Input CSV or TSV")
    p_ncbi.add_argument("--out_file", required=True, help="Output CSV")
    p_ncbi.add_argument(
        "--column_name",
        default=None,
        help="Column to translate (auto-detected if omitted)",
    )
    p_ncbi.add_argument(
        "--out_column_name",
        default="gtdb_translated",
        help="Name for the translated column (default: gtdb_translated)",
    )
    p_ncbi.add_argument(
        "--multi_sep",
        default=None,
        help="Separator for multiple independent entries per cell "
             "(default: |). With --full_lineage, a cell may hold several "
             "lineages joined by --multi_sep, each internally delimited "
             "by --lineage_sep.",
    )
    p_ncbi.add_argument(
        "--lineage_sep",
        default=None,
        help="Separator within a lineage string (default: auto-detected, "
             "usually ';')",
    )
    p_ncbi.add_argument(
        "--full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as full lineages (default: auto-detected)",
    )
    p_ncbi.add_argument(
        "--from_taxids",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as NCBI tax IDs (default: auto-detected)",
    )
    p_ncbi.add_argument(
        "--genus_fallback",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Fall back to genus-level when exact species match fails (default: disabled)",
    )
    _add_output_args(p_ncbi)
    _add_report_args(p_ncbi)
    _add_bundle_args(p_ncbi)

    p_fwd = subparsers.add_parser(
        "forward",
        help="Forward-translate old GTDB names to the current release.",
    )
    p_fwd.add_argument("--in_file", required=True, help="Input CSV or TSV")
    p_fwd.add_argument("--out_file", required=True, help="Output CSV")
    p_fwd.add_argument(
        "--column_name",
        default=None,
        help="Column containing old GTDB names to forward-translate "
             "(auto-detected if omitted)",
    )
    p_fwd.add_argument(
        "--out_column_name",
        default="gtdb_forwarded",
        help="Name for the output column (default: gtdb_forwarded)",
    )
    p_fwd.add_argument(
        "--multi_sep",
        default=None,
        help="Separator for multiple independent entries per cell "
             "(default: none -- each cell holds one value). Never "
             "auto-detected; set it explicitly if your cells hold "
             "several entries.",
    )
    p_fwd.add_argument(
        "--lineage_sep",
        default=None,
        help="Separator between ranks within a full-lineage entry "
             "(default: auto-detected, usually ';'; only used with "
             "--full_lineage)",
    )
    p_fwd.add_argument(
        "--full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as full GTDB lineages, e.g. "
             "'d__Bacteria;p__Firmicutes;...' or the bare "
             "'Bacteria;Firmicutes;...' (default: auto-detected)",
    )
    _add_output_args(p_fwd)
    # Forward translation walks a rename graph rather than a vote tally,
    # so there is neither purity to report nor a threshold to apply.
    p_fwd.set_defaults(report_purity=False, min_purity=0.0)
    _add_bundle_args(p_fwd)

    p_silva = subparsers.add_parser(
        "silva",
        help="Batch-translate SILVA taxonomy to GTDB.",
    )
    p_silva.add_argument("--in_file", required=True, help="Input CSV or TSV")
    p_silva.add_argument("--out_file", required=True, help="Output CSV")
    p_silva.add_argument(
        "--column_name",
        default=None,
        help="Column to translate (auto-detected if omitted)",
    )
    p_silva.add_argument(
        "--out_column_name",
        default="gtdb_translated",
        help="Name for the translated column (default: gtdb_translated)",
    )
    p_silva.add_argument(
        "--multi_sep",
        default=None,
        help="Separator for multiple independent lineages per cell "
             "(default: none -- each cell is one lineage)",
    )
    p_silva.add_argument(
        "--lineage_sep",
        default=None,
        help="Separator between ranks within a lineage "
             "(default: auto-detected, usually ';')",
    )
    p_silva.add_argument(
        "--full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as full SILVA lineages rather than single "
             "taxon names (default: auto-detected)",
    )
    p_silva.add_argument(
        "--genus_fallback",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Fall back to binomial/genus lookups when exact matching "
             "fails (default: disabled)",
    )
    _add_output_args(p_silva)
    _add_report_args(p_silva)
    _add_bundle_args(p_silva)

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "build":
        _build(args)
    elif args.command == "ncbi":
        _ncbi(args)
    elif args.command == "silva":
        _silva(args)
    elif args.command == "forward":
        _forward(args)


if __name__ == "__main__":
    main()
