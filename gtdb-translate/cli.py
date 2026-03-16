"""Command-line interface for gtdb-translate.

Subcommands
-----------
build     Build a translation bundle from raw GTDB + NCBI files.
translate Batch-translate a column in a CSV/TSV file.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def _build(args: argparse.Namespace) -> None:
    from .ncbi import NCBITranslator

    translator = NCBITranslator(version=args.version)
    translator.build(
        metadata_paths=args.metadata,
        names_dmp_path=args.names_dmp,
        changelog_path=args.changelog,
    )
    out = args.output or f"gtdb_translate_{args.version}.msgpack.zst"
    translator.save(out)
    print(f"Bundle saved to {out}")
    print(f"  NCBI mappings:        {len(translator.ncbi_name_to_gtdb):,}")
    print(f"  NCBI ID → scientific: {len(translator.ncbi_id_to_scientific):,}")
    print(f"  GTDB lineage entries: {len(translator.gtdb_name_to_lineage):,}")
    if translator.forward:
        print(f"  Forward translations: {len(translator.forward):,}")


def _translate(args: argparse.Namespace) -> None:
    import pandas as pd

    from .ncbi import NCBITranslator

    # Load translator
    if args.bundle:
        translator = NCBITranslator.load(args.bundle)
    else:
        translator = NCBITranslator.default(version=args.version)

    # Load input table
    sep_in = "," if args.in_file.endswith(".csv") else "\t"
    df = pd.read_csv(args.in_file, dtype=str, sep=sep_in)

    # Detect or validate column
    column_name = args.column_name
    if column_name is None:
        column_name = translator.detect_column(df, sep=args.sep)
        if column_name is None:
            print(
                "ERROR: Could not auto-detect the column to translate. "
                "Use --column_name to specify it.",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Auto-detected column: '{column_name}'")

    column = df[column_name]

    # Sanitize if full-lineage mode
    if args.full_lineage and not args.from_taxids:
        column = NCBITranslator.sanitize_lineages(
            column,
            lineage_sep=args.lineage_sep,
            check_merge_species=True,
            replace_symbols={"_": " "},
            check_remove_sk=True,
        )

    # Translate
    if args.from_taxids:
        df[args.out_column_name] = translator.translate_ids(
            column, sep=args.lineage_sep, full_lineage=args.full_lineage
        )
    else:
        df[args.out_column_name] = translator.translate(
            column, sep=args.lineage_sep, full_lineage=args.full_lineage
        )

    # Optional extras
    if args.full_lineage:
        df[column_name + "_lowest"] = [
            t.split(args.lineage_sep)[-1] for t in df[args.out_column_name]
        ]

    if args.output_full_lineage:
        df[column_name + "_lineage"] = [
            translator.gtdb_name_to_lineage.get("s__" + t, "") if t else ""
            for t in df[args.out_column_name]
        ]

    df.to_csv(args.out_file, index=False)
    print(f"Output written to {args.out_file}")


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

    # -- build ---------------------------------------------------------
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
        required=True,
        help="Path to NCBI names.dmp",
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

    # -- translate -----------------------------------------------------
    p_trans = subparsers.add_parser(
        "translate",
        help="Batch-translate a column in a CSV/TSV file.",
    )
    p_trans.add_argument("--in_file", required=True, help="Input CSV or TSV")
    p_trans.add_argument("--out_file", required=True, help="Output CSV")
    p_trans.add_argument(
        "--column_name",
        default=None,
        help="Column to translate (auto-detected if omitted)",
    )
    p_trans.add_argument(
        "--out_column_name",
        default="gtdb_translated",
        help="Name for the translated column (default: gtdb_translated)",
    )
    p_trans.add_argument(
        "--sep",
        default="|",
        help="Separator for multiple names per cell (default: |)",
    )
    p_trans.add_argument(
        "--lineage_sep",
        default=";",
        help="Separator within a lineage string (default: ;)",
    )
    p_trans.add_argument(
        "--full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as full lineages",
    )
    p_trans.add_argument(
        "--output_full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Add a column with the full GTDB lineage",
    )
    p_trans.add_argument(
        "--from_taxids",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as NCBI tax IDs",
    )
    p_trans.add_argument(
        "--bundle",
        default=None,
        help="Path to a local .msgpack.zst bundle (downloads latest if omitted)",
    )
    p_trans.add_argument(
        "--version",
        default=None,
        help="GTDB version to download (default: latest release)",
    )

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "build":
        _build(args)
    elif args.command == "translate":
        _translate(args)


if __name__ == "__main__":
    main()
