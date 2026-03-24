"""Command-line interface for gtdb-translate.

Subcommands
-----------
build    Build a translation bundle from raw GTDB + NCBI files.
ncbi     Batch-translate NCBI/SILVA names or tax IDs to GTDB.
forward  Forward-translate old GTDB names to the current release.
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
        print(f"  Forward translations: {translator.forward}")


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


# -- ncbi subcommand (was: translate) -------------------------------------

def _ncbi(args: argparse.Namespace) -> None:
    import pandas as pd

    from .ncbi import NCBITranslator

    translator = _load_translator(args)
    df = _read_table(args.in_file)

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

    # Sanitize SILVA lineages (opt-in)
    if args.from_silva and not args.from_taxids:
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
            column, sep=args.lineage_sep, full_lineage=args.full_lineage,
            genus_fallback=args.genus_fallback,
        )

    # Optional extras
    if args.full_lineage:
        df[column_name + "_lowest"] = [
            t.split(args.lineage_sep)[-1] for t in df[args.out_column_name]
        ]

    if args.output_full_lineage:
        def lookup_lineage_single(name):
            if not name or name == "no_translation":
                return ""
            for prefix in ("s__", "g__", "f__", "o__", "c__", "p__", "d__"):
                result = translator.gtdb_name_to_lineage.get(prefix + name)
                if result:
                    return result
            return ""

        def lookup_lineage(entry):
            if not entry or entry == "no_translation":
                return ""
            parts = entry.split(args.sep)
            return args.sep.join(lookup_lineage_single(p.strip()) for p in parts)

        df[args.out_column_name + "_lineage"] = [
            lookup_lineage(t) for t in df[args.out_column_name]
        ]

    if args.empty_on_fail:
        df[args.out_column_name] = df[args.out_column_name].replace("no_translation", "")

    df.to_csv(args.out_file, index=False)
    print(f"Output written to {args.out_file}")

def _forward(args: argparse.Namespace) -> None:
    import pandas as pd

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

    if column_name not in df.columns:
        print(
            f"ERROR: Column '{column_name}' not found in {args.in_file}. "
            f"Available: {', '.join(df.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    NO_TRANS = "no_translation"
    fwd = translator.forward
    lineage_dict = translator.gtdb_name_to_lineage

    # Prefixes to try, from most specific to least
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

        name = name.strip()

        # 1. Already in current GTDB?
        for prefix, rank in RANK_PREFIXES:
            if (prefix + name) in lineage_dict:
                return name

        # 2. Try forward mapping
        # Species-level (contains a space)
        if " " in name:
            mapped = fwd.translate(name)
            if mapped != name and ("s__" + mapped) in lineage_dict:
                return mapped
        # Try all rank dicts
        for prefix, rank in RANK_PREFIXES[1:]:  # skip species
            mapped = fwd.translate_rank(name, rank)
            if mapped != name and (prefix + mapped) in lineage_dict:
                return mapped

        return NO_TRANS

    sep = args.sep

    def forward_entry(entry: str) -> str:
        """Forward-translate an entry that may contain multiple names."""
        if not isinstance(entry, str) or not entry.strip():
            return NO_TRANS
        parts = entry.split(sep)
        if args.full_lineage:
            translated = []
            for part in parts:
                part = part.strip()
                if not part:
                    translated.append(NO_TRANS)
                    continue
                result = fwd.translate_lineage(part, lineage_dict)
                translated.append(result if result is not None else NO_TRANS)
        else:
            translated = [forward_one(p.strip()) for p in parts]
        if all(t == NO_TRANS for t in translated):
            return NO_TRANS
        return sep.join(translated)

    df[args.out_column_name] = [forward_entry(v) for v in df[column_name]]

    # Count results
    n_total = len(df)
    n_translated = (df[args.out_column_name] != NO_TRANS).sum()
    n_unchanged = (df[args.out_column_name] == df[column_name]).sum()
    n_mapped = n_translated - n_unchanged
    print(f"Forward mapping: {n_total} entries")
    print(f"  Already current: {n_unchanged}")
    print(f"  Forward-mapped:  {n_mapped}")
    print(f"  No translation:  {n_total - n_translated}")

    # Optional lineage column (skip if --full_lineage, output is already a lineage)
    if args.output_full_lineage and not args.full_lineage:
        def get_lineage_single(name):
            if name == NO_TRANS:
                return ""
            for prefix, _ in RANK_PREFIXES:
                result = lineage_dict.get(prefix + name)
                if result:
                    return result
            return ""

        def get_lineage(entry):
            if entry == NO_TRANS:
                return ""
            parts = entry.split(sep)
            return sep.join(get_lineage_single(p.strip()) for p in parts)

        df[args.out_column_name + "_lineage"] = [
            get_lineage(t) for t in df[args.out_column_name]
        ]

    if args.empty_on_fail:
        df[args.out_column_name] = df[args.out_column_name].replace("no_translation", "")

    df.to_csv(args.out_file, index=False)
    print(f"Output written to {args.out_file}")


# -- CLI entry point -------------------------------------------------------

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

    # -- ncbi ----------------------------------------------------------
    p_ncbi = subparsers.add_parser(
        "ncbi",
        help="Batch-translate NCBI/SILVA names or tax IDs to GTDB.",
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
        "--sep",
        default="|",
        help="Separator for multiple names per cell (default: |)",
    )
    p_ncbi.add_argument(
        "--lineage_sep",
        default=";",
        help="Separator within a lineage string (default: ;)",
    )
    p_ncbi.add_argument(
        "--full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as full lineages",
    )
    p_ncbi.add_argument(
        "--output_full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Add a column with the full GTDB lineage",
    )
    p_ncbi.add_argument(
        "--from_taxids",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as NCBI tax IDs",
    )
    p_ncbi.add_argument(
        "--from_silva",
        action=argparse.BooleanOptionalAction,
        help="Sanitize input as SILVA lineages (remove sk__, replace _ with space, merge genus into species)",
    )
    p_ncbi.add_argument(
        "--genus_fallback",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Fall back to genus-level when exact species match fails (default: disabled)",
    )
    p_ncbi.add_argument(
        "--empty_on_fail",
        action="store_true",
        help="Output empty string instead of 'no_translation' for failed lookups",
    )
    _add_bundle_args(p_ncbi)

    # -- forward -------------------------------------------------------
    p_fwd = subparsers.add_parser(
        "forward",
        help="Forward-translate old GTDB names to the current release.",
    )
    p_fwd.add_argument("--in_file", required=True, help="Input CSV or TSV")
    p_fwd.add_argument("--out_file", required=True, help="Output CSV")
    p_fwd.add_argument(
        "--column_name",
        required=True,
        help="Column containing old GTDB names to forward-translate",
    )
    p_fwd.add_argument(
        "--out_column_name",
        default="gtdb_forwarded",
        help="Name for the output column (default: gtdb_forwarded)",
    )
    p_fwd.add_argument(
        "--sep",
        default="|",
        help="Separator for multiple names per cell (default: |)",
    )
    p_fwd.add_argument(
        "--full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Treat entries as full GTDB lineages (e.g. d__Bacteria;p__Firmicutes;...)",
    )
    p_fwd.add_argument(
        "--output_full_lineage",
        action=argparse.BooleanOptionalAction,
        help="Add a column with the full current GTDB lineage",
    )
    p_fwd.add_argument(
        "--empty_on_fail",
        action="store_true",
        help="Output empty string instead of 'no_translation' for failed lookups",
    )
    _add_bundle_args(p_fwd)

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "build":
        _build(args)
    elif args.command == "ncbi":
        _ncbi(args)
    elif args.command == "forward":
        _forward(args)


if __name__ == "__main__":
    main()
