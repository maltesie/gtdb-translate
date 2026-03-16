# gtdb-translate

Translate taxonomy names from NCBI/SILVA to GTDB, and forward-translate
renamed GTDB species across releases.

Pre-built translation dictionaries are downloaded automatically from the
[latest release](https://github.com/maltesie/gtdb-translate/releases)
on first use — no manual setup required.

## Installation

```bash
pip install git+https://github.com/maltesie/gtdb-translate.git
```

## Translating a table

The most common use case: you have a CSV or TSV with a column of NCBI
taxonomy names (or tax IDs, or full lineages) and want to add a column
with the corresponding GTDB names.

```bash
gtdb-translate translate \
    --in_file my_data.tsv \
    --out_file my_data_translated.csv \
    --column_name taxonomy
```

On first run, the tool automatically downloads the translation bundle
for the latest GTDB release and caches it in `~/.cache/gtdb_translate/`.

If you omit `--column_name`, the tool will try to auto-detect which
column contains translatable names by sampling a few rows against the
dictionary.

### Options

```
--column_name       Column to translate (auto-detected if omitted)
--out_column_name   Name for the output column (default: gtdb_translated)
--sep               Separator for multiple names per cell (default: |)
--lineage_sep       Separator within a lineage string (default: ;)
--full_lineage      Treat entries as full lineages
--output_full_lineage  Add a column with the full GTDB lineage
--from_taxids       Treat entries as NCBI tax IDs instead of names
--version           GTDB release to target (default: latest release)
--bundle            Path to a local bundle file (skips download)
```

### Examples

Translate a column of species names:

```bash
gtdb-translate translate \
    --in_file samples.csv \
    --out_file samples_gtdb.csv \
    --column_name species
```

Translate NCBI tax IDs with full lineage output:

```bash
gtdb-translate translate \
    --in_file otus.tsv \
    --out_file otus_gtdb.csv \
    --column_name tax_id \
    --from_taxids \
    --full_lineage
```

Translate full NCBI lineages (e.g. from SILVA or QIIME):

```bash
gtdb-translate translate \
    --in_file silva_table.csv \
    --out_file silva_gtdb.csv \
    --column_name lineage \
    --full_lineage \
    --lineage_sep ";"
```

Pin to a specific GTDB release:

```bash
gtdb-translate translate \
    --in_file data.tsv \
    --out_file data_gtdb.csv \
    --column_name taxonomy \
    --version r226
```

## Building your own bundle

If you want to build translation dictionaries yourself (e.g. for a
newer GTDB release before an official bundle is published), you need
the GTDB metadata TSVs, NCBI `names.dmp`, and optionally the
[gtdb-taxdump](https://github.com/shenwei356/gtdb-taxdump) changelog
for forward translation of renamed species:

```bash
gtdb-translate build \
    --metadata bac120_metadata_r226.tsv ar53_metadata_r226.tsv \
    --names_dmp names.dmp \
    --changelog gtdb-taxid-changelog.csv \
    --version r226 \
    -o gtdb_translate_r226.msgpack.zst
```

The resulting `.msgpack.zst` bundle can be used locally with
`--bundle path/to/file` or uploaded as a GitHub release for others to
download automatically.

## Bundle format

Bundles are serialized with msgpack + zstandard for fast loading and
compact size. Legacy `translation_dicts_rXXX.json.gz` files (from
earlier versions of this tool) are also supported.

---

## Python API

All functionality is also available as a Python package for use in
scripts and pipelines.

### Quick start

```python
from gtdb_translate import NCBITranslator

# Auto-downloads the latest bundle on first use (~/.cache/gtdb_translate/)
t = NCBITranslator.default()

# Or pin to a specific version
t = NCBITranslator.default(version="r226")

# Translate NCBI names → GTDB
t.translate(["Escherichia coli", "Staphylococcus aureus"])

# Translate NCBI tax IDs → GTDB
t.translate_ids(["562", "1280"])

# Forward-translate renamed GTDB species (included in the bundle)
if t.forward:
    t.forward.translate("Lactobacillus oldname")
```

### Building from Python

```python
from gtdb_translate import NCBITranslator

t = NCBITranslator(version="r226")
t.build(
    metadata_paths=["bac120_metadata_r226.tsv", "ar53_metadata_r226.tsv"],
    names_dmp_path="names.dmp",
    changelog_path="gtdb-taxid-changelog.csv",
)
t.save("gtdb_translate_r226.msgpack.zst")
```

### Loading a local bundle

```python
t = NCBITranslator.load("gtdb_translate_r226.msgpack.zst")
```

### Using components independently

```python
from gtdb_translate import GTDBTaxonomy, ForwardTranslator

# GTDB taxonomy index
tax = GTDBTaxonomy.from_tsv("bac120_taxonomy_r226.tsv")
tax.get_rank("Bacillus subtilis", "phylum")  # → "Bacillota"

# Forward translator (standalone, without NCBI dicts)
fwd = ForwardTranslator()
fwd.build("gtdb-taxid-changelog.csv")
fwd.save("forward.json")

fwd = ForwardTranslator.load("forward.json")
fwd.translate("Lactobacillus oldname")
```
