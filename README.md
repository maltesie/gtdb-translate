# gtdb-translate

Translate taxonomy names from NCBI/SILVA to GTDB, and forward-translate
renamed GTDB names across releases.

Pre-built translation dictionaries are downloaded automatically from the
[latest release](https://github.com/maltesie/gtdb-translate/releases)
on first use — no manual setup required.

## Installation

```bash
pip install git+https://github.com/maltesie/gtdb-translate.git
```

## Translating NCBI/SILVA names to GTDB

```bash
gtdb-translate ncbi \
    --in_file my_data.tsv \
    --out_file my_data_translated.csv \
    --column_name taxonomy
```

On first run, the tool automatically downloads the translation bundle
for the latest GTDB release and caches it in `~/.cache/gtdb_translate/`.

If you omit `--column_name`, the tool will try to auto-detect which
column contains translatable names.

### Options

```
--column_name       Column to translate (auto-detected if omitted)
--out_column_name   Name for the output column (default: gtdb_translated)
--sep               Separator for multiple names per cell (default: |)
--lineage_sep       Separator within a lineage string (default: ;)
--full_lineage      Treat entries as full lineages
--output_full_lineage  Add a column with the full GTDB lineage
--from_taxids       Treat entries as NCBI tax IDs instead of names
--from_silva        Sanitize SILVA lineages before translation
--genus_fallback    Fall back to genus when species match fails (default: off)
--version           GTDB release to translate against (default: latest)
--bundle            Path to a local bundle file (skips download)
```

### Examples

Translate a column of species names:

```bash
gtdb-translate ncbi \
    --in_file samples.csv \
    --out_file samples_gtdb.csv \
    --column_name species
```

Translate NCBI tax IDs with full lineage output:

```bash
gtdb-translate ncbi \
    --in_file otus.tsv \
    --out_file otus_gtdb.csv \
    --column_name tax_id \
    --from_taxids \
    --full_lineage
```

Translate SILVA lineages:

```bash
gtdb-translate ncbi \
    --in_file silva_table.csv \
    --out_file silva_gtdb.csv \
    --column_name lineage \
    --from_silva \
    --full_lineage \
    --lineage_sep ";"
```

## Forward-translating old GTDB names

If you have a table with GTDB names from an older release:

```bash
gtdb-translate forward \
    --in_file old_data.csv \
    --out_file updated_data.csv \
    --column_name host_species \
    --output_full_lineage
```

This checks each name against the current GTDB. Names already present
are kept as-is. Outdated names are forward-mapped to their current
equivalents. Names that cannot be resolved get `no_translation`.

Works at all taxonomic ranks — species, genus, family, order, class,
and phylum names are all forward-translated.

## Building your own bundle

```bash
gtdb-translate build \
    --metadata bac120_metadata_r226.tsv ar53_metadata_r226.tsv \
    --names_dmp names.dmp \
    --changelog gtdb-taxid-changelog.csv \
    --version r226 \
    -o gtdb_translate_r226.msgpack.zst
```

## Bundle format

Bundles are serialized with msgpack + zstandard for fast loading and
compact size. Legacy `translation_dicts_rXXX.json.gz` files are also
supported by `NCBITranslator.load()`.

---

## Python API

### Quick start

```python
from gtdb_translate import NCBITranslator

# Auto-downloads the latest bundle
t = NCBITranslator.default()

# Translate NCBI names → GTDB
t.translate(["Escherichia coli", "Staphylococcus aureus"])

# Translate NCBI tax IDs → GTDB
t.translate_ids(["562", "1280"])

# Forward-translate old GTDB species names
if t.forward:
    t.forward.translate("Bacillus_C megaterium")

    # Forward-translate at a specific rank
    t.forward.translate_rank("Firmicutes", "phylum")

    # Forward-translate a full lineage (bottom-up, validates against current GTDB)
    t.forward.translate_lineage(
        "d__Bacteria;p__Firmicutes;g__Bacillus_C;s__Bacillus_C megaterium",
        t.gtdb_name_to_lineage,
    )
```

### Building from Python

```python
t = NCBITranslator(version="r226")
t.build(
    metadata_paths=["bac120_metadata_r226.tsv", "ar53_metadata_r226.tsv"],
    names_dmp_path="names.dmp",
    changelog_path="gtdb-taxid-changelog.csv",
)
t.save("gtdb_translate_r226.msgpack.zst")
```

### Using components independently

```python
from gtdb_translate import GTDBTaxonomy, ForwardTranslator

# GTDB taxonomy index
tax = GTDBTaxonomy.from_tsv("bac120_taxonomy_r226.tsv")
tax.get_rank("Bacillus subtilis", "phylum")  # → "Bacillota"

# Forward translator (standalone)
fwd = ForwardTranslator()
fwd.build("gtdb-taxid-changelog.csv")
fwd.save("forward.json")

fwd = ForwardTranslator.load("forward.json")
fwd.translate("Bacillus_C megaterium")
fwd.translate_rank("Firmicutes", "phylum")
```
