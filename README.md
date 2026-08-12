# gtdb-translate

Translate taxonomy names from NCBI or SILVA to GTDB, and forward-translate
renamed GTDB names across releases.

Pre-built translation bundles are downloaded automatically from the
[latest release](https://github.com/maltesie/gtdb-translate/releases)
on first use and cached in `~/.cache/gtdb_translate/`.

## Installation

```bash
pip install git+https://github.com/maltesie/gtdb-translate.git
```

## Subcommands

| Command | Input |
| --- | --- |
| `ncbi` | NCBI names, NCBI lineages, or NCBI tax IDs |
| `silva` | SILVA taxon names or SILVA lineages |
| `forward` | GTDB names or lineages from an older release |
| `build` | Raw GTDB metadata + NCBI `names.dmp` (produces a bundle) |

## Output columns

The main output column always holds the full GTDB lineage, whatever shape
the input had. Failed lookups get `no_translation`.

| Column | Written | Contents |
| --- | --- | --- |
| `<out_column_name>` | always | Full GTDB lineage, e.g. `d__Bacteria;…;s__Escherichia coli` |
| `<out_column_name>_lowest` | `--output_lowest_rank` | Most specific translated rank, no rank prefix, e.g. `Escherichia coli` |
| `<out_column_name>_purity` | `--report_purity` (on by default) | Winning share of the vote, `0`–`1`; `1.0` is unanimous. Empty where no votes back the mapping |
| `<out_column_name>_support` | `--report_purity` (on by default) | Total votes behind the mapping. Empty where no votes back the mapping |

`--report_purity` and `--min_purity` are not available for `forward`, which
walks a rename graph rather than a vote tally.

Translations whose purity falls below `--min_purity` (default `0.5`) become
`no_translation`; in lineage mode the walk continues to the next rank up
instead. Mappings with no votes behind them — a GTDB taxon that maps to
itself — are never filtered.

## Auto-detection

All three subcommands auto-detect the column to translate when
`--column_name` is omitted, scoring each string column by how many of its
sampled values are covered by the relevant dictionary. Only one column
can be translated in one go, pass `--column_name` if you have multiple.

They also auto-detect the input format of the chosen column and print what
they inferred. A column counts as holding lineages when its taxa nest —
each resolvable rank an ancestor of the next — in more than half its cells;
the printed percentage is that share. This works for bare lineages and for
every rank-prefix scheme (`d__`, `sk__`, `k__`, `D_0__`). Any inferred
setting can be overridden explicitly:

| Detected | Subcommands | Overridden by |
| --- | --- | --- |
| Values are full lineages rather than single names | `ncbi`, `silva`, `forward` | `--full_lineage` / `--no-full_lineage` |
| Separator between ranks within a lineage | `ncbi`, `silva`, `forward` | `--lineage_sep` |
| Values are NCBI tax IDs rather than names | `ncbi` | `--from_taxids` / `--no-from_taxids` |

`--multi_sep` is never auto-detected. It defaults to none — one value per
cell — and must be set explicitly if your cells hold several entries. When
set, the rank separator is chosen from the remaining candidates.

---

## `ncbi`

```bash
gtdb-translate ncbi \
    --in_file taxa.csv \
    --out_file taxa_gtdb.csv
```

Translate NCBI tax IDs, adding the lowest-rank column:

```bash
gtdb-translate ncbi \
    --in_file taxa.tsv \
    --out_file taxa_gtdb.csv \
    --column_name tax_id \  # define column name explicitly
    --from_taxids \         # treat input as NCBI tax IDs
    --output_lowest_rank    # add column with lowest rank of translation
```

### Arguments

```
--in_file              Input CSV or TSV (required)
--out_file             Output CSV (required)
--column_name          Column to translate (auto-detected if omitted)
--out_column_name      Name for the output column (default: gtdb_translated)
--multi_sep            Separator for multiple entries per cell
                       (default: none — one value per cell)
--lineage_sep          Separator between ranks within a lineage
                       (default: auto-detected, usually ;)
--full_lineage         Treat entries as full lineages (default: auto-detected)
--from_taxids          Treat entries as NCBI tax IDs (default: auto-detected)
--genus_fallback       Fall back to binomial/genus lookups when the exact
                       match fails (default: disabled)
--output_lowest_rank   Add the <out_column>_lowest column (default: disabled)
--min_purity           Reject translations below this vote purity
                       (default: 0.5; pass 0 to disable)
--report_purity        Add the _purity and _support columns (default: enabled)
--empty_on_fail        Write an empty string instead of 'no_translation'
--bundle               Path to a local .msgpack.zst bundle (skips download)
--version              GTDB release to translate against (default: latest)
--force                Re-download the bundle even if cached locally
```

---

## `silva`

```bash
gtdb-translate silva \
    --in_file taxa.csv \
    --out_file taxa_gtdb.csv \
    --column_name silva_taxonomy    # define column name explicitly
```

Accepts bare SILVA paths (`Bacteria;Bacillota;…`) and prefixed ones as
emitted by QIIME/DADA2 (`d__Bacteria;p__Bacillota;…`), as well as single
taxon names.

### Arguments

```
--in_file              Input CSV or TSV (required)
--out_file             Output CSV (required)
--column_name          Column to translate (auto-detected if omitted)
--out_column_name      Name for the output column (default: gtdb_translated)
--multi_sep            Separator for multiple entries per cell
                       (default: none — one value per cell)
--lineage_sep          Separator between ranks within a lineage
                       (default: auto-detected, usually ;)
--full_lineage         Treat entries as full SILVA lineages rather than single
                       taxon names (default: auto-detected)
--genus_fallback       Fall back to binomial/genus lookups when the exact
                       match fails (default: disabled)
--output_lowest_rank   Add the <out_column>_lowest column (default: disabled)
--min_purity           Reject translations below this vote purity
                       (default: 0.5; pass 0 to disable)
--report_purity        Add the _purity and _support columns (default: enabled)
--empty_on_fail        Write an empty string instead of 'no_translation'
--bundle               Path to a local .msgpack.zst bundle (skips download)
--version              GTDB release to translate against (default: latest)
--force                Re-download the bundle even if cached locally
```

### Note on `ncbi` vs `silva`

The two subcommands draw on different evidence: the NCBI dictionary is
built from each genome's curated NCBI lineage, the SILVA dictionary from
SILVA classifications of recovered rRNA genes. The same input name can
therefore translate differently, and carry different `_purity` and
`_support` values, depending on which subcommand you use. Use the one
matching the source of your taxonomy.

---

## `forward`

Update a table carrying GTDB names from an older release. Names already
present in the current GTDB are kept; outdated ones are forward-mapped;
unresolvable ones get `no_translation`. Works at every rank from species
to phylum.

```bash
gtdb-translate forward \
    --in_file old_gtdb.csv \
    --out_file updated_gtdb.csv \
    --column_name gtdb_taxa \       # define column name explicitly
```

### Arguments

```
--in_file              Input CSV or TSV (required)
--out_file             Output CSV (required)
--column_name          Column to forward-translate
                       (auto-detected if omitted)
--out_column_name      Name for the output column (default: gtdb_forwarded)
--multi_sep            Separator for multiple entries per cell
                       (default: none — one value per cell)
--lineage_sep          Separator between ranks within a lineage
                       (default: auto-detected, usually ;)
--full_lineage         Treat entries as full GTDB lineages, prefixed or bare
                       (default: auto-detected)
--output_lowest_rank   Add the <out_column>_lowest column (default: disabled)
--empty_on_fail        Write an empty string instead of 'no_translation'
--bundle               Path to a local .msgpack.zst bundle (skips download)
--version              GTDB release to translate against (default: latest)
--force                Re-download the bundle even if cached locally
```

---

## `build`

```bash
gtdb-translate build \
    --metadata bac120_metadata_r226.tsv ar53_metadata_r226.tsv \
    --names_dmp names.dmp \
    --changelog gtdb-taxid-changelog.csv \
    --version r226 \
    -o gtdb_translate_r226.msgpack.zst
```

Metadata files must be uncompressed. After building, the bundle is scored
against the metadata it came from and a report is printed.

### Arguments

```
--metadata             GTDB metadata TSV(s) (required)
--names_dmp            NCBI names.dmp. Without it, synonym expansion and
                       tax-ID translation are unavailable
--changelog            gtdb-taxid-changelog.csv, enabling forward translation
--silva_columns        Metadata columns to pool SILVA votes from
                       (default: ssu_silva_taxonomy lsu_silva_23s_taxonomy)
--no_silva             Skip the SILVA dictionary entirely
--self_test            Score the finished bundle (default: enabled)
--min_purity           Purity threshold the self-test scores at (default: 0.5)
--self_test_limit      Score only this many genomes
--version              GTDB version label (default: r226)
--output, -o           Output path
                       (default: gtdb_translate_<version>.msgpack.zst)
```

---

## Python API

### Quick start

```python
from gtdb_translate import NCBITranslator, SILVATranslator

t = NCBITranslator.default()          # auto-downloads the latest bundle

t.translate(["Escherichia coli", "Staphylococcus aureus"])
t.translate_ids(["562", "1280"])

# purity and vote count alongside the translation
translations, purity, support = t.translate(
    ["Escherichia coli"], with_support=True, min_purity=0.5
)

# vote statistics for a single name
t.support_for("Escherichia coli")     # → [votes, purity]
```

### SILVA

`SILVATranslator.from_ncbi()` shares an already-loaded bundle; use
`.default()` or `.load()` to read one directly.

```python
s = SILVATranslator.from_ncbi(t)

s.translate([
    "Bacteria;Bacillota;Clostridia;Lachnospirales;"
    "Lachnospiraceae;Lachnospiraceae NK4A136 group"
])

# one lineage, with vote statistics
lineage, support = s.translate_lineage("Bacteria;Pseudomonadota;…")

# one token
s.lookup_token("Escherichia-Shigella", rank="genus")
```

### Forward translation

```python
if t.forward:
    t.forward.translate("Bacillus_C megaterium")
    t.forward.translate_rank("Firmicutes", "phylum")
    t.forward.translate_lineage(
        "d__Bacteria;p__Firmicutes;g__Bacillus_C;s__Bacillus_C megaterium",
        t.gtdb_name_to_lineage,
    )
```

### Building from Python

```python
t = NCBITranslator.build(
    metadata_paths=["bac120_metadata_r226.tsv", "ar53_metadata_r226.tsv"],
    names_dmp_path="names.dmp",
    changelog_path="gtdb-taxid-changelog.csv",
    version="r226",
)
t.save("gtdb_translate_r226.msgpack.zst")
```

### Using components independently

```python
from gtdb_translate import GTDBTaxonomy, ForwardTranslator

tax = GTDBTaxonomy.from_tsv("bac120_taxonomy_r226.tsv")
tax.get_rank("Bacillus subtilis", "phylum")   # → "Bacillota"

fwd = ForwardTranslator()
fwd.build("gtdb-taxid-changelog.csv")
fwd.save("forward.json")

fwd = ForwardTranslator.load("forward.json")
fwd.translate("Bacillus_C megaterium")
fwd.translate_rank("Firmicutes", "phylum")
```
