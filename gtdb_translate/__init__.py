"""gtdb_translate — Translate taxonomy names across GTDB releases and from NCBI.

Repository: https://github.com/maltesie/gtdb-translate

Quick start
-----------

.. code-block:: python

    from gtdb_translate import NCBITranslator

    # Option A: auto-download latest pre-built bundle
    t = NCBITranslator.default()

    # Option B: pin to a specific GTDB version
    t = NCBITranslator.default(version="r226")

    # Option C: load a local bundle
    t = NCBITranslator.load("gtdb_translate_r226.msgpack.zst")

    # Translate NCBI names
    t.translate(["Escherichia coli", "Staphylococcus aureus"])

    # Translate NCBI tax IDs
    t.translate_ids(["562", "1280"])

    # Forward-translate a renamed GTDB species (if changelog was included)
    if t.forward:
        t.forward.translate("Lactobacillus oldname")

    # Translate SILVA lineages, sharing the same loaded bundle
    from gtdb_translate import SILVATranslator

    s = SILVATranslator.from_ncbi(t)
    s.translate(["Bacteria;Bacillota;Clostridia;Lachnospirales;"
                 "Lachnospiraceae;Lachnospiraceae NK4A136 group"])
"""

from .forward import ForwardTranslator
from .ncbi import NCBITranslator
from .silva import SILVATranslator
from .taxonomy import GTDBTaxonomy

__all__ = [
    "GTDBTaxonomy",
    "ForwardTranslator",
    "NCBITranslator",
    "SILVATranslator",
]
