# GenomeDiff — Hilbert Curve Visualizer

Visualise pairwise mitochondrial genome alignments as a space-filling Hilbert curve.

```
fetch_align.py  →  <pair_id>.json  →  visualizer.html
     Part 1                                Part 2
```

The two parts are intentionally decoupled — re-run the renderer any number of
times without re-fetching data from NCBI.

---

## Quick start

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Fetch & align (default: woolly mammoth vs African elephant)
python fetch_align.py

# 3. Open the visualizer in a browser, drop in the generated JSON
open visualizer.html        # macOS
xdg-open visualizer.html    # Linux
start visualizer.html       # Windows
```

The JSON file can be tens of megabytes; modern browsers handle it fine.

---

## Swapping species

Edit **`species_config.json`** — no Python changes needed.

1. Add a new object to the `"pairs"` array (copy an existing entry).
2. Set `"id"`, `"label"`, `"species_a"` and `"species_b"`.
3. Use NCBI RefSeq accession numbers (e.g. `NC_012920` for human mtDNA).
4. Run with `--pair <your_new_id>`.

A human vs chimpanzee example is included in the config as `_swap_example`.

```bash
# Copy the example into the pairs array, then:
python fetch_align.py --pair human_vs_chimp --out human_chimp.json
```

Common mitochondrial accessions:

| Species            | Accession  |
|--------------------|------------|
| Woolly Mammoth     | NC_007596  |
| African Elephant   | NC_000934  |
| Human              | NC_012920  |
| Chimpanzee         | NC_001643  |
| Neanderthal        | NC_011137  |
| Gray Wolf          | NC_002008  |
| Domestic Dog       | NC_002008  |
| House Cat          | NC_001700  |
| Horse              | NC_001640  |

---

## CLI reference — fetch_align.py

```
usage: fetch_align.py [-h] [--pair PAIR_ID] [--config CONFIG] [--out OUTPUT] [--email EMAIL]

Options:
  --pair    Species pair ID from species_config.json  (default: mammoth_vs_elephant)
  --config  Path to species config file               (default: species_config.json)
  --out     Output JSON path                          (default: <pair_id>.json)
  --email   Your email for NCBI Entrez (required by NCBI policy)
```

---

## Visualizer controls

| Control | Description |
|---------|-------------|
| **Hilbert order** | Resolution: order 7 = 128×128 grid; order 9 = 512×512 |
| **Canvas size** | Pixel dimensions; 2048+ for print-quality export |
| **Colors** | Match / mismatch / gap colours, fully customisable |
| **Re-render** | Apply changed settings |
| **Export PNG** | Downloads the canvas at its full pixel resolution |

Hover over the canvas to inspect individual alignment positions.

---

## JSON format

```jsonc
{
  "meta": {
    "pair_id": "mammoth_vs_elephant",
    "label": "Woolly Mammoth vs African Elephant",
    "species_a": { "name": "…", "accession": "…", "length": 16842 },
    "species_b": { "name": "…", "accession": "…", "length": 16820 },
    "alignment_length": 17100,
    "stats": { "match": 15200, "mismatch": 1400, "gap": 500, "identity_pct": 88.9 }
  },
  "colors": { "match": "#2196F3", "mismatch": "#FF5722", "gap": "#9E9E9E" },
  "positions": [
    { "pos": 0, "a": "G", "b": "G", "type": "match" },
    { "pos": 1, "a": "T", "b": "C", "type": "mismatch" },
    …
  ]
}
```
