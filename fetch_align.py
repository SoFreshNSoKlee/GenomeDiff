#!/usr/bin/env python3
"""
fetch_align.py — Part 1 of GenomeDiff Visualizer

Fetches two mitochondrial genomes from NCBI, aligns them with BioPython,
and writes a JSON file consumed by the HTML visualizer (Part 2).

Usage:
    python fetch_align.py [--pair PAIR_ID] [--config CONFIG] [--out OUTPUT]
    python fetch_align.py --pair mammoth_vs_african_elephant
    python fetch_align.py --pair dodo_vs_nicobar_pigeon --out dodo_pigeon.json
    python fetch_align.py --pair thylacine_vs_tasmanian_devil

Available pairs (see species_config.json):
    mammoth_vs_african_elephant   Woolly Mammoth vs African Elephant
    mammoth_vs_asian_elephant     Woolly Mammoth vs Asian Elephant
    dodo_vs_nicobar_pigeon        Dodo vs Nicobar Pigeon
    thylacine_vs_dunnart          Thylacine vs Fat-tailed Dunnart
    thylacine_vs_tasmanian_devil  Thylacine vs Tasmanian Devil

Species with no 'accession' in config are resolved automatically via
Entrez esearch — the longest matching sequence is chosen.
"""

import argparse
import json
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO, Align
from Bio.Align import PairwiseAligner


# ── NCBI requires a contact email ────────────────────────────────────────────
Entrez.email = "genomediff@example.com"


def resolve_accession(species: dict, retries: int = 3) -> str:
    """
    Return the accession to use for this species entry.
    Uses 'accession' if present; otherwise runs an Entrez esearch on
    'search_term' to find the longest complete mitochondrial genome.
    """
    if species.get("accession"):
        return species["accession"]

    search_term = species.get("search_term")
    if not search_term:
        raise ValueError(
            f"Species '{species.get('name')}' has neither 'accession' nor 'search_term'."
        )

    print(f"  No accession for {species['name']} — searching: {search_term!r}", flush=True)
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=8, sort="relevance")
            result = Entrez.read(handle)
            handle.close()
            ids = result.get("IdList", [])
            if not ids:
                raise RuntimeError("No results found for search term.")
            # Fetch summaries and pick the longest sequence (most complete genome)
            handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
            summaries = Entrez.read(handle)
            handle.close()
            best = max(summaries, key=lambda s: int(s.get("Length", 0)))
            acc = best["AccessionVersion"]
            print(f"  Resolved: {acc}  ({best.get('Title', '')[:72]})")
            return acc
        except Exception as exc:
            wait = 2 ** attempt
            print(f"  Search attempt {attempt + 1}/{retries} failed: {exc}")
            if attempt < retries - 1:
                time.sleep(wait)
    raise RuntimeError(
        f"Could not resolve accession for '{species.get('name')}' via Entrez search."
    )


def fetch_sequence(accession: str, retries: int = 3) -> str:
    """Fetch a nucleotide sequence from NCBI by accession number."""
    for attempt in range(retries):
        try:
            print(f"  Fetching {accession} from NCBI...", end=" ", flush=True)
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text",
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
            seq = str(record.seq).upper()
            print(f"OK ({len(seq):,} bp)")
            return seq
        except Exception as exc:
            wait = 2 ** attempt
            print(f"FAILED (attempt {attempt + 1}/{retries}): {exc}")
            if attempt < retries - 1:
                print(f"  Retrying in {wait}s...")
                time.sleep(wait)
    raise RuntimeError(f"Could not fetch {accession} after {retries} attempts.")


def align_sequences(seq_a: str, seq_b: str) -> tuple[str, str]:
    """
    Global pairwise alignment via BioPython PairwiseAligner.
    Returns (aligned_a, aligned_b) with gap characters inserted.
    """
    print("  Aligning sequences (global, nucleotide)...", end=" ", flush=True)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq_a, seq_b)
    best = next(iter(alignments))          # highest-scoring alignment

    # Extract the two aligned strings with gaps
    aligned_a, aligned_b = best[0], best[1]
    aligned_a = str(aligned_a)
    aligned_b = str(aligned_b)

    total = len(aligned_a)
    matches = sum(
        1 for a, b in zip(aligned_a, aligned_b) if a == b and a != "-"
    )
    pct = 100.0 * matches / total
    print(f"OK  (alignment length {total:,}, identity {pct:.1f}%)")
    return aligned_a, aligned_b


def build_positions(aligned_a: str, aligned_b: str) -> list[dict]:
    """
    Walk the alignment and return one record per position:
      { "pos": int, "a": char, "b": char, "type": "match"|"mismatch"|"gap" }
    """
    positions = []
    for i, (a, b) in enumerate(zip(aligned_a, aligned_b)):
        if a == "-" or b == "-":
            kind = "gap"
        elif a == b:
            kind = "match"
        else:
            kind = "mismatch"
        positions.append({"pos": i, "a": a, "b": b, "type": kind})
    return positions


def load_pair(config_path: str, pair_id: str) -> dict:
    """Load a species-pair definition from the config file."""
    with open(config_path) as f:
        config = json.load(f)

    for pair in config.get("pairs", []):
        if pair["id"] == pair_id:
            return pair

    available = [p["id"] for p in config.get("pairs", [])]
    raise ValueError(
        f"Pair '{pair_id}' not found in {config_path}. "
        f"Available: {available}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch, align, and export genome diff data to JSON."
    )
    parser.add_argument(
        "--pair",
        default="mammoth_vs_african_elephant",
        help="Species pair ID from species_config.json (default: mammoth_vs_african_elephant)",
    )
    parser.add_argument(
        "--config",
        default="species_config.json",
        help="Path to species config file (default: species_config.json)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output JSON path (default: <pair_id>.json)",
    )
    parser.add_argument(
        "--email",
        default=None,
        help="Email address for NCBI Entrez (required by NCBI policy)",
    )
    args = parser.parse_args()

    if args.email:
        Entrez.email = args.email

    out_path = args.out or f"{args.pair}.json"

    print(f"\n=== GenomeDiff: {args.pair} ===")
    print(f"Config : {args.config}")
    print(f"Output : {out_path}\n")

    # ── 1. Load species pair definition ──────────────────────────────────────
    pair = load_pair(args.config, args.pair)
    sp_a = pair["species_a"]
    sp_b = pair["species_b"]

    print(f"Species A : {sp_a['name']}  [{sp_a.get('accession', 'via search')}]")
    print(f"Species B : {sp_b['name']}  [{sp_b.get('accession', 'via search')}]")
    print()

    # ── 2. Resolve accessions (direct or via esearch) then fetch ─────────────
    acc_a = resolve_accession(sp_a)
    acc_b = resolve_accession(sp_b)
    seq_a = fetch_sequence(acc_a)
    seq_b = fetch_sequence(acc_b)

    # ── 3. Align ──────────────────────────────────────────────────────────────
    aligned_a, aligned_b = align_sequences(seq_a, seq_b)

    # ── 4. Build position list ────────────────────────────────────────────────
    print("  Building position list...", end=" ", flush=True)
    positions = build_positions(aligned_a, aligned_b)
    n_match    = sum(1 for p in positions if p["type"] == "match")
    n_mismatch = sum(1 for p in positions if p["type"] == "mismatch")
    n_gap      = sum(1 for p in positions if p["type"] == "gap")
    print(f"OK  ({n_match:,} match / {n_mismatch:,} mismatch / {n_gap:,} gap)")

    # ── 5. Write JSON ─────────────────────────────────────────────────────────
    output = {
        "meta": {
            "pair_id": pair["id"],
            "label": pair["label"],
            "species_a": {
                "name": sp_a["name"],
                "accession": acc_a,
                "length": len(seq_a),
            },
            "species_b": {
                "name": sp_b["name"],
                "accession": acc_b,
                "length": len(seq_b),
            },
            "alignment_length": len(aligned_a),
            "stats": {
                "match": n_match,
                "mismatch": n_mismatch,
                "gap": n_gap,
                "identity_pct": round(100.0 * n_match / len(aligned_a), 2),
            },
        },
        "colors": {
            "match":    sp_a.get("color_match",    "#2196F3"),
            "mismatch": sp_a.get("color_mismatch", "#FF5722"),
            "gap":      sp_a.get("color_gap",      "#9E9E9E"),
        },
        "positions": positions,
    }

    Path(out_path).write_text(json.dumps(output, separators=(",", ":")))
    size_kb = Path(out_path).stat().st_size / 1024
    print(f"\nWrote {out_path}  ({size_kb:.1f} KB)")
    print("Done. Run the HTML visualizer to render the Hilbert curve.")


if __name__ == "__main__":
    main()
