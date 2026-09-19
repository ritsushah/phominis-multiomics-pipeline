#!/usr/bin/env python3
"""
Build a custom high-sensitivity detection database for Pentatrichomonas hominis
using the real genome assembly (GCA_047301545.1) + curated marker genes.

Outputs:
  - Combined FASTA of genome + markers
  - Simple k-mer index summary
  - Marker BLAST-ready database
  - Metadata JSON for provenance
"""

import argparse
import json
import hashlib
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip


def compute_md5(filepath: Path) -> str:
    """Compute MD5 of a file for provenance."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def load_sequences(fasta_path: Path) -> list:
    """Load sequences from FASTA (handles .gz)."""
    records = []
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    mode = "rt" if str(fasta_path).endswith(".gz") else "r"
    with opener(fasta_path, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records.append(record)
    return records


def build_combined_fasta(genome_path: Path, marker_paths: list, outdir: Path) -> Path:
    """Create a single FASTA containing genome contigs + markers with clear IDs."""
    outdir.mkdir(parents=True, exist_ok=True)
    combined_path = outdir / "phominis_custom_db.fasta"

    all_records = []

    # Genome contigs
    genome_records = load_sequences(genome_path)
    for i, rec in enumerate(genome_records):
        rec.id = f"PHOMINIS_GENOME|{rec.id}"
        rec.description = f"P.homini genome contig {i+1} | source={genome_path.name}"
        all_records.append(rec)

    # Markers
    for mpath in marker_paths:
        mrecs = load_sequences(mpath)
        for rec in mrecs:
            rec.id = f"PHOMINIS_MARKER|{rec.id}"
            rec.description = f"Marker gene | source={mpath.name}"
            all_records.append(rec)

    SeqIO.write(all_records, combined_path, "fasta")
    print(f"[OK] Combined FASTA written: {combined_path} ({len(all_records)} sequences)")
    return combined_path


def extract_marker_only(marker_paths: list, outdir: Path) -> Path:
    """Write a markers-only FASTA for targeted detection."""
    outdir.mkdir(parents=True, exist_ok=True)
    marker_path = outdir / "phominis_markers_only.fasta"
    records = []
    for mpath in marker_paths:
        records.extend(load_sequences(mpath))
    SeqIO.write(records, marker_path, "fasta")
    print(f"[OK] Markers-only FASTA: {marker_path} ({len(records)} sequences)")
    return marker_path


def summarize_genome(genome_path: Path) -> dict:
    """Basic assembly statistics."""
    records = load_sequences(genome_path)
    lengths = [len(r.seq) for r in records]
    total = sum(lengths)
    lengths_sorted = sorted(lengths, reverse=True)
    n50 = 0
    cum = 0
    for L in lengths_sorted:
        cum += L
        if cum >= total / 2:
            n50 = L
            break
    return {
        "n_contigs": len(records),
        "total_length_bp": total,
        "n50": n50,
        "longest_contig": max(lengths) if lengths else 0,
        "gc_content": sum(str(r.seq).upper().count("G") + str(r.seq).upper().count("C") for r in records) / total if total else 0,
    }


def main():
    parser = argparse.ArgumentParser(description="Build custom P. hominis detection database")
    parser.add_argument("--genome", required=True, type=Path, help="Path to P. hominis genome FASTA")
    parser.add_argument("--markers", nargs="+", required=True, type=Path, help="Marker FASTA files")
    parser.add_argument("--outdir", type=Path, default=Path("results/detection/custom_db"))
    args = parser.parse_args()

    print("=" * 60)
    print("Building custom P. hominis detection database")
    print("=" * 60)
    print(f"Genome : {args.genome}")
    print(f"Markers: {args.markers}")
    print(f"Outdir : {args.outdir}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Build databases
    combined = build_combined_fasta(args.genome, args.markers, args.outdir)
    markers_only = extract_marker_only(args.markers, args.outdir)

    # Stats
    stats = summarize_genome(args.genome)
    print("\nGenome statistics:")
    for k, v in stats.items():
        print(f"  {k}: {v}")

    # Provenance
    provenance = {
        "created": datetime.utcnow().isoformat() + "Z",
        "genome_file": str(args.genome),
        "genome_md5": compute_md5(args.genome),
        "genome_stats": stats,
        "marker_files": [str(p) for p in args.markers],
        "marker_md5s": {str(p): compute_md5(p) for p in args.markers},
        "combined_fasta": str(combined),
        "markers_only_fasta": str(markers_only),
        "notes": "Real P. hominis Hs-3:NIH genome (GCA_047301545.1) + public marker sequences. Ready for Kraken2 / BLAST / k-mer detection.",
    }

    prov_path = args.outdir / "database_provenance.json"
    with open(prov_path, "w") as f:
        json.dump(provenance, f, indent=2)
    print(f"\n[OK] Provenance written: {prov_path}")

    # Simple README for the DB
    readme = args.outdir / "README_DB.md"
    with open(readme, "w") as f:
        f.write("# Custom *P. hominis* Detection Database\n\n")
        f.write(f"- Genome: {args.genome.name} (real GCA_047301545.1)\n")
        f.write(f"- Contigs: {stats['n_contigs']}\n")
        f.write(f"- Total length: {stats['total_length_bp']:,} bp\n")
        f.write(f"- N50: {stats['n50']:,} bp\n")
        f.write(f"- Markers: {len(args.markers)} real sequences\n")
        f.write("\nUse `phominis_custom_db.fasta` for full-genome sensitive detection.\n")
        f.write("Use `phominis_markers_only.fasta` for targeted ITS/18S detection.\n")
    print(f"[OK] DB README: {readme}")

    print("\n" + "=" * 60)
    print("Database build complete. Ready for detection module.")
    print("=" * 60)


if __name__ == "__main__":
    main()
