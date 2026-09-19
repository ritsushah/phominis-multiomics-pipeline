#!/usr/bin/env python3
"""
Detect Pentatrichomonas hominis in metagenomic reads using marker genes
and simple k-mer matching against the real custom database.

This is a lightweight, pure-Python demonstration module that works without
external aligners. For production scale, replace the core matching with
Kraken2 / Bowtie2 / minimap2 against the same database.
"""

import argparse
import json
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from datetime import datetime


def load_kmers(fasta_path: Path, k: int = 31) -> set:
    """Extract unique k-mers from a FASTA (genome or markers)."""
    kmers = set()
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    mode = "rt" if str(fasta_path).endswith(".gz") else "r"
    with opener(fasta_path, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper().replace("N", "")
            for i in range(len(seq) - k + 1):
                kmers.add(seq[i:i+k])
                # also reverse complement for double-stranded matching
                rc = str(Seq(seq[i:i+k]).reverse_complement())
                kmers.add(rc)
    return kmers


def load_markers(marker_fasta: Path) -> dict:
    """Load marker sequences as dict id -> sequence."""
    markers = {}
    with open(marker_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            markers[record.id] = str(record.seq).upper()
    return markers


def count_marker_hits(fastq_paths: list, markers: dict, min_identity: float = 0.9) -> dict:
    """
    Simple exact / near-exact substring search for markers in reads.
    For real production use Bowtie2 or BLAST; this demonstrates the concept.
    """
    hits = defaultdict(int)
    total_reads = 0
    opener = gzip.open if any(str(p).endswith(".gz") for p in fastq_paths) else open

    for fq in fastq_paths:
        mode = "rt" if str(fq).endswith(".gz") else "r"
        with opener(fq, mode) as handle:
            # Assume interleaved or single; simple line-based for demo
            for i, line in enumerate(handle):
                if i % 4 == 1:  # sequence line
                    total_reads += 1
                    read = line.strip().upper()
                    for mid, mseq in markers.items():
                        # exact substring (highly specific for short markers)
                        if mseq in read or str(Seq(mseq).reverse_complement()) in read:
                            hits[mid] += 1
                        # also check short windows for partial matches
                        elif len(mseq) > 50:
                            for j in range(0, len(mseq) - 40, 20):
                                window = mseq[j:j+40]
                                if window in read:
                                    hits[mid] += 1
                                    break
    return dict(hits), total_reads


def kmer_coverage(fastq_paths: list, ref_kmers: set, k: int = 31, max_reads: int = 100000) -> dict:
    """Estimate fraction of reference k-mers observed in the sample (lightweight)."""
    observed = set()
    total_reads = 0
    opener = gzip.open if any(str(p).endswith(".gz") for p in fastq_paths) else open

    for fq in fastq_paths:
        mode = "rt" if str(fq).endswith(".gz") else "r"
        with opener(fq, mode) as handle:
            for i, line in enumerate(handle):
                if i % 4 == 1:
                    total_reads += 1
                    if total_reads > max_reads:
                        break
                    read = line.strip().upper()
                    for j in range(0, len(read) - k + 1, k // 2):  # stride for speed
                        kmer = read[j:j+k]
                        if kmer in ref_kmers:
                            observed.add(kmer)
            if total_reads > max_reads:
                break

    coverage = len(observed) / len(ref_kmers) if ref_kmers else 0.0
    return {
        "observed_kmers": len(observed),
        "total_ref_kmers": len(ref_kmers),
        "coverage_fraction": coverage,
        "reads_examined": total_reads,
    }


def main():
    parser = argparse.ArgumentParser(description="Detect P. hominis in metagenomic FASTQs")
    parser.add_argument("--fastq", nargs="+", required=True, type=Path, help="Input FASTQ(s)")
    parser.add_argument("--db", type=Path, required=True, help="Path to custom_db directory")
    parser.add_argument("--outdir", type=Path, default=Path("results/detection"))
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--max_reads_for_kmer", type=int, default=50000)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("P. hominis Detection Module (real custom database)")
    print("=" * 60)

    markers_fasta = args.db / "phominis_markers_only.fasta"
    combined_fasta = args.db / "phominis_custom_db.fasta"

    if not markers_fasta.exists():
        raise FileNotFoundError(f"Markers FASTA not found: {markers_fasta}. Run build_custom_db.py first.")

    # Load markers
    markers = load_markers(markers_fasta)
    print(f"Loaded {len(markers)} marker sequences")

    # Marker hits
    print("\nScanning for marker gene hits...")
    marker_hits, total_reads = count_marker_hits(args.fastq, markers)
    print(f"Total reads examined: {total_reads:,}")
    print("Marker hits:")
    for mid, count in marker_hits.items():
        print(f"  {mid}: {count}")

    # Lightweight k-mer coverage (using markers for speed; full genome is large)
    print("\nEstimating k-mer coverage (markers)...")
    marker_kmers = load_kmers(markers_fasta, k=args.k)
    kmer_stats = kmer_coverage(args.fastq, marker_kmers, k=args.k, max_reads=args.max_reads_for_kmer)
    print(f"  Observed k-mers: {kmer_stats['observed_kmers']} / {kmer_stats['total_ref_kmers']}")
    print(f"  Coverage fraction: {kmer_stats['coverage_fraction']:.4f}")

    # Decision logic
    total_marker_hits = sum(marker_hits.values())
    detected = total_marker_hits >= 1 or kmer_stats["coverage_fraction"] > 0.01
    confidence = "HIGH" if total_marker_hits >= 5 else ("MEDIUM" if total_marker_hits >= 1 else "LOW")

    result = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "fastq_files": [str(p) for p in args.fastq],
        "total_reads_examined": total_reads,
        "marker_hits": marker_hits,
        "total_marker_hits": total_marker_hits,
        "kmer_stats": kmer_stats,
        "detected": detected,
        "confidence": confidence,
        "notes": "Lightweight demonstration detector. For production use Kraken2/Bowtie2 against the same custom DB.",
    }

    out_json = args.outdir / "detection_result.json"
    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n[OK] Result written: {out_json}")

    # Human-readable report
    report = args.outdir / "detection_report.txt"
    with open(report, "w") as f:
        f.write("P. hominis Detection Report\n")
        f.write("=" * 40 + "\n")
        f.write(f"Detected: {detected}\n")
        f.write(f"Confidence: {confidence}\n")
        f.write(f"Total marker hits: {total_marker_hits}\n")
        f.write(f"K-mer coverage (markers): {kmer_stats['coverage_fraction']:.4f}\n")
        f.write(f"Reads examined: {total_reads:,}\n")
        f.write("\nMarker breakdown:\n")
        for mid, count in marker_hits.items():
            f.write(f"  {mid}: {count}\n")
    print(f"[OK] Report written: {report}")

    print("\n" + "=" * 60)
    print(f"Detection complete. Detected = {detected} ({confidence})")
    print("=" * 60)


if __name__ == "__main__":
    main()
