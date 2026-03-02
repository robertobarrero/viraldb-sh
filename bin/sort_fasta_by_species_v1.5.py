#!/usr/bin/env python3

import argparse
from pathlib import Path
from collections import defaultdict
import re
import sys
from typing import Iterable, Tuple, Optional

SPECIES_RE = re.compile(r"\|species=([^|]+)")
LENGTH_RE = re.compile(r"\|length=([0-9]+)")


def parse_fasta(fasta_path: Path) -> Iterable[Tuple[str, str]]:
    """Yield (header, sequence) tuples from FASTA."""
    header = None
    seq_chunks = []

    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_chunks)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header:
            yield header, "".join(seq_chunks)


def iter_merged_fastas(f1: Path, f2: Optional[Path]) -> Iterable[Tuple[str, str]]:
    yield from parse_fasta(f1)
    if f2 is not None:
        yield from parse_fasta(f2)


def extract_species(header: str) -> str:
    m = SPECIES_RE.search(header)
    return m.group(1) if m else "NA"


def extract_length(header: str, seq: str) -> int:
    m = LENGTH_RE.search(header)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            pass
    return len(seq)


def default_merged_name(f1: Path, f2: Path) -> str:
    return f"merged_{f1.stem}_and_{f2.stem}.fasta"


def main():
    ap = argparse.ArgumentParser(
        description="Sort FASTA by species and filter by length and N content. Optionally merge two FASTA inputs."
    )
    ap.add_argument("--fasta", required=True, help="Input FASTA #1")
    ap.add_argument("--fasta2", default=None, help="Optional input FASTA #2")
    ap.add_argument("--out_sorted_fasta", required=True, help="Output FASTA (full path allowed)")

    ap.add_argument(
        "--out_merged_fasta",
        default=None,
        help="Optional output path for merged FASTA. "
             "If omitted, a default name based on input FASTA filenames is used.",
    )

    ap.add_argument("--min_len", type=int, default=200, help="Minimum sequence length (default: 200)")
    ap.add_argument("--max_Ns", type=int, default=None, help="Maximum allowed Ns (absolute)")
    ap.add_argument("--max_Ns_fraction", type=float, default=3.0,
                    help="Maximum percent Ns allowed (default: 3.0)")

    args = ap.parse_args()

    fasta1 = Path(args.fasta)
    fasta2 = Path(args.fasta2) if args.fasta2 else None
    fasta_out = Path(args.out_sorted_fasta)

    if not fasta1.exists():
        sys.exit(f"[ERROR] Input FASTA not found: {fasta1}")
    if fasta2 and not fasta2.exists():
        sys.exit(f"[ERROR] Input FASTA2 not found: {fasta2}")

    fasta_out.parent.mkdir(parents=True, exist_ok=True)

    removed_short = fasta_out.with_name("removed_short.tsv")
    removed_ns = fasta_out.with_name("removed_Ns.tsv")

    merged_fasta_path = None
    if fasta2:
        if args.out_merged_fasta:
            merged_fasta_path = Path(args.out_merged_fasta)
        else:
            merged_fasta_path = fasta_out.parent / default_merged_name(fasta1, fasta2)

    by_species = defaultdict(list)

    total = kept = dropped_short = dropped_ns = 0

    merged_out = merged_fasta_path.open("w", encoding="utf-8") if merged_fasta_path else None

    try:
        with removed_short.open("w") as rs, removed_ns.open("w") as rn:
            rs.write("header\tlength\n")
            rn.write("header\tlength\tNs\tNs_percent\treason\n")

            for header, seq in iter_merged_fastas(fasta1, fasta2):
                total += 1
                seq = seq.upper()

                if merged_out:
                    merged_out.write(f"{header}\n{seq}\n")

                seqlen = extract_length(header, seq)

                if seqlen < args.min_len:
                    rs.write(f"{header}\t{seqlen}\n")
                    dropped_short += 1
                    continue

                Ns = seq.count("N")
                Ns_pct = (Ns / seqlen) * 100 if seqlen else 0.0

                if args.max_Ns is not None and Ns > args.max_Ns:
                    rn.write(f"{header}\t{seqlen}\t{Ns}\t{Ns_pct:.2f}\tNs_absolute\n")
                    dropped_ns += 1
                    continue

                if Ns_pct > args.max_Ns_fraction:
                    rn.write(f"{header}\t{seqlen}\t{Ns}\t{Ns_pct:.2f}\tNs_fraction\n")
                    dropped_ns += 1
                    continue

                species = extract_species(header)
                by_species[species].append((header, seq))
                kept += 1
    finally:
        if merged_out:
            merged_out.close()

    with fasta_out.open("w") as out:
        for species in sorted(by_species):
            for header, seq in by_species[species]:
                out.write(f"{header}\n{seq}\n")

    print("[OK] FASTA processing complete")
    print(f"  Input FASTA #1     : {fasta1}")
    if fasta2:
        print(f"  Input FASTA #2     : {fasta2}")
    print(f"  Output FASTA       : {fasta_out}")
    if merged_fasta_path:
        print(f"  Merged FASTA       : {merged_fasta_path}")
    print(f"  Min length         : {args.min_len}")
    print(f"  Max Ns (absolute)  : {args.max_Ns}")
    print(f"  Max Ns fraction %  : {args.max_Ns_fraction}")
    print(f"  Total read         : {total}")
    print(f"  Kept               : {kept}")
    print(f"  Dropped (short)    : {dropped_short}")
    print(f"  Dropped (Ns)       : {dropped_ns}")
    print(f"  Logs:")
    print(f"    - {removed_short}")
    print(f"    - {removed_ns}")


if __name__ == "__main__":
    main()
