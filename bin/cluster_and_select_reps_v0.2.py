#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# -----------------------------
# FASTA parsing
# -----------------------------

@dataclass
class Record:
    accession: str
    header: str          # full header line INCLUDING leading ">"
    seq: str
    meta: Dict[str, str]

def fasta_iter(path: Path) -> Iterable[Tuple[str, str]]:
    header: Optional[str] = None
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts)

def parse_header_kv(header: str) -> Tuple[str, Dict[str, str]]:
    """
    Header format expected:
      >OQ993360|species=Pepino_mosaic_virus|strain=DSMZ_PV-0632|...|lineage=k__...;...;f__X;g__Y;s__Z|...
    Returns:
      accession (string before first '|')
      meta dict of key->value
    """
    h = header[1:].strip()
    parts = h.split("|")
    accession = parts[0].strip()
    meta: Dict[str, str] = {}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            meta[k.strip()] = v.strip()
    return accession, meta

def parse_lineage_rank(lineage: str, rank_prefix: str) -> str:
    """
    lineage example:
      k__Viruses;p__Kitrinoviricota;c__Alsuviricetes;o__Tymovirales;f__Alphaflexiviridae;g__Potexvirus_pepini;s__Pepino_mosaic_virus
    rank_prefix: 'f__' or 'g__' or 's__'
    """
    if not lineage or lineage == "NA":
        return "NA"
    for token in lineage.split(";"):
        token = token.strip()
        if token.startswith(rank_prefix):
            return token[len(rank_prefix):] or "NA"
    return "NA"

def load_records(fasta: Path) -> List[Record]:
    recs: List[Record] = []
    for h, seq in fasta_iter(fasta):
        acc, meta = parse_header_kv(h)
        recs.append(Record(accession=acc, header=h, seq=seq.upper(), meta=meta))
    return recs

# -----------------------------
# Clustering helpers
# -----------------------------

def which_or_die(cmd: str) -> None:
    if shutil.which(cmd) is None:
        raise SystemExit(f"[ERROR] Required executable not found on PATH: {cmd}")

def run(cmd: List[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise SystemExit(
            f"[ERROR] Command failed ({p.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )

def normalize_identity(x: float) -> float:
    """
    Accept either:
      1.0 / 0.995 / 0.99
    or:
      100 / 99.5 / 99
    Convert to 0-1 scale.
    """
    if x > 1.0:
        x = x / 100.0
    if x <= 0 or x > 1.0:
        raise ValueError(f"Identity must be in (0,1] or (0,100]; got {x}")
    return x

def pick_word_size(c: float) -> int:
    """
    cd-hit-est typical word sizes:
      c >= 0.99 -> n=10
      c >= 0.95 -> n=8
      c >= 0.90 -> n=7
      else -> n=6 (not used here)
    """
    if c >= 0.99:
        return 10
    if c >= 0.95:
        return 8
    if c >= 0.90:
        return 7
    return 6

def seq_type_priority(seq_type: str) -> int:
    """
    Higher priority = chosen earlier as representative.
    """
    s = (seq_type or "NA").lower()
    if s == "complete_genome":
        return 3
    if s == "complete_segment":
        return 2
    if s == "partial":
        return 1
    return 0

def record_sort_key(r: Record) -> Tuple[int, int]:
    st = r.meta.get("sequence", "NA")
    pri = seq_type_priority(st)
    return (pri, len(r.seq))

def group_key(r: Record, mode: str) -> str:
    mode = mode.lower()
    if mode == "all":
        return "ALL"
    if mode == "species":
        return r.meta.get("species", "NA") or "NA"
    if mode == "genus":
        lin = r.meta.get("lineage", "NA")
        g = parse_lineage_rank(lin, "g__")
        return g or "NA"
    if mode == "family":
        lin = r.meta.get("lineage", "NA")
        f = parse_lineage_rank(lin, "f__")
        return f or "NA"
    raise ValueError(f"Unknown group mode: {mode}")

def write_fasta(records: List[Record], out_fa: Path) -> None:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w", encoding="utf-8") as out:
        for r in records:
            out.write(r.header + "\n")
            out.write(r.seq + "\n")

# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Cluster viral/viroid sequences with CD-HIT-EST and select representatives, optionally within taxonomic groups."
    )
    ap.add_argument("--fasta", required=True, help="Input FASTA with enriched headers")
    ap.add_argument("--outdir", default="clusters_out", help="Output directory")
    ap.add_argument("--group", choices=["species", "genus", "family", "all"], default="species",
                    help="Clustering scope (default: species)")
    ap.add_argument("--identity", nargs="+", type=float, default=[1.0, 0.995, 0.99],
                    help="One or more identity thresholds (e.g. 1.0 0.995 0.99 OR 100 99.5 99). Default: 1.0 0.995 0.99")
    ap.add_argument("--threads", type=int, default=8, help="Threads for cd-hit-est (default 8)")
    ap.add_argument("--memory", type=int, default=0,
                    help="Memory limit for cd-hit-est in MB (-M). 0 = unlimited (default 0)")
    ap.add_argument("--min_length", type=int, default=0,
                    help="Optional minimum sequence length filter (default 0 = keep all)")
    ap.add_argument("--verbose", action="store_true", help="Verbose logging")
    args = ap.parse_args()

    fasta = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    which_or_die("cd-hit-est")

    # Normalize identities
    identities = [normalize_identity(x) for x in args.identity]

    # Load
    recs = load_records(fasta)
    if args.min_length > 0:
        recs = [r for r in recs if len(r.seq) >= args.min_length]

    if args.verbose:
        print(f"[INFO] Loaded records: {len(recs)}", flush=True)

    # Group
    groups: Dict[str, List[Record]] = {}
    for r in recs:
        k = group_key(r, args.group)
        groups.setdefault(k, []).append(r)

    if args.verbose:
        print(f"[INFO] Groups ({args.group}): {len(groups)}", flush=True)

    summary_tsv = outdir / "summary.tsv"
    with summary_tsv.open("w", encoding="utf-8") as sm:
        sm.write("group\tidentity\tinput_n\tclusters_n\treps_fasta\tclstr\n")

        for c in identities:
            reps_all: List[Path] = []

            for gi, (gname, grecs) in enumerate(groups.items(), start=1):
                # Sort so CD-HIT picks the "best" rep (first sequence in a cluster is representative)
                grecs_sorted = sorted(grecs, key=record_sort_key, reverse=True)

                # Write group fasta
                safe_g = re.sub(r"[^A-Za-z0-9_.-]+", "_", gname)[:120]
                grp_dir = outdir / f"{args.group}_groups" / safe_g / f"c{c:.6f}"
                grp_dir.mkdir(parents=True, exist_ok=True)

                in_fa = grp_dir / "input.fasta"
                write_fasta(grecs_sorted, in_fa)

                out_fa = grp_dir / "reps.fasta"
                out_clstr = grp_dir / "reps.fasta.clstr"

                n = pick_word_size(c)

                cmd = [
                    "cd-hit-est",
                    "-i", str(in_fa),
                    "-o", str(out_fa),
                    "-c", f"{c:.6f}",
                    "-n", str(n),
                    "-T", str(args.threads),
                    "-M", str(args.memory),
                    "-d", "0",
                ]

                if args.verbose:
                    print(f"[INFO] ({gi}/{len(groups)}) cd-hit-est group={gname} c={c} n={n} input={len(grecs_sorted)}", flush=True)

                run(cmd)

                reps_all.append(out_fa)

                # Count clusters from .clstr
                clusters_n = 0
                if out_clstr.exists():
                    with out_clstr.open("r", encoding="utf-8", errors="replace") as fh:
                        for line in fh:
                            if line.startswith(">Cluster"):
                                clusters_n += 1

                sm.write(
                    f"{gname}\t{c:.6f}\t{len(grecs_sorted)}\t{clusters_n}\t{out_fa}\t{out_clstr}\n"
                )

            # Merge all reps across groups for this identity
            merged_reps = outdir / f"representatives__{args.group}__c{c:.6f}.fasta"
            with merged_reps.open("w", encoding="utf-8") as out:
                for p in reps_all:
                    with p.open("r", encoding="utf-8", errors="replace") as fh:
                        shutil.copyfileobj(fh, out)

            if args.verbose:
                print(f"[OK] Merged reps for c={c}: {merged_reps}", flush=True)

    print(f"[OK] Summary: {summary_tsv}")
    print(f"[OK] Outputs in: {outdir}")
    print("[TIP] Count reps with: grep -c '^>' representatives__*.fasta")

if __name__ == "__main__":
    main()
