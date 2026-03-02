#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


# -------------------------
# FASTA parsing
# -------------------------
def fasta_iter(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def header_id(header: str) -> str:
    # CD-HIT stores headers truncated in .clstr display, but the ID is the first token after '>'
    # Your headers: >OQ993360|species=...
    h = header[1:].strip()
    return h.split()[0].split("|")[0]


# -------------------------
# Header field parsing
# -------------------------
PAIR_RE = re.compile(r"\|([^=|]+)=([^|]*)")

def parse_fields(header: str) -> Dict[str, str]:
    # Parses |key=value|key=value ... into a dict
    d: Dict[str, str] = {}
    for m in PAIR_RE.finditer(header):
        k = (m.group(1) or "").strip()
        v = (m.group(2) or "").strip()
        if k:
            d[k] = v if v else "NA"
    return d


def is_complete_like(fields: Dict[str, str]) -> bool:
    seqtype = (fields.get("sequence", "NA") or "NA").lower()
    return seqtype in {"complete_genome", "complete_segment"}


def has_strain(fields: Dict[str, str]) -> bool:
    s = (fields.get("strain", "NA") or "NA").strip()
    return s != "" and s.upper() != "NA"


def accession_looks_refseq(acc: str) -> bool:
    # crude but useful: RefSeq accessions often have prefixes like NC_, NR_, NM_, NP_, etc.
    return bool(re.match(r"^(NC|NR|NM|NP|NG|NT|NW|XM|XR|XP)_\d+", acc))


def is_refseq(fields: Dict[str, str], acc: str) -> bool:
    src = (fields.get("source", "NA") or "NA").lower()
    if src in {"refseq"}:
        return True
    # Some pipelines store sourceDatabase separately; keep fallback
    if "refseq" in src:
        return True
    return accession_looks_refseq(acc)


def get_tax_label(fields: Dict[str, str]) -> str:
    # prefer species, then genus, else NA
    sp = fields.get("species", "NA") or "NA"
    if sp != "NA":
        return sp
    # genus is embedded in lineage as g__Something; if present in lineage
    lineage = fields.get("lineage", "NA") or "NA"
    m = re.search(r"g__([^;|]+)", lineage)
    if m:
        g = m.group(1).strip()
        return g if g else "NA"
    return "NA"


# -------------------------
# CD-HIT .clstr parsing
# -------------------------
@dataclass
class Cluster:
    cid: str
    members: List[str]          # sequence IDs (accession w/out extra fields)
    cdhit_rep: str              # seq ID marked with '*'


CLSTR_HEADER_RE = re.compile(r"^>Cluster\s+(\d+)\s*$")
CLSTR_LINE_RE = re.compile(r"^\s*\d+\s+\d+nt,\s+>([^\.]+)\.\.\.\s*(\*|at\s+.+)?\s*$")

def parse_cdhit_clstr(clstr_path: Path) -> List[Cluster]:
    clusters: List[Cluster] = []
    current_id: Optional[str] = None
    current_members: List[str] = []
    current_rep: Optional[str] = None

    with clstr_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = CLSTR_HEADER_RE.match(line)
            if m:
                # flush previous
                if current_id is not None:
                    if not current_rep and current_members:
                        current_rep = current_members[0]
                    clusters.append(Cluster(cid=current_id, members=current_members, cdhit_rep=current_rep or "NA"))
                current_id = m.group(1)
                current_members = []
                current_rep = None
                continue

            if not line.strip():
                continue

            # Example line:
            # 0   23635nt, >NC_031463|species=C... *
            # But CD-HIT truncates; we only get the prefix before "..."
            # We'll take token before first '|' if present.
            if ">" in line and "..." in line:
                try:
                    after_gt = line.split(">", 1)[1]
                    prefix = after_gt.split("...", 1)[0]  # e.g. "NC_031463|species=C"
                    seqid = prefix.split("|", 1)[0].strip()
                    if seqid:
                        current_members.append(seqid)
                        if line.strip().endswith("*"):
                            current_rep = seqid
                except Exception:
                    continue

    # flush last
    if current_id is not None:
        if not current_rep and current_members:
            current_rep = current_members[0]
        clusters.append(Cluster(cid=current_id, members=current_members, cdhit_rep=current_rep or "NA"))

    return clusters


# -------------------------
# Phase2 selection
# -------------------------
@dataclass
class Candidate:
    seqid: str
    header: str
    seq: str
    length: int
    fields: Dict[str, str]
    refseq: bool
    complete_like: bool
    strain_ok: bool
    tax_label: str


def choose_tax_label(policy: str, members: List[Candidate]) -> Tuple[str, bool, List[str]]:
    labels = [c.tax_label for c in members if c.tax_label and c.tax_label != "NA"]
    if not labels:
        return "NA", False, []

    counts = Counter(labels)
    most_common = counts.most_common()
    top_label, top_n = most_common[0]
    total = sum(counts.values())

    # majority: >50%
    if policy == "majority":
        if top_n > total / 2:
            return top_label, False, sorted(set(labels))
        # no majority -> fall back to refseq_wins then longest
        # We'll handle in rep selection; still report mixed
        return top_label, True, sorted(set(labels))

    # refseq_wins: prefer label among refseq members if it exists, else most_common
    if policy == "refseq_wins":
        refseq_labels = [c.tax_label for c in members if c.refseq and c.tax_label != "NA"]
        if refseq_labels:
            rc = Counter(refseq_labels).most_common()
            return rc[0][0], len(set(refseq_labels)) > 1, sorted(set(labels))
        return top_label, len(set(labels)) > 1, sorted(set(labels))

    # longest: just report most common label, mixed if >1
    return top_label, len(set(labels)) > 1, sorted(set(labels))


def phase2_pick_rep(
    cluster: Cluster,
    cand_map: Dict[str, Candidate],
    refseq_min_len_frac: float,
    policy: str,
) -> Tuple[str, str]:
    """
    Returns: (chosen_seqid, reason)
    """
    members: List[Candidate] = [cand_map[sid] for sid in cluster.members if sid in cand_map]
    if not members:
        return cluster.cdhit_rep, "fallback_no_members"

    # CD-HIT rep length reference
    cdhit_rep_len = cand_map[cluster.cdhit_rep].length if cluster.cdhit_rep in cand_map else max(c.length for c in members)
    min_len_for_refseq = int(refseq_min_len_frac * cdhit_rep_len)

    chosen_label, mixed_flag, label_set = choose_tax_label(policy, members)

    # Candidate ranking
    def score(c: Candidate) -> Tuple[int, int, int, int, int]:
        # Higher better.
        # 1) tax label agreement (policy-driven) — prefer those matching chosen_label when not NA
        tax_ok = 1 if (chosen_label != "NA" and c.tax_label == chosen_label) else 0
        # 2) refseq preference but only if long enough
        ref_ok = 1 if (c.refseq and c.length >= min_len_for_refseq) else 0
        # 3) complete-like
        comp_ok = 1 if c.complete_like else 0
        # 4) strain
        strain_ok = 1 if c.strain_ok else 0
        # 5) length (tie-break)
        return (tax_ok, ref_ok, comp_ok, strain_ok, c.length)

    members_sorted = sorted(members, key=score, reverse=True)
    best = members_sorted[0]

    reason_bits = []
    reason_bits.append(f"policy={policy}")
    reason_bits.append(f"chosen_label={chosen_label}")
    if mixed_flag:
        reason_bits.append("mixed_tax_labels=1")
    if best.refseq and best.length >= min_len_for_refseq:
        reason_bits.append("refseq=1")
    if best.complete_like:
        reason_bits.append("complete_like=1")
    if best.strain_ok:
        reason_bits.append("strain=1")
    reason_bits.append(f"len={best.length}")
    reason = ";".join(reason_bits)

    return best.seqid, reason


# -------------------------
# Main
# -------------------------
def find_threshold_dirs(reps_all_dir: Path) -> List[Tuple[str, Path]]:
    # expects: reps_all/all_groups/ALL/c0.990000 etc
    base = reps_all_dir / "all_groups" / "ALL"
    if not base.exists():
        raise SystemExit(f"[ERROR] Expected folder not found: {base}")

    out: List[Tuple[str, Path]] = []
    for p in sorted(base.glob("c*")):
        if p.is_dir():
            label = p.name  # e.g. c0.990000
            out.append((label, p))
    if not out:
        raise SystemExit(f"[ERROR] No threshold directories found under: {base}")
    return out


def main():
    ap = argparse.ArgumentParser(description="Phase 2 representative selection from CD-HIT clusters (reps_all).")
    ap.add_argument("--reps_all_dir", required=True, help="Path to reps_all directory.")
    ap.add_argument("--policy", choices=["majority", "refseq_wins", "longest"], default="majority",
                    help="Cluster taxonomy policy (default majority).")
    ap.add_argument("--refseq_min_len_frac", type=float, default=0.90,
                    help="RefSeq candidate must be >= this fraction of CD-HIT rep length (default 0.90).")
    ap.add_argument("--outdir", default=None,
                    help="Output directory for phase2 files. Default: parent of reps_all_dir.")
    ap.add_argument("--prefix",default="phase2",
                    help="Prefix for output files (default: phase2). E.g. viralDB_phase2 → viralDB_phase2_representatives__all__c*.fasta")    
    ap.add_argument("--verbose", action="store_true", help="Verbose logging.")
    args = ap.parse_args()

    reps_all_dir = Path(args.reps_all_dir).resolve()
    if not reps_all_dir.exists():
        raise SystemExit(f"[ERROR] reps_all_dir not found: {reps_all_dir}")

    outdir = Path(args.outdir).resolve() if args.outdir else reps_all_dir.parent
    outdir.mkdir(parents=True, exist_ok=True)

    thr_dirs = find_threshold_dirs(reps_all_dir)

    for thr_label, thr_dir in thr_dirs:
        input_fa = thr_dir / "input.fasta"
        clstr = thr_dir / "reps.fasta.clstr"

        if not input_fa.exists():
            print(f"[WARN] Skipping {thr_label}: missing {input_fa}", file=sys.stderr)
            continue
        if not clstr.exists():
            print(f"[WARN] Skipping {thr_label}: missing {clstr}", file=sys.stderr)
            continue

        if args.verbose:
            print(f"[INFO] Processing {thr_label}: {thr_dir}", flush=True)

        # Load sequences from input.fasta
        cand_map: Dict[str, Candidate] = {}
        for h, s in fasta_iter(input_fa):
            sid = header_id(h)
            fields = parse_fields(h)
            seq = s.upper()
            L = len(seq)
            cand_map[sid] = Candidate(
                seqid=sid,
                header=h,
                seq=seq,
                length=L,
                fields=fields,
                refseq=is_refseq(fields, sid),
                complete_like=is_complete_like(fields),
                strain_ok=has_strain(fields),
                tax_label=get_tax_label(fields),
            )

        clusters = parse_cdhit_clstr(clstr)

        # Outputs
        prefix = args.prefix
        out_fa = outdir / f"{prefix}_representatives__all__{thr_label}.fasta"
        out_tsv = outdir / f"{prefix}_summary__all__{thr_label}.tsv"        

        # Select reps and write
        selected_ids: List[str] = []
        with out_tsv.open("w", encoding="utf-8") as tsv:
            tsv.write(
                "threshold\tcluster_id\tcdhit_rep\tphase2_rep\tphase2_reason\t"
                "cluster_size\tcdhit_rep_len\tphase2_rep_len\t"
                "species_set_size\tspecies_set\tmixed_species_flag\n"
            )

            for cl in clusters:
                members = [sid for sid in cl.members if sid in cand_map]
                if not members:
                    continue

                chosen_id, reason = phase2_pick_rep(
                    cluster=cl,
                    cand_map=cand_map,
                    refseq_min_len_frac=args.refseq_min_len_frac,
                    policy=args.policy,
                )

                # cluster label stats based on species field (not genus fallback)
                species_labels = []
                for sid in members:
                    sp = (cand_map[sid].fields.get("species", "NA") or "NA")
                    if sp != "NA":
                        species_labels.append(sp)
                species_set = sorted(set(species_labels))
                mixed_flag = 1 if len(species_set) > 1 else 0

                cdhit_len = cand_map[cl.cdhit_rep].length if cl.cdhit_rep in cand_map else cand_map[members[0]].length
                chosen_len = cand_map[chosen_id].length if chosen_id in cand_map else 0

                tsv.write(
                    f"{thr_label}\t{cl.cid}\t{cl.cdhit_rep}\t{chosen_id}\t{reason}\t"
                    f"{len(members)}\t{cdhit_len}\t{chosen_len}\t"
                    f"{len(species_set)}\t{';'.join(species_set) if species_set else 'NA'}\t{mixed_flag}\n"
                )

                selected_ids.append(chosen_id)

        # Write FASTA of selected reps (de-dup while preserving order)
        seen = set()
        with out_fa.open("w", encoding="utf-8") as fa:
            for sid in selected_ids:
                if sid in seen:
                    continue
                seen.add(sid)
                c = cand_map.get(sid)
                if not c:
                    continue
                fa.write(c.header + "\n")
                fa.write(c.seq + "\n")

        print(f"[OK] {thr_label}: phase2 reps FASTA -> {out_fa}")
        print(f"[OK] {thr_label}: phase2 summary TSV -> {out_tsv}")

    print("[DONE] Phase2 representative selection complete.")


if __name__ == "__main__":
    main()
