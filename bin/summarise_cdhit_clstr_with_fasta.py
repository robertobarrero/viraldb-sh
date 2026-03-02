#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


# CD-HIT member line has: "... >ID... *"  or "... >ID... at +/99.00%"
CDHIT_ID_RE = re.compile(r">\s*([^,\s]+)")
# Your enriched FASTA uses pipe-delimited key=value fields after accession
PAIR_RE = re.compile(r"=", 1)


def fasta_iter(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def parse_header_kv(header: str, seq_len: int) -> Tuple[str, Dict[str, str]]:
    """
    Header example:
      OQ993360|species=Pepino_mosaic_virus|...|lineage=k__...;...;s__...
    Returns accession (no version) + dict fields.
    """
    parts = header.split("|")
    acc = parts[0].strip()
    acc_nover = acc.split(".", 1)[0]  # drop .1 if present

    meta: Dict[str, str] = {}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            meta[k.strip()] = v.strip()

    meta.setdefault("length", str(seq_len))
    return acc_nover, meta


def load_fasta_metadata(fasta: Path) -> Dict[str, Dict[str, str]]:
    md: Dict[str, Dict[str, str]] = {}
    for h, seq in fasta_iter(fasta):
        acc, meta = parse_header_kv(h, len(seq))
        md[acc] = meta
    return md


def safe(x: Optional[str]) -> str:
    v = (x or "").strip()
    return v if v else "NA"


def lineage_get(lineage: str, prefix: str) -> str:
    if not lineage or lineage == "NA":
        return "NA"
    for chunk in lineage.split(";"):
        chunk = chunk.strip()
        if chunk.startswith(prefix):
            return safe(chunk[len(prefix):])
    return "NA"


def parse_cdhit_clstr(clstr: Path) -> List[Tuple[int, str, bool]]:
    """
    Returns list of (cluster_id, member_id, is_representative).
    member_id is accession (no version) as best-effort.
    """
    out: List[Tuple[int, str, bool]] = []
    cluster_id: Optional[int] = None

    with clstr.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">Cluster"):
                # >Cluster 0
                toks = line.split()
                cluster_id = int(toks[1])
                continue

            if cluster_id is None:
                continue

            m = CDHIT_ID_RE.search(line)
            if not m:
                continue

            ident = m.group(1).strip()
            # ident can be "NC_031463|species=C..." in clstr (truncated after that)
            # We only need the first token before '|' and before possible '.1'
            ident = ident.split("|", 1)[0]
            ident = ident.split(".", 1)[0]

            is_rep = line.endswith("*")
            out.append((cluster_id, ident, is_rep))

    return out


def main():
    ap = argparse.ArgumentParser(
        description="Summarise CD-HIT .clstr by looking up full metadata in enriched FASTA (avoids truncated headers)."
    )
    ap.add_argument("--clstr", required=True, help="CD-HIT output .clstr file")
    ap.add_argument("--fasta", required=True, help="Enriched FASTA used as CD-HIT input (full headers)")
    ap.add_argument("--out_summary", required=True, help="Output TSV: one row per cluster")
    ap.add_argument("--out_mixed", required=True, help="Output TSV: members for clusters mixing species labels")
    ap.add_argument("--out_flat", default=None, help="Optional TSV: identity_label/group/rep/member flattened table")
    ap.add_argument("--identity_label", default="c1000", help="Label to write into outputs (e.g. c1000, c995, c990)")
    ap.add_argument("--group", default="ALL", help="Group label to write into outputs (e.g. species/genus/ALL)")
    ap.add_argument("--treat_na_as_label", action="store_true",
                    help="If set, NA counts as a species label when detecting mixed clusters (default: NA ignored).")
    args = ap.parse_args()

    clstr = Path(args.clstr)
    fasta = Path(args.fasta)
    out_summary = Path(args.out_summary)
    out_mixed = Path(args.out_mixed)
    out_flat = Path(args.out_flat) if args.out_flat else None

    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_mixed.parent.mkdir(parents=True, exist_ok=True)
    if out_flat:
        out_flat.parent.mkdir(parents=True, exist_ok=True)

    md = load_fasta_metadata(fasta)

    triples = parse_cdhit_clstr(clstr)
    if not triples:
        raise SystemExit(f"[ERROR] No members parsed from: {clstr}")

    # cluster_id -> list of members
    members_by_cluster: Dict[int, List[str]] = defaultdict(list)
    rep_by_cluster: Dict[int, str] = {}

    for cid, mem, is_rep in triples:
        members_by_cluster[cid].append(mem)
        if is_rep:
            rep_by_cluster[cid] = mem

    # Some clstrs might not have '*' for weird cases; fallback to first member
    for cid, mems in members_by_cluster.items():
        if cid not in rep_by_cluster and mems:
            rep_by_cluster[cid] = mems[0]

    # write outputs
    with out_summary.open("w", encoding="utf-8") as s_out, out_mixed.open("w", encoding="utf-8") as m_out:
        s_out.write("\t".join([
            "identity_label", "group", "cluster_id",
            "representative", "n_members",
            "rep_species", "rep_genus", "rep_family", "rep_sequence", "rep_length",
            "n_species_labels", "species_labels", "mixed_species",
            "n_missing_metadata"
        ]) + "\n")

        m_out.write("\t".join([
            "identity_label", "group", "cluster_id", "representative", "member",
            "member_species", "member_genus", "member_family",
            "member_sequence", "member_length", "member_lineage"
        ]) + "\n")

        if out_flat:
            with out_flat.open("w", encoding="utf-8") as f_out:
                f_out.write("\t".join(["identity_label", "group", "representative", "member"]) + "\n")
                for cid, mems in members_by_cluster.items():
                    rep = rep_by_cluster[cid]
                    for mem in mems:
                        f_out.write(f"{args.identity_label}\t{args.group}\t{rep}\t{mem}\n")

        for cid in sorted(members_by_cluster.keys()):
            rep = rep_by_cluster[cid]
            mems = members_by_cluster[cid]

            rep_meta = md.get(rep, {})
            rep_lineage = safe(rep_meta.get("lineage", "NA"))
            rep_species = safe(rep_meta.get("species", "NA"))
            rep_genus = lineage_get(rep_lineage, "g__")
            rep_family = lineage_get(rep_lineage, "f__")
            rep_seqtype = safe(rep_meta.get("sequence", "NA"))
            rep_len = safe(rep_meta.get("length", "NA"))

            species_set = set()
            missing = 0
            per_member_rows = []

            for mem in mems:
                meta = md.get(mem)
                if not meta:
                    missing += 1
                    lin = "NA"
                    sp = "NA"
                    g = "NA"
                    f = "NA"
                    seqtype = "NA"
                    length = "NA"
                else:
                    lin = safe(meta.get("lineage", "NA"))
                    sp = safe(meta.get("species", "NA"))
                    g = lineage_get(lin, "g__")
                    f = lineage_get(lin, "f__")
                    seqtype = safe(meta.get("sequence", "NA"))
                    length = safe(meta.get("length", "NA"))

                if args.treat_na_as_label or sp != "NA":
                    species_set.add(sp)

                per_member_rows.append((mem, sp, g, f, seqtype, length, lin))

            n_species = len(species_set)
            mixed = "YES" if n_species > 1 else "NO"
            species_labels = ",".join(sorted(species_set)) if species_set else "NA"

            if mixed == "YES":
                for mem, sp, g, f, seqtype, length, lin in per_member_rows:
                    m_out.write("\t".join([
                        args.identity_label, args.group, str(cid), rep, mem,
                        sp, g, f, seqtype, length, lin
                    ]) + "\n")

            s_out.write("\t".join([
                args.identity_label, args.group, str(cid),
                rep, str(len(mems)),
                rep_species, rep_genus, rep_family, rep_seqtype, rep_len,
                str(n_species), species_labels, mixed,
                str(missing)
            ]) + "\n")

    print(f"[OK] Summary: {out_summary}")
    print(f"[OK] Mixed:   {out_mixed}")
    if out_flat:
        print(f"[OK] Flat:    {out_flat}")


if __name__ == "__main__":
    main()
