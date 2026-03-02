#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# Fixed order as requested (do not change order unless you intentionally version headers)
FIELDS_ORDER = [
    "species", "strain", "status", "source", "sequence", "lineage",
    "host", "country", "paper", "length", "prev_db"
]

BAD_CHARS = re.compile(r"[|=]")

def sanitize(val: str) -> str:
    val = (val or "").strip()
    if not val:
        return "NA"
    # keep underscores for stability in downstream parsing
    val = val.replace(" ", "_")
    val = BAD_CHARS.sub("_", val)
    return val

def fasta_iter(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
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

def parse_accession_and_version(raw_header: str) -> Tuple[str, str]:
    """
    Input examples:
      >OQ993360.1 Pepino mosaic virus ...
    Returns:
      ('OQ993360', 'OQ993360.1')
    """
    h = raw_header[1:].strip()
    tok = h.split()[0].strip("|")
    if re.match(r"^[A-Z]{1,4}_?\d+\.\d+$", tok):
        return tok.split(".")[0], tok
    if re.match(r"^[A-Z]{1,4}_?\d+$", tok):
        return tok, "NA"
    # fallback
    if "." in tok and tok.split(".")[-1].isdigit():
        return tok.split(".")[0], tok
    return tok, "NA"

def normalize_sequence_type(completeness: Optional[str], segment: Optional[str]) -> str:
    """
    Datasets gives completeness like PARTIAL/COMPLETE plus may have 'segment' field.
    """
    comp = (completeness or "").strip().upper()
    seg = (segment or "").strip()
    if comp == "COMPLETE":
        return "complete_segment" if seg else "complete_genome"
    if comp == "PARTIAL":
        return "partial"
    return "NA"

def lineage_list_to_kpc_ofgs(lineage: List[dict], organism_name: Optional[str]) -> str:
    """
    Convert datasets 'virus.lineage' list of dicts [{name,taxId}, ...] to:
      k__...;p__...;c__...;o__...;f__...;g__...;s__...
    Strategy:
      - k__: always 'Viruses' if present in list else NA
      - f__: last entry whose name endswith 'viridae' (case-insensitive) if present
      - g__: choose the entry immediately before organism_name (often genus)
      - o__/c__/p__: best-effort via suffix heuristics
    """
    names = [x.get("name", "") for x in lineage if isinstance(x, dict)]
    names_clean = [n.strip() for n in names if n and isinstance(n, str)]

    k = "NA"
    if any(n == "Viruses" for n in names_clean):
        k = "Viruses"

    # family heuristic
    f = "NA"
    fams = [n for n in names_clean if n.lower().endswith("viridae")]
    if fams:
        f = fams[-1]

    # species: prefer organism_name; else last lineage entry
    s = organism_name.strip() if organism_name else (names_clean[-1] if names_clean else "NA")

    # genus heuristic: entry just before organism_name in lineage
    g = "NA"
    if organism_name and organism_name in names_clean:
        idx = names_clean.index(organism_name)
        if idx > 0:
            g = names_clean[idx - 1]

    # order heuristic: endswith 'virales'
    o = "NA"
    orders = [n for n in names_clean if n.lower().endswith("virales")]
    if orders:
        o = orders[-1]

    # class heuristic: endswith 'cetes' (common in virus taxonomy)
    c = "NA"
    classes = [n for n in names_clean if n.lower().endswith("cetes")]
    if classes:
        c = classes[-1]

    # phylum heuristic: endswith 'viricota' (common in virus taxonomy)
    p = "NA"
    phyla = [n for n in names_clean if n.lower().endswith("viricota")]
    if phyla:
        p = phyla[-1]

    return (
        f"k__{sanitize(k)};"
        f"p__{sanitize(p)};"
        f"c__{sanitize(c)};"
        f"o__{sanitize(o)};"
        f"f__{sanitize(f)};"
        f"g__{sanitize(g)};"
        f"s__{sanitize(s)}"
    )

def build_header(accession_nover: str, meta: Dict[str, str]) -> str:
    parts = [f">{accession_nover}"]
    for k in FIELDS_ORDER:
        parts.append(f"{k}={sanitize(meta.get(k, 'NA'))}")
    return "|".join(parts)

def load_datasets_jsonl(report_path: Path) -> Dict[str, dict]:
    """
    Key by accession version (e.g., EF599580.1). Also store without version as fallback.
    """
    m: Dict[str, dict] = {}
    with report_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue
            acc = rec.get("accession")
            if not acc:
                continue
            m[acc] = rec
            # also map accession without version for fallback
            if isinstance(acc, str) and "." in acc and acc.split(".")[-1].isdigit():
                m[acc.split(".")[0]] = rec
    return m

def infer_sequence_from_defline(desc: str) -> str:
    t = (desc or "").lower()
    if "complete genome" in t:
        return "complete_genome"
    if "complete sequence" in t:
        return "complete_genome"
    if "complete cds" in t:
        return "CDS"
    if "partial" in t:
        return "partial"
    if "segment" in t and "complete" in t:
        return "complete_segment"
    return "NA"

def parse_species_strain_from_defline(raw_header: str) -> Tuple[str, str, str, str]:
    """
    Parse species + strain + sequence + original_description from a FASTA header
    when JSONL record is missing.

    Example:
      >KX013550.1 Citrus viroid VII clone LD3, complete genome
    -> species=Citrus_viroid_VII
       strain=LD3
       sequence=complete_genome
    """
    text = raw_header[1:].strip()
    # Remove accession token
    parts = text.split(maxsplit=1)
    desc = parts[1] if len(parts) > 1 else ""
    seq_type = infer_sequence_from_defline(desc)

    # Normalize punctuation
    desc_clean = desc.replace(";", ",")
    low = desc_clean.lower()

    # Identify isolate/clone/strain/variant label
    strain = "NA"
    for key in [" clone ", " isolate ", " strain ", " variant "]:
        if key in low:
            idx = low.index(key)
            after = desc_clean[idx + len(key):]
            strain = after.split(",")[0].strip()
            # species is everything before the key (trim to first comma as well)
            species = desc_clean[:idx].split(",")[0].strip()
            return species, strain, seq_type, desc

    # No explicit key: species before comma, strain NA
    species = desc_clean.split(",")[0].strip() if desc_clean else "NA"
    return species, "NA", seq_type, desc

def main():
    ap = argparse.ArgumentParser(
        description="Enrich FASTA headers using merged_data_report.jsonl (virus.lineage) and fallback parsing from FASTA deflines when JSONL is missing."
    )
    ap.add_argument("--fasta", required=True, help="Input FASTA")
    ap.add_argument("--datasets_report", required=True, help="Input JSONL (merged_data_report.jsonl)")
    ap.add_argument("--out_fasta", required=True, help="Output enriched FASTA")
    ap.add_argument("--out_metadata", required=True, help="Output metadata TSV")
    ap.add_argument("--source", default="NCBI_Virus", help="Header source= (default NCBI_Virus)")
    ap.add_argument("--status", default="NA", help="Header status= (default NA)")
    ap.add_argument("--prev_db", default="NA", help="Header prev_db= (default NA)")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    report_path = Path(args.datasets_report)

    ds = load_datasets_jsonl(report_path)

    out_fa = Path(args.out_fasta)
    out_md = Path(args.out_metadata)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)

    with out_fa.open("w", encoding="utf-8") as fa_out, out_md.open("w", encoding="utf-8") as md_out:
        # Add original_defline column so you never lose information
        md_out.write(
            "accession\taccession_version\ttaxid\torganismName\tsourceDatabase\toriginal_defline\t"
            + "\t".join(FIELDS_ORDER) + "\n"
        )

        seen = set()

        for raw_h, seq in fasta_iter(fasta_path):
            acc_nover, acc_ver = parse_accession_and_version(raw_h)
            if not acc_nover:
                continue
            if acc_nover in seen:
                continue
            seen.add(acc_nover)

            rec = ds.get(acc_ver) or ds.get(acc_nover)

            # Extract original defline (everything after accession token in the raw header)
            raw_text = raw_h[1:].strip()
            parts = raw_text.split(maxsplit=1)
            original_defline = parts[1] if len(parts) > 1 else ""

            # Defaults
            meta = {
                "species": "NA",
                "strain": "NA",
                "status": args.status,
                "source": args.source,
                "sequence": "NA",
                "lineage": "NA",
                "host": "NA",
                "country": "NA",
                "paper": "NA",
                "length": str(len(seq)),
                "prev_db": args.prev_db,
            }

            taxid = "NA"
            org = "NA"
            srcdb = "NA"

            if rec:
                srcdb = rec.get("sourceDatabase", "NA") or "NA"

                # virus object
                virus = rec.get("virus", {}) if isinstance(rec.get("virus", {}), dict) else {}
                taxid = virus.get("taxId", "NA") or "NA"
                org = virus.get("organismName", "NA") or "NA"
                if org != "NA":
                    meta["species"] = org

                # isolate name
                iso = rec.get("isolate", {}) if isinstance(rec.get("isolate", {}), dict) else {}
                meta["strain"] = iso.get("name", "NA") or meta["strain"]

                # host name if present (note: host can be nested object with organismName)
                host = rec.get("host")
                if isinstance(host, dict):
                    meta["host"] = host.get("organismName", "NA") or meta["host"]
                elif isinstance(host, str):
                    meta["host"] = host

                # country / location
                loc = rec.get("location", {}) if isinstance(rec.get("location", {}), dict) else {}
                meta["country"] = loc.get("geographicLocation", "NA") or meta["country"]

                # sequence type from completeness + segment
                seqtype = normalize_sequence_type(rec.get("completeness"), rec.get("segment"))
                if seqtype != "NA":
                    meta["sequence"] = seqtype

                # lineage from virus.lineage
                lineage = virus.get("lineage", [])
                if isinstance(lineage, list) and lineage:
                    meta["lineage"] = lineage_list_to_kpc_ofgs(lineage, org)

            # Fallback parsing from FASTA header if JSONL missing (or fields still NA)
            # This retains info like: "Citrus viroid VII clone LD3, complete genome"
            if (not rec) or meta["species"] == "NA" or meta["strain"] == "NA" or meta["sequence"] == "NA":
                sp, st, seqtype2, defline2 = parse_species_strain_from_defline(raw_h)
                if meta["species"] == "NA" and sp != "NA":
                    meta["species"] = sp
                if meta["strain"] == "NA" and st != "NA":
                    meta["strain"] = st
                if meta["sequence"] == "NA" and seqtype2 != "NA":
                    meta["sequence"] = seqtype2
                # if the original defline was empty, fill from parsed defline
                if not original_defline and defline2:
                    original_defline = defline2

            # Write FASTA (one-line sequence)
            fa_out.write(build_header(acc_nover, meta) + "\n")
            fa_out.write(seq.upper() + "\n")

            md_out.write(
                f"{acc_nover}\t{acc_ver}\t{taxid}\t{sanitize(org)}\t{sanitize(srcdb)}\t{sanitize(original_defline)}\t"
                + "\t".join(sanitize(meta[k]) for k in FIELDS_ORDER) + "\n"
            )

    print(f"[OK] Enriched FASTA: {out_fa}")
    print(f"[OK] Metadata TSV:  {out_md}")

if __name__ == "__main__":
    main()
