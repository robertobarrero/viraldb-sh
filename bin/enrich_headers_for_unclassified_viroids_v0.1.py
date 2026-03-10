#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


FIELDS_ORDER = [
    "species",
    "virus_name",
    "strain",
    "genus",
    "family",
    "source",
    "sequence",
    "lineage",
    "host",
    "country",
    "collectionDate",
    "releaseDate",
    "length",
    "taxId",
]

BAD_CHARS = re.compile(r"[|=]")


def sanitize(val: str) -> str:
    val = str(val or "").strip()
    if not val:
        return "NA"
    val = val.replace(" ", "_")
    val = BAD_CHARS.sub("_", val)
    return val


def clean_date(val: str) -> str:
    """
    Convert e.g. 2023-12-11T00:00:00Z -> 2023-12-11
    Keep partial dates unchanged.
    """
    val = str(val or "").strip()
    if not val:
        return "NA"

    m = re.match(r"^(\d{4}-\d{2}-\d{2})T", val)
    if m:
        return m.group(1)

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
    Examples:
      >OQ993360.1 ...
      >NC_028271 ...
    """
    h = raw_header[1:].strip()
    tok = h.split()[0].strip("|")

    if re.match(r"^[A-Z]{1,4}_?\d+\.\d+$", tok):
        return tok.split(".")[0], tok
    if re.match(r"^[A-Z]{1,4}_?\d+$", tok):
        return tok, "NA"

    if "." in tok and tok.split(".")[-1].isdigit():
        return tok.split(".")[0], tok

    return tok, "NA"


def normalize_sequence_type(completeness: Optional[str], segment: Optional[str]) -> str:
    comp = (completeness or "").strip().upper()
    seg = (segment or "").strip()

    if comp == "COMPLETE":
        return "complete_segment" if seg else "complete_genome"
    if comp == "PARTIAL":
        return "partial"

    return "NA"


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
    Example:
      >KX013550.1 Citrus viroid VII clone LD3, complete genome
    -> species=Citrus viroid VII
       strain=LD3
       sequence=complete_genome
    """
    text = raw_header[1:].strip()
    parts = text.split(maxsplit=1)
    desc = parts[1] if len(parts) > 1 else ""
    seq_type = infer_sequence_from_defline(desc)

    desc_clean = desc.replace(";", ",")
    low = desc_clean.lower()

    strain = "NA"
    for key in [" clone ", " isolate ", " strain ", " variant "]:
        if key in low:
            idx = low.index(key)
            after = desc_clean[idx + len(key):]
            strain = after.split(",")[0].strip()
            species = desc_clean[:idx].split(",")[0].strip()
            return species, strain, seq_type, desc

    species = desc_clean.split(",")[0].strip() if desc_clean else "NA"
    return species, "NA", seq_type, desc


def load_datasets_jsonl(report_path: Path) -> Dict[str, dict]:
    """
    Key by accession version and accession without version.
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
            if isinstance(acc, str) and "." in acc and acc.split(".")[-1].isdigit():
                m[acc.split(".")[0]] = rec
    return m


def extract_taxid(rec: dict) -> str:
    """
    Prefer virus.taxId, else last lineage taxId.
    """
    virus = rec.get("virus", {})
    if isinstance(virus, dict):
        taxid = virus.get("taxId")
        if taxid not in (None, "", "NA"):
            return str(taxid)

        lineage = virus.get("lineage", [])
        if isinstance(lineage, list) and lineage:
            last = lineage[-1]
            if isinstance(last, dict):
                last_taxid = last.get("taxId")
                if last_taxid not in (None, "", "NA"):
                    return str(last_taxid)

    return "NA"


def extract_collection_date(rec: dict) -> str:
    val = rec.get("collectionDate")
    if val:
        return clean_date(val)

    biosample = rec.get("biosample", {})
    if isinstance(biosample, dict):
        val = biosample.get("collectionDate")
        if val:
            return clean_date(val)

    return "NA"


def extract_release_date(rec: dict) -> str:
    return clean_date(rec.get("releaseDate", "NA"))


def lineage_list_to_ranks(lineage: List[dict], organism_name: Optional[str]) -> Tuple[str, str, str, str, str, str]:
    """
    Convert datasets lineage to a virus-style lineage:
      k__...;p__...;c__...;o__...;f__...;g__...;s__...;virusName__...

    For viroids / viroid-like RNAs, taxonomy can be sparse or unusual.
    Practical rules:
      - k__Viruses if present else NA
      - p__/c__/o__ via suffix heuristics where available
      - f__ from last 'viridae' or 'viroidae'
      - g__ from the entry immediately before organism_name, if informative
      - s__:
          * if entry immediately before organism_name exists and looks species-like, use it
          * else use organism_name
      - virusName__ always organism_name
    """
    names = [x.get("name", "") for x in lineage if isinstance(x, dict)]
    names = [n.strip() for n in names if n and isinstance(n, str)]

    k = "Viruses" if any(n == "Viruses" for n in names) else "NA"

    p = "NA"
    c = "NA"
    o = "NA"
    f = "NA"
    g = "NA"
    s = "NA"
    virus_name = organism_name.strip() if organism_name else "NA"

    phyla = [n for n in names if n.lower().endswith("viricota")]
    if phyla:
        p = phyla[-1]

    classes = [n for n in names if n.lower().endswith("cetes")]
    if classes:
        c = classes[-1]

    orders = [n for n in names if n.lower().endswith("virales")]
    if orders:
        o = orders[-1]

    families = [n for n in names if n.lower().endswith("viridae") or n.lower().endswith("viroidae")]
    if families:
        f = families[-1]

    if organism_name and organism_name in names:
        idx = names.index(organism_name)

        if idx > 0:
            prev = names[idx - 1]

            # Do not treat higher ranks as genus
            prev_low = prev.lower()
            is_higher = (
                prev == "Viruses"
                or prev_low.endswith("viridae")
                or prev_low.endswith("viroidae")
                or prev_low.endswith("virales")
                or prev_low.endswith("viricota")
                or prev_low.endswith("cetes")
            )

            if not is_higher:
                g = prev

        # Species rank:
        # If there is an informative node before organism_name and it is NOT the same as family/order/etc,
        # use it as s__. Otherwise use organism_name.
        if idx > 0:
            prev = names[idx - 1]
            prev_low = prev.lower()
            is_higher = (
                prev == "Viruses"
                or prev_low.endswith("viridae")
                or prev_low.endswith("viroidae")
                or prev_low.endswith("virales")
                or prev_low.endswith("viricota")
                or prev_low.endswith("cetes")
            )
            if not is_higher:
                s = prev
            else:
                s = organism_name
        else:
            s = organism_name
    else:
        s = organism_name if organism_name else "NA"

    lineage_str = (
        f"k__{sanitize(k)};"
        f"p__{sanitize(p)};"
        f"c__{sanitize(c)};"
        f"o__{sanitize(o)};"
        f"f__{sanitize(f)};"
        f"g__{sanitize(g)};"
        f"s__{sanitize(s)};"
        f"virusName__{sanitize(virus_name)}"
    )

    return lineage_str, g, f, s, virus_name, k


def build_header(accession_nover: str, meta: Dict[str, str]) -> str:
    return ">" + accession_nover + "\t|" + "|".join(
        f"{k}={sanitize(meta.get(k, 'NA'))}" for k in FIELDS_ORDER
    )


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Enrich viroid FASTA headers so they match the NCBI virus header structure, "
            "using datasets JSONL plus FASTA defline fallback."
        )
    )
    ap.add_argument("--fasta", required=True, help="Input FASTA")
    ap.add_argument("--datasets_report", required=True, help="Input JSONL")
    ap.add_argument("--out_fasta", required=True, help="Output enriched FASTA")
    ap.add_argument("--out_metadata", required=True, help="Output metadata TSV")
    ap.add_argument("--source", default="ViroidDB", help="Header source=")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    report_path = Path(args.datasets_report)

    ds = load_datasets_jsonl(report_path)

    out_fa = Path(args.out_fasta)
    out_md = Path(args.out_metadata)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)

    with out_fa.open("w", encoding="utf-8") as fa_out, out_md.open("w", encoding="utf-8") as md_out:
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

            raw_text = raw_h[1:].strip()
            parts = raw_text.split(maxsplit=1)
            original_defline = parts[1] if len(parts) > 1 else ""

            meta = {
                "species": "NA",
                "virus_name": "NA",
                "strain": "NA",
                "genus": "NA",
                "family": "NA",
                "source": args.source,
                "sequence": "NA",
                "lineage": "NA",
                "host": "NA",
                "country": "NA",
                "collectionDate": "NA",
                "releaseDate": "NA",
                "length": str(len(seq)),
                "taxId": "NA",
            }

            org = "NA"
            srcdb = "NA"

            if rec:
                srcdb = rec.get("sourceDatabase", "NA") or "NA"

                virus = rec.get("virus", {}) if isinstance(rec.get("virus", {}), dict) else {}
                org = virus.get("organismName", "NA") or "NA"

                if org != "NA":
                    meta["species"] = org
                    meta["virus_name"] = org

                meta["taxId"] = extract_taxid(rec)

                iso = rec.get("isolate", {}) if isinstance(rec.get("isolate", {}), dict) else {}
                meta["strain"] = iso.get("name", "NA") or meta["strain"]

                host = rec.get("host")
                if isinstance(host, dict):
                    meta["host"] = host.get("organismName", "NA") or meta["host"]
                elif isinstance(host, str):
                    meta["host"] = host

                loc = rec.get("location", {}) if isinstance(rec.get("location", {}), dict) else {}
                meta["country"] = loc.get("geographicLocation", "NA") or meta["country"]

                seqtype = normalize_sequence_type(rec.get("completeness"), rec.get("segment"))
                if seqtype != "NA":
                    meta["sequence"] = seqtype

                meta["collectionDate"] = extract_collection_date(rec)
                meta["releaseDate"] = extract_release_date(rec)

                lineage = virus.get("lineage", [])
                if isinstance(lineage, list) and lineage:
                    lineage_str, genus, family, species_rank, virus_name, _k = lineage_list_to_ranks(lineage, org)
                    meta["lineage"] = lineage_str
                    meta["genus"] = genus
                    meta["family"] = family
                    if meta["virus_name"] == "NA":
                        meta["virus_name"] = virus_name

            # FASTA defline fallback
            if (not rec) or meta["species"] == "NA" or meta["strain"] == "NA" or meta["sequence"] == "NA":
                sp, st, seqtype2, defline2 = parse_species_strain_from_defline(raw_h)

                if meta["species"] == "NA" and sp != "NA":
                    meta["species"] = sp
                if meta["virus_name"] == "NA" and sp != "NA":
                    meta["virus_name"] = sp
                if meta["strain"] == "NA" and st != "NA":
                    meta["strain"] = st
                if meta["sequence"] == "NA" and seqtype2 != "NA":
                    meta["sequence"] = seqtype2
                if not original_defline and defline2:
                    original_defline = defline2

            # If lineage is still NA, create a minimal lineage with virusName
            if meta["lineage"] == "NA":
                meta["lineage"] = (
                    f"k__Viruses;"
                    f"p__NA;"
                    f"c__NA;"
                    f"o__NA;"
                    f"f__{sanitize(meta['family'])};"
                    f"g__{sanitize(meta['genus'])};"
                    f"s__{sanitize(meta['species'])};"
                    f"virusName__{sanitize(meta['virus_name'])}"
                )

            fa_out.write(build_header(acc_nover, meta) + "\n")
            fa_out.write(seq.upper() + "\n")

            md_out.write(
                f"{acc_nover}\t{acc_ver}\t{sanitize(meta['taxId'])}\t{sanitize(org)}\t{sanitize(srcdb)}\t{sanitize(original_defline)}\t"
                + "\t".join(sanitize(meta[k]) for k in FIELDS_ORDER) + "\n"
            )

    print(f"[OK] Enriched FASTA: {out_fa}")
    print(f"[OK] Metadata TSV:  {out_md}")


if __name__ == "__main__":
    main()
