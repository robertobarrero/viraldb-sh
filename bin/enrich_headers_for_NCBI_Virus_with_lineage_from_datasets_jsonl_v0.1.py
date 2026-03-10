#!/usr/bin/env python3

from __future__ import annotations
import argparse
import json
import re
from pathlib import Path
from typing import Iterable, List, Tuple


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


def sanitize(val) -> str:
    val = str(val or "").strip()
    if not val:
        return "NA"
    val = val.replace(" ", "_")
    val = BAD_CHARS.sub("_", val)
    return val


def clean_date(val) -> str:
    """
    Convert values like 2023-12-11T00:00:00Z -> 2023-12-11
    Keep partial dates such as 2023 or 2023-08 unchanged.
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

    with path.open() as fh:
        for line in fh:
            line = line.rstrip()

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)

                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())

        if header is not None:
            yield header, "".join(seq_parts)


def parse_accession_and_version(header: str):
    tok = header[1:].split()[0]
    if "." in tok:
        return tok.split(".")[0], tok
    return tok, "NA"


def load_jsonl(path: Path):
    """
    Load JSONL records keyed by accession with and without version.
    """
    records = {}

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            rec = json.loads(line)
            acc = rec.get("accession")
            if not acc:
                continue

            records[acc] = rec
            records[acc.split(".")[0]] = rec

    return records


def normalize_sequence_type(completeness, segment):
    comp = (completeness or "").upper()

    if comp == "COMPLETE":
        if segment:
            return "complete_segment"
        return "complete_genome"

    if comp == "PARTIAL":
        return "partial"

    return "NA"


def extract_taxid(rec):
    """
    Prefer virus.taxId.
    If missing, fall back to the last lineage element taxId.
    """
    virus = rec.get("virus", {})

    taxid = virus.get("taxId")
    if taxid not in (None, ""):
        return str(taxid)

    lineage = virus.get("lineage", [])
    if isinstance(lineage, list) and lineage:
        last_item = lineage[-1]
        if isinstance(last_item, dict):
            last_taxid = last_item.get("taxId")
            if last_taxid not in (None, ""):
                return str(last_taxid)

    return "NA"

def extract_collection_date(rec):
    """
    Try common locations for collection date.
    Priority:
      1. top-level collectionDate
      2. isolate.collectionDate
      3. biosample.collectionDate
    """
    val = rec.get("collectionDate")
    if val not in (None, "", "NA"):
        return clean_date(val)

    isolate = rec.get("isolate", {})
    if isinstance(isolate, dict):
        val = isolate.get("collectionDate")
        if val not in (None, "", "NA"):
            return clean_date(val)

    biosample = rec.get("biosample", {})
    if isinstance(biosample, dict):
        val = biosample.get("collectionDate")
        if val not in (None, "", "NA"):
            return clean_date(val)

    return "NA"

def is_higher_rank_name(name: str, family: str, order: str, phylum: str, class_name: str) -> bool:
    if not name:
        return True

    higher = {
        "Viruses",
        "Riboviria",
        family,
        order,
        phylum,
        class_name,
    }

    if name in higher:
        return True

    nl = name.lower()

    if nl.endswith("viridae"):
        return True
    if nl.endswith("virales"):
        return True
    if nl.endswith("viricota"):
        return True
    if nl.endswith("cetes"):
        return True

    return False


def is_unclassified_name(name: str) -> bool:
    if not name or name == "NA":
        return False
    return name.lower().startswith("unclassified_")


def looks_like_species_for_genus(name: str, genus: str) -> bool:
    """
    Heuristics for a taxonomic species node.
    Examples:
      Potexvirus citriflavivenae
      Emaravirus actinidiae
      Potexvirus sp.
    """
    if not name or name == "NA":
        return False

    name_clean = name.replace("_", " ").strip().lower()

    if re.search(r"\bsp\.$", name_clean, flags=re.IGNORECASE):
        return True

    if genus and genus != "NA":
        genus_clean = genus.replace("_", " ").strip().lower()
        if name_clean.startswith(genus_clean + " "):
            return True

    return False


def get_last_index(names: List[str], organism: str):
    idx = None
    for i, n in enumerate(names):
        if n == organism:
            idx = i
    return idx


def infer_genus_species_and_virus_name(lineage_names: List[str], organism: str, family: str, order: str, phylum: str, class_name: str):
    """
    Rules implemented to match the requested behavior:

    1. genus + taxonomic species + virus/common name
       Potexvirus -> Potexvirus citriflavivenae -> Citrus yellow vein clearing virus
       => g__Potexvirus ; s__Potexvirus_citriflavivenae ; virusName__Citrus_yellow_vein_clearing_virus

    2. only virus/common name available
       ... -> Alternanthera mosaic virus
       => g__NA ; s__NA ; virusName__Alternanthera_mosaic_virus

    3. unclassified genus + species placeholder
       ... -> unclassified Potexvirus -> Potexvirus sp.
       => keep genus, keep species, virusName same as organism

    4. unclassified genus/family + species-like organism
       ... -> unclassified Alphaflexiviridae -> Pistachio potex-like virus
       => keep genus, keep species same as organism, virusName same as organism
    """
    genus = "NA"
    species_rank = "NA"
    virus_name = organism if organism else "NA"

    org_idx = get_last_index(lineage_names, organism)

    if org_idx is None:
        return genus, species_rank, virus_name

    prev1 = lineage_names[org_idx - 1] if org_idx - 1 >= 0 else None
    prev2 = lineage_names[org_idx - 2] if org_idx - 2 >= 0 else None

    prev1_is_higher = is_higher_rank_name(prev1, family, order, phylum, class_name) if prev1 else True
    prev2_is_higher = is_higher_rank_name(prev2, family, order, phylum, class_name) if prev2 else True

    # Case A: genus + species + virus name
    if prev1 and prev2 and not prev1_is_higher and not prev2_is_higher:
        if looks_like_species_for_genus(prev1, prev2):
            genus = prev2
            species_rank = prev1
            return genus, species_rank, virus_name

    # Case B: single informative predecessor
    if prev1 and not prev1_is_higher:
        # unclassified genus case: keep genus and decide species
        if is_unclassified_name(prev1):
            genus = prev1

            # If organism is already species-like or placeholder-like, keep it as species
            if organism and (
                re.search(r"\bsp\.$", organism.replace("_", " "), flags=re.IGNORECASE)
                or "virus" in organism.lower()
            ):
                species_rank = organism
            else:
                species_rank = "NA"

            return genus, species_rank, virus_name

        # ordinary genus only before organism
        genus = prev1

        # if organism looks like a species attached to that genus, use it
        if looks_like_species_for_genus(organism, genus):
            species_rank = organism
        else:
            species_rank = "NA"

        return genus, species_rank, virus_name

    # Case C: nothing informative before organism -> only virusName available
    genus = "NA"
    species_rank = "NA"
    return genus, species_rank, virus_name


def lineage_to_string_and_ranks(lineage, organism):
    names = [x.get("name") for x in lineage if isinstance(x, dict) and x.get("name")]

    family = "NA"
    order = "NA"
    phylum = "NA"
    class_name = "NA"

    for n in names:
        nl = n.lower()

        if nl.endswith("viridae"):
            family = n

        if nl.endswith("virales"):
            order = n

        if nl.endswith("viricota"):
            phylum = n

        if nl.endswith("cetes"):
            class_name = n

    genus, species_rank, virus_name = infer_genus_species_and_virus_name(
        lineage_names=names,
        organism=organism,
        family=family,
        order=order,
        phylum=phylum,
        class_name=class_name,
    )

    lineage_str = (
        f"k__Viruses;"
        f"p__{sanitize(phylum)};"
        f"c__{sanitize(class_name)};"
        f"o__{sanitize(order)};"
        f"f__{sanitize(family)};"
        f"g__{sanitize(genus)};"
        f"s__{sanitize(species_rank)};"
        f"virusName__{sanitize(virus_name)}"
    )

    return lineage_str, family, genus, species_rank, virus_name


def build_header(acc, meta):
    return ">" + acc + "\t|" + "|".join(
        f"{k}={sanitize(meta[k])}" for k in FIELDS_ORDER
    )


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Enrich FASTA headers using NCBI datasets JSONL report, including "
            "corrected lineage genus/species/virusName parsing, taxId, and collectionDate."
        )
    )
    ap.add_argument("--fasta", required=True, help="Input FASTA")
    ap.add_argument("--datasets_report", required=True, help="NCBI datasets JSONL report")
    ap.add_argument("--out_fasta", required=True, help="Output enriched FASTA")
    ap.add_argument("--source", default="NCBI_Virus", help="Value to use for source field")
    args = ap.parse_args()

    fasta = Path(args.fasta)
    report = Path(args.datasets_report)

    records = load_jsonl(report)

    with open(args.out_fasta, "w") as out:
        for header, seq in fasta_iter(fasta):
            acc_nover, acc_ver = parse_accession_and_version(header)
            rec = records.get(acc_ver) or records.get(acc_nover)

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

            if rec:
                virus = rec.get("virus", {})
                organism = virus.get("organismName", "NA")
                lineage = virus.get("lineage", [])

                lineage_str, family, genus, species_rank, virus_name = lineage_to_string_and_ranks(
                    lineage, organism
                )

                # Keep species field as organism/common name to match your requested header style
                meta["species"] = species_rank if species_rank != "NA" else organism
                meta["virus_name"] = virus_name
                meta["genus"] = genus
                meta["family"] = family
                meta["lineage"] = lineage_str
                meta["taxId"] = extract_taxid(rec)

                iso = rec.get("isolate", {})
                if isinstance(iso, dict):
                    meta["strain"] = iso.get("name", "NA")

                host = rec.get("host")
                if isinstance(host, dict):
                    meta["host"] = host.get("organismName", "NA")

                loc = rec.get("location", {})
                if isinstance(loc, dict):
                    meta["country"] = loc.get("geographicLocation", "NA")

                meta["sequence"] = normalize_sequence_type(
                    rec.get("completeness"),
                    rec.get("segment")
                )

                meta["collectionDate"] = extract_collection_date(rec)
                meta["releaseDate"] = clean_date(rec.get("releaseDate", "NA"))

            out.write(build_header(acc_nover, meta) + "\n")
            out.write(seq.upper() + "\n")

    print("[OK] Enriched FASTA written to", args.out_fasta)


if __name__ == "__main__":
    main()
