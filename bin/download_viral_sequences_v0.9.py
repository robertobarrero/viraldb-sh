#!/usr/bin/env python3
"""
download_viral_sequences_v0.9.py

Download viral/viroid sequences from one or more sources, depending on flags.

Sources:
  --refseq_viral : NCBI Datasets virus genome package (RefSeq only)
  --ncbi_viral   : NCBI Datasets virus genome package (GenBank+RefSeq via NCBI Virus)
  --ncbi_nt      : NCBI nucleotide via Entrez Direct (EDirect) using a user query
  --viroid_db    : ViroidDB snapshot (GitHub zip fallback)

Family mode:
  --viral_families <tsv/txt> : one family name per line (comments with # allowed).
                               If provided, downloads refseq/ncbi viral BY FAMILY and
                               skips families that fail after retries.

Key outputs (family mode):
  - Always writes merged files at TOP-LEVEL outdir:
      <outdir>/<LABEL>__ALL_FAMILIES.fasta
      <outdir>/<LABEL>__ALL_FAMILIES.data_report.jsonl  (if --merge_reports)

  - Also keeps per-family files under:
      <outdir>/families/<FAMILY>/<LABEL>.fasta
      <outdir>/families/<FAMILY>/<LABEL>_data_report.jsonl

Notes:
- We do NOT restrict to "complete genome" at download time; Datasets returns both COMPLETE/PARTIAL.
  You can flag/choose representatives later using the data_report.jsonl completeness field.
"""

from __future__ import annotations

import argparse
import os
import random
import shutil
import subprocess
import time
import zipfile
from pathlib import Path
from typing import List, Optional, Tuple


VIRUSES_TAXID = "10239"  # NCBI Taxonomy: Viruses


def log(msg: str, verbose: bool = True) -> None:
    if verbose:
        print(msg, flush=True)


def which_or_die(cmd: str, hint: str = "") -> None:
    if shutil.which(cmd) is None:
        msg = f"[ERROR] Required executable not found on PATH: {cmd}"
        if hint:
            msg += f"\n        Hint: {hint}"
        raise SystemExit(msg)


def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if p.returncode != 0:
        raise SystemExit(
            f"[ERROR] Command failed ({p.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )


def safe_unzip(zip_path: Path, extract_dir: Path) -> None:
    extract_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(extract_dir)


def write_one_line_fasta(in_fa: Path, out_fa: Path) -> None:
    """Convert multi-line sequences to single-line per record; keep headers as-is."""
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with in_fa.open("r", encoding="utf-8", errors="replace") as fin, out_fa.open(
        "w", encoding="utf-8"
    ) as fout:
        header = None
        seq_chunks: List[str] = []
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    fout.write(header + "\n")
                    fout.write("".join(seq_chunks).upper() + "\n")
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            fout.write(header + "\n")
            fout.write("".join(seq_chunks).upper() + "\n")


def concat_fastas(fastas: List[Path], out_fa: Path) -> None:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w", encoding="utf-8") as fout:
        for fa in fastas:
            with fa.open("r", encoding="utf-8", errors="replace") as fin:
                shutil.copyfileobj(fin, fout)


def _pick_largest_fasta(paths: List[Path]) -> Optional[Path]:
    paths = [p for p in paths if p.exists() and p.is_file()]
    if not paths:
        return None
    return sorted(paths, key=lambda p: p.stat().st_size, reverse=True)[0]


def _find_data_report_jsonl(extract_dir: Path) -> Optional[Path]:
    canonical = extract_dir / "ncbi_dataset" / "data" / "data_report.jsonl"
    if canonical.exists():
        return canonical
    hits = list(extract_dir.rglob("data_report.jsonl"))
    return hits[0] if hits else None


def run_datasets_download(
    cmd: List[str],
    out_zip: Path,
    retries: int,
    no_progressbar: bool,
    verbose: bool,
    base_backoffs: Optional[List[int]] = None,
) -> None:
    """
    Run an NCBI datasets download command with retries + exponential backoff.

    - Optional: --no-progressbar to avoid ANSI escape spam in logs
    - Writes to *.partial then renames on success
    - Not reliably resumable; on retry, starts fresh.
    """
    which_or_die(
        "datasets",
        hint="Install NCBI Datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/",
    )

    out_zip = Path(out_zip)
    out_zip.parent.mkdir(parents=True, exist_ok=True)
    tmp_zip = out_zip.with_suffix(out_zip.suffix + ".partial")

    # Always start clean; datasets downloads are not reliably resumable
    if tmp_zip.exists():
        tmp_zip.unlink()

    backoffs = base_backoffs or [30, 60, 120, 240, 480]

    def _ensure_filename_and_progress(original: List[str], filename: Path) -> List[str]:
        cmd2: List[str] = []
        i = 0
        replaced_filename = False
        while i < len(original):
            if original[i] == "--filename" and i + 1 < len(original):
                cmd2.extend(["--filename", str(filename)])
                i += 2
                replaced_filename = True
            else:
                cmd2.append(original[i])
                i += 1
        if not replaced_filename:
            cmd2.extend(["--filename", str(filename)])

        if no_progressbar and "--no-progressbar" not in cmd2:
            cmd2.append("--no-progressbar")
        return cmd2

    for attempt in range(1, retries + 1):
        cmd2 = _ensure_filename_and_progress(cmd, tmp_zip)

        log(f"[INFO] datasets attempt {attempt}/{retries}: {' '.join(cmd2)}", verbose)
        p = subprocess.run(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if p.returncode == 0 and tmp_zip.exists() and tmp_zip.stat().st_size > 0:
            tmp_zip.replace(out_zip)
            log(f"[OK] Downloaded: {out_zip}", verbose)
            return

        msg = (p.stderr or p.stdout or "").strip()
        log(f"[WARN] datasets failed (attempt {attempt}): {msg}", verbose)

        if attempt < retries:
            sleep_s = backoffs[min(attempt - 1, len(backoffs) - 1)]
            sleep_s = int(sleep_s * (0.8 + 0.4 * random.random()))
            log(f"[INFO] Sleeping {sleep_s}s then retrying...", verbose)
            time.sleep(sleep_s)

            if tmp_zip.exists():
                tmp_zip.unlink()

    raise RuntimeError(f"datasets download failed after {retries} attempts: {out_zip}")


def datasets_download_virus_taxon(
    outdir: Path,
    taxon: str,
    refseq_only: bool,
    label: str,
    copy_data_report: bool,
    ncbi_api_key: Optional[str],
    retries: int,
    no_progressbar: bool,
    verbose: bool,
) -> Tuple[Path, Optional[Path]]:
    """
    Download virus genomes for a given taxon name or taxid,
    extract genome FASTA, and write one-line FASTA.

    Returns: (fasta_path, copied_data_report_path_or_None)
    """
    pkg_zip = outdir / f"{label}.zip"
    extract_dir = outdir / f"{label}_pkg"
    extract_dir.mkdir(parents=True, exist_ok=True)

    cmd = ["datasets", "download", "virus", "genome", "taxon", str(taxon)]
    if refseq_only:
        cmd.append("--refseq")
    cmd.extend(["--include", "genome"])
    if ncbi_api_key:
        cmd.extend(["--api-key", ncbi_api_key])
    cmd.extend(["--filename", str(pkg_zip)])

    run_datasets_download(
        cmd=cmd,
        out_zip=pkg_zip,
        retries=retries,
        no_progressbar=no_progressbar,
        verbose=verbose,
    )

    safe_unzip(pkg_zip, extract_dir)

    copied_report = None
    if copy_data_report:
        report = _find_data_report_jsonl(extract_dir)
        if report:
            copied_report = outdir / f"{label}_data_report.jsonl"
            shutil.copy2(report, copied_report)
            log(f"[OK] Copied data_report.jsonl -> {copied_report}", verbose)
        else:
            log(f"[WARN] Could not find data_report.jsonl inside {extract_dir}", verbose)

    preferred_root = extract_dir / "ncbi_dataset" / "data"
    preferred_candidates: List[Path] = []
    if preferred_root.exists():
        preferred_candidates = (
            list(preferred_root.rglob("*.fna"))
            + list(preferred_root.rglob("*.fasta"))
            + list(preferred_root.rglob("*.fa"))
        )

    all_candidates = (
        list(extract_dir.rglob("*.fna"))
        + list(extract_dir.rglob("*.fasta"))
        + list(extract_dir.rglob("*.fa"))
    )

    best = _pick_largest_fasta(preferred_candidates) or _pick_largest_fasta(all_candidates)
    if not best:
        raise RuntimeError(f"Could not locate FASTA inside datasets package: {pkg_zip}")

    out_fa = outdir / f"{label}.fasta"
    write_one_line_fasta(best, out_fa)
    return out_fa, copied_report


def edirect_download_nt(outdir: Path, query: str, max_records: int, label: str) -> Path:
    which_or_die("esearch", hint="Install NCBI EDirect: https://www.ncbi.nlm.nih.gov/books/NBK179288/")
    which_or_die("efetch", hint="Install NCBI EDirect: https://www.ncbi.nlm.nih.gov/books/NBK179288/")

    out_raw = outdir / f"{label}.raw.fasta"
    out_fa = outdir / f"{label}.fasta"

    bash_cmd = (
        f'esearch -db nucleotide -query "{query}" '
        f"| efetch -format fasta "
        f'| head -n {max_records * 200} '
    )
    run(["bash", "-lc", f"{bash_cmd} > {out_raw}"])

    write_one_line_fasta(out_raw, out_fa)
    return out_fa


def viroiddb_download(outdir: Path, label: str, tag: str = "2021-06-06") -> Path:
    which_or_die("curl", hint="Need curl to download ViroidDB zip snapshot.")
    zip_url = f"https://github.com/Benjamin-Lee/viroiddb/archive/refs/tags/{tag}.zip"
    zip_path = outdir / f"{label}_{tag}.zip"
    extract_dir = outdir / f"{label}_{tag}_src"
    extract_dir.mkdir(parents=True, exist_ok=True)

    run(["curl", "-L", "-o", str(zip_path), zip_url])
    safe_unzip(zip_path, extract_dir)

    candidates = (
        list(extract_dir.rglob("*.fasta"))
        + list(extract_dir.rglob("*.fa"))
        + list(extract_dir.rglob("*.fna"))
    )
    if not candidates:
        raise SystemExit("[ERROR] Could not find FASTA inside ViroidDB snapshot.")

    best = _pick_largest_fasta(candidates)
    if not best:
        raise SystemExit("[ERROR] ViroidDB FASTA candidates were found but none were readable.")

    out_fa = outdir / f"{label}.fasta"
    write_one_line_fasta(best, out_fa)
    return out_fa


def copy_unclassified_fasta_from_viroiddb(outdir: Path, viroiddb_tag: str) -> Optional[Path]:
    src = outdir / f"viroiddb_{viroiddb_tag}_src" / f"viroiddb-{viroiddb_tag}" / "db" / "unclassified.fasta"
    if not src.exists():
        hits = list((outdir / f"viroiddb_{viroiddb_tag}_src").rglob("unclassified.fasta"))
        if hits:
            src = hits[0]
        else:
            print(f"[WARN] Could not find ViroidDB unclassified.fasta at expected location: {src}")
            return None

    dest = outdir / "unclassified.fasta"
    write_one_line_fasta(src, dest)
    print(f"[OK] Copied unclassified.fasta -> {dest}")
    return dest


def maybe_sleep(seconds: int, reason: str, verbose: bool) -> None:
    if seconds <= 0:
        return
    log(f"[INFO] Sleeping {seconds}s ({reason})...", verbose)
    time.sleep(seconds)


def read_family_list(path: Path) -> List[str]:
    fams: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fams.append(line.split("\t")[0].strip())
    seen = set()
    out: List[str] = []
    for f in fams:
        if f and f not in seen:
            out.append(f)
            seen.add(f)
    return out


def merge_jsonl(inputs: List[Path], out_jsonl: Path) -> None:
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with out_jsonl.open("w", encoding="utf-8") as outj:
        for rp in inputs:
            with rp.open("r", encoding="utf-8", errors="replace") as inj:
                shutil.copyfileobj(inj, outj)


def maybe_copy_aliases(
    label: str,
    merged_fa: Optional[Path],
    merged_jsonl: Optional[Path],
    outdir: Path,
    also_write_refseq_alias: bool,
    verbose: bool,
) -> None:
    """
    If you ran label=ncbi_viral in family mode, optionally also write:
      refseq_viral__ALL_FAMILIES.fasta
      refseq_viral__ALL_FAMILIES.data_report.jsonl
    at top-level outdir, as plain copies, for backwards compatibility / convenience.
    """
    if not also_write_refseq_alias:
        return
    if label != "ncbi_viral":
        return

    if merged_fa and merged_fa.exists():
        alias_fa = outdir / "refseq_viral__ALL_FAMILIES.fasta"
        shutil.copy2(merged_fa, alias_fa)
        log(f"[OK] Alias copy -> {alias_fa}", verbose)

    if merged_jsonl and merged_jsonl.exists():
        alias_j = outdir / "refseq_viral__ALL_FAMILIES.data_report.jsonl"
        shutil.copy2(merged_jsonl, alias_j)
        log(f"[OK] Alias copy -> {alias_j}", verbose)


def main():
    ap = argparse.ArgumentParser(description="Download viral/viroid sequences from selected sources.")
    ap.add_argument("-o", "--outdir", default="downloads_viral_db", help="Output directory")

    ap.add_argument("--refseq_viral", action="store_true",
                    help="Download viral genomes via NCBI Datasets (RefSeq only)")
    ap.add_argument("--ncbi_viral", action="store_true",
                    help="Download viral genomes via NCBI Datasets (GenBank+RefSeq)")

    ap.add_argument("--viral_families", default=None,
                    help="TSV/TXT with one family name per line. If set, refseq/ncbi viral downloads run per-family and failures are skipped.")

    ap.add_argument("--ncbi_nt", action="store_true",
                    help="Download sequences from NCBI nucleotide using EDirect (requires --nt_query)")
    ap.add_argument("--nt_query", default=None,
                    help='EDirect query for nucleotide DB, e.g. \'txid10239[Organism] AND "complete genome"[Title]\'')
    ap.add_argument("--nt_max_records", type=int, default=50000,
                    help="Safety cap for nt downloads (default 50k; adjust carefully)")

    ap.add_argument("--viroid_db", action="store_true",
                    help="Download ViroidDB (GitHub snapshot fallback)")
    ap.add_argument("--viroiddb_tag", default="2021-06-06",
                    help="ViroidDB GitHub tag to download (default 2021-06-06)")

    ap.add_argument("--out_combined", default="combined_downloads.fasta",
                    help="Combined FASTA filename (top-level downloads only)")

    ap.add_argument(
        "--merge_reports",
        action="store_true",
        default=True,
        help="When downloading by family, merge per-family data_report.jsonl into a single JSONL per label (default: on).",
    )
    ap.add_argument(
        "--no-merge_reports",
        dest="merge_reports",
        action="store_false",
        help="Disable merging per-family data_report.jsonl files.",
    )

    # Reliability knobs
    ap.add_argument("--datasets_retries", type=int, default=5,
                    help="Retries for NCBI Datasets downloads (default 5)")
    ap.add_argument("--ncbi_api_key", default=None,
                    help="NCBI API key for datasets downloads (optional). If not provided, uses env NCBI_API_KEY if set.")
    ap.add_argument("--sleep_between_downloads", type=int, default=10,
                    help="Seconds to sleep between download tasks (default 10; set 0 to disable)")

    ap.add_argument("--no-progressbar", dest="no_progressbar", action="store_true", default=True,
                    help="Pass --no-progressbar to datasets (default ON).")
    ap.add_argument("--progressbar", dest="no_progressbar", action="store_false",
                    help="Allow datasets progress bar (may emit ANSI escapes in logs).")

    ap.add_argument("--verbose", action="store_true", help="Verbose logging (recommended for family downloads).")

    # ✅ CHANGED DEFAULT: alias OFF unless explicitly requested
    ap.add_argument("--also_write_refseq_alias", action="store_true", default=False,
                    help="When label is ncbi_viral in family mode, also write alias copies named "
                         "refseq_viral__ALL_FAMILIES.fasta / .data_report.jsonl at top-level outdir (default OFF).")
    ap.add_argument("--no_refseq_alias", dest="also_write_refseq_alias", action="store_false",
                    help="Disable writing refseq_viral__ALL_FAMILIES.* alias copies.")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ncbi_api_key = args.ncbi_api_key or os.environ.get("NCBI_API_KEY")

    family_list: Optional[List[str]] = None
    if args.viral_families:
        family_list = read_family_list(Path(args.viral_families))
        if not family_list:
            raise SystemExit(f"[ERROR] --viral_families provided but no family names found: {args.viral_families}")
        log(f"[INFO] Loaded {len(family_list)} families from {args.viral_families}", args.verbose)

    skipped_path = outdir / "skipped_families.tsv"
    skipped: List[Tuple[str, str]] = []  # (family, reason)

    fastas_top: List[Path] = []
    downloaded_any = False

    # ---- RefSeq viral ----
    if args.refseq_viral:
        label = "refseq_viral"
        if family_list:
            fam_fastas: List[Path] = []
            fam_reports: List[Path] = []

            for i, fam in enumerate(family_list, start=1):
                fam_dir = outdir / "families" / fam
                fam_dir.mkdir(parents=True, exist_ok=True)

                log(f"[INFO] ({i}/{len(family_list)}) RefSeq family download: {fam}", args.verbose)

                try:
                    fa, rep = datasets_download_virus_taxon(
                        outdir=fam_dir,
                        taxon=fam,
                        refseq_only=True,
                        label=label,
                        copy_data_report=True,
                        ncbi_api_key=ncbi_api_key,
                        retries=args.datasets_retries,
                        no_progressbar=args.no_progressbar,
                        verbose=args.verbose,
                    )
                    fam_fastas.append(fa)
                    if rep:
                        fam_reports.append(rep)
                    downloaded_any = True
                except Exception as e:
                    skipped.append((fam, f"refseq_viral_failed: {str(e)}"))
                    log(f"[WARN] Skipping family '{fam}' after {args.datasets_retries} retries (RefSeq).", True)

                maybe_sleep(args.sleep_between_downloads, f"between family downloads (RefSeq): {fam}", args.verbose)

            if fam_fastas:
                merged_fa_top = outdir / f"{label}__ALL_FAMILIES.fasta"
                concat_fastas(fam_fastas, merged_fa_top)
                log(f"[OK] Merged {label} FASTA -> {merged_fa_top}", True)

            if args.merge_reports and fam_reports:
                merged_rep_top = outdir / f"{label}__ALL_FAMILIES.data_report.jsonl"
                merge_jsonl(fam_reports, merged_rep_top)
                log(f"[OK] Merged {label} JSONL -> {merged_rep_top}", True)

        else:
            fa, _rep = datasets_download_virus_taxon(
                outdir=outdir,
                taxon=VIRUSES_TAXID,
                refseq_only=True,
                label=label,
                copy_data_report=True,
                ncbi_api_key=ncbi_api_key,
                retries=args.datasets_retries,
                no_progressbar=args.no_progressbar,
                verbose=args.verbose,
            )
            fastas_top.append(fa)
            downloaded_any = True
            log(f"[OK] refseq_viral -> {fa}", True)
            maybe_sleep(args.sleep_between_downloads, "between downloads", args.verbose)

    # ---- NCBI viral (GenBank+RefSeq) ----
    if args.ncbi_viral:
        label = "ncbi_viral"
        if family_list:
            fam_fastas: List[Path] = []
            fam_reports: List[Path] = []

            for i, fam in enumerate(family_list, start=1):
                fam_dir = outdir / "families" / fam
                fam_dir.mkdir(parents=True, exist_ok=True)

                log(f"[INFO] ({i}/{len(family_list)}) NCBI_Virus family download: {fam}", args.verbose)

                try:
                    fa, rep = datasets_download_virus_taxon(
                        outdir=fam_dir,
                        taxon=fam,
                        refseq_only=False,
                        label=label,
                        copy_data_report=True,
                        ncbi_api_key=ncbi_api_key,
                        retries=args.datasets_retries,
                        no_progressbar=args.no_progressbar,
                        verbose=args.verbose,
                    )
                    fam_fastas.append(fa)
                    if rep:
                        fam_reports.append(rep)
                    downloaded_any = True
                except Exception as e:
                    skipped.append((fam, f"ncbi_viral_failed: {str(e)}"))
                    log(f"[WARN] Skipping family '{fam}' after {args.datasets_retries} retries (NCBI).", True)

                maybe_sleep(args.sleep_between_downloads, f"between family downloads (NCBI): {fam}", args.verbose)

            merged_fa_top: Optional[Path] = None
            merged_rep_top: Optional[Path] = None

            if fam_fastas:
                merged_fa_top = outdir / f"{label}__ALL_FAMILIES.fasta"
                concat_fastas(fam_fastas, merged_fa_top)
                log(f"[OK] Merged {label} FASTA -> {merged_fa_top}", True)

            if args.merge_reports and fam_reports:
                merged_rep_top = outdir / f"{label}__ALL_FAMILIES.data_report.jsonl"
                merge_jsonl(fam_reports, merged_rep_top)
                log(f"[OK] Merged {label} JSONL -> {merged_rep_top}", True)

            # ✅ Alias now OFF by default
            maybe_copy_aliases(
                label=label,
                merged_fa=merged_fa_top,
                merged_jsonl=merged_rep_top,
                outdir=outdir,
                also_write_refseq_alias=args.also_write_refseq_alias,
                verbose=True,
            )

        else:
            fa, _rep = datasets_download_virus_taxon(
                outdir=outdir,
                taxon=VIRUSES_TAXID,
                refseq_only=False,
                label=label,
                copy_data_report=True,
                ncbi_api_key=ncbi_api_key,
                retries=args.datasets_retries,
                no_progressbar=args.no_progressbar,
                verbose=args.verbose,
            )
            fastas_top.append(fa)
            downloaded_any = True
            log(f"[OK] ncbi_viral -> {fa}", True)
            maybe_sleep(args.sleep_between_downloads, "between downloads", args.verbose)

    # ---- NCBI nt via edirect ----
    if args.ncbi_nt:
        if not args.nt_query:
            raise SystemExit("[ERROR] --ncbi_nt requires --nt_query (nt is huge; be explicit).")
        fa = edirect_download_nt(outdir, query=args.nt_query, max_records=args.nt_max_records, label="ncbi_nt")
        fastas_top.append(fa)
        downloaded_any = True
        log(f"[OK] ncbi_nt -> {fa}", True)
        maybe_sleep(args.sleep_between_downloads, "between downloads", args.verbose)

    # ---- ViroidDB ----
    if args.viroid_db:
        fa = viroiddb_download(outdir, label="viroiddb", tag=args.viroiddb_tag)
        fastas_top.append(fa)
        downloaded_any = True
        log(f"[OK] viroid_db -> {fa}", True)
        _ = copy_unclassified_fasta_from_viroiddb(outdir, viroiddb_tag=args.viroiddb_tag)
        maybe_sleep(args.sleep_between_downloads, "between downloads", args.verbose)

    # ---- Write skipped families ----
    if skipped:
        with skipped_path.open("w", encoding="utf-8") as fh:
            fh.write("family\treason\n")
            for fam, reason in skipped:
                fh.write(f"{fam}\t{reason}\n")
        log(f"[WARN] Skipped {len(skipped)} families. See: {skipped_path}", True)

    if not downloaded_any and not fastas_top:
        raise SystemExit("[ERROR] No sources selected or all downloads failed.")

    # ---- Combined top-level FASTA ----
    if len(fastas_top) >= 2:
        combined = outdir / args.out_combined
        concat_fastas(fastas_top, combined)
        log(f"[OK] Combined FASTA (top-level) -> {combined}", True)
    elif len(fastas_top) == 1:
        log(f"[INFO] Only one top-level FASTA produced ({fastas_top[0].name}); not writing {args.out_combined}.", True)


if __name__ == "__main__":
    main()
