#!/usr/bin/env python3
"""Run paired-end Whippet quantification once per biological sample."""

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--julia", required=True, type=Path)
    parser.add_argument("--whippet-quant", required=True, type=Path)
    parser.add_argument("--index", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=6)
    return parser.parse_args()


def read_pairs(path):
    pairs = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected two tab-delimited fields"
                )
            pairs.append((fields[0], fields[1]))
    return pairs


def gzip_valid(path):
    path = Path(path)
    if not path.is_file():
        return False
    return subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode == 0


def quantify(pair, args):
    r1, r2 = pair
    fq1 = args.input_dir / f"{r1}.fastq.gz"
    fq2 = args.input_dir / f"{r2}.fastq.gz"

    if not fq1.is_file() or not fq2.is_file():
        return "FAILED", r1, f"Missing input: {fq1} or {fq2}"

    output_prefix = args.output_dir / r1
    psi_file = Path(f"{output_prefix}.psi.gz")
    log_file = args.log_dir / f"{r1}.log"

    if gzip_valid(psi_file):
        return "SKIP", r1, "existing valid PSI"

    cmd = [
        str(args.julia),
        str(args.whippet_quant),
        str(fq1),
        str(fq2),
        "-x", str(args.index),
        "-o", str(output_prefix),
    ]

    with log_file.open("w") as log:
        result = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
        )

    if result.returncode != 0:
        return "FAILED", r1, f"exit code {result.returncode}"
    if not gzip_valid(psi_file):
        return "FAILED", r1, "missing or invalid .psi.gz"

    return "DONE", r1, ""


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    if not args.index.is_file():
        raise FileNotFoundError(f"Whippet index not found: {args.index}")

    pairs = read_pairs(args.pairs)
    print(f"Biological samples: {len(pairs)}", flush=True)

    failed = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(quantify, pair, args)
            for pair in pairs
        ]

        for future in as_completed(futures):
            status, sample, message = future.result()
            suffix = f" ({message})" if message else ""
            print(f"{status} {sample}{suffix}", flush=True)
            if status == "FAILED":
                failed.append(sample)

    if failed:
        raise SystemExit("FAILED SAMPLES: " + ",".join(sorted(failed)))

    print(
        f"ALL {len(pairs)} PAIRED-END WHIPPET QUANTIFICATIONS COMPLETE",
        flush=True,
    )


if __name__ == "__main__":
    main()
