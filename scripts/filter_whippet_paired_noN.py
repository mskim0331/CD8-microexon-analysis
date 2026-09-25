#!/usr/bin/env python3
"""Remove paired-end FASTQ records containing ambiguous N bases.

If either mate contains an N in its sequence, both mates are discarded so
record synchronization is preserved for downstream paired-end Whippet
quantification.
"""

import argparse
import gzip
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True, type=Path,
                        help="Two-column TSV containing R1 and R2 sample labels.")
    parser.add_argument("--input-dir", required=True, type=Path,
                        help="Directory containing <sample>.fastq.gz inputs.")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Directory for pair-filtered FASTQ files.")
    parser.add_argument("--workers", type=int, default=6,
                        help="Number of biological samples processed concurrently.")
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


def filter_pair(task):
    r1_name, r2_name, input_dir, output_dir = task

    r1_in = input_dir / f"{r1_name}.fastq.gz"
    r2_in = input_dir / f"{r2_name}.fastq.gz"
    r1_out = output_dir / f"{r1_name}.fastq.gz"
    r2_out = output_dir / f"{r2_name}.fastq.gz"
    r1_tmp = output_dir / f".{r1_name}.fastq.gz.tmp"
    r2_tmp = output_dir / f".{r2_name}.fastq.gz.tmp"

    if not r1_in.is_file() or not r2_in.is_file():
        raise FileNotFoundError(f"Missing input FASTQ pair: {r1_in}, {r2_in}")

    total = removed = kept = 0

    try:
        with gzip.open(r1_in, "rt") as f1,              gzip.open(r2_in, "rt") as f2,              gzip.open(r1_tmp, "wt") as o1,              gzip.open(r2_tmp, "wt") as o2:

            while True:
                h1 = f1.readline()
                h2 = f2.readline()

                if not h1 and not h2:
                    break
                if not h1 or not h2:
                    raise RuntimeError(
                        f"R1/R2 length mismatch: {r1_name}, {r2_name}"
                    )

                s1, s2 = f1.readline(), f2.readline()
                p1, p2 = f1.readline(), f2.readline()
                q1, q2 = f1.readline(), f2.readline()

                if not all([s1, s2, p1, p2, q1, q2]):
                    raise RuntimeError(
                        f"Incomplete FASTQ record: {r1_name}, {r2_name}"
                    )

                total += 1

                if "N" in s1.upper() or "N" in s2.upper():
                    removed += 1
                    continue

                o1.write(h1 + s1 + p1 + q1)
                o2.write(h2 + s2 + p2 + q2)
                kept += 1

        os.replace(r1_tmp, r1_out)
        os.replace(r2_tmp, r2_out)
    except Exception:
        r1_tmp.unlink(missing_ok=True)
        r2_tmp.unlink(missing_ok=True)
        raise

    return r1_name, r2_name, total, removed, kept


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs = read_pairs(args.pairs)

    print(f"Biological samples: {len(pairs)}", flush=True)

    tasks = [
        (r1, r2, args.input_dir, args.output_dir)
        for r1, r2 in pairs
    ]

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(filter_pair, task) for task in tasks]
        for future in as_completed(futures):
            r1, r2, total, removed, kept = future.result()
            print(
                f"DONE {r1} + {r2} | "
                f"total_pairs={total} removed_pairs={removed} kept_pairs={kept}",
                flush=True,
            )


if __name__ == "__main__":
    main()
