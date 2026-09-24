import gzip
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE = Path(".")
INPUT_DIR = BASE / "FASTQ"
OUTPUT_DIR = BASE / "FASTQ_whippet_paired_clean"
PAIR_FILE = BASE / "Paired_samples.tsv"

OUTPUT_DIR.mkdir(exist_ok=True)

def filter_pair(pair):
    r1_name, r2_name = pair

    r1_in = INPUT_DIR / f"{r1_name}.fastq.gz"
    r2_in = INPUT_DIR / f"{r2_name}.fastq.gz"

    r1_out = OUTPUT_DIR / f"{r1_name}.fastq.gz"
    r2_out = OUTPUT_DIR / f"{r2_name}.fastq.gz"

    r1_tmp = OUTPUT_DIR / f".{r1_name}.fastq.gz.tmp"
    r2_tmp = OUTPUT_DIR / f".{r2_name}.fastq.gz.tmp"

    total = 0
    removed = 0
    kept = 0

    with gzip.open(r1_in, "rt") as f1, \
         gzip.open(r2_in, "rt") as f2, \
         gzip.open(r1_tmp, "wt") as o1, \
         gzip.open(r2_tmp, "wt") as o2:

        while True:
            h1 = f1.readline()
            h2 = f2.readline()

            if not h1 and not h2:
                break

            if not h1 or not h2:
                raise RuntimeError(
                    f"R1/R2 length mismatch: {r1_name}, {r2_name}"
                )

            s1 = f1.readline()
            s2 = f2.readline()
            p1 = f1.readline()
            p2 = f2.readline()
            q1 = f1.readline()
            q2 = f2.readline()

            if not all([s1, s2, p1, p2, q1, q2]):
                raise RuntimeError(
                    f"Incomplete FASTQ record: {r1_name}, {r2_name}"
                )

            total += 1

            # 둘 중 하나라도 N이 있으면 pair 전체 제거
            if "N" in s1.upper() or "N" in s2.upper():
                removed += 1
                continue

            o1.write(h1)
            o1.write(s1)
            o1.write(p1)
            o1.write(q1)

            o2.write(h2)
            o2.write(s2)
            o2.write(p2)
            o2.write(q2)

            kept += 1

    os.replace(r1_tmp, r1_out)
    os.replace(r2_tmp, r2_out)

    return r1_name, r2_name, total, removed, kept


pairs = []

with open(PAIR_FILE) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        cols = line.split("\t")

        if len(cols) != 2:
            raise ValueError(f"Invalid Paired_samples.tsv line: {line}")

        pairs.append((cols[0], cols[1]))

print(f"Biological samples: {len(pairs)}")
print(f"FASTQ pairs to process: {len(pairs)}")

with ProcessPoolExecutor(max_workers=6) as ex:
    futures = [ex.submit(filter_pair, p) for p in pairs]

    for fut in as_completed(futures):
        r1, r2, total, removed, kept = fut.result()
        print(
            f"DONE {r1} + {r2} | "
            f"total_pairs={total} "
            f"removed_pairs={removed} "
            f"kept_pairs={kept}",
            flush=True
        )
