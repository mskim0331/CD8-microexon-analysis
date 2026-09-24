import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = Path(".")
PAIR_FILE = BASE / "Paired_samples.tsv"
INPUT_DIR = BASE / "FASTQ_whippet_paired_clean"

OUTPUT_DIR = BASE / "Whippet_Paired" / "Quant"
LOG_DIR = BASE / "Whippet_Paired" / "logs"

JULIA = "/home/pth1612/anaconda3/envs/smakemake2/bin/julia"
WHIPPET = "/home/pth1612/Whippet.jl/bin/whippet-quant.jl"
INDEX = BASE / "Whippet" / "Index" / "whippet.jls"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

pairs = []

with open(PAIR_FILE) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        r1, r2 = line.split("\t")
        pairs.append((r1, r2))

print("Biological samples:", len(pairs), flush=True)


def gzip_valid(path):
    if not Path(path).exists():
        return False

    return subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    ).returncode == 0


def run_whippet(pair):
    r1, r2 = pair

    fq1 = INPUT_DIR / (r1 + ".fastq.gz")
    fq2 = INPUT_DIR / (r2 + ".fastq.gz")

    output_prefix = OUTPUT_DIR / r1
    psi_file = Path(str(output_prefix) + ".psi.gz")
    log_file = LOG_DIR / (r1 + ".log")

    if gzip_valid(psi_file):
        return "SKIP", r1

    cmd = [
        JULIA,
        WHIPPET,
        str(fq1),
        str(fq2),
        "-x", str(INDEX),
        "-o", str(output_prefix)
    ]

    with open(log_file, "w") as log:
        result = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT
        )

    if result.returncode != 0:
        return "FAILED", r1

    if not gzip_valid(psi_file):
        return "FAILED", r1

    return "DONE", r1


failed = []

with ThreadPoolExecutor(max_workers=6) as ex:
    futures = [ex.submit(run_whippet, pair) for pair in pairs]

    for future in as_completed(futures):
        status, sample = future.result()

        print(status, sample, flush=True)

        if status == "FAILED":
            failed.append(sample)

if failed:
    print("FAILED SAMPLES:", ",".join(failed), flush=True)
    raise SystemExit(1)

print("ALL 24 PAIRED-END WHIPPET QUANTIFICATIONS COMPLETE", flush=True)
