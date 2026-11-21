"""Check the archived research inputs and results from the repository root."""
import csv
import hashlib
import json
from pathlib import Path


def verify():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "results/run.json").read_text(encoding="utf8"))
    if not manifest["complete"]:
        raise ValueError("Incomplete computation")
    for name, expected in manifest["sha256"].items():
        actual = hashlib.sha256((root/name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch for {name}")
    for name, expected in [("diffusion", 84), ("sampling", 1200)]:
        with (root / "results" / f"{name}.csv").open(encoding="utf8", newline="") as f:
            count = sum(1 for _ in csv.DictReader(f))
        if count != expected:
            raise ValueError(f"Wrong row count for {name}: {count}")
    print("All archived input and output digests match. 84 cases and 1,200 sample replicates.")


if __name__ == "__main__":
    verify()
