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
    directory = root / "results/trajectories"
    trajectories = json.loads((directory / "manifest.json").read_text(encoding="utf8"))
    for folder, group in [(root, "inputs"), (directory, "outputs")]:
        for name, expected in trajectories[group].items():
            if hashlib.sha256((folder / name).read_bytes()).hexdigest() != expected:
                raise ValueError(f"Trajectory SHA-256 mismatch for {name}")
    for name in ("harmonic", "passive"):
        with (directory / f"{name}.csv").open(encoding="utf8", newline="") as f:
            rows = list(csv.DictReader(f))
        data = json.loads((directory / f"{name}.json").read_text(encoding="utf8"))
        keys = {f"{r['order']}/{r['block']}/{r['lag']}" for r in rows}
        if len(rows) != 60 or len(keys) != 60 or keys != set(data["counts"]):
            raise ValueError(f"Missing endpoint scores or directional matrices for {name}")
        for r in rows:
            counts = data["counts"][f"{r['order']}/{r['block']}/{r['lag']}"]
            for matrix, multiplier in [("one", 1), ("two", 2)]:
                expected = int(r["frames"]) - multiplier * int(r["lag"])
                if sum(map(sum, counts[matrix])) != expected:
                    raise ValueError(f"Endpoint count mismatch in {name}, {matrix}")
    print("Verified 84 grid cases, 1,200 simulations, and two measured records (120 endpoint rows).")


if __name__ == "__main__":
    verify()
