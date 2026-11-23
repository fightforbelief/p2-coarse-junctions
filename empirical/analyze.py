"""Run from the repository root with python -m empirical.analyze."""
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import h5py
import numpy as np
import scipy
from scipy.io import loadmat

from empirical.trajectory import endpoint_score, excursion_table, excursion_test, labels_for_trace, marked_visit_test


def sha(path, algo="sha256"):
    with path.open("rb") as f:
        return hashlib.file_digest(f, algo).hexdigest()


def read_harmonic(path):
    if sha(path) != "992deadeaf290153c45d7ff2bbb974d520a32587b1fa4500dbedf18d16ea2059":
        raise ValueError("Wrong DeepCalib recording")
    x = loadmat(path, variable_names=["x"])["x"]
    if x.shape != (1, 50000):
        raise ValueError("Expected the published one-dimensional recording")
    return x.T, {"source": "DeepCalib bea01cc, Data_Example1b.mat",
                 "window": "complete published x record"}


def analyze(name, path, values, metadata, out):
    labels, split = labels_for_trace(values)
    d = 2 * values.shape[1]
    metadata.update(split, frames=len(values), axes=values.shape[1], sha256=sha(path),
                    numpy=np.__version__, scipy=scipy.__version__, h5py=h5py.__version__,
                    python=platform.python_version(), platform=platform.system())
    raw, rows, excursions = {}, [], {}
    for control, labs in [("chronological", labels),
                          ("shuffled", np.random.default_rng(20261123).permutation(labels))]:
        table = excursion_table(labs, d)
        excursions[control] = excursion_test(table, 20261123)
        excursions[control]["marked_visits"] = marked_visit_test(labs, d)
        pieces = [("pooled", labs)]
        if control == "chronological":
            pieces.extend((str(i+1), v) for i,v in enumerate(np.array_split(labs, 8)))
        for block, sequence in pieces:
            for lag in (1, 4, 16, 64, 256, 1024):
                score, counts = endpoint_score(sequence, d, lag)
                key = f"{control}/{block}/{lag}"
                raw[key] = counts
                row = dict(order=control, block=block, lag=lag)
                for label, value in score.items():
                    if isinstance(value, list):
                        row.update({f"{label}{i+1}": x for i,x in enumerate(value)})
                    else:
                        row[label] = value
                rows.append(row)
    out.mkdir(parents=True, exist_ok=True)
    with (out / f"{name}.csv").open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    data = dict(metadata=metadata, excursions=excursions, counts=raw)
    (out / f"{name}.json").write_text(json.dumps(data, sort_keys=True, indent=1)+"\n",
                                      encoding="utf8", newline="\n")


def main():
    from empirical.passive import read_passive
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--passive", required=True, type=Path)
    cli.add_argument("--harmonic", required=True, type=Path)
    cli.add_argument("--out", type=Path, default=Path("results/trajectories"))
    args = cli.parse_args()
    for name, path, reader in [("harmonic", args.harmonic, read_harmonic),
                               ("passive", args.passive, read_passive)]:
        values, metadata = reader(path)
        analyze(name, path, values, metadata, args.out)
        print(f"{name}: {len(values)} frames, {values.shape[1]} axes")
    inputs = [Path("junction.py"), *sorted(Path("empirical").glob("*.py")),
              Path("empirical/protocol.md"), Path("empirical/requirements.txt")]
    outputs = [args.out / f"{name}.{ext}" for name in ("passive", "harmonic")
               for ext in ("json", "csv")]
    manifest = {"inputs": {p.as_posix(): sha(p) for p in inputs},
                "outputs": {p.name: sha(p) for p in outputs}}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n",
                                          encoding="utf8", newline="\n")


if __name__ == "__main__":
    main()
