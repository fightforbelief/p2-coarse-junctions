"""Extract the passive calibration segment from the pinned C-Trap file."""
import hashlib

import h5py
import numpy as np


def read_passive(path):
    with path.open("rb") as f:
        if hashlib.file_digest(f, "md5").hexdigest() != "cdc5f55039e4ff52d9c78eebf11f33c5":
            raise ValueError("Wrong passive calibration recording")
    columns = []
    clock = None
    with h5py.File(path) as h:
        for axis in ("Force 1x", "Force 1y"):
            channel = h[f"Force HF/{axis}"]
            calibration = h[f"Calibration/6/{axis}"].attrs
            t0 = int(channel.attrs["Start time (ns)"])
            rate = float(channel.attrs["Sample rate (Hz)"])
            dt = round(1e9/rate)
            begin = int(calibration["Start time (ns)"])
            end = int(calibration["Stop time (ns)"])
            current = (rate, t0, begin, end)
            if clock is not None and clock != current:
                raise ValueError("Axes must be synchronous")
            clock = current
            if dt*rate != 1e9 or (begin-t0) % dt or (end-t0) % dt:
                raise ValueError("Calibration clock is inconsistent with samples")
            i, j = (begin-t0)//dt, (end-t0)//dt
            if not 0 <= i < j <= len(channel):
                raise ValueError("Missing calibration frames")
            columns.append(channel[i:j])
    return np.column_stack(columns), dict(source="Zenodo 7729824, passive_calibration.h5",
                                         window="Calibration/6", source_begin=i, source_end=j,
                                         start_ns=begin, stop_ns=end, rate_hz=rate)
