"""Endpoint measurements and a regenerative one-state comparison."""
from itertools import combinations

import numpy as np
from scipy.stats import random_table

from junction import envelope


def labels_for_trace(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or x.shape[1] not in (1, 2) or not np.isfinite(x).all():
        raise ValueError("Expected finite, time-ordered one- or two-axis observations")
    train = len(x) // 5
    center, scale = x[:train].mean(axis=0), x[:train].std(axis=0)
    if not (scale > 0).all():
        raise ValueError("Constant calibration coordinate")
    held = (x[train:] - center) / scale
    dim = x.shape[1]
    labels = np.full(len(held), 2*dim + 1, dtype=np.int8)
    labels[np.max(np.abs(held), axis=1) <= .75] = 0
    for axis in range(dim):
        side = True if dim == 1 else np.abs(held[:, 1-axis]) <= .375
        labels[(held[:, axis] < -.75) & side] = 1 + 2*axis
        labels[(held[:, axis] > .75) & side] = 2 + 2*axis
    return labels, dict(training_frames=train, held_out=len(held),
                        center=center.tolist(), scale=scale.tolist())


def joint_counts(labels, d, lag):
    size = d + 2
    if not 0 < lag < len(labels):
        raise ValueError("Lag out of range")
    pair_ids = size * labels[:-lag].astype(np.int64) + labels[lag:]
    return np.bincount(pair_ids, minlength=size**2).reshape(size, size)


def endpoint_score(labels, d, lag):
    one = joint_counts(labels, d, lag)
    two = joint_counts(labels, d, 2*lag)
    j = (one + one.T) / (2 * one.sum())
    j2 = (two + two.T) / (2 * two.sum())
    masses = j.sum(axis=1)
    flows = j[0, 1:d+1]
    pairs = list(combinations(range(1, d+1), 2))
    observed = sum(j2[i, l] for i, l in pairs)
    # A label with no observations contributes zero joint mass.
    keep = masses > 0
    ck = (j[:, keep] / masses[keep]) @ j[keep, :]
    exact = [envelope(flows, masses[0], k)[0] for k in range(1, d)]
    relaxed = [envelope(flows, masses[0], k, curvature=True)[0] for k in range(1, d)]
    summary = dict(frames=len(labels), hub=masses[0], observed=observed,
                   ck=sum(ck[i, l] for i, l in pairs), flows=flows.tolist(),
                   exact=exact, relaxed=relaxed,
                   plugin_lower=1+sum(bool(f-observed > 1e-12) for f in exact),
                   one_step_asymmetry=float(np.max(np.abs(one-one.T))/one.sum()),
                   two_step_asymmetry=float(np.max(np.abs(two-two.T))/two.sum()))
    return summary, dict(one=one.tolist(), two=two.tolist())


def excursion_table(labels, d):
    """Pair completed port-to-port excursions in disjoint chronological pairs."""
    is_hub = labels == 0
    departures = np.flatnonzero(is_hub[:-1] & ~is_hub[1:]) + 1
    returns = np.flatnonzero(~is_hub[:-1] & is_hub[1:])
    following = np.searchsorted(returns, departures, side="left")
    valid = following < len(returns)
    start = labels[departures[valid]]
    end = labels[returns[following[valid]]]
    both_ports = (start <= d) & (end <= d)
    start, end = start[both_ports], end[both_ports]
    n = len(start)//2
    cell = d * (end[:2*n:2].astype(int)-1) + start[1:2*n:2] - 1
    return np.bincount(cell, minlength=d*d).reshape(d, d)


def excursion_test(table, seed):
    size = int(table.sum())
    if size == 0:
        raise ValueError("No complete excursion pairs")
    margins = table.sum(axis=1), table.sum(axis=0)
    matches = int(np.trace(table))
    row = np.flatnonzero(margins[0])
    col = np.flatnonzero(margins[1])
    if min(len(row), len(col)) == 1:
        null_matches = np.full(9999, matches)
    else:
        draws = random_table.rvs(margins[0][row], margins[1][col], size=9999,
                                 method="patefield", random_state=np.random.default_rng(seed))
        null_matches = draws[:, row[:, None] == col[None, :]].sum(axis=1)
    return dict(pairs=size, table=table.tolist(), matches=matches,
                expected=float(np.dot(margins[0].astype(float), margins[1])/size),
                p=(1+int(np.sum(null_matches >= matches)))/10000)


def marked_visit_test(labels, d):
    """Simultaneous anytime-valid tests of port-specific exit probabilities."""
    in_hub = labels == 0
    entries = np.flatnonzero(~in_hub[:-1] & in_hub[1:])  # last exterior frame
    exits = np.flatnonzero(in_hub[:-1] & ~in_hub[1:]) + 1  # first exterior frame
    following = np.searchsorted(exits, entries, side="right")
    complete = following < len(exits)
    before = labels[entries[complete]]
    after = labels[exits[following[complete]]]
    alternatives = np.array([.55, .65, .75, .90])
    rows = []
    for port in range(1, d+1):
        n = int(np.count_nonzero(before == port))
        hits = int(np.count_nonzero((before == port) & (after == port)))
        logs = hits * np.log(2*alternatives) + (n-hits)*np.log(2*(1-alternatives))
        log_e = float(np.logaddexp.reduce(logs) - np.log(len(logs)))
        rows.append(dict(visits=n, same=hits, log_e=log_e, rejected=bool(log_e >= np.log(120))))
    return dict(ports=rows, lower=max(1,sum(r["rejected"] for r in rows)))
