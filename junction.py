"""Local hidden-state bounds from reversible coarse transition masses."""
from itertools import combinations
from math import comb, fsum, isfinite, sqrt


def partitions(items, groups):
    """Unlabelled partitions into at most groups nonempty cells."""
    blocks = []

    def visit(i):
        if i == items:
            yield tuple(tuple(b) for b in blocks)
            return
        for block in blocks:
            block.append(i)
            yield from visit(i + 1)
            block.pop()
        if len(blocks) < groups:
            blocks.append([i])
            yield from visit(i + 1)
            blocks.pop()

    yield from visit(0)


def optimal_masses(sizes, costs, total):
    """Minimize sum(cost / mass) over masses >= sizes, by an active set."""
    mass = list(sizes)
    active = {i for i, b in enumerate(costs) if b > 0}
    if not active:
        mass[0] += total - fsum(sizes)
        return mass
    while active:
        budget = total - fsum(m for i, m in enumerate(mass) if i not in active)
        denom = fsum(sqrt(costs[i]) for i in active)
        proposed = {i: budget * sqrt(costs[i]) / denom for i in active}
        fixed = {i for i in active if proposed[i] < sizes[i]}
        if not fixed:
            for i in active:
                mass[i] = proposed[i]
            break
        active -= fixed
    return mass


def envelope(flows, total, k, weights=None, *, curvature=False):
    """Exact curvature bound or exact finite-lag bound for metric pair costs.

    Returns the bound, a minimizing port partition, and the hub masses.
    At curvature level zero masses denote limits. Small numbers of ports only.
    """
    d = len(flows)
    if d < 1 or not isinstance(k, int) or k < 1 or not isfinite(total) or total <= 0:
        raise ValueError("invalid dimensions or hub mass")
    if any(not isfinite(a) or a < 0 for a in flows):
        raise ValueError("flows must be finite and nonnegative")
    if not curvature and fsum(flows) > total:
        raise ValueError("transition masses exceed the hub mass")
    if weights is None:
        weights = [[float(i != j) for j in range(d)] for i in range(d)]
    if len(weights) != d or any(len(row) != d for row in weights):
        raise ValueError("wrong cost matrix shape")
    for i in range(d):
        for j in range(d):
            w = weights[i][j]
            if not isfinite(w) or w < 0 or w != weights[j][i] or (i == j and w != 0):
                raise ValueError("expected symmetric, nonnegative, zero-diagonal costs")
            if not curvature and any(w > weights[i][h] + weights[h][j] for h in range(d)):
                raise ValueError("finite-lag partition formula requires metric costs")
    best = None
    for blocks in partitions(d, min(k, d)):
        sizes = [0 if curvature else fsum(flows[i] for i in block) for block in blocks]
        costs = [fsum(weights[i][j] * flows[i] * flows[j] for i, j in combinations(block, 2))
                 for block in blocks]
        masses = optimal_masses(sizes, costs, total)
        score = fsum(b / m for b, m in zip(costs, masses) if b)
        if best is None or score < best[0]:
            best = score, blocks, tuple(masses)
    return best


def equal_port_bound(d, k, a, total):
    if not isinstance(d, int) or not isinstance(k, int) or d < 1 or k < 1:
        raise ValueError("state and port counts must be positive integers")
    if not isfinite(a) or not isfinite(total) or total <= 0 or a < 0 or d*a > total:
        raise ValueError("equal port masses must fit in the hub")
    if k >= d or a == 0:
        return 0.0
    return comb(d - k + 1, 2) * a * a / (total - (k - 1) * a)
