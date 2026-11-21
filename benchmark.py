"""Compute the reflecting heat benchmark and independent endpoint experiments."""
import csv
import hashlib
import json
import math
import platform
import random
import sys
from functools import lru_cache
from itertools import combinations
from pathlib import Path

from junction import envelope

ROOT = Path(__file__).resolve().parent
LAGS = [1 / 2048, 1 / 512, 1 / 128, 1 / 32, 1 / 8, 1 / 2, 2]
TOL = 5e-12


def regions(n, width):
    u, v = n + (n - width) // 2, n + (n + width) // 2
    r = n // 3
    return [(n, 2*n, n, 2*n), (n-r, n, u, v), (2*n, 2*n+r, u, v),
            (u, v, n-r, n), (u, v, 2*n, 2*n+r)]


@lru_cache(None)
def modes(n, left, right):
    size = 3*n
    return [math.sqrt((1 + (j != 0)) / size) * math.fsum(
        math.cos(math.pi*j*(x+0.5)/size) for x in range(left, right)) for j in range(size)]


@lru_cache(None)
def heat_1d(n, lag, a, b, c, d):
    eig = (4*n*n * math.sin(math.pi*j/(6*n))**2 for j in range(3*n))
    return math.fsum(x*y*math.exp(-lag*z) for x, y, z in zip(modes(n, a, b), modes(n, c, d), eig))/n


def joint(n, width, dy, t):
    rects = regions(n, width)
    mass = [(b-a)*(d-c)/(n*n) for a, b, c, d in rects]
    rows = []
    for a, b, c, d in rects:
        row = [heat_1d(n, t, a, b, e, f) * heat_1d(n, t*dy, c, d, g, h)
               for e, f, g, h in rects]
        rows.append(row)
    for m, row in zip(mass, rows):
        row.append(m - math.fsum(row))
    rest = 9 - math.fsum(mass)
    tail = [row[-1] for row in rows]
    rows.append(tail + [rest - math.fsum(tail)])
    return mass + [rest], rows


def records():
    for n in [9, 18, 36]:
        for width in [n//3, n]:
            for dy in [1, 0.25]:
                for lag in LAGS:
                    mass, one = joint(n, width, dy, lag)
                    _, two = joint(n, width, dy, 2*lag)
                    flows = one[0][1:5]
                    y = math.fsum(two[i][j] for i, j in combinations(range(1, 5), 2))
                    ck = math.fsum(one[i][b]*one[b][j]/mass[b]
                                   for i, j in combinations(range(1, 5), 2) for b in range(6))
                    exact = [envelope(flows, mass[0], k)[0] for k in [1, 2, 3]]
                    relaxed = [envelope(flows, mass[0], k, curvature=True)[0] for k in [1, 2, 3]]
                    yield dict(density=n, width=width/n, vertical_diffusivity=dy, lag=lag,
                               west=flows[0], east=flows[1], south=flows[2], north=flows[3],
                               observed=y, ck=ck, **{f'bound{k+1}': b for k, b in enumerate(exact)},
                               **{f'relaxed{k+1}': b for k, b in enumerate(relaxed)},
                               certified=1+sum(b-y>TOL for b in exact),
                               relaxed_certified=1+sum(b-y>TOL for b in relaxed))


def binomial_limit(count, n, lower, delta):
    """Invert the one-sided Bernoulli KL-Chernoff tail bound."""
    q = count/n
    left, right = (0.0, q) if lower else (q, 1.0)
    target = math.log(1/delta)/n
    for _ in range(52):
        p = (left+right)/2
        if p == 0 or p == 1:
            break
        kl = (q*math.log(q/p) if q else 0.0) + ((1-q)*math.log((1-q)/(1-p)) if q<1 else 0.0)
        if (lower and kl>target) or (not lower and kl<target):
            left = p
        else:
            right = p
    return left if lower else right


def endpoint_counts(rng, n, probs):
    remaining, available = n, 1.0
    counts = []
    for p in probs:
        k = rng.binomialvariate(remaining, p/available)
        counts.append(k)
        remaining -= k
        available -= p
    return counts


def sampling():
    rng = random.Random(20261121)
    mass, one = joint(18, 6, 1, 1/128)
    _, two = joint(18, 6, 1, 1/64)
    flows = one[0][1:5]
    area = sum(mass[1:5])
    heat_y = sum(two[i][j] for i, j in combinations(range(1, 5), 2))
    null_y = envelope(flows, 1, 3)[0]
    for label, y in [('heat', heat_y), ('three_phase', null_y)]:
        for n in [1000, 10000, 100000]:
            for rep in range(200):
                counts = endpoint_counts(rng, n, flows)
                hit = rng.binomialvariate(n, max(0, 2*y/area))
                lower = [binomial_limit(x, n, True, .01) for x in counts]
                upper = area * binomial_limit(hit, n, False, .01)/2
                bounds = [envelope(lower, 1, k)[0] for k in [1, 2, 3]]
                yield dict(case=label, trials=n, replicate=rep, west=counts[0], east=counts[1],
                           south=counts[2], north=counts[3], hits=hit, upper=upper,
                           bound1=bounds[0], bound2=bounds[1], bound3=bounds[2],
                           certified=1 + sum(b > upper for b in bounds))


def write_csv(path, data):
    it = iter(data)
    first = next(it)
    with path.open('w', encoding='utf8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=first.keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerow(first)
        writer.writerows(it)


def main():
    out = ROOT / 'results'
    out.mkdir(exist_ok=True)
    write_csv(out/'diffusion.csv', records())
    write_csv(out/'sampling.csv', sampling())
    inputs = ['benchmark.py', 'junction.py', 'tests/test_junction.py', 'protocol.md', 'execute.ps1']
    names = inputs + ['results/diffusion.csv', 'results/sampling.csv']
    manifest = {'python':sys.version, 'platform':platform.platform(), 'executable':str(Path(sys.executable).relative_to(ROOT)) if Path(sys.executable).is_relative_to(ROOT) else Path(sys.executable).name,
                'command':r'.venv\Scripts\python.exe benchmark.py', 'complete':True,
                'sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names}}
    (out/'run.json').write_text(json.dumps(manifest, indent=2)+'\n',encoding='utf8', newline='\n')
    print('Wrote diffusion and sampling results.')


if __name__ == '__main__':
    main()
