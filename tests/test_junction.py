import math
import unittest
from itertools import combinations, product

from benchmark import binomial_limit, joint, regions
from junction import envelope, equal_port_bound


class Bounds(unittest.TestCase):
    def test_balanced_formula_including_saturation(self):
        for d in range(2, 8):
            for k in range(1, d+1):
                for total in [d, d+0.1, 2*d, 100*d]:
                    self.assertAlmostEqual(envelope([1]*d, total, k)[0],
                                           equal_port_bound(d, k, 1, total))
                    self.assertAlmostEqual(envelope([1]*d, total, k, curvature=True)[0],
                                           math.comb(d-k+1, 2)/total)

    def test_exhaustive_fractional_allocations(self):
        # Six units of mass in a hub. Exhaust all two-phase assignments
        # in quanta of 1/2, including allocation of the unobserved slack.
        a = [1, 1.5, 2]
        w = [[0, 1, 3], [1, 0, 2], [3, 2, 0]]
        best = math.inf
        for parts in product(*(range(int(2*x)+1) for x in a+[1.5])):
            first = [x/2 for x in parts]
            second = [x-y for x, y in zip(a+[1.5], first)]
            score = sum(w[i][j]*v[i]*v[j]/sum(v)
                        for v in [first, second] if sum(v)
                        for i, j in combinations(range(3), 2))
            best = min(best, score)
        self.assertAlmostEqual(best, envelope(a, 6, 2, w)[0])

    def test_sharp_reversible_kernel(self):
        a = [.14, .22, .19, .09, .12]
        value, blocks, mass = envelope(a, 1, 3)
        # Individual ports have unit mass, hub states have fitted masses.
        d = len(a)
        stat = [1]*d + list(mass)
        flux = [[0.] * len(stat) for _ in stat]
        for h, group in enumerate(blocks, d):
            for i in group:
                flux[i][h] = flux[h][i] = a[i]
        for i, pi in enumerate(stat):
            flux[i][i] = pi - sum(flux[i])
            self.assertGreaterEqual(flux[i][i], -1e-15)
        observed = sum(flux[i][x]*flux[x][j]/stat[x]
                       for i, j in combinations(range(d), 2) for x in range(len(stat)))
        self.assertAlmostEqual(value, observed)

    def test_approximate_fit_radius(self):
        delta = (math.sqrt(1.25)-1.1)/2
        for increment in [-1e-4, 0, 1e-4]:
            radius = delta + increment
            gap = envelope([.2-radius]*4, 1, 3)[0] - (.05+radius)
            if increment:
                self.assertLess(gap*increment, 0)
            else:
                self.assertAlmostEqual(gap, 0)

    def test_invalid_flows(self):
        with self.assertRaises(ValueError):
            envelope([.6,.5], 1, 2)
        with self.assertRaises(ValueError):
            equal_port_bound(2, 1, -.1, 1)

    def test_odd_cycle_zero_threshold(self):
        w = [[int((i-j) % 5 in [1, 4]) for j in range(5)] for i in range(5)]
        self.assertAlmostEqual(envelope([1]*5, 1, 2, w, curvature=True)[0], 1)
        self.assertEqual(envelope([1]*5, 1, 3, w, curvature=True)[0], 0)
        with self.assertRaises(ValueError):
            envelope([.1]*5, 1, 2, w)

    def test_zero_flow_and_endpoint_ci(self):
        self.assertEqual(envelope([0,.1,.3], .4, 2)[0], 0)
        self.assertAlmostEqual(binomial_limit(0, 100, False, .01), 1-.01**.01)
        self.assertAlmostEqual(binomial_limit(100, 100, True, .01), .01**.01)


class Diffusion(unittest.TestCase):
    def test_uniformization_independent_of_cosine_solver(self):
        n, dy, t = 9, .25, 1/128
        side = 3*n
        rects = regions(n, 3)
        sets = [[y*side+x for x in range(a,b) for y in range(c,d)] for a,b,c,d in rects]
        _, exact = joint(n, 3, dy, t)
        mu = 2*(1+dy)*n*n*t
        px, py = 1/(2*(1+dy)), dy/(2*(1+dy))
        for index, start in enumerate(sets):
            v = [0.] * side**2
            for i in start:
                v[i] = 1/n**2
            p = math.exp(-mu)
            ans = [p*x for x in v]
            for order in range(1, 65):
                new = [0.] * len(v)
                for i, x in enumerate(v):
                    col, row = i % side, i // side
                    for j, weight in [(i-1 if col else i, px), (i+1 if col+1<side else i, px),
                                      (i-side if row else i, py), (i+side if row+1<side else i, py)]:
                        new[j] += weight*x
                v = new
                p *= mu/order
                for i, x in enumerate(v):
                    ans[i] += p*x
            for target, nodes in enumerate(sets):
                self.assertAlmostEqual(sum(ans[i] for i in nodes), exact[index][target], places=12)


if __name__ == '__main__':
    unittest.main()
