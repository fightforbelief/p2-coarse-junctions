import unittest

import numpy as np

from empirical.trajectory import endpoint_score, excursion_table, excursion_test, joint_counts, marked_visit_test


class TraceChecks(unittest.TestCase):
    def test_counts_keep_direction_and_lag(self):
        labels = np.array([0, 1, 2, 1, 0, 0, 2])
        c = joint_counts(labels, 2, 2)
        self.assertEqual(c.sum(), 5)
        self.assertEqual(c[0, 2], 2)
        self.assertEqual(c[2, 0], 1)
        self.assertEqual(c[1, 1], 1)

    def test_boundary_and_exterior_excursions(self):
        # Complete marked excursions are 1->2, 2->2, 1->1, 2->1.
        # Drop initial/final partial paths, self holds and a 3->1 excursion.
        lab = np.array([2,0,0,1,2,0,2,0,3,1,0,1,0,2,1,0,1])
        self.assertEqual(excursion_table(lab, 2).tolist(), [[0, 1], [0, 1]])
        self.assertEqual(excursion_test(np.array([[0,1],[0,1]]), 1)["p"], 1)

    def test_detects_memory_in_excursion_order(self):
        out = excursion_test(np.diag([20,20]), 2)
        self.assertEqual(out["matches"], 40)
        self.assertAlmostEqual(out["expected"], 20)
        self.assertLess(out["p"], .01)

    def test_marked_visits_include_exterior_failures(self):
        # Last visit is censored; port 1 has one same-port exit and one
        # exterior exit, port 2 has one exit to port 1.
        out = marked_visit_test(np.array([1,0,0,1,0,3,2,0,1,0]), 2)
        self.assertEqual([(p["visits"], p["same"]) for p in out["ports"]], [(2,1),(1,0)])
        self.assertEqual(out["lower"], 1)
        full = np.concatenate([np.tile([i,0,i], 40) for i in range(1,5)])
        self.assertEqual(marked_visit_test(full, 4)["lower"], 4)

    def test_unknown_hub_mass_and_unused_label(self):
        # Rotations of the cycle give equal endpoint occupancies.
        lab = np.tile(np.array([0, 1, 0, 2]), 100)
        out, _ = endpoint_score(lab, 2, 1)
        self.assertAlmostEqual(out["hub"], .5)
        self.assertAlmostEqual(out["exact"][0], out["flows"][0]*out["flows"][1]/.5)


if __name__ == "__main__":
    unittest.main()
