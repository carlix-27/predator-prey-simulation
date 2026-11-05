import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'docs'))
from euler_cosner import resolver_euler
from rk2_cosner import resolver_rk2

class Integradores(unittest.TestCase):
    def test_intervalo_y_validacion(self):
        for solver in (resolver_euler, resolver_rk2):
            f = lambda t, y, p: np.ones(3)
            t, y = solver(f, [1, 1, 1], 0, 1, .3, {})
            self.assertAlmostEqual(t[-1], 1)
            np.testing.assert_allclose(y[-1], [2, 2, 2])
            for y0, end, h in [([1, 1, 1], 1, 0), ([1, 1, 1], -1, .1), ([1, -1, 1], 1, .1), ([1, np.nan, 1], 1, .1), ([1, 1, 1], 1, np.inf)]:
                with self.assertRaises(ValueError): solver(f, y0, 0, end, h, {})
            for bad in (np.inf, -100):
                with self.assertRaises(ValueError): solver(lambda t, y, p: np.full(3, bad), [1, 1, 1], 0, 1, .1, {})

    def test_orden(self):
        f = lambda t, y, p: -y
        for solver, expected in ((resolver_euler, 2), (resolver_rk2, 4)):
            errors = []
            for h in (.1, .05):
                _, y = solver(f, [1, 1, 1], 0, 1, h, {})
                errors.append(abs(y[-1, 0] - np.exp(-1)))
            self.assertGreater(errors[0] / errors[1], expected * .8)
            self.assertLess(errors[0] / errors[1], expected * 1.2)

if __name__ == '__main__': unittest.main()
