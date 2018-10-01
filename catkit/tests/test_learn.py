import unittest
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
from sklearn.model_selection import train_test_split
from catkit.learn import online_learning
from catkit.pawprint.db import FingerprintDB
from catkit.learn import optimizer
import numpy as np
from glob import glob
import pkg_resources
import os
np.random.seed(42)

# Load fingerprint test data
data_path = pkg_resources.resource_filename(
    'catkit', 'tests/data/fingerprints.db')


class TestLearn(unittest.TestCase):
    """Test features of the learn module."""

    def tearDown(self):
        """Clean up any plotted images."""
        for f in glob('*.png'):
            os.unlink(f)

    def test_learning(self):
        """General testing for catkit.learn"""
        with FingerprintDB(data_path) as fp:
            X = fp.get_fingerprints(params=np.arange(20) + 1)
            y = fp.get_fingerprints(params=['Ef']).T[0]

        X0, X1, y0, y1 = train_test_split(X, y, test_size=0.6, shuffle=True)

        kernel = DotProduct() + WhiteKernel()
        gp = GaussianProcessRegressor(
            kernel=kernel, optimizer=optimizer,
            n_restarts_optimizer=0, alpha=0)
        gp.fit(X0, y0)

        # This is an ugly way to define default properties
        # The 3rd argument is for a certain number of global
        # optimization steps using the basinhopping algorithm
        optimizer.__defaults__ = (True, 'L-BFGS-B', 3)

        gp = GaussianProcessRegressor(
            kernel=kernel, optimizer=optimizer,
            n_restarts_optimizer=0, alpha=0)
        gp.fit(X0, y0)

        samples = online_learning(
            X, y, [0, 1, 2], factors=[1, 1], nsteps=3, plot=True)

        test_array = np.array([0, 1, 2, 24, 15, 23])
        np.testing.assert_allclose(samples, test_array)


if __name__ == '__main__':
    unittest.main()
