import catgen import utils
import numpy as np
import os


def test_classifier():
    with utils.cd('temp'):
        assert(os.getcwd().split('/')[-1] == 'temp')

    # Test rmean
    test_range = utils.rmean(np.arange(4))
    assert(test_range == [1.5])

    test_range = utils.rmean(np.arange(7))
    assert(test_range.tolist() == [2.0, 3.0, 4.0])

    # Test trilaterate
    points = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0]])

    intersect = utils.trilaterate([points[0]], [1])
    test_point = np.array([1, 0, 1])
    np.testing.assert_allclose(test_point, intersect)

    intersect = utils.trilaterate(points[:-1], [1, 1])
    test_point = np.array([0.5, 0.5, 0.70710678])
    np.testing.assert_allclose(test_point, intersect)

    intersect = utils.trilaterate(points, [2, 2, 2])
    test_point = np.array([0.5, 0.5, -1.87082869])
    np.testing.assert_allclose(test_point, intersect)


if __name__ == "__main__":
    test_classifier()
