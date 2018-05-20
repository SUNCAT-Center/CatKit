from scipy.optimize import minimize
from scipy.optimize import basinhopping
import matplotlib.pyplot as plt
from sklearn.gaussian_process.kernels import WhiteKernel, DotProduct
from sklearn.metrics import mean_absolute_error
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn import preprocessing
import numpy as np


def optimizer(obj_func,
              initial_theta,
              bounds,
              gradient=True,
              minimizer='L-BFGS-B',
              hopping=0,
              **kwargs):
    """Substitute optimizer in scikit-learn Gaussian Process function.

    Note 'L-BFGS-B' is equivalent to the standard optimizer used in
    scikit-learn. This function allows for more direct control over the
    arguments. https://docs.scipy.org/doc/scipy/reference/optimize.html

    Parameters
    ----------
    obj_func : function
        scikit-learn objective function.
    initial_theta : array (n,)
        Hyperparameters to be optimized against.
    bounds : list of tuples (n, 2)
        Lower and upper bounds for each hyper parameter.
    gradient : bool
        Include the gradient for the optimization function.
    minimizer : str
        A scipy minimization method.
    hopping : int
        Perform a number of basin hopping steps.

    Returns
    -------
    theta_opt : list (n,)
        Optimized hyperparameters.
    func_min : float
        Value of the minimized objective function.
    """
    margs = {
        'method': minimizer,
        'args': (gradient, ),
        'jac': gradient,
        'bounds': bounds,
    }

    if hopping:
        m = basinhopping(
            obj_func,
            initial_theta,
            minimizer_kwargs=margs,
            niter=hopping,
            T=1e2,
            stepsize=2,
        )
    else:
        m = minimize(obj_func, initial_theta, **margs)
    theta_opt = m.x
    func_min = m.fun

    return theta_opt, func_min


def online_learning(X, y, samples, factors=[1.0, 1.0], nsteps=40, plot=False):
    """A simple utility for performing online learning. The main
    components required are a regression method and a scoring
    technique.

    Currently, the scoring methodology and regressor are baked in.
    These need to be made modular.

    Minimum 3 samples are required for 3 fold cross validation.
    """
    ids = np.arange(len(y))

    kernel = DotProduct() + WhiteKernel()
    regressor = GaussianProcessRegressor(
        kernel=kernel, n_restarts_optimizer=5, alpha=0)

    step = 0
    while step < nsteps:
        X0 = X[samples]
        y0 = y[samples]

        regressor.fit(X0, y0)
        yp, ys = regressor.predict(X, return_std=True)

        # Provides some form of normalization.
        # Multiples denote relative importance
        yp_scale = preprocessing.scale(yp) * factors[0]
        ys_scale = preprocessing.scale(ys) * factors[1]
        score = ys_scale - yp_scale
        srt = np.argsort(score)[::-1]

        for s in srt:
            if s not in samples:
                samples = np.concatenate([samples, [s]])
                break

        if plot:
            mae = np.round(mean_absolute_error(yp, y), 3)
            n = len(samples)

            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(ids, y, 'o', zorder=0)
            ax.errorbar(ids, yp, yerr=ys, fmt='o', zorder=1)
            ax.plot(samples, y[samples], 'o', zorder=3)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            ax.text(xlim[0] / 9.0, ylim[0] / 9.0, mae)
            plt.tight_layout()
            plt.savefig('./online-learning-RBF-{}.png'.format(n))
            plt.close()

        step += 1

    return samples
