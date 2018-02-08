from scipy.optimize import minimize, basinhopping


def optimizer(
        obj_func,
        initial_theta,
        bounds,
        gradient=True,
        minimizer='L-BFGS-B',
        hopping=0,
        **kwargs
):
    """Substitute optimizer in scikit-learn Gaussian Process function.

    Note 'L-BFGS-B' is equivalent to the standard optimizer used in
    scikit-learn. This function allows for more direct control over the
    arguments. https://docs.scipy.org/doc/scipy/reference/optimize.html

    Parameters:
    -----------
    obj_func: function
      scikit-learn objective function.
    initial_theta: array (n,)
      Hyperparameters to be optimized against.
    bounds: list of tuples (n, 2)
      lower and upper bounds for each hyper parameter.
    gradient: bool
      Include the gradient for the optimization function.
    minimizer: str
      a scipy minimization method.
    hopping: int
      Perform a number of basin hopping steps.

    Returns:
    --------
    theta_opt: list (n,)
      Optimized hyperparameters.
    func_min: float
      Value of the minimized objective function.
    """
    margs = {
        'method': minimizer,
        'args': (gradient,),
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
        m = minimize(
            obj_func,
            initial_theta,
            **margs
        )
    theta_opt = m.x
    func_min = m.fun

    return theta_opt, func_min
