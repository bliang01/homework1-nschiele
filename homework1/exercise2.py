def gradient_step(xk, df, sigma):
    """Returns a step x_k+1 in the gradient descent method

        Parameters
    ----------
    xf: float
    df : function
    sigma : float

    Returns
    -------
    float
        The next step of a gradient descent
    """
    return xk - sigma * df(xk)

def gradient_descent(f, df, x, sigma=0.5, epsilon=1e-8):
    """Returns a minima of `f` using the Gradient Descent method.

    A local minima, x*, is such that `f(x*) <= f(x)` for all `x` near `x*`.
    This function returns a local minima which is accurate to within `epsilon`.

    `gradient_descent` raises a ValueError if `sigma` is not strictly between
    zero and one.

    Parameters
    ----------
    f : function
    df : function
    x : float
    sigma : float
    epsilon : float

    Returns
    -------
    float
        The local min of f
    """
    if sigma <= 0 or sigma >= 1:
        raise ValueError('Negative Values not allowed for Collatz')

    xk = x+1;
    xk1 = x
    while abs(xk1 - xk) > epsilon:
        xk = xk1;
        xk1 = gradient_step(xk, df, sigma)

    return xk1