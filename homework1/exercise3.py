# Hint: you should only need the following functions from numpy and scipy
from numpy import diag, tril, triu, dot, array
from numpy.linalg import norm
from scipy.linalg import solve_triangular

def decompose(A):
    """Decomposes a matrix into its compenent parts

    The component parts: 
        L = lower triangular matrix with 0s on the diagonal
        U = upper triangular matrix with 0 on the diagonal
        D = Diagnol matrix of A

    Parameters
    ----------
    A : a numpy matrix

    Returns
    -------
    D, L, U: numpy matricies
    """
    L = tril(A, -1)
    U = triu(A, 1)
    D = A - U - L
    return [D, L, U]

def is_sdd(A):
    """Checks if a matrix is strictly diagonally dominant.

    Parameters
    ----------
    A : a numpy matrix

    Returns
    -------
    int 1 if the matrix is strictly diagonally dominant
    int 0 if the matrix is not strictly diagonally dominant
    """
    for r in range(0, len(A[:])):
        if(abs(A[r,r]) <= sum(abs(A[r,:])) - abs(A[r,r])):
            return 0
    return 1

def jacobi_step(D, L, U, b, xk):
    """Conducts a single step of the Jacobi Method

    The Jacobi Step: x_(k+1)[i] = (b - sum j from 1 to n of (U + L)[i,j] * x_k)/D[i,i]

    Parameters
    ----------
    D: a numpy n x n matrix containing the diagonals
    L: a numpy n x n matrix containing the lower triangular matrix
    U: a numpy n x n matrix containing the upper triangular matrix
    b: an n-length vector which is the solution to Ax = b
    xk = an n-length vector which is the current x in Ax = b

    Returns
    -------
    xk1: an n-length numpy list containing the next step in the Jacobi iteration
    """
    T = L + U
    xk1 = [0]*len(xk);

    for r in range(0, len(L[:])):
        xk1[r] = (b[r] - dot(T[r], xk))/D[r,r]
    return array(xk1)

def jacobi_iteration(A, b, x0, epsilon=1e-8):
    """Conducts the entire Jacobi sequence

    The Jacobi Method: while 2-norm of x_(k+1) - x_k is greater than the tolerance
        x_(k+1)[i] = (b - sum j from 1 to n of (U + L)[i,j] * x_k)/D[i,i]


    Parameters
    ----------
    A: a numpy n x n matrix as a part of Ax=b
    b: an n-length vector which is the solution to Ax = b
    x0 = an n-length vector which contains the inital guess to Ax=b
    epsilon = a keyword argument float which is the tolerance Jacobi operates against

    Returns
    -------
    xk1: an n-length numpy list containing the solution to Ax = b
    """
    D, L, U = decompose(A)
    xk = x0
    xk1 = jacobi_step(D, L, U, b, xk)
    print type(xk), type(xk1)

    while norm(xk1 - xk, 2) > epsilon:
        xk = xk1
        xk1 = jacobi_step(D, L, U, b, xk)
        print norm(xk1 - xk, 2)
    return xk1

def gauss_seidel_step(D, L, U, b, xk):
    """Conducts a single step of the Jacobi Method 

    The Gauss-Seidel Step: x_(k+1)[i] = (b - sum j from 1 to n of U[i,j] * x_k - sum j from 1 to n of L[i,j] * x_(k+1) )/D[i,i]

    Parameters
    ----------
    D: a numpy n x n matrix containing the diagonals
    L: a numpy n x n matrix containing the lower triangular matrix
    U: a numpy n x n matrix containing the upper triangular matrix
    b: an n-length vector which is the solution to Ax = b
    xk = an n-length vector which is the current x in Ax = b

    Returns
    -------
    xk1: an n-length numpy list containing the next step in the Gauss-Seidel iteration
    """
    xk1 = [0]*len(xk);

    for r in range(0, len(L[:])):
        xk1[r] = (b[r] - dot(L[r], xk1) - dot(U[r], xk))/D[r,r]
    return array(xk1)

def gauss_seidel_iteration(A, b, x0, epsilon=1e-8):
    """Conducts the entire Gauss-Seidel sequence

    The Gauss-Seidel Method: while 2-norm of x_(k+1) - x_k is greater than the tolerance
        x_(k+1)[i] = (b - sum j from 1 to n of U[i,j] * x_k - sum j from 1 to n of L[i,j] * x_(k+1) )/D[i,i]


    Parameters
    ----------
    A: a numpy n x n matrix as a part of Ax=b
    b: an n-length vector which is the solution to Ax = b
    x0 = an n-length vector which contains the inital guess to Ax=b
    epsilon = a keyword argument float which is the tolerance Gauss-Seidel operates against

    Returns
    -------
    xk1: an n-length numpy list containing the solution to Ax = b
    """
    D, L, U = decompose(A)
    xk = x0
    xk1 = gauss_seidel_step(D, L, U, b, xk)
    print type(xk), type(xk1)
    while norm(xk1 - xk, 2) > epsilon:
        xk = xk1
        xk1 = gauss_seidel_step(D, L, U, b, xk)
        print norm(xk1 - xk, 2)
    return xk1
