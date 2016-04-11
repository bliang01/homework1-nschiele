# the documenation has been written for you in this  exercise

def collatz_step(n):
    """Returns the result of the Collatz function.

    The Collatz function C : N -> N is used in `collatz` to generate collatz
    sequences. Raises an error if n < 1.

    Parameters
    ----------
    n : int

    Returns
    Testing whether git will work with SMC
    -------
    int
        The result of C(n).

    """
    if n < 0:
        raise ValueError('Negative Values not allowed for Collatz')
    if n == 1:
        return 1
    if n % 2 == 0:
        return n/2
    return 3*n + 1


def collatz(n):    
    """Returns the Collatz sequence beginning with `n`.

    It is conjectured that Collatz sequences all end with `1`. Calls
    `collatz_step` at each iteration.

    Parameters
    ----------
    n : int

    Returns
    -------
    sequence : list
        A Collatz sequence.

    """
    cList = [n]
    while(n > 1):
        n = collatz_step(cList[-1])
        cList.append(n)
    return cList
