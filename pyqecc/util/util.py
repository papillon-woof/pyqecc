import numpy as np

X_OR_Z = 2


def arr2int(arr):
    results = 0
    lengths = len(arr)
    for i in range(lengths):
        results += 2 ** (lengths - i - 1) * arr[i]
    return results


def int2arr(i, k):
    if type(i) == list:
        ValueError("i is list!")
    results = []
    for ind in range(k):
        results.append(i >> (k - ind - 1) & 1)
    return np.array(results, dtype="i1")


def symplex_binary_inner_product(a, b):
    n = a.T.shape[0] // 2
    z = np.zeros((n, n), dtype="i1")
    i = np.identity(n, dtype="i1")
    Lam = np.c_[np.r_[z, i], np.r_[i, z]]
    return np.mod(np.dot(np.dot(a, Lam), b), 2)


def gaussjordan(X, change=False):
    """Compute the binary row reduced echelon form of X.
    Parameters
    ----------
    X: array (m, n)
    change : boolean (default, False). If True returns the inverse transform
    Returns
    -------
    if `change` == 'True':
        A: array (m, n). row reduced form of X.
        P: tranformations applied to the identity
    else:
        A: array (m, n). row reduced form of X.
    """
    A = np.copy(X)
    m, n = A.shape

    if change:
        P = np.identity(m).astype(int)

    pivot_old = -1
    for j in range(n):
        filtre_down = A[pivot_old + 1 : m, j]
        pivot = np.argmax(filtre_down) + pivot_old + 1

        if A[pivot, j]:
            pivot_old += 1
            if pivot_old != pivot:
                aux = np.copy(A[pivot, :])
                A[pivot, :] = A[pivot_old, :]
                A[pivot_old, :] = aux
                if change:
                    aux = np.copy(P[pivot, :])
                    P[pivot, :] = P[pivot_old, :]
                    P[pivot_old, :] = aux

            for i in range(m):
                if i != pivot_old and A[i, j]:
                    if change:
                        P[i, :] = abs(P[i, :] - P[pivot_old, :])
                    A[i, :] = abs(A[i, :] - A[pivot_old, :])

        if pivot_old == m - 1:
            break

    if change:
        return A, P
    return A


def any2arr(i, k):
    if type(i) == int:
        return int2arr(i, k)
    if type(i) == list:
        return np.array(i, dtype="i1")
    if (
        type(i) == np.int64
        or type(i) == np.int32
        or type(i) == np.int16
        or type(i) == np.int8
    ):
        return int2arr(i, k)
    return i


def bitwise_to_blockwise_error_probability(bitwise_error_probabilityrobability):
    """
    Input: np.array[qubit_number][4]
    Output: np.array[2 ** (2 * qubit_number)]
    """
    n = bitwise_error_probabilityrobability.shape[0]
    blockwise_error_probabilityrobability = np.ones(2 ** (2 * n))
    for ind in range(2 ** (2 * n)):
        ind_list = int2arr(ind, 2 * n)
        for ei in range(n):
            i = (ind_list[ei] << 1) + ind_list[ei + n]
            blockwise_error_probabilityrobability[
                ind
            ] *= bitwise_error_probabilityrobability[ei][i]
    return blockwise_error_probabilityrobability


def blockwise_to_bitwise_error_probability(blockwise_error_probabilityrobability):
    """
    Input: np.array[2 ** (2 * qubit_number)]
    Output: np.array[qubit_number][4]
    """
    n = int(np.log2(blockwise_error_probabilityrobability.shape[0])) // 2
    bitwise_error_probabilityrobability = np.zeros((n, 4))
    for ind in range(len(blockwise_error_probabilityrobability)):
        ind_list = int2arr(ind, 2 * n)
        for ei in range(n):
            i = (ind_list[ei] << 1) + ind_list[ei + n]
            bitwise_error_probabilityrobability[ei][
                i
            ] += blockwise_error_probabilityrobability[ind]
    return bitwise_error_probabilityrobability


def pishift(a, div_max=100000):
    if a >= 0:
        for i in range(div_max):
            if a > np.sqrt(np.pi):
                a -= 2 * np.sqrt(np.pi)
            else:
                return a
            if a < 0:
                return a
    else:
        for i in range(div_max):
            if a < np.sqrt(np.pi):
                a += 2 * np.sqrt(np.pi)
            else:
                return a
            if a >= 0:
                return a
    raise AssertionError("Div_max is proceeded.")


pishifts = np.frompyfunc(pishift, 1, 1)
