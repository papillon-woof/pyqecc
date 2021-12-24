import numpy as np

def arr2int(arr):
    results=0
    lengths = len(arr)
    for i in range(lengths):
        results+=2**(lengths-i-1)*arr[i]
    return results
def int2arr(i,k):
    if type(i)==list:
        ValueError("i is list!")
    results = []
    for ind in range(k):
        results.append(i>>(k-ind-1)&1)
    return np.array(results)
def symplex_binary_inner_product(a,b):
    n = a.T.shape[0]//2
    z = np.zeros((n,n),dtype='i1')
    i = np.identity(n,dtype='i1')
    Lam = np.c_[np.r_[z,i],np.r_[i,z]]
    return np.mod(np.dot(np.dot(a,Lam),b),2)

def gaussjordan(X, change=0):
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
        filtre_down = A[pivot_old+1:m, j]
        pivot = np.argmax(filtre_down)+pivot_old+1

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
                        P[i, :] = abs(P[i, :]-P[pivot_old, :])
                    A[i, :] = abs(A[i, :]-A[pivot_old, :])

        if pivot_old == m-1:
            break

    if change:
        return A, P
    return A
