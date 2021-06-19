import numpy as np

def symplex_binary_inner_product(a,b):
    print(a,b)
    n = len(a)//2
    print(n)
    print((sum(a[:n]*b[n:])%2),sum(a[n:]*b[:n])%2,a[n:])
    return (sum(a[:n]*b[n:])%2)^(sum(a[n:]*b[:n])%2)
print(symplex_binary_inner_product(np.zeros([0,0,1,0]),np.zeros([1,0,0,0])))
