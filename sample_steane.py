import numpy as np
from src import *

myst = FIVE();print(myst)

print(myst.get_syndrome(myst._L[0]))
print(myst.get_syndrome(myst._L[1]))
#exit()
def BLER(myqec,MONTE=10000):
    result = []
    d =10
    input_BER = [1/(2**d)*(2**i) for i in range(d)]
    input_BER = [1/4]
    n = myqec.n
    print(n)
    for p in input_BER:
        ble = 0
        for monte in range(1,MONTE+1):
            e = depolarizing_noise(n,p)
            #e = np.array([0,0,0,0,0,0,0,0,0,1,0,0,0,0],dtype="i1")
            #ee = myst.hard_decode(e)
            L = myqec.ML_decode(np.array([1-p,p/3,p/3,p/3]),myqec.get_T(myqec.get_syndrome(e)))
            #print(e,myqec.get_syndrome(e))
            #print(ee,myqec.get_syndrome(ee))
            #print(e^ee)
            #exit()
            #print("Monte:",monte," p:",p," BLER:",ble/monte,ee^e)
            for i in range(len(myqec.L)):
                if not myst.in_S(e^ee):
                    ble+=1
                    break
        print("誤り確率:",ble/MONTE)
        result.append(ble/MONTE)
    print(input_BER)
    print(result)
BLER(FIVE())
