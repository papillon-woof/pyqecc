import numpy as np
from src import *

myst = FIVE();print(myst)

#exit()
def BLER(myqec,MONTE=100):
    result = []
    d =10
    input_BER = [1/(2**d)*(2**i) for i in range(d)]
    #input_BER = [1/3]
    n = myqec.n
    for p in input_BER:
        ble = 0
        myqec.set_P(np.array([1-p,p/3,p/3,p/3]));
        for monte in range(1,MONTE+1):
            #e = depolarizing_noise(n,p)
            e = np.array([0,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0,0])
            syndrome = myqec.get_syndrome(e)
            ee = myqec.decode(syndrome,mode='ML');
            #print(myqec.get_syndrome(e),myqec.get_syndrome(ee))
            #ee = myqec.decode(syndrome);
            print(e,myqec.get_syndrome(e),myqec.get_syndrome(ee))
            #print(e)
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
