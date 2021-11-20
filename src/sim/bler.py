import numpy as np
from ..channel import *
def BLER(myqec,monte=1000,err_stop=1000):
    result = []
    d =10
    input_BER = [0.5,0.4,0.3,0.2,0.1,0.05,0.01,0.05,0.01,0.005,0.001,0.0005,0.0001]
    #input_BER = [1/3]
    n = myqec.n
    for p in input_BER:
        ble = 0
        myqec.set_P(np.array([1-p,p/3,p/3,p/3]))
        for mc in range(1,monte+1):
            #e = one_bit_channel(n,1)
            e = depolarizing_noise(n,p)
            syndrome = myqec.get_syndrome(e)
            ee = myqec.decode(syndrome)
            #print(e,ee)
            for i in range(len(myqec.L)):
                if not myqec.in_S(e^ee):
                    ble+=1
                    break
            print("p",p,"誤り確率:",mc,ble/mc)
            if ble>err_stop:
                break
        result.append(ble/mc)
        if ble/mc==0:
            break
    print(input_BER)
    print(result)
    return input_BER,result
