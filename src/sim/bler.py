import numpy as np
def BLER(myqec,monte=10000):
    result = []
    d =10
    input_BER = [1/(2**d)*(2**i) for i in range(d)]
    #input_BER = [1/3]
    n = myqec.n
    for p in input_BER:
        ble = 0
        myqec.set_P(np.array([1-p,p/3,p/3,p/3]));
        for mc in range(1,monte+1):
            e = depolarizing_noise(n,p)
            syndrome = myqec.get_syndrome(e)
            ee = myqec.decode(syndrome,mode='ML');
            for i in range(len(myqec.L)):
                if not myqec.in_S(e^ee):
                    ble+=1
                    break
            print("誤り確率:",mc,ble/mc)
        result.append(ble/monte)
    print(input_BER)
    print(result)
    return input_BER,result
