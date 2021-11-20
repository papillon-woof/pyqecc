import numpy as np
from ..channel import *
def BLER(myQECC,MONTE=1000,ERR_STOP=1000,PHYSICAL_ERROR_PROB=[0.5,0.4,0.3,0.2,0.1,0.05,0.01,0.05,0.01,0.005,0.001,0.0005,0.0001],DEBUG=False):
    RESULTS = {}
    RESULTS['LOGICAL_ERROR_PROB'] = []
    RESULTS['PHYSICAL_ERROR_PROB'] = PHYSICAL_ERROR_PROB
    n = myQECC.n
    for p in PHYSICAL_ERROR_PROB:
        ble = 0
        myQECC.set_P(np.array([1-p,p/3,p/3,p/3]))
        for mc in range(1,MONTE+1):
            E = depolarizing_noise(n,p)
            syndrome = myQECC.get_syndrome(E)
            EE = myQECC.decode(syndrome)
            for i in range(len(myQECC.L)):
                if not myQECC.in_S(E^EE):
                    ble+=1
                    break
            print("p",p,"誤り確率:",mc,ble/mc)
            if ble>ERR_STOP:
                break
        RESULTS['BLER'].append(ble/mc)
        if ble/mc==0:
            break
    if DEBUG:
        print('PHYSICAL_ERROR_PROB',PHYSICAL_ERROR_PROB)
        print('LOGICAL_ERROR_PROB',RESULTS['LOGICAL_ERROR_PROB'])
    return RESULTS
