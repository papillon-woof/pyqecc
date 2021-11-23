import numpy as np
from ..channel import *
import os
from datetime import datetime
def dec_sim(myQECC,MONTE=1000,ERR_STOP=500,DEPOLARIZING_PROB=[0.01,0.005,0.001,0.0005,0.0001],CHANNEL=None,DEBUG=False,LABEL=['DEPOLARIZING_PROB','PHYSICAL_ERROR_PROB','LOGICAL_ERROR_PROB'],LOG_OUTPUT=True):
    RESULTS = {}
    RESULTS['LOGICAL_ERROR_PROB'] = []
    RESULTS['PHYSICAL_ERROR_PROB'] = []
    RESULTS['DEPOLARIZING_PROB'] = DEPOLARIZING_PROB
    #RESULTS['PHYSICAL_ERROR_PROB'] = PHYSICAL_ERROR_PROB
    n = myQECC.n
    phisical_error_count = 0
    for p in DEPOLARIZING_PROB:
        ble = 0
        myQECC.set_P(np.array([1-p,p/3,p/3,p/3]))
        for mc in range(1,MONTE+1):
            E = depolarizing_noise(n,p)
            #E = one_bit_noise(n,p)
            phisical_error_count += sum(E)
            syndrome = myQECC.get_syndrome(E)
            EE = myQECC.decode(syndrome)
            for i in range(len(myQECC.L)):
                if not myQECC.in_S(E^EE):
                    ble+=1
                    break
            #print("DEPOLARIZING_ERROR_PROB",p,"PHYSICAL_ERROR_PROB:",phisical_error_count/(mc*myQECC.n)," LOGICAL_ERROR_PROB:",mc,ble/mc)
            print("DEPOLARIZING_ERROR_PROB",p,"PHYSICAL_ERROR_PROB:",p," LOGICAL_ERROR_PROB:",mc,ble/mc)
            if ble>ERR_STOP:
                break
        RESULTS['LOGICAL_ERROR_PROB'].append(ble/mc)
        RESULTS['PHYSICAL_ERROR_PROB'].append(p)
        if ble/mc==0:
            break
    if DEBUG:
        print("DEPOLARIZING_ERROR_PROB",DEPOLARIZING_PROB)
        print('PHYSICAL_ERROR_PROB',p)
        print('LOGICAL_ERROR_PROB',RESULTS['LOGICAL_ERROR_PROB'])
    if LOG_OUTPUT:
        DEC_DATA_DIR = "./dec_data"
        if not os.path.exists(DEC_DATA_DIR):
            # ディレクトリが存在しない場合、ディレクトリを作成する(Windows 用?)
            os.makedirs(DEC_DATA_DIR)
        #データ書き込み
        #データ名
        path_w = myQECC.name+"_"+str(myQECC.n)+"_"+str(myQECC.k)+"_monte_"+str(MONTE)+"_"+datetime.now().strftime('%Y%m%d%H%M%S')+".csv"
        with open('dec_data/'+path_w, mode='w') as f:
            # paramを書き込む．
            for label in LABEL:
                if label in RESULTS.keys():
                    f.write(str(label)+","+str(RESULTS[label])+"\n")
    return RESULTS
