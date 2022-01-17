import numpy as np
from ..channel import *
import os
from datetime import datetime


def dec_sim(
    myQECC,
    MONTE=1000,
    ERR_STOP=1000,
    PROB=[0.1, 0.01, 0.001, 0.0001],
    CHANNEL_MODEL="DEPOLARIZING",
    DEBUG=False,
    LABEL=["DEPOLARIZING_PROB", "PHYSICAL_ERROR_PROB", "LOGICAL_ERROR_PROB"],
    LOG_OUTPUT=True,
    LOG_OUTPUT_SPAN=100,
):
    RESULTS = {}
    RESULTS["LOGICAL_ERROR_PROB"] = []
    RESULTS["PHYSICAL_ERROR_PROB"] = []
    RESULTS["DEPOLARIZING_PROB"] = PROB
    n = myQECC.n
    for p in PROB:
        ble = 0
        myQECC.set_error_probability(np.array([1 - p, p / 3, p / 3, p / 3]), iid=True)
        for mc in range(1, MONTE + 1):
            E = channel(n, p, CHANNEL_MODEL=CHANNEL_MODEL)
            syndrome = myQECC.get_syndrome(E)
            EE = myQECC.decode(syndrome)["LT"]
            # print(E,EE)
            if not myQECC.in_S(E ^ EE):
                ble += 1
            if not mc % LOG_OUTPUT_SPAN:
                print(
                    "MONTE",
                    mc,
                    "BLE",
                    ble,
                    "DEPOLARIZING_ERROR_PROB",
                    p,
                    " LOGICAL_ERROR_PROB:",
                    mc,
                    ble / mc,
                )
            if ble > ERR_STOP:
                break
        RESULTS["LOGICAL_ERROR_PROB"].append(ble / mc)
        RESULTS["PHYSICAL_ERROR_PROB"].append(p)
        if ble / mc == 0:
            break

    if LOG_OUTPUT:
        DEC_DATA_DIR = "./dec_data"
        if not os.path.exists(DEC_DATA_DIR):
            os.makedirs(DEC_DATA_DIR)
        path_w = (
            myQECC.name
            + "_"
            + str(myQECC.n)
            + "_"
            + str(myQECC.k)
            + "_monte_"
            + str(MONTE)
            + "_"
            + datetime.now().strftime("%Y%m%d%H%M%S")
            + ".csv"
        )
        with open("dec_data/" + path_w, mode="w") as f:
            for label in LABEL:
                if label in RESULTS.keys():
                    f.write(str(label) + "," + str(RESULTS[label]) + "\n")
    if DEBUG:
        print(RESULTS)
    return RESULTS
