import os
from datetime import datetime
import numpy as np
from ..channel import *



def dec_sim(
    myQECC,
    MONTE=1000,
    ERR_STOP=1000,
    channel_instance = None,
    DEBUG=False,
    LABEL=["CHANNEL_PARAMATER", "LOGICAL_ERROR_PROB"],
    LOG_OUTPUT=True,
    LOG_OUTPUT_SPAN=100, 
):
    RESULTS = {}
    RESULTS["LOGICAL_ERROR_PROB"] = []
    
    n = myQECC.n
    if channel_instance is None:
        channel_instance = DepolarizingChannel(p = [0.1, 0.01, 0.001, 0.0001])
    RESULTS["CHANNEL_PARAMETER"] = channel_instance.channel_parameter
    for ind in range(len(channel_instance.channel_parameter["PARAM_SET"])):
        ble = 0
        myQECC.set_channel_param(channel_instance.channel_parameter["PARAM_SET"][ind])
        for mc in range(1, MONTE + 1):
            channel_output = channel_instance.channel(n,ind)
            syndrome = myQECC.get_syndrome(channel_output)
            if not myQECC.in_S(channel_output["E"] ^ myQECC.decode(syndrome)["LT"]):
                ble += 1
            if not mc % LOG_OUTPUT_SPAN:
                print(
                    "MONTE",
                    mc,
                    "BLOCK_ERROR",
                    ble,
                    channel_instance.channel_parameter["PARAM_SET"][ind],
                    "LOGICAL_ERROR_PROB:",
                    mc,
                    ble / mc,
                )
            if ble > ERR_STOP:
                break
        RESULTS["LOGICAL_ERROR_PROB"].append(ble / mc)
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
