import os
from datetime import datetime
import numpy as np
from ..channel import *


def dec_sim(
    myQECC,
    MONTE=1000,
    ERR_STOP=1000,
    channel_instance=None,
    DEBUG=False,
    LOG_OUTPUT=True,
    LOG_OUTPUT_SPAN=100,
):
    RESULTS = {}

    n = myQECC.n
    if channel_instance is None:
        channel_instance = DepolarizingChannel(n, p=[0.1, 0.01, 0.001, 0.0001])
    for k, v in channel_instance.channel_parameter.items():
        RESULTS[k] = v
    for ind in range(channel_instance.param_len):
        ble = 0
        myQECC.set_channel_param(channel_instance.get_param(ind))
        for mc in range(1, MONTE + 1):
            channel_output = channel_instance.channel(ind=ind)
            syndrome = myQECC.get_syndrome(channel_output)
            if not myQECC.in_S(channel_output["E"] ^ myQECC.decode(syndrome)["LT"]):
                ble += 1
            if not mc % LOG_OUTPUT_SPAN:
                text = "MONTE: " + str(mc) + " BLOCK_ERROR: " + str(ble) + " "
                for k, v in channel_instance.channel_parameter.items():
                    text += k + ": " + str(v[ind]) + " "
                text += "LOGICAL_ERROR_PROB: " + str(ble / mc)
                print(text)

            if ble > ERR_STOP:
                break
        if "LOGICAL_ERROR_PROB" not in RESULTS:
            RESULTS["LOGICAL_ERROR_PROB"] = []
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
            for label in RESULTS.keys():
                f.write(str(label) + ",")
            f.write("\n")
            for i in range(len(RESULTS["LOGICAL_ERROR_PROB"])):
                for label in RESULTS.keys():
                    val = RESULTS[label][i]
                    f.write(str(val) + ",")
                f.write("\n")
    if DEBUG:
        print(RESULTS)
    return RESULTS
