import numpy as np
import pytest
from pyqecc.qecc import *

TEST_BOSONIC_BASIC = {
    "TEST_CASE": 4,
    "CODE_INSTANCE": [
        BitFlipCode(),
        PhaseFlipCode(),
        SteaneCode(),
        FiveCode(),
        ConcCode(
            [PhaseFlipCode(), ParaCode([BitFlipCode(), BitFlipCode(), BitFlipCode()])]
        ),
    ],
}

TEST_BOSONIC_CALC_LLR = {
    "TEST_CASE": 1,
    "CODE_INSTANCE": BitFlipCode(),
    "PARAM": (0, 1, ((np.sqrt(np.pi)) ** 2) / (2 * 1**2)),
}


def test_bosonic_basic():
    # Preparing test
    for i in range(TEST_BOSONIC_BASIC["TEST_CASE"]):
        assert GKP(TEST_BOSONIC_BASIC["CODE_INSTANCE"][i])

    for i in range(TEST_BOSONIC_BASIC["TEST_CASE"]):
        my_gkp = GKP(TEST_BOSONIC_BASIC["CODE_INSTANCE"][i])
        channel_output = {"E": np.zeros(2 * my_gkp.n, dtype="i1")}
        assert sum(my_gkp.get_syndrome(channel_output)) == 0
        channel_output["shift_error"] = np.ones(2 * my_gkp.n)
        assert sum(my_gkp.get_syndrome(channel_output)[0]) == 0


def test_bosonic_calc_llr():
    my_gkp = GKP(TEST_BOSONIC_CALC_LLR["CODE_INSTANCE"])
    for i in TEST_BOSONIC_CALC_LLR["PARAM"]:
        assert (
            my_gkp.calc_llr(
                TEST_BOSONIC_CALC_LLR["PARAM"][0], TEST_BOSONIC_CALC_LLR["PARAM"][1]
            )
            == TEST_BOSONIC_CALC_LLR["PARAM"][2]
        )
