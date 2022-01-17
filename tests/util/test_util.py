import numpy as np
import pytest
from pyqecc.util import *

DIFF_THRESHOLD = 1e-4
TEST_DATA_SYMPLEX_BINARY_INNER_PRODUCT = {
    "NUM_OF_CASE": 4,
    "RESULT": [1, 0, 0, 1],
    "TEST_CASE": np.array(
        [
            [[0, 0, 1, 0], [1, 0, 0, 0]],
            [[1, 0, 0, 0], [1, 0, 0, 0]],
            [[1, 0, 1, 0], [1, 0, 1, 0]],
            [[1, 1, 1, 0], [1, 0, 0, 0]],
        ],
        dtype="i1",
    ),
}

TEST_ARR_2_INT = {
    "NUM_OF_CASE": 16,
    "TEST_CASE": np.array(
        [
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [0, 0, 1, 1],
            [0, 1, 0, 0],
            [0, 1, 0, 1],
            [0, 1, 1, 0],
            [0, 1, 1, 1],
            [1, 0, 0, 0],
            [1, 0, 0, 1],
            [1, 0, 1, 0],
            [1, 0, 1, 1],
            [1, 1, 0, 0],
            [1, 1, 0, 1],
            [1, 1, 1, 0],
            [1, 1, 1, 1],
        ],
        dtype="i1",
    ),
}

TEST_BITWISE_TO_BLOCKWISE_PROBABILITY = {
    "NUM_OF_CASE": 6,
    "TEST_CASE": [1, 2, 3, 4, 5, 6, 7],
}


def test_symplex_binary_inner_product():
    for i in range(TEST_DATA_SYMPLEX_BINARY_INNER_PRODUCT["NUM_OF_CASE"]):
        assert TEST_DATA_SYMPLEX_BINARY_INNER_PRODUCT["RESULT"][
            i
        ] == symplex_binary_inner_product(
            TEST_DATA_SYMPLEX_BINARY_INNER_PRODUCT["TEST_CASE"][i][0],
            TEST_DATA_SYMPLEX_BINARY_INNER_PRODUCT["TEST_CASE"][i][1],
        )


def test_arr2int():
    for i in range(TEST_ARR_2_INT["NUM_OF_CASE"]):
        assert i == arr2int(TEST_ARR_2_INT["TEST_CASE"][i])


def test_int2arr():
    m = 5
    for i in range(2 ** m):
        assert i == arr2int(int2arr(i, m))


def test_bitwise_to_blockwise_probability():
    for n in TEST_BITWISE_TO_BLOCKWISE_PROBABILITY["TEST_CASE"]:
        p = np.random.rand(4 * n).reshape(n, 4)
        for i in range(n):
            p[i] /= sum(p[i])
        assert DIFF_THRESHOLD > np.sum(
            np.abs(
                blockwise_to_bitwise_error_probability(
                    bitwise_to_blockwise_error_probability(p)
                )
                - p
            )
        )
