<div align="center">
<img src="https://user-images.githubusercontent.com/72004949/148188473-22ea4600-d1d1-46b6-814b-0d3414af5750.png" alt="pyqecc" title="pyqecc">
</div>

[![Documentation Status](https://readthedocs.org/projects/pyqecc/badge/?version=latest)](https://pyqecc.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads](https://pepy.tech/badge/pyqecc)](https://pepy.tech/project/pyqecc)
[![Downloads](https://pepy.tech/badge/pyqecc/month)](https://pepy.tech/project/pyqecc)
[![Downloads](https://pepy.tech/badge/pyqecc/week)](https://pepy.tech/project/pyqecc)
# Overview
PyQecc mainly provide quantum error correction code (QECC) simulator.
- [Installation](https://pyqecc.readthedocs.io/en/latest/?) (This page)
- [Quick start](https://pyqecc.readthedocs.io/en/latest/?) (This page)
- [Features](https://pyqecc.readthedocs.io/en/latest/features.html)
- [Source code](https://github.com/papillon-woof/pyqecc)
- [PyPI](https://pypi.org/project/pyqecc/0.0.1/)
# Installation

```
pip install pyqecc
```
In some cases, installation required setting for `--proxy`, `--user` or `sudo`. PyQecc is written by python3.

# Quick start
We explane the tutorial usage. 
Please prepare the new `.py` file (e.g. `test.py`). Please copy and paste following code
```python
from pyqecc import FiveCode, dec_sim

my_code = FiveCode(mode="ML")
print(my_code)
dec_sim(my_code)
```
The steps of evalutation for decoding performance are `Import the PyQecc.`, `Create the instance for QECC.`, `Prepare the decoding simulator`, `Start the decoding simulation`, and `Confirm the decoding result`.

## 1. Import the PyQecc.
```python
from pyqecc import FiveCode, dec_sim
```
## 2. Create the instance for QECC.
For example, we prepare the 5-qubit code.
```python
my_code = FiveCode()
```
We confirm the information for QECC `my_code` by
```python
print(my_code)
```
```
NAME               : FIVE_CODE
PHYSICAL QUBITS (n): 5
LOGICAL QUBITS (k) : 1
CODE RATE (R = k/n): 0.2
DECODING_MODE      : ML
```

## 3. Prepare the decoding simulator.
```python
dec_sim(my_code)
```
default settings:
- depolarizing channel
- 1000 codeward
- maximum likelihood decoding.

## 4. Start the decoding simulation. 
```
python test.py
```
Please wait patiently. 
## 5. Confirm the decoding results.
In `/dec_data`, PyQecc generates the simulation results.

```console
MONTE: 100 BLOCK_ERROR: 5 p: 0.1 LOGICAL_ERROR_PROB: 0.05
MONTE: 200 BLOCK_ERROR: 19 p: 0.1 LOGICAL_ERROR_PROB: 0.095
MONTE: 300 BLOCK_ERROR: 24 p: 0.1 LOGICAL_ERROR_PROB: 0.08
MONTE: 400 BLOCK_ERROR: 29 p: 0.1 LOGICAL_ERROR_PROB: 0.0725
MONTE: 500 BLOCK_ERROR: 35 p: 0.1 LOGICAL_ERROR_PROB: 0.07
MONTE: 600 BLOCK_ERROR: 44 p: 0.1 LOGICAL_ERROR_PROB: 0.07333333333333333
MONTE: 700 BLOCK_ERROR: 53 p: 0.1 LOGICAL_ERROR_PROB: 0.07571428571428572
MONTE: 800 BLOCK_ERROR: 63 p: 0.1 LOGICAL_ERROR_PROB: 0.07875
MONTE: 900 BLOCK_ERROR: 71 p: 0.1 LOGICAL_ERROR_PROB: 0.07888888888888888
...
```

directory structure
```
├── test.py
└── dec_data (Folder)
```

`dec_data/FIVE_CODE_5_1_monte_1000_20220220212203.csv`
```
p,LOGICAL_ERROR_PROB,
0.1,0.08,
0.01,0.0,
```

# Features
See the detail for [features](features.md)

## Stabilizer Code
- 5-qubit code
- 7-qubit code (STEANE code)
- bit flip code
- phase flip code
- 9-qubit shor code (concatenated bit and phase flip code.)
- concatenated code

## decoder
- syndrome decoding
- maximum likelihood (ML) decoding
- belief propagation decoding (concatenated code only)
- decoding with analog informatuon [3] (GKP qubit only)

## Decoding simulation
- block error rate

## Channel Model
- depolarizing channel
- pauli channel
- bit flip channel
- quantum gaussian channel

# Simulation example
Concatenated 5-qubit codes (concatenation for 1, 2, and 3) [2, Fig. 1].  
![image](https://user-images.githubusercontent.com/72004949/148180717-3c523204-3acc-48c6-a736-503b14dece4e.png)
```python
#Source code
from pyqecc import FiveCode, dec_sim, ConcCode, ParaCode, DepolarizingChannel

NUM_OF_CONCATENATE = 3
for num_of_concatenate in range(1, NUM_OF_CONCATENATE + 1):
    conc_code = [FiveCode()]
    for i in range(1, num_of_concatenate):
        conc_code += [ParaCode([FiveCode() for i in range(5**i)])]
    my_code = ConcCode(conc_code)
    print(my_code)
    dec_sim(
        my_code,
        channel_instance=DepolarizingChannel(
            my_code.n, p=[0.13, 0.15, 0.17, 0.18, 0.1885, 0.19]
        ),
        MONTE=5000,
    )
```

# Future works
- surface code
- color code
- quantum LDPC code
- quantum polar code

# References
[1] Nielsen, Michael A., and Isaac Chuang. "Quantum computation and quantum information." (2002): 558-559.

[2] Poulin, David. "Optimal and efficient decoding of concatenated quantum block codes." Physical Review A 74.5 (2006): 052333.

[3] Kosuke Fukui, Akihisa Tomita, and Atsushi Okamoto Phys. Rev. Lett. 119, 180507 – Published 3 November 2017
