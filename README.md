<div align="center">
<img src="https://user-images.githubusercontent.com/72004949/148188473-22ea4600-d1d1-46b6-814b-0d3414af5750.png" alt="pyqecc" title="pyqecc">
</div>

# Overview
PyQecc mainly provide quantum error correction code (QECC) simulator.
- [Installation](https://pyqecc.readthedocs.io/en/latest/?) (This page)
- [Quick start](https://pyqecc.readthedocs.io/en/latest/?) (This page)
- [Features](https://pyqecc.readthedocs.io/en/latest/features.html)
- [Source code](https://github.com/papillon-woof/pyqecc)
# Installation

```
pip install pyqecc
```
In some cases, installation required setting for `--proxy`, `--user` or `sudo`. PyQecc is written by python3

# Quick start
We explane the tutorial usage. 
Please prepare the `.py` file (e.g. `test.py`). Please copy and paste following code
```python
from pyqecc import *
my_code = FIVE()
print(my_code)
dec_sim(my_code)
```
The steps of evalutation for decoding performance are `Import the PyQecc.`, `Create the instance for QECC.`, `Prepare the decoding simulator`, `Start the decoding simulation`, and `Confirm the decoding result`.

## 1. Import the PyQecc.
```python
from pyqec import *
```
## 2. Create the instance for QECC.
For example, we prepare the 5-qubit code.
```python
my_code = FIVE()
```
We confirm the information for QECC `my_code` by
```python
print(my_code)
```
```
NAME            :FIVE_CODE
N               : 5
K               : 1
R               : 0.2
DECODING_MODE   : ML_LUT
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

```
[decoding result]
```

directory structure
```
├── test.py
└── dec_data (Folder)
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

## Decoding simulation
- block error rate

## Channel Model
- depolarizing channel

# Simulation example
Concatenated 5-qubit codes (concatenation for 1, 2, and 3) [2, Fig. 1].  
![image](https://user-images.githubusercontent.com/72004949/148180717-3c523204-3acc-48c6-a736-503b14dece4e.png)
```python
#Source code
from pyqecc import *
NUM_OF_CONCATENATE = 3
for num_of_concatenate in range(1,NUM_OF_CONCATENATE+1):
    conc_code = [FiveCode()]
    for i in range(1,num_of_concatenate):
        conc_code += [ParaCode([FiveCode() for i in range(5 ** i)])]
    my_code = ConcCode(conc_code)
    print(my_code)
    dec_sim(my_code,PROB=[0.13, 0.15, 0.17, 0.18, 0.1885, 0.19],MONTE=5000)
```

# Future works
- quantum LDPC code
- quantum polar code
- surface code
- pauli channel
- amplitude damping channel

# References
[1] Nielsen, Michael A., and Isaac Chuang. "Quantum computation and quantum information." (2002): 558-559.

[2] Poulin, David. "Optimal and efficient decoding of concatenated quantum block codes." Physical Review A 74.5 (2006): 052333.
