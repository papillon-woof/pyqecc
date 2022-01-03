<div align="center">
<img src="https://user-images.githubusercontent.com/72004949/147950760-7b073f0b-1efb-4d65-a8db-b347c0e115e0.png" alt="pyqec" title="pyqec">
</div>

**Warning: Pre-alpha version for PyQec is provided. It may have a large bugs.**

PyQec mainly provide quantum error correction code (QECC) simulator. 

# Installation

```
pip install pyqec
```
In some cases, installation required setting for `--proxy`, `--user` or `sudo`. PyQec is written by python3

# Usage
We explane the tutorial usage. 
Please prepare the `.py` file (e.g. `test.py`). Please copy and paste following code
```python
from pyqec import *
my_code = FIVE()
print(my_code)
dec_sim(my_code)
```
The steps of evalutation for decoding performance are `Import the PyQec.`, `Create the instance for QECC.`, `Prepare the decoding simulator`, `Start the decoding simulation`, and `Confirm the decoding result`.

## 1. Import the PyQec.
```python
from pyqec import *
```
## 2. Create the instance for QECC. For example, we prepare the 5-qubit code.
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

## 3. Prepare the decoding simulator with `my_code`
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
## 5. Confirm the decoding result
In `/decdata`, PyQec generates the simulation results.

```
[decoding result]
```

# Features
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

# Future work
- quantum LDPC code
- quantum polar code
- surface code
- (Future) Pauli channel
- (Future) Amplitude Damping channel
