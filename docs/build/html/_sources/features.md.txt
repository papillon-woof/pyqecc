# Feature

# `pyqecc.sim.dec_sim`
```python
def dec_sim(
    myQECC,
    MONTE=1000,
    ERR_STOP=1000,
    channel_instance=None,
    DEBUG=False,
    LOG_OUTPUT=True,
    LOG_OUTPUT_SPAN=100,
):
```

This function `pyqecc.sim.dec_sim` is numerical simulation for quantum error correction code.   
The `dec_sim` requires arguments for
- `myQECC`: `CODE` class instance (e.g. `SC`, `BitFlipCode`, `GKP`, ...)
and has the optional arguments for 
- `channel_instance`: `Channel` class instance (e.g. `DepolarizingChannel`,...) (default: `DepolarizingChannel(n, p=[0.1, 0.01, 0.001, 0.0001])`)
- `MONTE`: number of monte carlo (default: `1000`)
- `ERR_STOP`: Threshold of error count for stopping simulation (default: `1000`)
- `DEBUG`: Debug mode (default: `False`)
- `LOG_OUTPUT`: Whether log output console (default: `True`)
- `LOG_OUTPUT_SPAN`: log output span of monte calro (default: `100`)

# `pyqecc.qecc.toy.SteaneCode()`
Return the instance of steane code

# `pyqecc.qecc.toy.FiveCode()`
Return the instance of 5-qubits code

# `pyqecc.qecc.toy.ShorCode()`
Return the instance of shor code

# `pyqecc.qecc.toy.BitFlipCode()`
Return the instance of bit flip code

# `pyqecc.qecc.toy.PhaseFlipCode()`
Return the instance of phase flip code

# `pyqecc.qecc.stabilizer.SC(CODE)`
Stabilizer Code class

**Member**  

`n :int`  
Code length (Number of physical qubit)  
`k :int`  
Information length (Number of logical qubit)  
`nk :int`  
n - k  
`R :float`  
Code rate (R:=k/n)  
`H :np.array((nk,2*n),dtype=np.int8`)`  
Parity check matrix. For example `bit flip()` is defined by
```python
np.array([
    [0,0,0,1,1,0]
    [0,0,0,0,1,1]
])
```
the first `n` bits and others correspond to Pauli `X` and `Z` in Stabilizer respectively.

`T :np.array (int) or None`  
  Recovery operator.  
`L :np.array (int) or None`  
  Logical operator.  
`blockwise_error_probability :np.array(2**(2*n)) or None`  
  blockwise channel error probability．  
`bitwise_error_probability :np.array((n,4)) or None`  
  bitwise channel error probability．  
`iid :boolean`
　Whether channel error probability is iid   

**Constructer**  
Setup for the stabilizer need to `n`,`k` and `H`. But, in most cases, you need to set `T` and `L`.
The `P` is channel probability. if `BIT_WISE` is true, it calculates the bitwise error probability. In contrast, if `BIT_WISE` is false, it calculates the block wise probability.
```python
__init__(
        self,
        n,
        k,
        H = None,
        T = None,
        L = None,
        P = None,
        mode = 'HD',
        BITWISE = True,
        ):
```

**Method**  
- `set_channel_param(self,error_probability,BITWISE=True,iid=True,OUTPUT_LOG=False):`: Stores channel probability distribution information `error_probability`. `OUTPUT_LOG` can confirm whether `error_probability` is a varid value.
- `get_syndrome(self,e)`: It receives the binary representation vector `e` of the error`E` belonging to the Pauli group `G` and calculates the syndrome`syndrome`. The `PROB` is list of depolarizing probability.
- `get_T(self,idx)`: Given the index or array `idx`, it returns the corresponding error operator`T`.
- `get_L(self,idx)`: Given the index or array `idx`, it returns the corresponding error operator`L`.
- `get_S(self,idx)`: Given the index or array `idx`, it returns the corresponding error operator`S`.
- `in_S(self,b)`: Given the binary representation of the Pauli group `b`, it is determined whether the element` b` is included in the fixed subgroup.
- `decode(self,syndrome,mode='hd')`: Syndrome (Hard decision) decoding or maximum likelihood decoding is performed according to `mode`, and the recovery operator `T` (or `LT`) and likelihood probability are returned in` dict` type.
  - mode=`HD`: syndrome decoding
  - mode=`ML_LUT`: ML decoding (using lookup table)
  - mode=`ML`: ML decoding (calculate logical error probability `L` and `L=argmax(logical_error_probability)`, `L^T` is returned)

# `pyqecc.qecc.ParaCode(SC)`
`code_instances`: The code instances (Write the `[Code1,Code2,...]`.)

# `pyqecc.qecc.ConcCode(SC)`
Concatenated code class.

**constructor**  
```python
def __init__(self,code_instances,P=None,BITWISE = True,mode="BP"):
```

**Member**  
`code_instances`: The code instances (Write the `[Code1,Code2,...]`.)

**example: `[FiveCode(),Para([FiveCode(),FiveCode(),FiveCode(),FiveCode(),FiveCode()])]` is concatenated 5-qubit code, which is a code in which one qubit is encoded into five qubits by "five-qubit code", are encoded  into twenty five qubits by five "five-qubit code"s, one for each.** 

**method**  
`calc_grobal_H(self)`: calculation global parity check matrix

# `pyqecc.qecc.GKP(CODE)`
GKP Qubit class
`
def __init__(
        self,
        code_instance,
        sigma=0.1,
        mode="SYNDROME",
    ):
`
The `GKP` class has the arguments for 
- `code_instance`: code instance for `CODE` class
- `sigma=0.1`: Standard deviation for Quantum Gaussian Channel
- `mode="SYNDROME"`: decoding mode

# `pyqecc.qecc.abstruct.CODE()`
```python
class CODE(metaclass=ABCMeta):
```
The all of code class inherits this class.
**Constructor**  
```python
def __init__(self,**param):
```
`n :int`  
code length (number of physical qubit)  
`k :int`  
information length (number of logical qubit)  
`R :float`  
code rate (R:=k/n)  
