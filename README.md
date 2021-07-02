# What is pyQec
pyQec is quantum error correction code (QECC) library. This library mainly provide a decoding algorithm for QECC.  
The QECC is needed to overcome the quantum noise arising from time evolution, which cannot avoid in principle, the imperfections quantum circuit and channel.   

# Contents
## Quantum Error Correction Code
### basic code and toys
- Stabilizer Code
- Toy code (5-q code, 7-q(STEANE) code)

### Construction of QECC from classical code
- (Future) quantum LDPC code
- (Future) quantum polar code
- concatenated code

### Construction QECC from classical physics
- Topological code (Toric, colar)

## decoder
- syndrome decoding (stabilizer)
- renormalize group decoding (toric code)
- SC decoder (SC decoder)

## Channel Model
- depolarizing channel
- (Future) Pauli channel
