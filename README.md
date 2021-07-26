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

# Explanation Class and Function
## `pyQec.fec.abstruct.CODE()`
`pyQec.fec.CODE`は，量子誤り訂正符号クラスの抽象クラスです．下記のメンバおよびメソッドを持ちます．
### メソッド
- `__init__(self)`: コンストラクタ．メンバの初期値を設定する．
- `__str__(self)`: パラメータ標準出力用関数．メンバを出力する．
### メンバ
- enc_circ: 符号化回路(@property)．(未実装)
- `n`: 符号長(@property)
- `k`: 情報長(@property)
- `R`: 符号化率(@property)
- `name`: クラス名(@property)
-
## `pyQec.fec.stabilizer.SC(CODE)`
`pyQec.fec.SC(CODE)`は，Stabilizer符号のクラスです．クラス`CODE`を継承します．
### メソッド
- `set_P(self,P,iid=True)`: 通信路確率分布情報`P`を格納します．
- `get_syndrome(self,e)`: 誤り`e`が与えられた際のシンドロームを返します．
- `get_T(self,idx)`:インデックス`idx`が与えられると，対応する誤り演算子`T`を返します．
- `get_L(self,idx)`:インデックス`idx`が与えられると，対応する論理演算子`L`を返します．
- `get_S(self,idx)`:インデックス`idx`が与えられると，対応する固定部分群の元`S`を返します．
- `in_S(self,b)`: パウリ群の元の二値表現`b`が与えられると，その元`b`が固定部分群に含まれるか判定します．
- `hard_decode(self,syndrome)`: シンドローム`syndrome`が与えられると，硬判定復号(シンドローム復号)を行い，訂正のための演算子を返します．
- `ML_decode(self,syndrome)`: シンドローム`syndrome`が与えられると，メンバ`P`を元に最尤復号を行います．
- `decode(self,syndrome,mode='hd')`: シンドローム`syndrome`が与えられると，`mode`に合わせて硬判定復号もしくは最尤復号を行います．

### メンバ
- `P`: 通信路情報(@property)
- `L`: 論理演算子(@property)
- `T`: 誤り演算子(@property)
- `H`: 検査行列(@property)

## `pyQec.fec.toy.STEANE()`
STEANE符号[7,1]Stabilizer符号インスタンスを返します．

## `pyQec.fec.toy.FIVE()`
FIVE符号の[5,1]Stabilizer符号インスタンスを返します．

## `pyQec.util.util`
### 関数
- `symplex_binary_inner_product(a,b)`: 二値ベクトルのsympletics内積$\braket{a,b}_{symp}$を返します．
- `gaussjordan(X)`: `X`の上三角行列表現を返します．
-
