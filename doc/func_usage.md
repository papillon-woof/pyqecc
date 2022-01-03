# 基本的な使用方法

### 1. 符号インスタンスを自由に作る．  
### 2. 復号シミュレーションにインスタンスを渡す．

# 仕様
## 符号クラス
### abstruct.py
抽象クラス．
```
class CODE(metaclass=ABCMeta):
```

``` 
def __init__(self,**param):
```
コンストラクタ．符号長，情報長，符号化率が少なくとも求められる．

### Stabilizer.py: 
```
class SC(CODE)
```
Stabilizer符号が定義されている．
```
def __init__(self,n,k,H='random',T=None,L=None,P=None,iid=True)
```

`n :int`  
  符号長  
`k :int`  
　情報長  
`R :float`  
  符号化率  
`H :np.array or string`  
  検査行列自身，あるいは検査行列を指定する文字列．現状未実装  
`T :np.array (int) or None`  
  シンドロームに対応する回復演算子．シンドロームと一対一対応する必要がある．  
`P :np.array (float) or None`  
  通信路情報．  
`iid :boolean`
　確率分布が独立同一分布かどうか．  

```
def get_syndrome(self,e):
```
パウリ群`G`に属する誤り`E`の二値表現ベクトル`e`を受け取り，シンドローム`s`を計算します．
`e: np.array(2*n,dtype='i8')`
  符号長`n`としたときの`E`の二値表現ベクトル  

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
