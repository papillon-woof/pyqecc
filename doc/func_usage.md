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
