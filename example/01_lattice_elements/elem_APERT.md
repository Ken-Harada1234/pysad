# elem_APERT (開口)
種別: elem。フォルダ: 01_lattice_elements。ラティスは `FODO_ring.sad` と同じ光学のリング (定義はファイル内、ORIGIN の次に `APERT AP` を 1 個置いた。RF なし)。エネルギー 3 GeV。

## 何を見る例か
`APERT AP = (AX=0.01 AY=0.005)`: **追跡 (TrackParticles) だけで効く**楕円の開口。粒子が $\frac{(x-DX)^2}{AX^2}+\frac{(y-DY)^2}{AY^2}>1$ の外に出ると失われ、`flag` が 0 になる (`DX1 DX2 DY1 DY2` を使う矩形も指定できる。SADHelp の APERT の項)。
x=1 mm (内側) と x=20 mm (外側。AX=10 mm より外) の 2 粒子を 1 周追跡して flag を見る。

## スクリプト
```
APERT AP = (AX=0.01 AY=0.005);
LINE RING = (ORIGIN AP 24*CELL);
...
beam = {1, {{0.001,0.02}, {0,0}, {0.001,0.001}, {0,0}, {0,0}, {0,0}, {1,1}}};
r = TrackParticles[beam, LINE["LENGTH"], 1, 1];      (* 1 周 *)
Print["1st call: flag = ", r[[2,7]], "  x = ", r[[2,1]]];
```
`TrackParticles[{ls, {x,px,y,py,z,dp,flag}}, 到達位置, 周回の始め, 終わり]`。粒子ごとの列の並びは 6 座標 + flag。K64 は最初の `TrackParticles` では最初から開口の外の粒子を落とさない癖があるため (docs/11 の A52)、同じ呼び出しを 2 回している。

## 結果の読み方
`flag` が 1 なら生存、0 なら失われた粒子。失われた粒子の `x` は落ちた時点の値で止まる。

| 呼び出し | gs の flag | pysad の flag |
|---|---|---|
| 1 回目 | flag = [1.0, 0.0]、x = [0.001715, 0.02] | flag = [1.0, 1.0]、x = [0.001715, 0.037768] |
| 2 回目 | flag = [1.0, 0.0]、x = [0.001715, 0.02] | flag = [1.0, 1.0]、x = [0.001715, 0.037768] |

## 物理・解析式
なし (開口の判定は上の楕円の式のとおり)。

## gs と pysad の比較
**△ 不一致 (pysad の未対応)。** gs は x=20 mm の粒子を flag=0 にする (内側の粒子は 1 周後 x=1.7 mm)。pysad は 2 粒子とも flag=1 のまま追跡し、外側の粒子の x は 37.8 mm まで進む (開口なしの追跡と同じ)。
`TYPE AP` も違う: gs は `APERT AP =(AX =.01 AY =.005 )`、pysad は `DRIFT AP =(L =0 )` (APERT を長さ 0 のドリフトとして読んでいる)。原因はコードまでは追っていない。pysad の `ffscmd.py`/`parser.py` には APERT の語はあるので、キーワード (AX AY) の取り込みか TrackParticles 側のどこかと思われる。

## 未確認・限界
- DX DY (開口の中心)、矩形 (DX1 DX2 DY1 DY2)、DA (動的口径) での使い方 (動的口径は 05_tracking の例で扱う) は未確認。
- LOSSAMPL・LOSSDZ (特殊変数の外側の開口) とは別のもの。
