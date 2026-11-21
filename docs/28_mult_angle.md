# MULT (ANGLE≠0、曲げつき多極)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `tmulta.f` (追跡 `tmulta`、Twiss・Emittance[] の `tmultae`)。[`pysad/multa.py`](../pysad/multa.py) (追跡)、[`pysad/multae.py`](../pysad/multae.py) (COD 写像と行列)。

```
MULT M1 = (L=1 ANGLE=0.05 K1=0.1 K2=0.5 E1=0.5 E2=0.5);
```

曲げ磁石 (ANGLE) の上に、K0 (双極の過不足)、K1、K2 … と SK0 … (スキュー) の多極をのせる要素。`ndiv` 個のスライスに割り、
「半スライスの曲げ → 曲率つきの多極キック → 全スライスの曲げ → キック → … → 半スライスの曲げ」の順に進む。
曲げは BEND と同じ `tbend` (K1=0 は tbendcore、K1≠0 は tbendi) をスライスの引数で呼ぶ。キックは `gkn` の表 (`gkninit`) を使う解析式。

## 使われる所

| 計算 | 扱い |
| --- | --- |
| 追跡 (track6n) | `K_MULTA` (`multa.k_multa`)。放射 (RAD)、揺らぎ (FLUC)、TRPT の ndivrad ×10 も BEND と同じ |
| Twiss (CALC) | `lattice._el_cod_step` → `multae.lattice_step`。COD 写像は解析、行列は各段の解析行列の積 |
| Emittance[]・CALC6D・Touschek | `emit6d._tmulte` → `multae.emit` (誤差なしの BEND の `_tbende_core` を段ごとに呼ぶ) |
| 幾何 | BEND と同じ (ANGLE を使う) |

## 検証 (gs との比較、`pysad/tests/test_mult_angle.py` 21 本)

- 追跡 (11 通り、20 粒子、3 周): 相対 0 (K1≠0/=0、K0、スキュー SK2 SK3、E1 E2、フリンジ FB1 FB2 FRINGE=3、
  DX DY ROTATE DROTATE、CHI2、DISFRIN、RAD、RAD+FLUC (乱数も一致)、ANGLE=0.3 の大角度 (ndiv 大))。
- Twiss (CALC): チューン 1e-15、β・α 1e-13。
- Emittance[] (12 セルのリング、MULT が曲げ磁石): ε 1e-11、チューン 1e-13 (K1、K2、K3、K4、SK2、SK3、E1、E2、K0、DX DY ROTATE)。

## 限界・未確認

- 未対応 (NotImplementedError): BZ≠0 (K64 自体が abort)、VOLT (RF)、F1、TAPER、スキュー双極 SK0 (K64 の ndiv 計算が `fact(-1)` を読む。結果が定まらない)、L=0。
- 追跡の光子 (PHOTONS) は未対応 (ffscmd が SADError で拒否)。偏極 (CALPOL・POL) のスピンは、ANGLE=0 の MULT は docs/36 のとおり gs と一致し、ANGLE≠0 の MULT (`tmultae`) も gs と一致した (2026-11-16 夜: 偏極時間 14 桁、スピンチューン 1e-14。平面リングなので偏極ベクトルは 0,1,0 で、縦成分のある偏極ベクトルは確かめていない。TrackParticles のスピン追跡での ANGLE≠0 の MULT は未確認)。
- Twiss・Emittance[] で FB1/FB2 (tblfre) と CHI2 は未対応 (NotImplementedError)。追跡は対応。
- DROTATE≠0 の Emittance[] は ε で 5e-6 ずれる (チューン 1e-9)。追跡は一致。tchge の DROTATE の行列が数値微分のため (BEND と同じ)。
- K0 + DX DY + ROTATE を同時に与えたリングの Emittance[] は gs と一致しない (δ の閉軌道が 4e-3 と 5e-4 に分かれる)。水平チューンが 1e-3 で整数に近く、
  閉軌道の探索が不安定なためとみて追っていない。2 つずつの組み合わせはすべて一致する。
- tmultae の tradke は m ループの中にあり 22 回/スライス呼ばれる。結果に影響しない (1 回でも同じ) が、そのまま写してある (`multae.TRADKE_CALLS`)。
- 追跡と tmultae は K64 の書き方どおり違う: tmultae は ndiv に放射の項がなく、キックを (1+δ) で割らない。
