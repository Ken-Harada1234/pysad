# flag_JITTER
種別: フラグ / フォルダ: 05_tracking / ラティス: `00_lattice/FODO_ring.sad` (3 GeV) / フラグ: `JITTER`、`NORFSW NORAD NOFLUC`

## 何を見る例か
`JITTER` (反対語 `QUIET`、`NOJITTER`、`NOQUIET`) は、TRPT での入射ビーム生成で、重心のゆらぎを許すスイッチ (SADHelp: 「JITTER or NOQUIET allows jitter of the center-of-mass of the incoming beam in the case of TRPT. Otherwise the center-of-mass statistically fluctuates...」)。K64 のソースでは、`ttinit.f` の初期粒子の生成 (MAIN レベルの `TRACK`、`MEASURE`、`TRCOD` 用) でしか参照されない。`TrackParticles` に渡す粒子は利用者が与えるので、影響しない。この例では、4 通りのフラグ指定で同じ入力を追跡し、結果が変わらないこと (`TrackParticles` では効かないこと) を確かめる。

## スクリプト
```
FFS[{"JITTER", "NOJITTER", "QUIET", "NOQUIET"}[[k]]];   ! Do の中のフラグ切り替えは FFS[...]
r = TrackParticles[{1, beam}, 1, 1, 10][[2]];            ! 100 粒子 10 周
Print["case ", k, ": sum x=", Plus@@r[[1]], ...];
```
フラグは `Do` の中では `FFS["..."]` で切り替える (そのまま `JITTER;` と書いても `Do` の中ではコマンドとして実行されない)。

## 結果の読み方
4 ケースとも $\sum x=4.48582049\times10^{-5}$、$\sum y=-6.0993114\times10^{-4}$、$\sum p_x=7.92596428\times10^{-6}$ で、`max difference between cases` は 0。

## 物理・解析式
なし。重心のゆらぎとは、粒子数 $N$ が有限の乱数サンプルでは、重心が $\sigma/\sqrt N$ ほど揺れることを指す。`NOJITTER` は、その平均を引いて重心を 0 にそろえる (`ttinit.f` の `sx, sy, ...` の平均を引く処理)。

## gs と pysad の比較
○ 一致 (全 6 行)。pysad では JITTER は値を保持するだけ (`pysad/ffscmd.py` の `_FLAGS_DEFAULT`、docs/33 で「実装 (値を保持するだけ)」)。MAIN レベルの `TRACK`/`MEASURE` は K64 でも動かない (docs/10 §8) ため、効果を比べられる経路が無い。

## 未確認・限界
- JITTER が効く経路 (ttinit の初期分布生成) は、gs でも pysad でも試していない。この例は「`TrackParticles` に影響しない」ことの確認だけ。
- 3 GeV、100 粒子、10 周のみ。
