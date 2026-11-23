# fn_CSRInit と fn_CSRTrack
種別: 関数 / フォルダ: 06_collective / ラティス: なし (縦の 1 周写像を直接与える)。参考: `docs/15_csr.md`

## 何を見る例か
CSR (コヒーレント放射) の組み込み関数のうち、`CSRInit` (ダクト断面の格子上で静電場を解く) と `CSRTrack` (ウェイクのキックを入れたマクロ粒子の縦方向追跡) の使い方。`CSRTrack` の縦の線形写像を、ウェイクなしの極限で式と照合する。

## スクリプト
```
r = CSRInit[{{0.5e-3, 0.5e-3}, 4, 1e-3, {4, 6, 8, 8, 8, 8, 6, 4}}];   ! {{dx,dy}, 中心列, rho, pipe}
z1 = newz[]; dp1 = newd[];                          ! 決定的な 2000 粒子 (z=3e-4 sin(0.37 i), dp=1e-3 cos(0.23 i))
zl = Table[0, {i, nr}];                              ! ウェイクなし (インピーダンス 0)
r = CSRTrack[z1, dp1, zl, {dz, 4, 0, 0.05, -0.3, 0, 0, 0, 0}, 6];
         ! {dz, ndiv, inten, r56, r65, frf, phirf, damp, sigmae}, nturn
```
- `CSRTrack[z, dp, zl, {dz, ndiv, inten, r56, r65, frf, phirf, damp, sigmae}, nturn]`: `zl` は長さ 2 のべきの実 FFT 形式のインピーダンス、`dz` はビン幅、`ndiv` は 1 周の分割数、`inten` は強度、`r56`、`r65` は 1 周の縦の行列要素 (m、1/m)。戻りは各周の $\{\sigma_z,\sigma_\delta\}$ を並べた長さ $2\,n_{turn}$ のリスト。
- `z`、`dp` は**その場で書き換わる** (gs も pysad も)。この例の 2 行目の出力 (`in-place`) で確認している。
- `zl` は人工の値で、物理的な CSR インピーダンスではない (`CSR.n` が作る ZL を使う場合は docs/15 の 2.5)。

## 結果の読み方
- `CSRInit: Dimensions={2, 8, 8}`: 戻りは 2 成分 (電場の 2 成分) で各 8×8 の格子 (pipe の列 8 本、高さは列ごとに 4〜8)。
- `CSRTrack zl=0, r56=0.05 r65=-0.3`: 各周の $(\sigma_z,\sigma_\delta)$。
- `inten=...  last (sigz, sigdp)`: 人工 ZL で `inten` を 0、1e-9、1e-8 にした 6 周後の値。

## 物理・解析式
ウェイクなし (inten=0、減衰 0) では、1 周は縦の線形写像: 位置 $z$ と運動量 $\delta$ が $z\to z+r_{56}\delta$、$\delta\to\delta+r_{65}z$ の組み合わせ (ここでは ndiv=4 の分割で、$r_{56}/n$ のドリフトと $r_{65}/n$ のキックを交互に置く対称分割)。この写像の rms の発展 (Python で計算):

| 周 | $\sigma_z$ (式) | $\sigma_z$ (gs = pysad) | $\sigma_\delta$ (式) | $\sigma_\delta$ (gs = pysad) |
|---|---|---|---|---|
| 初期 | 2.122046e-4 | 2.122046e-4 | 7.071168e-4 | 7.071168e-4 |
| 1 | 2.137810e-4 | 2.137717e-4 | 7.042615e-4 | 7.042655e-4 |
| 2 | 2.178913e-4 | 2.178566e-4 | 6.966606e-4 | 6.966857e-4 |
| 6 | 2.503055e-4 | 2.500725e-4 | 6.278989e-4 | 6.281577e-4 |

分割の式は 4e-5 (1 周) から 9e-4 (6 周) の差で一致する。差は、K64 の ndiv 分割の細部 (ドリフトの置き方) を正確に再現していないためと考える (`csrtrack.f` の読みは未確認)。1 周の縦振動は、$\omega^2=-r_{56}r_{65}=0.015$ ($\omega=0.1225$)、$\cos\mu=\cos\omega\simeq1+r_{56}r_{65}/2=0.9925$ で、シンクロトロンチューンは約 $0.0195$ /周。
ウェイクを入れると (inten>0)、$\sigma_z$ は $+2.4\times10^{-9}$ (inten=1e-9)、$+2.45\times10^{-8}$ (1e-8) 変わる。変化は inten にほぼ比例する (人工の ZL が小さいため微小)。

## gs と pysad の比較
○ 一致 (全 6 行が印字桁まで同じ)。

## 未確認・限界
- ZL は人工値。実際のダクトの CSR インピーダンスは `CSR.n` で作るが、約 7 分かかる (pysad) ので例にしていない。
- ndiv 分割の正確な式は未確認 (上の差)。`CSRHaissin`、`CSRSetupOY`、`CSROYMatrix`、`CSRStability` は動かしていない (`CSROYMatrix` は K64 の gs では動かない: docs/15 の B60)。
- 乱数 (`damp`、`sigmae`) を使う設定は見ていない。
