# flag_TRPT
フラグ・05_tracking・ラティス `00_lattice/linac_line.sad` (新規: 5 セルの簡易 FODO 直線加速器、入口 100 MeV、ACC 10 台 (各 L=1 m、10 MV、2.856 GHz))・`TRPT; INS; RFSW; NORAD; NOFLUC`

## 何を見る例か
`TRPT` (輸送ライン) では、空洞 (CAVI) で基準運動量 (設計運動量) が上がる。リング (NOTRPT、既定) では空洞はエネルギーを与えないドリフト扱い。加速すると断熱減衰 ($x'\propto1/p$) でエミッタンスが $p_0/p$ に比例して小さくなる。CAVI の位相 `PHI` を変えて、加速量と減衰を確かめる。

## スクリプト
```
RFSW; NORAD; NOFLUC;
USE LINAC;
TRPT;                  ! 輸送ライン: 空洞で基準運動量が上がる
INS; CALC;
phis = {0, -0.25*Pi, -0.5*Pi};
Do[ Element["PHI","ACC"] = phis[[k]]; CALC;
    r = TrackParticles[{1, {{1e-4},{0},{0},{0},{0},{0},{1}}}]; ...
    inv = gf x^2 + 2 af x px + bf px^2 ; ratio = inv/inv0 ,  {k,1,3}];
```
- `MARK ORIGIN = (AX=0 BX=5 AY=0 BY=5 EMITX=1E-8 ...)` が入口の Twiss。TRPT ではラインの入口条件から追跡と光学を行う。
- CAVI の位相: $\Delta E=-eV\sin(\mathrm{PHI})$ (SADHelp)。電子では `PHI=−π/2` で最大の加速。

## 結果の読み方
各 `PHI` の行: 追跡した 1 粒子の終端の $x$、$p_x$ (入口の基準運動量で割った値)、$dp$ (その位置の基準運動量に対する値)。次の行: 終端の Twiss ($\beta_x$、$\alpha_x$)、入口の不変量 inv0 = $\gamma_0x_0^2=2\times10^{-9}$、終端の不変量 invf、比 `ratio`。

## 物理・解析式
基準運動量の増加: 1 台あたり $\Delta E=V|\sin\phi|$ なので 10 台で $\Delta E_{tot}=100\,\mathrm{MeV}\times|\sin\phi|$。$p_0=100$ MeV/c、$E_0=\sqrt{p_0^2+m^2}$ から $p_f=\sqrt{(E_0+\Delta E_{tot})^2-m^2}$。
断熱減衰: 幾何エミッタンス $\varepsilon\propto1/p$ なので、同じ粒子の不変量は $J_f/J_0=p_0/p_f$。

| PHI | $\Delta E_{tot}$ [MeV] | $p_f$ [MeV/c] | 式の $p_0/p_f$ | pysad の ratio | gs の ratio |
|---|---|---|---|---|---|
| 0 | 0 | 100.000 | 1.0000 | 1.00003 | 1.0812 (*) |
| −π/4 | 70.71 | 170.711 | 0.58578 | 0.58678 | 0.6055 (*) |
| −π/2 | 100.00 | 200.001 | 0.49999 | 0.50175 | 0.50175 |

(*) の gs の値は、終端の Twiss が更新されないために生じた見かけの差 (下記)。粒子の座標 $x, p_x, dp$ そのものは gs と pysad が 8 桁以上一致している。
pysad の ratio は式に 0.1〜0.2% 近い (残りは RF の横集束 (ponderomotive)・エッジ項と、$10^{-4}$ の有限振幅によると考えられるが、確認はしていない)。

## gs と pysad の比較
△ 粒子の座標は一致 (表示の全桁、3 つの PHI すべて)。差が出たのは `Element["PHI","ACC"] = ...` の後の `CALC` の Twiss (BXend、AXend): pysad は PHI の変更ごとに終端の Twiss (加速で変わる) が更新され、gs は最初の PHI (−π/2) のときの値 ($\beta_x=19.153$) のまま変わらない。
gs の TRPT の光学は、設計運動量の表 (`gammab`) を `USE`・最初の `CALC` で作り、`Element` で CAVI の PHI を変えても作り直さないとみられる (`USE` を再実行しても変わらないことを確認)。pysad は毎回作り直す。どちらが SAD の意図かは未確認 (gs の動作が本家)。

## 未確認・限界
- $E_0=100$ MeV、5 セル、振幅 $10^{-4}$ m の 1 粒子だけ。RAD は OFF。`PHI=+π/2` (減速) では運動量がほぼ 0 になり、gs は NaN、pysad は `math domain error` で止まる (比較しない)。
- 加速のあるライン (TRPT) の `Emittance[]` や wake は見ていない (wake は [flag_LWAKE_TWAKE_WAKEOPT](../06_collective/flag_LWAKE_TWAKE_WAKEOPT.md))。
- 式との差 (0.2%) の原因は未確認。
