# 04_emittance_radiation 一覧
ラティスは `00_lattice/FODO_ring_V4MV.sad` (FODO_ring の RF 電圧だけ 4 MV に変更)。`RINGQ` を使う (RADINT の NaN を避けるため)。式の検算は `formula_check.py` → `formula_check.txt`。

| ファイル | 説明 | gs と pysad |
|---|---|---|
| FFS_EMITTANCE | EMIT の報告と放射積分の式による検算 (U0, τ, εx, σδ, νs, σz) | ○ |
| FFS_RADINT | RADINT の I1〜I5 と自前の積分 | ○ |
| fn_Emittance | Emittance[] の規則 | ○ |
| flag_RFSW_RAD_RADCOD_FLUC | 4 つのフラグの効き方 | △ (Emittance[] は一致、FLUC の追跡は乱数列が違う) |
| flag_RADTAPER | RADTAPER の前後 | ○ (k64 モードで gs と一致。v2 の既定では εx が 2.4% 違う: Jx の式の違い B75) |
| flag_CALC6D | CALC6D の縦の光学 | ○ |
| fn_SynchroBetaEmittance | シンクロ・ベータ結合の emitx (平面リング) | ○ (結合のあるラティスは gs が落ちる) |
| [flag_RADCOD_closure](flag_RADCOD_closure.sad) | RADTAPER・RADCOD の閉軌道が、追跡 (RAD) の 1 周で閉じるか | △ (gs は閉じない: x 4e-7 m、RADTAPER ありで 2e-6 m。v2 は閉じる (1e-9 m)。k64 は gs と同じ) |
| [flag_RADTAPER_optics](flag_RADTAPER_optics.sad) | RADTAPER の閉軌道 (収束)、CALC6D の RADTAPER、CalculateOptics の TRPT の加速、加速区間の分散 | 軌道 −2.1732e-7 (gs は既定の CODCONV で −2.1977e-7、CODCONV=1e-4 で −2.1732e-7。CalculateOptics = CALC) |
| [flag_NORFSW_orbit](flag_NORFSW_orbit.sad) | `NORFSW; RAD` の Emittance[] と閉軌道のずれ (12 セル FODO、3 GeV、QF に DX = 1e-4)。K64 の tcod は RFSW OFF で x px y py だけ閉じ、δ・z は動かさない | ○ (k64 モードは gs と一致: εx = 2.70000645e-7 (gs 2.70000647e-7)。v2 は追跡の放射の式が違うので少し違う) |
| [var_MASS_CHARGE_proton_ring](var_MASS_CHARGE_proton_ring.sad) | TrackParticles・Emittance[]・CALC6D の `MASS`・`CHARGE` (30 GeV/c の陽子) | ○ (Qs、1 周・100 周の追跡が一致 (表示 6 桁)。U0 は v2 が 6.89e-4 eV (古典式)、gs・k64 が 1.265 eV: K64 の放射は電子の古典半径を使う) |
| [var_MASS_proton_beta_p1p5](var_MASS_proton_beta_p1p5.sad) | 陽子 (MASS = 938.272088 MeV) の放射の β の因子 (B112)。Emittance[] の EnergyLossU0 を古典式 U0 = (2/3) r_p m c² β³ γ⁴ Σ θ²/L と比べる (p = 1.5 GeV/c) | ○ (k64 モードは gs と全桁一致。v2 の U0 は古典式の 1/β³ 倍 (限界)) |
| [var_MASS_proton_beta_p3](var_MASS_proton_beta_p3.sad) | 同 (p = 3 GeV/c) | 同上 |
| [var_MASS_proton_beta_p30](var_MASS_proton_beta_p30.sad) | 同 (p = 30 GeV/c) | 同上 |
| [elem_BEND_F1_emittance](elem_BEND_F1_emittance.sad) | Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE (K1≠0 の複合機能も) | ○ (k64 モードは全桁一致 (表示 10 桁)。v2 は端の角の行路の補正で U0 が 1e-3 違う) |
| [FFS_RADINT_I4V_coupled](FFS_RADINT_I4V_coupled.sad) | 結合のあるリング (45 度傾いた四極 4 台) の RADINT の I4V。I4U + I4V = 物理的な I4 (1e-5)、DU・DV・Je と自前の 6x6 行列 + 1 次の摂動論との一致 (3e-3) | v2 の出力と独立な式の照合 (gs の出力あり) |
| [FFS_RADINT_k1_bend_vdisp](FFS_RADINT_k1_bend_vdisp.sad) | K1≠0 の BEND + 結合による Y 方向の分散の RADINT (B107) | △ (k64 モードは gs と全桁一致 (I4U + I4V = 0.420949)。v2 は独立な I4 (0.414392) と 1e-7) |
| [FFS_RADINT_quad_offset](FFS_RADINT_quad_offset.sad) | 据え付けのずれた四極の RADINT (B107b) | △ (k64 モードは gs と全桁一致 (ΔJx = +0.0028)。v2 は独立な積分・追跡・Emittance[] の ΔJx = +0.0078) |
| [FFS_RADINT_pf](FFS_RADINT_pf.sad) | PF リング (2.5 GeV、187 m) の RADINT・Emittance[]・EMITDIV と Polarization。自前の積分 (3x3 行列) との比較で I1〜I5 が 1e-6 以下、τp が 2e-5 (Polarization = 0.92376) | v2 の出力と独立な式の照合 (gs の出力あり) |
| [fn_SynchroBetaEmittance_coupled](fn_SynchroBetaEmittance_coupled.sad) | 結合のあるリング (24 セル FODO + skew 四極) の SynchroBetaEmittance | ○ (gs も走り、k64 モードは 7 桁以上一致。v2 は追跡の放射の式が違うので SB が 1.3% 違う) |
| [fn_SynchroBetaEmittance_dispersion_start](fn_SynchroBetaEmittance_dispersion_start.sad) | 始点を変えた 24 セル FODO の SynchroBetaEmittance (B128)。始点の分散 η ≠ 0 のとき K64 は Emittance[] の数倍〜40 倍を返す | △ (gs と k64 は SB が Emittance[] の 6.54 倍 (1e-7 で一致)。v2 は 1.0003〜1.0004 倍) |
| [fn_SynchroBetaEmittance_rings](fn_SynchroBetaEmittance_rings.py) | 大きなリング A (kek-ls_v03_12.lat)・PF・ASTRID2 の SB と Emittance[] の比 (結果は `.out.txt`、A の gs は `_ring_a.gs.txt`) | △ (A の k64 は 39.88 倍で gs と 2e-7 一致 (6.178294e-9)。v2 は 1.0084 倍) |
| [flag_EMITDIV_rect](flag_EMITDIV_rect.py) | Emittance[] の EMITDIV の収束 (矩形・扇形の FODO)。v2 の Emittance[] と RADINT の差の原因は EMITDIV の既定 (1) の離散化 (結果は `.out.txt`) | gs は使わない |
| [flag_RAD_tradke_exact](flag_RAD_tradke_exact.py) | 放射反作用の ddpx の符号 (B111)。放射反作用つきのローレンツ力の積分と、tradke・追跡 (tradkn) の比較 (結果は `.out.txt`) | v2 の tradke・追跡は厳密な積分と一致、k64 の tradke は既定の刻みで 9% ずれる。gs は走らない (Python の内部関数を使う) |
| [flag_RADCOD_derivative](flag_RADCOD_derivative.py) | RADCOD の解析 T + D と数値微分 (4 次の中心差分) の差 (結果は `.out.txt`) | v2 は BEND 1 台で 1.7e-10 (k64 は 5e-7)。gs は使わない |

上の表のうち、拡張子つきのリンクの行 (`.sad`・`.py`) は説明の `.md` を持たない。内容と結果の要点は表と、台本・スクリプトの先頭のコメントに書いた。`.py` の例は、リポジトリ直下から `py -3.12 example/04_emittance_radiation/NAME.py` で走らせる (出力は `NAME.out.txt`)。gs と同じ動作で走らせるときは `PYSAD_COMPAT=k64`。

## このフォルダで分かったこと
- RADINT の $J_x$ (0.9415) と EMIT の $J_x$ (0.9651) は 2.4% 違う。$I_4$ の端の項を 2 倍にすると EMIT に近づく (FFS_EMITTANCE.md)。原因: K64 の `Emittance[]`・EMIT が端の角の蹴りの行路を入れず $I_4$ の端の項が 2 倍になっていた (B75。RADINT が正しい)。v2 の既定では EMIT の $J_x$ は 0.9429 で RADINT に近い (gs と違う値になる)。`PYSAD_COMPAT=k64` では gs と同じ 0.9651。この節の他の例の pysad の出力 (`.pysad.txt`) は v2 のもの。
- 例 (`FFS_EMITTANCE`、`FFS_RADINT`、`fn_Emittance`、`flag_RFSW_RAD_RADCOD_FLUC`、`fn_SynchroBetaEmittance`) の判定は k64 モードでの gs との一致。v2 の既定では εx・減衰率・α などが B75・B77 の違いの分だけ gs と違う。
- `FODO_ring` の RF 1.5 MV は $U_0$ より小さく、バケットがない。
- 4 次元の CALC は RADTAPER を見ない (gs も同じ。試験 `test_v2_B1_optics.py` の `test_calc_4d_ignores_radtaper_like_k64`)。
