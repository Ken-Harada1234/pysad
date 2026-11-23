# v2.5.0 記録: 偏向磁石の端の Maxwell フリンジの COD の更新の因子 p (B324)

日付: 2026-11-23。試験 `pysad/tests/script_tests/test_v250_bend_maxwell_calc_p.py` (27 件)。gs の保存した出力は `pysad/tests/script_tests/scripts/v250_bend_maxwell_calc_p__ring_{b1,b2,b3,c1,c2,d1}.{sad,gs.txt}`。

## 1. 症状

`physics_guide/05_bend.md` 付録 B-1・`06_edge_fringe.md` §8.3 が、CALC の `tbedgemaxwell` (`tbfrie.f:190-229`) の COD の更新が、追跡・ヘルプ・端の磁場のローレンツ力の積分の式から因子 p = 1+δ だけ小さいと報告した (疑い。gs では未確認)。

- `tbfrie.f:209`: `cod(1)=cod(1)+bf*pr` (bf = b·y²/ρ/pv³)。正しくは Δx = bf·pr²。
- `tbfrie.f:211`: `cod(5)=cod(5)-bf*pxi`。正しくは Δz = −bf·pxi·pr。
- Δpy (`cod(4)`) は正しい。
- 同じ関数の解析行列 `trans(1,3) = af*pr**2/pvi2` ほかは p² の式の微分。COD の更新の数値微分と 1+δ 倍食い違い、更新のヤコビアンはシンプレクティックでない (δ = 0.1 で 1.7e-5。δ = 0 でも δ の列で 1e-6)。

pysad は k64・v2 とも K64 の直訳 (`lattice._tbedgemaxwell_cod`、numba の `elkern._tbedgemaxwell`) だった。

## 2. gs での確認

gs (SAD 1.4.4.2k64、WSL、`gsrun.sh` と同じ呼び出しで台本名を一意にして 1 本ずつ) で、検出しやすい条件を設計した。

- ラティス: BEND 4 台 (B1 = ANGLE 0.2、B2 = ANGLE −0.2、どちらも L = 1 m、ρ = ±5 m、FRINGE = 1、DISFRIN = 0、E1 = E2 = 0 または ±0.05)、QUAD QF (K1 = 0.9)・QD (K1 = −0.9) (L = 0.3 m、DY = +3 mm・−3 mm。垂直の閉軌道 y = 3.5〜8 mm を作る)、ドリフト 0.5 m。セル `ORG QF D B1 D QD D B2 D ENDM`。3 GeV/c、`NORFSW; NORAD; NOFLUC`。正負の曲率の磁石を置くのは、δ ≠ 0 でも水平の分散による x の閉軌道を小さく (x ≲ 3 mm) 保ち、δ の項を端の Maxwell に集めるため。
- 手順: `USE CELL; DP0 = δ; CELL; CALC;` (USE が DP0 を 0 に戻すので、USE のあとに DP0 を与える)。CALC の閉軌道 c = (Twiss["DX".."DDP", 1]) から `TrackParticles` を 1 周 (RFSW なし)。c の周りの ±1e-7 の 8 粒子から追跡の 4×4 M を作り、dc = (I − M)⁻¹ (T(c) − c) を求める。c + dc から追跡を 1 周した残りは 1e-14 以下 (追跡の閉軌道に収束)。dc は CALC の閉軌道と追跡の閉軌道の差。
- 疑いが正しければ、dc の px は δ·y²/(2ρ) 程度で δ に比例し、δ = 0 で 0。

| 条件 | gs の dc (x, px) | pysad k64 | pysad v2 (修正後) |
|---|---|---|---|
| δ = 0 (E = 0) | (−6.6e-14, −4.6e-13) | 同じ (1e-16 以下) | (−3.2e-14, 5.2e-14) |
| δ = 0.05 (E = 0) | (−1.16e-8, 8.46e-8) | 同じ (1e-16 以下) | (4.5e-10, 3.8e-10) |
| δ = 0.2 (E = 0) | (−1.02e-8, 3.80e-7) | 同じ | (3.5e-8, 3.2e-8) |
| δ = 0.2 (E1 = E2 = ±0.05) | (1.96e-9, 3.85e-7) | 同じ | (3.6e-8, 3.2e-8) |
| δ = −0.1 (E1 = E2 = ±0.05) | (3.71e-8, −1.49e-7) | 同じ | (1.3e-8, 1.4e-8) |
| δ = 0.1 (E = 0、DY = ±6 mm、y ≈ 8 mm) | (−8.42e-8, 6.90e-7) | 同じ | (7.6e-9, 5.3e-9) |

(gs の COD の x・y・px・py は pysad の k64 と 1e-15 以下で一致。保存した出力は `scripts/` の `.gs.txt`。)

読み: δ = 0 で dc が 1e-12 以下 (CALC の閉軌道と追跡の閉軌道が一致する条件) なのに、δ ≠ 0 で px に δ·y²/(2ρ) の大きさのずれ (4 つの端で足し合わさる。δ = 0.2・y = 4.9 mm・ρ = 5 m で δ·y²/(2ρ) = 4.8e-7 に対し 3.8e-7) が出る。符号は δ の符号に従う。疑いは正しい。さらに、pysad の k64 のまま COD の更新を p² 型にすると (v2)、dc の px は δ = 0.05 で 0.45%、δ = 0.2 で 8%、δ = 0.1・y = 8 mm で 0.8% に下がる。

v2 の dc の残り (δ = 0.2 で x・px とも 3e-8、y・py は 1e-8 台) は、端の Maxwell 以外の δ の項 (四極・本体・ほかの端) で、原因は特定していない (δ = 0 でも k64 の y に 3e-10 の差があり、QUAD の本体の積分の差と推定。未確認)。この調査の範囲外。

## 3. 修正 (v2 のみ)

| 場所 | 変更 |
|---|---|
| `pysad/lattice.py` `_tbedgemaxwell_cod` | v2: `cod[0] += bf*pr*pr`、`cod[4] -= bf*pxi*pr`。k64: K64 のまま (`bf*pr`、`bf*pxi`)。δ = 0 はビット同一 |
| `pysad/elkern.py` `_tbedgemaxwell(c, rhob, p2)` | `p2` (v2) のとき同じ式。`_edge_chain(..., p2)` に伝え、CALC の `_bend_core`・`_tsteee` が `edgefix` (= `compat.v2()`) を渡す |
| `pysad/track6n.py` `_edge_maxwell1`・`steer_edge_exact1` | v2 だけが使う追跡の厳密な端 (K1≠0 の BEND の v2 の端、ANGLE=0 の BEND の端) は常に p² 型 (追跡の `tbendcore` と同じ式にそろえた) |
| `pysad/track6.py` | `_lat._tbedge_cod` を呼ぶ v2 の端は `lattice._tbedgemaxwell_cod` の v2 の分岐を使う (変更なし) |
| `pysad/ampnf_k0.py` | 追跡の厳密な端の写像のコピー 2 か所 (`_edge_exact` の `maxwell`、`_maxwell_edge`) を p² 型に |

解析行列 (`_tbedgemaxwell_mat`、`elkern._tbedgemaxwell` の行列) は変更なし (もともと p² の式の微分)。v2 では更新の数値ヤコビアンが解析行列と 1e-9 以下で一致し、|MᵀSM − S| が 1e-9 以下。k64 は解析行列と 5e-7 以上食い違い、シンプレクティックでない (δ = 0 でも)。

## 4. 変わらないもの

- k64: gs と同じ (上の表)。δ = 0 は v2 も k64 もビット同一 (`test_delta_zero_same_in_both_modes_bitwise`)。
- DISFRIN: K1 = 0 の BEND の Maxwell フリンジを止めない (B134 の関連) 動作は変更なし。F1・FB1・FB2 の処理も触れていない。
- 追跡 (K1 = 0 の BEND の `tbendcore`、`k_bend0`) は変更なし (もとから p² 型)。

## 5. 未確認・限界

- 4 次元の閉軌道だけ gs と比べた。条件: 3 GeV/c、ρ = ±5 m、y = 3.5〜8 mm、δ = 0.05・0.2・0.1・−0.1、F1 = 0、FRINGE = 1・DISFRIN = 0、E1 = E2 = 0・±0.05、K1 = 0。
- Δz (`cod(5)`) は 4 次元の CALC の閉軌道に効かないので、gs と直接比べていない (ローレンツ力の積分・閉じた式・解析行列の数値微分で確認)。Emittance[]・CALC6D の 6 次元の量は未比較。
- F1 ≠ 0 (tbendef1 経由)、FRINGE = 2・3、K1 ≠ 0 の BEND、ANGLE = 0 の BEND (tsteee)、SOL 内は gs で比べていない。
- v2 の dc の残り (3e-8、δ = 0.2) は端の Maxwell 以外。
