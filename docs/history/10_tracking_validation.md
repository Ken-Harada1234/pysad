# 10_tracking 検証の詳細 (旧 08 §7 照合結果・§9.8 既知の差の原文)

本文 (docs/10_tracking.md) の「検証」には要点だけを置き、ベンチごとの比較表と gs の食い違いの切り分けの詳細をここに移した。

## 旧 §7 照合結果 (gs との比較表)

| ベンチ | 内容 | 結果 |
|---|---|---|
| Level 6 | 小リング 5 粒子 × 1/10/100 周、NORAD / RAD (NOFLUC) | ≤5e-11 (SAD の印字精度 1e-10) |
| Level 7 | K1≠0 BEND (tbendi)、誤差 DELX/DELY/DTHETA (tsolrot, TENT/TEXIT, tbshift, tbrot) | ≤5e-11 (z のみ 1.7e-10、trf0 の微差) |
| Level 8 | 四極 F1/FRINGE=3、長さのある空洞 L=0.5、MULT K2 | ≤5e-11 / ≤5.9e-10 / ≤1.3e-10 |
| Level 9 | ソレノイド領域 (BZ=2 T, D1+QD)、tsolqur の bzp の癖 | NORAD ≤5e-11、RAD ≤1.5e-10 |
| Level 10 | RAD+FLUC, SeedRandom[17] | ≤5e-11、追跡後の seed も一致 |
| Level 11 | MULT K0/K1/K2 + DX/DY/ROTATE (tsolqum の dv 項) | ≤5.9e-11 |
| Level 25 (a,b,c,d,e,f,h) | ROTATE: 垂直ステアリング (tsteer L>0 と L=0)、skew 四極 (L=0 の SQ は tthin(nord=4))、ROTATE 付き BEND (K1 あり)、6 極〜十極・MULT、組み合わせ | 20 周 NORAD ≤7e-16、RAD NOFLUC ≤4.4e-13 |
| Level 12 | 実ラティス 試験用のラティス (3563 要素、NLAT=3564) 5 粒子 × 1/20 周、NORAD/RAD/FLUC | ≤4.7e-11、乱数状態一致 (track6n も同一、track6 との差 4e-18) |
| Level 29e | FSHIFT=1e-6 の追跡 (実ラティス、NORAD/RAD/FLUC、5 粒子) | 座標 ≤1.9e-11 |
| Level 16j | 51 粒子 × 20 周 NORAD (DA の初期分布) | 5e-16 |
| Level 16g/16k (NPARA=1) | 51/32 粒子 RAD+FLUC、損失なし | 最終 seed まで完全一致 |
| Level 16l_j22/25/27 | 単粒子 (損失あり) を TrackParticles で 1,2,4,10 周 | 12 件すべて最終 seed 一致 |
| Level 16a | DynamicApertureSurvey NORAD 153 粒子 × 100 周 | SAD の地図 153 文字が全て一致 |
| Level 20 small | DynamicApertureSurvey RAD+FLUC 153 粒子 × 100 周 (実ラティス) | 地図・Score 63 一致。最終 seed は不一致 (カオス、下記)。再同期テストは 30 周で seed 全一致 |

**落ちる粒子の乱数消費** (詳細は [`docs/history/work_loss_rng.md`](history/work_loss_rng.md), [`docs/history/work_loss_rng2.md`](history/work_loss_rng2.md)): 粒子が落ちる周は、発散した粒子が tapert で除かれるまで光子の乱数を引き続ける (θ=π/2, anp≈94 なら約 100 個、NaN なら 1 個; §4) ので、その間の消費数を SAD と揃える必要がある。pysad が SAD に合わせている点は、`sqrt1`/`pxy2dpz`/`asinz`/`norm2`/`xsincos` の inf/NaN の出方 (§6.2)、`tapert` の la の数え方 (§2.3)、`rot1` が th2=0 のとき回さないこと (§3.2)、tradkfn の `max(gmin, g+dg)` が NaN で gmin を返すこと (`gfloor`, §4)。
判定は**再同期テスト**で行う: SAD のデバッグ版 (各要素の全粒子の状態と seed を印字) の状態から pysad の次の 1 要素を走らせて seed・座標を比べる。30 周 (小規模 DynamicApertureSurvey の 153 粒子、93 030 記録) のすべてで seed が SAD と一致する。
**自由に走らせた最終 seed は一致しない**: 落ちる粒子は発散中に光子を引き続け、その有無が 1e-11 の入力差 (丸め誤差の増幅) で逆になる。
SAD 自身も初期座標の相対 1e-16 の摂動で 1 周後の seed が変わる (12 回中 2 回、1e-10 では 12 回中 12 回) ので、この基準での一致は期待できない。

**未対応** (両モード): SOL 領域内・COORD・TCAVI・MAP のあるラティスの光子 (PHOTONS の PhotonList は BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND で対応: 9.7 (b))、BOUND=0 の SOL、wake のうち MULT の F1・FRINGE (23_wake.md)、BEAM/WIG 要素。k64 モードではさらに、SOL の DX/DY/DZ/CHI と F1、SOL 領域内の BEND・MULT・CAVI (§9)、空間電荷 (v2 モードは SPAC・PSPAC・WSPAC を入れる。9.7)。
(POL のスピン追跡、tmulta、COORD・TCAVI・UND・MAP、WAKE、BEAMBEAM は §9 で対応済み。LOSSMAP は両モードで対応。)


## 旧 §9.8 既知の差 (gs と合わない。いずれも `docs/history/work_unimpl.md`、スピンなしでも出る・未調査)

1. **SOL (BZ≠0) 領域内の四極 + RAD の追跡**: 保存済みの参照 `bench_level9.sad` (gs の出力、2026-10-14) は pysad と一致するが、現在の gs (2026-10-19 ビルド) は同じ入力で別の値を返す (`S R T1 P1` の x が −.0027867898 から −.0014539780)。NORAD は両方で一致。
   **保存済みの参照と現在の gs のどちらが正しいかは断定できないが、2026-11-16 夜の切り分け (単一粒子の 1 周追跡、2.5 GeV、RAD NOFLUC、emit のあと) で、pysad ではなく現在の gs が誤っている疑いが強くなった**。SOL の BZ を 0 → 0.01 T に変えると:
   (1) 領域内がドリフトだけ (`SB D1 SE`): gs = pysad (BZ=0・0.01 とも 9 桁以上で一致)。
   (2) 領域内が四極 (`SB QD SE`): BZ=0 は gs = pysad。BZ=0.01 で gs は x=6.33e-4、y=−5.6e-16 (横の結合なし)、pysad は x=4.04e-4、y=1.1e-8 (BZ に連続)。
   (3) `SB D1 QD SE`、BZ=0.01: gs の結果が四極なしの `SB D1 SE` と 8 桁で一致する (x=1.194038e-4、px=−5.363477e-4)。**gs は RAD ON かつ BZ≠0 の領域内で四極の収束・発散力を入れていないように振る舞う**。BZ=0.01 と 0.2 でも gs の x はほぼ同じ (2.95e-5、BZ=0 の −2.23e-4 とは大きく違う) ので BZ に連続でない。
   (4) 現在の gs の 2 つのビルド (10-19 の Linux6 と、11-14 に -O0 で作り直した `~/sad_dbg`) は同じ値 (x=2.954816161E-5、10 桁一致) で、最適化の違いではなく現在の K64 ソースから作った gs の振る舞い。保存済みの参照 (10-14) を出した 10-15 の古い実行ファイル (`obj/AMD64-Linux5/sad1.exe`) は symbol lookup error で今は動かず、同じ版かは確認できない。
   結論: pysad が壊れているのではなく、現在の gs が SOL 領域内の四極 + RAD を物理的に誤って扱う疑いが強い (pysad は物理と保存済みの参照に合う)。**K64 ソースのどこが原因かは未特定**。§7 の Level 9 の RAD の行 (≤1.5e-10) は保存済みの参照との比較である。

