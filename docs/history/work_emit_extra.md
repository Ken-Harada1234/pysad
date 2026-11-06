# 作業ログ: emit_extra.py (Emittance[] の OneTurnInformation / Orbit / Matrix の直訳) の解説と照合

日付: 2026-10-15 / 対象: `pysad/emit_extra.py` (317 行)。コードは未変更。
`bench_level21a_pysad.py` を実際に再実行し、保存済みの `bench_level21a_pysad.result` と時間行以外が完全に一致することを確認した
(emittance_rules 11.7 s / 保存時 16.2 s)。

## 1. 何をするモジュールか

`emit6d.py` が返すのは Emittances, MomentumSpread, BunchLength, Tunes, DampingRate, EquilibriumBeamMatrix, OrbitAtExit など一部だけだった。
`emit_extra.emittance_rules(lattice, momentum_eV, rfsw, radcod, one_turn, orbit, matrix, beam_matrices)` は Packages/emit.n の
`Emittance[opt]` が組み立てる規則一式を **SAD の規則名そのまま**のキーで dict に返す
(mode: Matrix=3 ⊃ Orbit=2 ⊃ OneTurnInformation=1 ⊃ 既定=0)。
`emit6d` の private 関数 (`_State`, `_tdrife0`, `_tbende`, `_tquade`, `_tcave`, `_tmulte`, `tecalc` …) と `calc6d` の `_closed_orbit`, `tfetwiss` を借用し、既存ファイルは変更しない。

処理の流れ: (1) tcod 相当の閉軌道 `_closed_orbit` → (2) `_turn_record` (tturne1 の plot=T 相当) で 1 周し、各要素入口で (累積 T, cod, beam) を記録、
(3) `tecalc` + `_normal_form` で params (安定性、tune、`cd`) を補完、(4) mode ごとに規則を詰める、
(5) beam_matrices=True なら temit.f:214-222 の「plot .and. postcal」の追加 1 周 (beam(1:21)=emitp から開始) で `BeamMatrices` を作る。

## 2. SAD 規則 ↔ Python キー ↔ 由来 の対応表

| SAD の規則名 | pysad の値 | 由来 (SAD 側) |
|---|---|---|
| Stable | Re log λ の絶対値 <1e-6 (3 モード) かつ閉軌道あり | tecalc の stab |
| Region | (1, nlat) | Region 指定は未対応 |
| Tunes | em['tunes'] (Im cd(4:6)/2π) | |
| EnergyLossU0 [eV], RfVoltageVc, EquilibriumPosition, MomentumCompaction, BucketHeight | em[...] (emit6d.tecalc) | |
| OrbitDilation | 0 (FSHIFT 未対応) | dleng = FSHIFT·circ |
| HarmonicNumber | wrfeff·c/ω0 (RF なしは 0) | tecalc |
| OrbitAtExit | 一周後の cod | |
| DampingRate | Re cd(1:3) | |
| Emittances / MomentumSpread / BunchLength | em[...] | |
| TuneShiftByRadiation | Im cd(1:3)/2π (rd = ri·D·r の 2×2 ブロックから) | tecalc |
| NominalSpinTune | h0·GSPIN (GSPIN=0.0011596521812800, tfstk.f:269) | |
| PolarizationTime / SpinTune / Polarization / Polarizations / PolarizationVector | 0 (偏極 calpol は未移植。SAD は params=0 で初期化) | |
| OrbitAtEntrance | codin | setiamat |
| OneTurnTransferMatrix | T = trans(:,1:6) (放射なしのシンプレクティック部分。emit6d の one_turn_matrix は T+D で別物) | setiamat |
| OneTurnDampingMatrix | D = trans(:,7:12) | |
| NormalCoordinates | ri (物理→正規モード) | |
| ExtendedTwissParameters | dict {AX, BX, …, ZPY} (28 個, TWISS_FUNS の順) = `tfetwiss(ri, cod, True)` | Packages/init.n FFS$TwissFuns |
| OneTurnExcitation | beamp (beam=0 から 1 周した量子励起、6×6) | tecalc |
| EquilibriumBeamMatrix | emitp (6×6) | |
| SpinPrecessionMatrix | 単位行列 (srotinit の初期値) | |
| ClosedOrbit | nlat 個の cod (各要素入口) | tturne1 の plot 分岐 |
| TransferMatrices | nlat 個の累積 T (T_1 = I) | 同上 |
| BeamMatrices (SAD の名前ではなく LINE["SIG",l]) | nlat 個の平衡ビーム行列 Σ(s) | temit.f:214-222、tftwiss.f:982- |

補助関数: `beam_size_at(e, l, key)` = `LINE["SIG"//key, l]` (key='' で 6×6、'X' で √Σxx、'XY' で Σxy)、`rules_to_sad_list(e, mode)` は SAD と同じ順
(`RULE_ORDER`) の [(名前, 値)] リスト。
未移植: SynchroBetaEmittance (temits.f)、Region 指定、InitialOrbit/InitialBeamMatrix、TRPT 分岐、IntrabeamExcitation (INTRA は touschek.py の領分)、
小数位置の LINE["SIG", l]。要素は Drift/Bend/Quad/薄レンズ多極 (SEXT 等)/Cavity/Multipole/長さ 0 (MARK, MONI) に対応し、それ以外は NotImplementedError。

## 3. 照合: bench_level21a

SAD 側 (`bench_level21a.sad`): pysad/tests/kek-ls_v03_12.lat、`CELL; RFSW; RAD; CALC6D; CALC;` のあと
`e=Emittance[OneTurnInformation->True, Matrix->True]` の各規則と、`ClosedOrbit`/`TransferMatrices` の 8 点 (l = 1, 2, 100, 500, 1000, 2000, 3000, 3564)、
`LINE["SIG",l]` (同 8 点)、`LINE["SIGX"/"SIGXY"/"SIGPX"/"SIGDP",100]`、`LINE["NAME",…]` を Print。
判定は「配列最大値に対する最大絶対差 (ノルム相対差)」、目標 1e-9。再実行結果:

| 項目 | ノルム相対差 | 備考 (要素最大相対差) |
|---|---|---|
| Stable, Region (1, 3564) | 一致 | |
| TUN (tunes) | 6.6e-11 | |
| U0V (U0, Vc, trf0, α, 膨張, bh, h) | 1.7e-16 | 1.1e-13 |
| OEX / OEN (閉軌道) | 1.7e-17 (絶対) | ほぼ 0 |
| DMP | 6.6e-13 | |
| EMT (エミッタンス) | 1.1e-12 | εy ≈1.4e-31 でも 4.9e-10 |
| SIG (σE, σz) | 1.6e-12 | |
| TSR | 1.3e-11 | 3.9e-9 |
| NSP, POL, SPM | 5e-16 / 0 / 0 | |
| OTM (T) | 1.7e-10 | 6.7e-10 |
| ODM (D) | 2.1e-10 (**最大**) | 3.8e-10 |
| NCO (ri) | 5.4e-12 | 7.0e-9 |
| ETK (28 名) | 順序一致 | |
| ETW (28 値) | 9.3e-12 | 結合・ y 由来の ≈1e-16 以下の項は雑音同士 (例: R1 pysad -1.3e-17 / SAD -2.1e-16)。BX 差 5.8e-11、BY 3.3e-11 |
| OTE / EBM | 1.1e-12 / 2.0e-12 | 9.1e-10 / 4.3e-10 |
| CO1…CO3564 | ≤1.7e-17 (絶対) | |
| TM1, TM2 | 0 | |
| TM100 / 500 / 1000 / 2000 / 3000 / 3564 | 3.9e-12 / 1.3e-11 / 5.4e-11 / 1.5e-10 / 2.0e-10 / 1.7e-10 | 周回とともに増える (誤差の蓄積) |
| SZ (LINE["SIG",l]) 8 点 + SZB 2 点 | 1.5e-12 〜 2.0e-12 | 要素ごとでも ≤1.4e-9 |
| SZX100 (σx, σxy, σpx, σδ) | 2.8e-13 | |
| 自己検査 Σ(nlat) vs emitp | 5.4e-14 | |
| NLAT / NAMES | {3564,3564,3564} 一致 / 名前の対応も一致 (SAD は MKS.1 のように通し番号付き) | |

**最大 2.126e-10 (ODM)、判定 OK (≤1e-9)**。前担当の報告「全項目 2.1e-10 で一致」と再現された。

注意点:
- ノルム相対差は零に近い成分 (y、結合) を配列最大値で割るため、それらの成分は雑音下限 (SAD 1e-21〜1e-16 / pysad 1e-17 程度) で一致を主張していない
  (平面ラティスなので物理的に 0)。
- TM は l が大きいほど差が増え 2.0e-10 まで達する。これは累積した 3000 要素分の丸めと、T の元の閉軌道 (RFSW の trf0) の 1e-10 級の差による。
- 偏極系 (PolarizationTime 等) は SAD も 0 を返す設定 (CALC6D のみ) の照合なので、calpol の正しさは検証されていない。

## 4. ffscmd.py への接続案 (今回は実装しない)

現状の `fn_emittance(args)` (739 行) は引数を無視して `run_emit` → 固定の 14 規則を返す。`OneTurnTransferMatrix` は `e['one_turn_matrix']` (= T+D で SAD の T とは別物) を入れており、
`Matrix`/`Orbit`/`OneTurnInformation` オプションも `LINE["SIG",l]` も未対応。

1. `fn_emittance(args)`: args から `OneTurnInformation->`, `Orbit->`, `Matrix->` の規則を読み (真偽、既定 False)、
   `emittance_rules(self.lat, self.momentum(), rfsw=flags['RFSW'], radcod=flags['RADCOD'], one_turn=…, orbit=…, matrix=…)` を呼ぶ。
   戻り dict は `rules_to_sad_list(e, mode)` で [(名前, 値)] に変換し、`Rule(S(名前), 値)` を並べる
   (ExtendedTwissParameters はさらに `Rule(S(AX), v)` のリストに。numpy 配列は `.tolist()`、真偽は 1.0/0.0)。
2. mode=0 (オプションなし) の既存の戻り値との整合: 既存 14 規則は emit_extra の mode=0 規則のサブセットなので、そのまま置き換えられる。
   ただし `Emittance[]` はグローバル (EMITX, SIGE …) を更新するので、`run_emit` を呼ぶ既存経路は残し、追加情報だけ `emittance_rules` で補う二段構えが安全
   (1 回 10 秒強なので重複計算が気になるなら `run_emit` を emit_extra に統合し `e['_emit']` を `self.emit` に使う)。
3. `LINE["SIG"//key, l]`: `beam_matrices=True` で `e['BeamMatrices']` を保持し (CALC6D 時に必要)、`beam_size_at(e, l, key)` を `LINE[…]` 実装から呼ぶ。
   これは CALC6D 時の `LINE["SIG",l]` / `Twiss["SIGX",l]` に相当するので、`run_emit` で CALC6D フラグが立っていたら自動で計算してキャッシュする。
4. 前提・制約: RFSW/RADCOD は emit と同じフラグ、RAD/FLUC は影響しない。未対応要素は NotImplementedError なので、`run_emit` 同様に捕えて
   `???-Emittance[]-closed orbit not found.` 経路に流すかエラーにするか決める必要がある。
5. 検証: 接続後は bench_level21a.sad をそのまま ffscmd 経由で流し、`bench_level21a_pysad.result` の最大 2.1e-10 と同等になることを確認する。

## 5. 判定

emit_extra.py は SAD の Emittance[OneTurnInformation->True, Matrix->True] の全規則と LINE["SIG",l] を最大 2.1e-10 (ノルム相対) で再現する。
偏極と SynchroBetaEmittance は未移植。次の手は上記の ffscmd.py 接続。
