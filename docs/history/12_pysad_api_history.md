# 12_pysad_api.md から取り除いた履歴の記述

12_pysad_api.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 12. pysad 公開 API リファレンス

整理前の 3 行目:

````text
最終更新: 2026-11-20 (v2.3.0 の現状のスナップショット。新しい API は §16。互換モード `pysad.compat` と v2 のモジュールは §16.6、v2.1.0 の `naff`・`specfn_cx`・`nbcache` は §16.7、v2.3.0 の周波数マップ・LMA・共鳴駆動項は §16.8)
````

## ## 16. 追加した API (§1〜15 に無いもの)

整理前の 407 行目:

````text
## 16. 追加した API (§1〜15 に無いもの)
````

## ### 16.5 限界・未確認

整理前の 474 行目:

````text
- 未翻訳: `wjfit/wfres/wfbin/wiord/wjdraw` (MAIN の TRACK の Fourier 解析。gs も SIGSEGV で参照が無い。チューンの抽出は v2.1.0 の `pysad.naff`: §16.7)。k64 モードでは MEASURE/TRCOD、SOL 領域内の MULT・BEND・CAVI、PHSROT、空間電荷、`Polarization` の有限値は動かない (v2 モードで実装。§16.6)。`measure.py` の `tlum/tluma`・MEASURE 系は gs と数値比較ができない (gs で MEASURE が動かない)。
````

## ### 16.6 互換モード `pysad.compat` と v2 で追加した API

整理前の 476 行目:

````text
### 16.6 互換モード `pysad.compat` と v2 で追加した API
````

整理前の 492 行目:

````text
**v2 で追加したモジュール** (`pysad/*.py`)
````

整理前の 502 行目:

````text
**v2 で追加・変更した主な API**
````

整理前の 511 行目:

````text
**使えない範囲・未確認**: k64 モードでは上の v2 だけの機能は動かないか無視される (それぞれの文書の「k64 モード」の欄)。`spch.py` の電子以外の粒子の式は 100 MeV の陽子で確かめた (`TrackParticles` は v2.1.0 から `MASS`・`CHARGE` を追跡に渡す。SAD スクリプトからの空間電荷の使用は未確認)。上のモジュールの公開関数名はソースの `def`・`class` から書いたもので、引数の細部は各ファイルの docstring を見ること。
````

## ### 16.7 v2.1.0 で追加したモジュール (`pysad/*.py`)

整理前の 513 行目:

````text
### 16.7 v2.1.0 で追加したモジュール (`pysad/*.py`)
````

整理前の 521 行目:

````text
| `specfn_cx2.py` | 複素の特殊関数の追加 (v2.2.0)。scipy の複素対応 (`exp1`・`expi`・`digamma`・`sici`・`fresnel`・`lambertw`・`hyp1f1`・`hyp2f1`)、整数次の ExpIntegralE の漸化式・級数、`polygamma` (Hurwitz ゼータ経由)、複素の次数の BesselJ の級数 (\|x\| ≤ 25)、Whittaker W (`specfn_cx.hyperu_c`)。`specfn.install` が複素の引数で切り替える (実数は従来のまま)。両モード |
````

整理前の 522 行目:

````text
| `spch_k64.py` | `K64Pspac` (v2.2.0。k64 モードの PSPAC。`SpaceCharge._kick2d` から呼ぶ)。K64 の形状関数: 横は CIC、縦は「ほぼ最近接」、電位は 5 点差分 (境界 φ = 0) の厳密解と双 3 次補間、線密度は 3 次スプライン。`PSPACNZ < 3` はエラー |
````

整理前の 525 行目:

````text
**v2.1.0 で追加・変更した主な API**
````

整理前の 535 行目:

````text
**v2.2.0 で追加・変更した主な API**
````

整理前の 538 行目:

````text
- `compat.set_bend_k1_def('log'|'lin'|None)`・`compat.bend_k1_def()`・`compat.bend_k1_lin()`・環境変数 `PYSAD_BEND_K1_DEF` (B160、v2.3.0 U2): K1≠0 の BEND の K1 の定義。'log' (既定) は K64 (中央面の磁場が B0 + (K1/h) ln(r/ρ0)、(1 + hx)∂b/∂x が一様)、'lin' は B' 一様 (elegant・MAD と同じ)。v2 だけ。k64 は常に 'log' (docs/11 A59、docs/47)。
````

整理前の 539 行目:

````text
- `compat.bend_dk_fix()`・`compat.set_bend_dk_fix(True|False|None)` (B161、v2.3.0 U2): CALC・`Emittance[]` の K1≠0 の BEND の解析行列の dk を kick の導関数に直す修正。v2 は既定で True、k64 は常に False。`set_bend_dk_fix` は測定・試験用の強制 (None で既定に戻す)。
````

整理前の 540 行目:

````text
- `spin.ST_EXACT`・`KINETIC`・`STEP_RATE`・`CHROMATIC` (v2 の Polarization。False で v2.1.0 の式)、`CHROM_KMAX`・`CHROM_NEL`。`spin.sideband_depolarization(..., C6=...)`、`spin._chromatic_matrix`。`spin_results` の `sdamp_rel` (厳密な緩和率)・`sdamp_rel_b2` (A2 の Σ r (b·n)²)・`kinetic`、不安定なリングの P = NaN (`depol_error`)。`track6._spin_rot_cs` (BEND の ROTATE のスピン)。
````

## ### 16.8 v2.3.0 で追加したモジュール (周波数マップ・LMA・共鳴駆動項・ADTS)

整理前の 546 行目:

````text
### 16.8 v2.3.0 で追加したモジュール (周波数マップ・LMA・共鳴駆動項・ADTS)
````

整理前の 555 行目:

````text
| `q1cmd.py` | 台本の関数 `FrequencyMap[位置, nturn, X->, Y->, DP->, Output->, Plot->]` (既定は RFSW・RAD・FLUC OFF、S2)、`MomentumAcceptance[位置, Cell->, Distance->, CheckEvery->, Verify->, TurnCheck->, Workers->, Output->, Plot->]` (`CheckEvery` の既定は 1、L2)、`TouschekLifetimeLMA[Formula->]` (既定は Piwinski 型、T2) | docs/48 §3.3・§4.3・§4.5 |
````

整理前の 558 行目:

````text
| `resdrive.py` | `resonance_terms(lat, terms=, method='analytic', nslice=8, delta=0, ref=None, rows=None, bend_cubic=True, bend_nslice=16, momentum_eV=3e9, order2=None, rdt4_options=None)` (引数の `None` は既定の動作。h_jklm、要素ごとの寄与。v2.3.0 R4: `bend_cubic=True` が既定で偏向磁石の 3 次の項の分 `.total_bend` を足す、`.total_sext` は 6 極のみ。δ≠0 では入れない)、`adts(lat, method=, include_octupole=True, delta=0, fourth=None, optics=None)` (6 極 2 次 + 8 極 1 次、`by_harmonic`・`per_element`。v2.3.0 V2: 既定は物理的な項だけ = `fourth='r3'`・`optics='calc'`。`fourth='r6'` は四極の積分器の台形則と偏向磁石の 4 次の項 (`.bend4`) を足す精密版、`optics='tracked'` は δ≠0 の光学を追跡の閉軌道から作る。docs/49 §11)、`chromaticity_numeric(lat, order=3, dmax=2e-3)`・`chromaticity_sextupole(lat)`・(R4) `chromaticity_natural(lat)` (自然クロマティシティの式 −(1/4π)∮(K1+h²)β + 端)・`chromaticity_bend(lat)` (偏向磁石の 3 次の項の分散による分)、`slice_thin`・`adts_tracking` (検証用)、v2.3.0 R2: `kinematic_adts`・`quad_fringe_adts` (`adts().kin`・`.fringe`・`.total`)。v2.3.0 R3: `adts(lat, bend_cubic=True, bend_nslice=16)` (偏向磁石の 3 次の項: `.bend`・`.cross`・`.total`、`bend_cubic=False` で R2 の値)、`bend_cubic_adts`。`bendnl.py`: `bend_decomposition(el, nslice)`・`bend_atoms`・`extract_map`・`Taylor2` (追跡の写像から偏向磁石の 3 次のリー生成関数を取り出す)。`rdttrack.py` (v2.3.0 R4、追跡から RDT を取り出す): `ring_map(elements, momentum, row0, amp=1e-5)`・`rdt_from_map(rm)` (1 周の写像の 2 次のテイラー展開から p=0 の h_jklm)・`delta_quadratic_total`・`delta2_orbit_h10002` (δ の項)・`tracked_rdt_table`・`periodic_points`・`predict_periodic_J` (3νx 共鳴線の近くの周期点)。`bendnl.py` (R4): `harmonic_amplitudes_shift`・`bend_atoms_rdt` (h への偏向磁石の寄与、δ の項は分散による軌道のずれの展開)。`adtsdirect.py`: `direct_adts(lat, momentum_eV, jmax, delta=0, order=2, nturn=1024)` (追跡 + NAFF + 多項式フィットの ADTS)・`tune_vs_delta`・`fit_tune_surface`・`fit_delta_curve`・`paraxial_tracker`。台本 `ResonanceDrivingTerms[]` (R4 で入力 `BendCubic->`・`BendSlices->` を追加)・`AmplitudeTuneShift[]` (R2 で `Kinematic->`・`QuadFringe->`・`Total->` を追加。R3 で入力 `BendCubic->`・`BendSlices->` と出力 `BendCubic->`・`SextBendCross->` を追加)・`ChromaticityTerms[]` (R4 で出力 `NaturalPart->`・`BendPart->` を追加) | docs/49 |
````

整理前の 561 行目:

````text
限界・未確認 (詳細は docs/48 §7・docs/49 §5・§11.8・docs/46): 周波数マップ・LMA の誤差つき・結合・FLUC は S2・L2 で調べた条件だけ (docs/46)。k64 モードの LMA、PF リング、SP8 (大きなリング B) は未確認。共鳴駆動項・ADTS は結合・垂直分散・端部・運動学的な項・6 極×8 極・10 極以上を入れない。
````

