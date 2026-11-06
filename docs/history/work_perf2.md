# 作業ログ: 性能 2 (誤差入り 6D emit と CodCorr の要素カーネル化) (2026-11-05)

残課題 1 位「性能の続き: 誤差入り 6D emit (print_all_error, print_cod_magnitude)、CodCorr 1 反復 46 s (BEND tchge・K1≠0 BEND・薄レンズの Jacobian と tturne の numba 化)」。
方針は 10/16 (work_performance.md) と同じく **結果を変えずに速くする**。
制約 (11/05 ユーザー指示): 同じマシンで FEL の計算が走っているので FEL 優先。pysad は 1 プロセスだけ、Windows の優先度 IDLE、
numba/BLAS のスレッド 1 で走らせる (`<work-perf>/lowprio.py` = run_sad を IDLE で、`<work-perf>/lowrun.py` = 任意のコマンドを IDLE で)。
以下の時間はすべてその条件 (FEL と同居、IDLE) で測ったもの。

## 1. プロファイル (作業前)

マイクロベンチ `<work-perf>/mb_emiterr_full.sad` (ESRF v03_12、3563 要素): `cell calc emit` (誤差なし) → print_all_error と同じ誤差
(DELX/DELY ae Q*,S*、DK fe B*,Q*,S*、DTHETA te B*,Q*,S*、seed 125188771) → `cell calc` → `cell calc emit`。

| 区間 | 作業前 | 主なホットスポット (cProfile) |
|---|---|---|
| emit (誤差なし) | 3.5 s | — |
| calc (誤差入り、4D) | 2〜3 s | `_el_cod_step` 0.6 s/周 (16 周): BEND の `_tchge_jac` (tbrote の数値微分 4.4 s/全体)、K1≠0 BEND の `_bend_k1_cod_matrix`、要素ごとの getattr |
| emit (誤差入り、6D) | 11.6 s | `tturne` 0.7 s/周 (14 周): piece ごとの Python (98k piece)、`_tbendebody`/`_tbendecorr` (numpy の小行列 4.5 s)、`_tchge_jac`、`ndarray.copy` 2.3M 回 |
| 合計 (wall) | 21.6 s | |

print_all_error は 1 実現に calc 1 回 + emit 3 回 (RMSWrite ×2 と誤差を戻した後) なので 1 実現 ≈ 20 s、300 実現で 100 分 (SAD は 37 s)。

## 2. 対策

| # | 対策 | 内容 | 効果 (mb_emiterr の wall) |
|---|---|---|---|
| 1 | tchge の解析ヤコビアン (`lattice._tbgrote_mat/_tbrote_jac/_tdrife1_mat`、`_tchge_jac` の書き直し) | SAD tchge.f の tbgrote/tbrote/tdrife1 の trans の直訳。(5,6) は SAD の式のまま (写像の微分ではない: dr=0 でも dpzi·(pr+1)/pzi が残る)。入口で th≠0 のとき回転が行列に入らない SAD の癖はそのまま | 21.6 → 17 s 程度 |
| 2 | K1≠0 BEND の本体を numba に (`pysad/perfkern_bend.py`: `bendecorr`, `bendebody`)、`_State.apply` を in-place の numba (`apply_tdb`) に | 演算順は Python 版と同じ (参照値 `pysad/tests/data_bend_kernels.npz` と 1e-12 で一致) | 17 → 13.1 s |
| 3 | **要素ごとのカーネル** `pysad/elkern.py` (numba): DRIFT (`k_drift`)、薄い多極 + tchge (`k_thin`)、QUAD (`k_quad`: tqfrie/tqlfre/tqente/tsolque 本体の写像と解析ヤコビアン、放射の分割、ずれ・ロール・K1<0 の回転)、BEND (`k_bend`: tchge 一般 → エッジ (tbfrie の解析行列) → 本体 (K1=0: tbendebody0、K1≠0: bendecorr/bendebody、ANGLE=0: tsteee/四極) → tchge)、tradke (放射) を 1 要素 1 回の numba 呼び出しに。6D (`emit6d.tturne`: 放射・dv 項・D/B・piece ごとの z 行の組み直し) と 4D (`lattice._el_cod_step`: 放射なし・dv なし) の両方から使う。解析ヤコビアンが作れない (pz²≤1e-12 等) ときは状態を書き戻さず Python 版に戻る | 13.1 → **2.8 s** (emit 誤差入り 11.6 → 0.9 s、calc 誤差入り 2.5 → 0.7 s) |
| 4 | `tturne` で cod が有限でなくなったら (tcod の反復が δ→−1・NaN に発散) 残りの要素を NaN で埋めて抜ける (`_nonfinite_state`) | NaN の cod ではカーネルが全要素で Python 版に戻り、QUAD 1 個あたり 0.09 s の数値微分を繰り返していた (CodCorr の v05_61 では emit 1 回 18 s、うち 17 s がこれ)。結果はどのみち NaN で、tcod は有限でない残差を「未収束」として次の経路 (trf0 固定の Newton 法) に進む | v05_61 の誤差入り emit 71 s → 5 s |
| 5 | `Random[]` を SAD の式に (`ffscmd.py`/`sadscript.py`: `SadRandom.uniform0`) | SAD の `Random` は RandomUniform0 (tfRandom_.c REG8) = (seed + 2³¹)/2³² の [0,1)。pysad は (0,1) の +0.5 の式を使っていたので、print_all_error の `Floor[Abs[Random[1000]*1e8]]*2-1` の seed 300 個のうち 4 個が 2 ずれ、その 4 実現が SAD と別の誤差になっていた (test_sadrandom.py) | (速度ではなく一致) |

細かい点:
- `_State` の u0/h1emit/dvemit は配列 `_S` (カーネルと共有) にし、属性は property。
- 4D の K1≠0 BEND は以前 `_State(cod, _Phys(1e9, 1))` (仮の運動量) の dv 項で z を進めていたが、カーネルは dv なし (SAD の qtwiss は dvemit=0)。
  cod の z が 1e-10 変わるだけで Twiss (x,px,y,py と分散) は変わらない (test_elkern の 4D 比較は z を 1e-6 で見る)。
- 6D のステアリング (ANGLE=0, K1=0 の BEND, `_tsteee`) は Python 版が数値微分だったのを SAD の解析行列 (tbfrie の trans、tsteee の trans1) にした (SAD と同じ)。
- 薄い多極のキック (x−iy)^k は CPython の complex ** int と同じ 2 進べき乗 (`_cpowi`) で計算する。numba の complex ** int は
  exp/log 経由で y=0 でも 1e-38 の虚部の雑音を作り、無結合のリングで x-y の結合に雑音が入って tecalc の RFSW OFF の εz (z 面の連立の条件数 ~1e6) が 4e-7 動いた。
- `_el_cod_step` の `_has_offset` (getattr 10 個) は設計軌道上のときだけ評価する。
- 回帰 26c と test_emit_report の比較は "Imag.tune:-0.0000000" と " 0.0000000" (−0 の丸め) を空白を無視して同じとみなすようにした。

## 3. 検証

- `pysad/tests/test_perf_kernels.py`: bendecorr/bendebody の参照値 (numba 化前の Python 版) との一致 1e-12、tbrote/tchge の解析ヤコビアンと数値微分の一致 (入口 th≠0 の癖も)、apply_tdb と numpy の一致。
- `pysad/tests/test_elkern.py`: 要素 11 種 (誤差入り) × cod 3 個 × radcod 有/無で 6D (cod/T/D/B/u0) が Python 版と 1e-10 (ステアリングは解析/数値の差で 1e-6)、4D (`_el_cod_step`) が 1e-10、ESRF 誤差入りの 1 周と Emittance が一致。
- `pysad/tests/test_sadrandom.py`: Random[] の式と seed 列。
- 既存: test_emit_jac (tiny ring の SAD 値との比較 1e-9 を含む) ほか pytest 全体、`QUICK=1 run_regression.py` (50 ベンチ)。

### SAD との照合 (<work-perf>/、比較は `cmp_pae.py`)

| スクリプト | pysad | SAD | 一致 |
|---|---|---|---|
| print_all_error (esrf v03_12、300 実現) | **完走 2343 s** (FEL と同居。3 実現のプロファイルでは 1 実現 ≈ 3.5 s) | 37 s | 300 実現のうち **292 実現が一致** (相対 4e-7 以下。うち 4 実現は Random[] の修正後に `pae4_full.sad` で確認)。残り 8 実現は不安定 (SAD は COD RMS 0・emittance 有限、pysad は COD 1e-7〜1e-12・emittance NaN): 不安定ラティスの扱いの差 (記録のみ) |
| print_cod_magnitude (esrf v03_12、ae 20 水準 × 1000 実現) | 2400 s で打ち切り (2388 実現、1 実現 ≈ 1.0 s) | 300 s で打ち切り (3178 実現) | 共通の **2388 実現すべて一致** (最悪 1.3e-5 は 7e-8 の値の雑音)。以前 (10/16) 「誤差の実現が違う」だったのは 11/04 の誤差の修正で解消 |
| correct_cod_v05_nomon の 1 実現 (`codcorr1_full.sad`: DA なし、40 反復) | **完走 86 s** (作業前は 1 反復 46 s、CodCorr 以降 900 s 超) | 270 s (全実現) | SAD と同じく誤差が COD に効かず (COD RMS 0、SAD は 4e-19)。補正量・DA は 10/16 の記録のとおり追わない |

プロファイル (3 実現の print_all_error、`<work-perf>/pae3.prof`): 14.5 s のうち 4D の閉軌道 (`_qcod_periodic` → `_one_turn_cod_matrix`、58 周) が 8.6 s。
1 周 0.15 s のほぼ全部が `_el_cod_step(_cached)` の要素ごとの Python 処理 (31 µs/要素 × 3563)。カーネル自体は数 µs なので、次に速くするなら
1 周分をまとめて numba に渡す (要素パラメータを配列に詰めて 1 回の呼び出しで 1 周) のが本筋。6D の `tturne` は 1 周 0.03 s。

回帰・テスト: pytest 243 passed (test_bend_feeddown/test_coord/test_sext_radiation/test_track6n_ext は従来どおり除外)。`QUICK=1 run_regression.py`: FAIL 0 (50 ベンチ、239 s)。

## 4. 変更したファイル

- 新規: `pysad/perfkern_bend.py`、`pysad/elkern.py`、`pysad/tests/test_perf_kernels.py`、`pysad/tests/test_elkern.py`、`pysad/tests/test_sadrandom.py`、`pysad/tests/data_bend_kernels.npz`。
- `pysad/lattice.py`: `_tdrife1_mat`/`_tbgrote_mat`/`_tbrote_jac`/`_tchge_jac`、`_bend_step_fast`/`_quad_step_fast`/`_thin_step_fast`、`_el_cod_step` の分岐、`_has_offset` の遅延評価。
- `pysad/emit6d.py`: `_State` (apply、_S)、`_tbendecorr`/`_tbendebody` (カーネル呼び出し)、`tturne` の分岐と `_nonfinite_state`、`_bende_fast`/`_quade_fast`/`_thine_fast`。
- `pysad/ffscmd.py`、`pysad/sadscript.py`: `Random[]` を `uniform0`。
- `pysad/tests/test_emit_report.py`、`cmp_emit_report.py`: −0 の表示差を無視。
- `<work-perf>/`: lowprio.py、lowrun.py、mb_emiterr*.sad、patch1.py/patch2.py (差し替えの記録)、gen_bend_ref.py、cmp_pae.py、codcorr1_full.sad、pae4_full.sad (seed 4 個の再実行)、
  長い実行の出力 (esrf.print_all_error.pysad.out、esrf.print_cod_magnitude.pysad.out、codcorr1.pysad.out) とログ (long_runs.log、chain*.log、pytest_full*.log、regression_quick*.log)。
  調査用の diag_*.py・*.prof・小さな出力は整理のときに削除した。
