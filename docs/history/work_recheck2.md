# 作業ログ: ユーザー実スクリプト 86 本の一括再実行 (B、2026-11-06)

ユーザー指示 (11/06): 誤差の補正やチューンサーベイは実現数・格子の点を間引いて、動作テストに無為に計算時間を掛けない。
試したコマンドはラティスと一体化し (前半ラティス、後半コマンド)、SAD の gs にそのまま渡せる `.sad` と pysad 側の `.py` を同じフォルダに置いて比較できるようにする。

## 1. 道具 (<work-userscripts>/)

- `build.py`: `<work-recheck>/manifest.txt` の (ラティス, スクリプト) から `<work-userscripts>/esrf/`, `<work-userscripts>/pf2/` に `<name>.sad` (一体化) と `<name>.py` (pysad で同じ .sad を実行) を作る。
  modules_*.n もコピー (esrf の modules_general.n は draw → Print の版)。間引きの規則 `RULES` は正規表現で、効いた規則を `.sad` の先頭コメントと `manifest.txt` に書く。
- 間引き (SAD/pysad 共通): NPARA=1、TkWait[] 行の削除と `FFS["draw …"]` → `Print["nodraw …"]` (WSL の SAD が止まる・落ちる)、DynamicApertureSurvey 5000/2000 周 → 100 周、
  DA の運動量点 dp/2 → 2dp、誤差の実現 300/1000/110 → 3/3/2、チューンサーベイの格子 54×54 → 3×3、DP0 走査 0.001 → 0.01、振幅走査の刻み 20 → 100、
  VOLT 走査 100 → 4 点、PBUNCH 走査 11 → 3 点。
- `run_sad_all.sh` (WSL の gs、1 本 300 s)、`run_py_all.sh` (py -3.12、最低優先度・1 スレッド、1 本 600 s)、`cmp_out.py` (数値の比較)、`status_table.py` + `notes.txt` (状態表)。
- 比較の方法 (`cmp_out.py`): 入力のエコー行・バナー・In[]/Out[] を除き、数値を `#` に置き換えた骨格を difflib で揃えて対応行の数値を比べる。
  |差| < 1e-9 と SAD の表示桁 1 ulp 以内は無視、相対差 1e-6 超を不一致と数える。`Momentum compact. alpha` の値は比べない (表示桁の差)。

## 2. 実行

- SAD: esrf 52 本 + pf2 34 本、全部 1〜50 s (間引き前に 300 s で打ち切られていた disp2/disp_da/disp_test/print_cod_magnitude は 47/39/178/4 s)。
  `TkWait[]` で止まる 2 本 (pf2 disp、esrf disp_L) は行を消して取り直した。
- pysad: pf2 34 本は 1〜38 s (10/16 の 164〜181 s は性能改善で 3〜4 s)。esrf は tune_survey_linear 107 s、chrom/matching3_01 が長い。
  600 s 打ち切りは amp_dep_ts_x_octa だけ (SAD は OF1B が無いので `General::invset` で Do を抜けて 1 s。pysad も C5 の実装で同じになる)。

## 3. 見つかって直したもの (pysad)

| 項目 | 現象 | 修正 |
|---|---|---|
| INTRA の EMIT 報告 | SAD は零電流の報告 → Touschek の表 (運動量アクセプタンスごとの寿命、Nx/Nz・Ny/Nz の寿命/100s) → `Parameters with intrabeam scattering` (IBS 込みの値) の 2 段。pysad は IBS 込みの値を 1 段で出していた (touscheck_r, touscheklife_r, disp_bunchlength_r, beamsize_r の骨格 28〜65 %) | `ffscmd._print_touschek`, `_print_emit_intra`, `print_emit(e0)` (temit.f 9103/9104/9105/9131)。4 本とも数値一致に |
| PBUNCH の既定値 | SAD は 1e10 (initb1.f)。pysad は未定義で 0 → 振幅の表が 999 | グローバル既定値に PBUNCH=1e10, NBUNCH=1, EMITXC/YC/ZC=0 |
| FFS[] のエラーで Do を抜ける (C5) | SAD は `?No variable.` 等 (` ???-FFS-Error-`) で `???FFS::error: Error in FFS command n in FFS["…"]` を出して最上位の文を打ち切る。`LINE[key,"NAME"]=v` の未定義要素は `???General::invset`。pysad は続行していた (disp_b, chrom_b, matching21, dx_sv, dy_sv, amp_dep_ts_x_octa) | `fn_ffs` (`_RecOut` でエラー行を検出)、`set_line`、`sadscript.run` の `... in 文`。<work-perf>/err_test.sad, err_test2.sad で SAD と照合 (DTHETA/disp/fit の未定義パターンと `???-Emittance[]` は打ち切らない) |
| qcod の溢れ | SAD は壊れた値で続行。pysad は直前の twiss が無い (VISIT 直後など) と例外で Do を抜けていた (amp_dep_ts_dp の dp=−0.04) | `lattice.nan_rows` で NaN の表にして続行 (`ffscmd.calc`, `match6._stale_rows`) |

SAD の仕様として記録: 11_sad_quirks A40 (INTRA の 2 段報告)、A41 (PBUNCH 既定)、C25 (FFS::error)、C26 (qcod 溢れ)、C27 (`LINE[key,"FAM"]=v` の読み出し)、C28 (USE/VISIT は定義から作り直す)。

## 4. 結果

状態表は `docs/12_user_scripts_status.md` §3 (`status_table.py` の出力)。残った不一致の分類は同 §4。

| 対象 | 本数 | 数値一致 | 一部一致 | 動くが不一致 | 例外 | 打ち切り | SAD 側無効 | 対象外 |
|---|---|---|---|---|---|---|---|---|
| esrf | 52 | 34 | 12 | 0 | 1 | 0 | 4 | 1 |
| pf2 | 34 | 28 | 6 | 0 | 0 | 0 | 0 | 0 |
| 合計 | 86 | 62 | 18 | 0 | 1 | 0 | 4 | 1 |

(10/16: 数値一致 31、一部一致 26、動くが不一致 3、例外 16、打ち切り 10。11/05 までの個別修正の後の推定: 48 / 23 / 2 / 6 / 7)

一部一致 18 本の内訳: fit の多解・FITP>1 (chrom, matching3_01, matching_v30_12_r, matching_v31_23_r, tune_sv_v30_24_1_r)、不安定域の縁・閉軌道探索の初期値
(amp_dep_ts, amp_dep_ts_x, amp_dep_ts_y, amp_dep_ts_dp, error_ring_r)、DA の乱数のカオス (disp2, disp_test, da_xy_r, correct_cod_9_16_r)、
誤差系の要調査 (correct_cod_nomon, correct_cod_nomon_dq, print_cod_magnitude)、小さな数値差 (amp_dep_ts_dz)。

## 5. 要調査 2 件の結論 (同日)

- **correct_cod_nomon / _dq**: 誤差 (DELX ae=5E-5, DK fe=5E-4, DTHETA te=1E-4、GCUT=1) を入れた後の `cell calc` は SAD も pysad も閉軌道が見つからない
  (`*****qcod---> Overflow & closed orbit not found` が 85 回ずつ)。SAD は最後の `cell calc` で qtwiss も溢れ (`***qtwiss---> Overflow (fb)`) Twiss が 0 になり、
  pysad は直前の値 (NX 47.15) を残す。誤差なしの応答行列・CodCorr の部分と実行の流れは一致。SAD の壊れた値の再現はしない (方針)。
- **print_cod_magnitude**: 外側のループは ae = 5E-6 〜 1E-4 (20 点、実現 3 つずつ)。ae ≤ 5.5E-5 の 11 点 (33 実現) のチューンのずれは SAD と一致。
  ae ≥ 6E-5 はラティスが不安定 (結果に NaN、チューンのずれ 0.3〜0.6) で、SAD の値 (qcod 溢れ後の停留値 0 を含む) と pysad の値が違う。誤差の実現 (乱数列) と
  DELX 0 の乱数消費 (<work-perf>/del0_test.sad: SAD と同じ) は一致している。
→ 両方とも「不安定ラティスの扱いの差」(docs/12 §4) に分類し、記録のみ。
