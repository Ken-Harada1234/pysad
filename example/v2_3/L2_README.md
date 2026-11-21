# example/v2_3: 担当 L2 の例 (局所運動量アクセプタンス (LMA) を信頼して使える状態にする、v2.3.0、2026-11-21)

記録は `docs/history/v23_L2.md`、解説は `docs/48_frequency_map_lma.md` §4.5、試験は `pysad/tests/test_v23_L2_lma.py`、コードは `pysad/lma.py`・`pysad/lmabucket.py`・`pysad/q1cmd.py`。
実行はリポジトリ直下から `PYTHONUTF8=1 py -3.12 example/v2_3/L2_xxx.py ...` (標準出力を `.out.txt` に保存した)。すべて v2 モード、RF・放射つき、LOSSAMPL = ±10 mm (例の値で実機の値ではない)。
使ったラティス: `example/astrid/lattice/ASTRID2_v01_01.sad` (ASTRID2、0.58 GeV、45.7 m、6 セル)、`pysad/tests/kek-ls_v03_12.lat` (KEK-LS、大きなリング A、3 GeV、603 m、20 セル)、`ASTRID2_v01_02_id.sad` の `RINGMPW`。
elegant は WSL (`~/elegant/env.sh`)、serial の elegant を位置ごと・種ごとに同時に走らせた (`L2_elegant_isr.py`、timeout 1800 s、終了後に `ps` で残りなし)。gs は使っていない。
並列は PC 全体の CPU 70% 未満に収めた (`workers` は位置ごと、`L2_fluc.py` は種ごと)。`.lma.txt` は `LMAResult.save` の表 (最後の列 `flags`)。

| スクリプト | 項目 | 出力 |
|---|---|---|
| `L2_common.py` | ラティスの読み込み、T1 の elegant の表の読み込み | — |
| `L2_check_every.py astrid2\|kekls` | 1. check_every 0・1・5 の δ±・時間・elegant との差 | `L2_check_every_<名前>.out.txt`、`_ce0/1/5.lma.txt` |
| `L2_limit_types.py astrid2\|kekls` | 2. 解析式・1 次元の写像・6 次元の追跡のバケットの限界、位置ごとの種類 (rf・aper) | `L2_limit_types_<名前>.out.txt`、`.table.txt`、`.lma.txt` |
| `L2_verify.py astrid2\|kekls [n]` | 3. `verify=n` で狭い損失の窓を数える | `L2_verify_<名前>.out.txt`、`.lma.txt` |
| `L2_turns.py` | 4. 大きなリング A の周回数 491・1000・1500・2000、1500 周の elegant との比較 (bucket 判定あり・なし) | `L2_turns_kekls.out.txt`、`L2_turns_kekls_n*.lma.txt` |
| `L2_turns_short.py` | 4. シンクロトロン周期の 0.25〜3 倍での収束 (ASTRID2・大きなリング A) | `L2_turns_short.out.txt` |
| `L2_cell.py astrid2\|kekls [workers] [tol]` | 5. 全セルの同じ位置の LMA (1 セルの結果を全周に延ばす誤差) | `L2_cell_<名前>.out.txt`、`L2_cell_astrid2_tol5e-6.out.txt` |
| `L2_fluc.py` / `L2_elegant_isr.py` | 6. FLUC つき (6 種、種ごとに並列) と elegant の ISR つき (6 種) | `L2_fluc_astrid2*.out.txt`・`.json`、`L2_elegant_isr_<名前>.lma.txt`・`.out.txt` |
| `L2_errors.py errors\|mpw` | 7. 誤差つき ASTRID2 (3 種)・MPW 版の LMA、周期性の警告 | `L2_errors_astrid2.out.txt`、`L2_errors_mpw.out.txt`、`L2_errors_seed*.lma.txt` |

## 結果の要約 (数値は `docs/history/v23_L2.md`)

- **`check_every` の既定を 1 (要素ごとの判定) に変えた**。SAD の tapert (15 要素ごと、0) は、アパーチャが限界の位置で δ を最大 1.5 倍大きく見積もる (ASTRID2 の位置 17: 1.05% 対 0.83%)。elegant との差は 1 で 3.5e-5・7.0e-5。時間は ASTRID2 で +7%、大きなリング A で ×2.4。
- LMA の値 = min(RF バケットの限界、アパーチャの限界)。ASTRID2 (±10 mm) の 25 位置はすべてアパーチャで決まる (T1 が正しい)。RF バケットの δ+ が δ_rf を 1〜2% 超えるのは、放射の減衰・α2・ベータトロン振幅による (1 次元の写像と 6 次元の追跡が 1e-4 以内で一致)。
- 狭い損失の窓: ASTRID2・大きなリング A とも 0 個 (`verify`)。周回数: シンクロトロン振動 2 周期で tol 以内、3 周期で収束、1 周期以下は危険 (警告)。
- 1 セルの結果を全周に使う誤差は tol・U0/E0 以下。誤差つき・MPW は周期性が崩れることを検出して警告し、全周で走る。FLUC の効果は 2 周期では二分の幅以下。
