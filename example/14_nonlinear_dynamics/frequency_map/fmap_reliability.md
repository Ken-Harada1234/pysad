# 周波数マップの信頼性の例

解説は `docs/26_frequency_map_lma.md` (§3.2・3.3・3.6・§10)、台本ごとのまとめは [`FMAP_summary.md`](FMAP_summary.md) §14、試験は `pysad/tests/test_v23_S2_fmapcheck.py`、コードは `pysad/fmapcheck.py`・`freqmap.py`・`fmapres.py` (`apply_alignment_errors`)。
実行はリポジトリ直下から `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_xxx.py ...` (標準出力をファイルへ落とすときは `PYTHONUTF8=1`)。結果の数値は `example/14_nonlinear_dynamics/frequency_map/fmap_check_data/*.jsonl` (1 行 = 1 つの周波数マップ。走った組は飛ばして追記する)。
共通の条件: v2、RFSW・RAD OFF、LOSSAMPL = 0.1 m、px = py = z = 0 から出発、NAFF の窓 Hann (p = 1)、規格化 = 標本の楕円、1024 周、**格子 21 × 16 = 336 点** (41 × 31 の格子の 1 つおき)。ASTRID2: 位置 2、x = 0.5〜18 mm、y = 0.5〜14 mm。KEK-LS: 位置 3、x = 0.2〜6 mm、y = 0.2〜8 mm。
使ったラティス: `example/astrid/lattice/ASTRID2_v01_01.sad`・`ASTRID2_v01_02_id.sad` (MPW)、`pysad/tests/kek-ls_v03_12.lat` (KEK-LS、大きなリング A、ESRF 型の HMBA)。
CPU は 4 コアまで (ProcessPoolExecutor で 1 つの周波数マップ = 1 プロセス、3 並列)。1 つの周波数マップ (336 点・1024 周) は ASTRID2 で 110〜150 s、KEK-LS で 90〜200 s (他の計算が同時に走っていた)。

| 項目 | ファイル | 内容・コマンド | 結果 |
|---|---|---|---|
| 1 誤差つき | `fmap_check_errors_scan.py` | `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_errors_scan.py astrid2 --seeds 1-8 --levels 0.3,1,3 --kinds all --jobs 3` (種類ごと: `--seeds 1-5 --levels 1,3 --kinds quad_dxy,sext_dxy,bend_rot,quad_rot,bend_dxy`、フィードダウンを切る: `--kinds all_cs`、KEK-LS: `kekls --seeds 1-5 --levels 0.03,0.1,0.3,1 --kinds all`、`kekls --seeds 1-3 --levels 0.1 --kinds quad_dxy,sext_dxy,oct_dxy,bend_rot,quad_rot,bend_dxy,all_cs`) | `errors_astrid2.jsonl`・`errors_kekls.jsonl` |
| 〃 集計 | `fmap_check_errors_summary.py` | `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_errors_summary.py astrid2` (または `kekls`) | `fmap_check_errors_summary_astrid2.out.txt`・`fmap_check_errors_summary_kekls.out.txt` |
| 〃 結合だけ | `fmap_check_coupling.py` | 1 台の QF3 のロールだけで |c−| = 0.0003〜0.01 | `coupling.jsonl` |
| 〃 KEK-LS の閉軌道 | `fmap_check_kekls_orbit.py` | 補正なしの誤差 ×0.03〜×3 × 種 1〜5 の閉軌道・安定性・結合の見積もり | `fmap_check_kekls_orbit.out.txt` |
| 2 MPW | `fmap_check_mpw.py` | `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_mpw.py --jobs 3`: MPW の強さ s = 0〜1、νy を ±0.01〜0.05 (四極 2 家族の Newton 法)、νx・νy を MPW なしに戻す | `mpw.jsonl`、`fmap_check_mpw.out.txt` |
| 3 RFSW・FLUC | `fmap_check_rfsw.py` | `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_rfsw.py --jobs 3`: base・rf・rfrad・rfradfluc (種 17・29)・radfluc | `rfsw.jsonl` |
| 〃 gs | `fmap_check_gs_rfsw.py` | `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_gs_rfsw.py` (gs は WSL、timeout 1700 s、終了後に sad1.exe の残りを確認): RFSW ON・RAD OFF の 6 点、512 周 | `fmap_check_gs_rfsw.{sad,gs.txt,pysad.txt,out.txt}` |
| 4 K1 の設定 | `fmap_check_k1opts.py` | `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_k1opts.py --jobs 2`: default・face・lin・face_lin・k64 | `k1opts.jsonl`、`fmap_check_k1opts.out.txt` |
| 5 周回数・窓・格子 | `fmap_check_nturn.py` | `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_nturn.py astrid2` (`kekls`、`--errors 1,1`): 1 回の追跡 4096 周から N = 256〜4096、窓 p = 0〜3。点ごとの傾き | `fmap_check_nturn_*.out.txt`、`fmap_check_data/nturn_*.npz` |
| 〃 格子 | `fmap_check_grid.py` | 保存済みの 41 × 31 の表から 1〜5 つおきの部分格子 | `fmap_check_grid.out.txt` |
| 7 gs・誤差つき | `fmap_check_gs_errors.py` | `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_gs_errors.py`: pysad の誤差の値を実体ごとの要素定義として gs の台本に書く (`LINE["DX", n] = 値` は同じ名前の全実体に効くので使えない)、6 点、512 周 | `fmap_check_gs_errors.{sad,gs.txt,pysad.txt,out.txt}` |

台本からの警告: `FrequencyMap[]` は既定で RFSW・RAD OFF (現在のフラグは使わない。`UseFlags->True` で使う)、警告は ` FrequencyMap: WARNING [コード] …` で出る (`test_script_frequency_map_safe_default`)。

再現の確認: `py -3.12 -m pytest pysad/tests/test_v23_S2_fmapcheck.py -q` (約 10 s。保存済みの出力を読む試験を含む)。gs の再走 (`fmap_check_gs_errors.py`・`fmap_check_gs_rfsw.py`、各 1 分以内) は保存済みの出力と一致する。
