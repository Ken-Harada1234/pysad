# 作業ログ: 性能 3 (4D 閉軌道の 1 周を numba 1 回に) (2026-11-05)

最終更新: 2026-11-16

残課題「性能の続き (任意): 4D 閉軌道の 1 周を numba 1 回に (誤差入りの cell calc emit 1 実現 3.5 s → 1 s 程度の見込み)」。
方針は work_perf2.md と同じく **結果を変えずに速くする**。実行条件も同じ (FEL と同居、IDLE 優先度、スレッド 1)。

## 1. 作業前の状態

work_perf2.md §3 のプロファイル: 誤差入りの cell calc emit 3 実現 14.5 s のうち 4D の閉軌道 (`_qcod_periodic` → `_one_turn_cod_matrix`、58 周) が 8.6 s。
要素カーネル (elkern) 自体は 1 要素数 µs で、1 周 0.04〜0.15 s のほぼ全部が `_el_cod_step(_cached)` の要素ごとの Python 処理
(getattr、キャッシュの照合、np.array の生成、6×6 の matmul)。誤差入りの 3563 要素のリングで測り直すと 1 周 41 ms (キャッシュが外れる Newton 反復中)、qcod 1 回 0.29 s。

## 2. 対策: `pysad/turnkern.py` (`turn4`) と `lattice._TurnPack`

- 要素のパラメータを配列に詰める (`lattice._pack_el`: 種別 kind、設計軌道で近道できないか has_off、パラメータ P (n,16)、線形要素の行列 Mlin とキック Klin、
  設計軌道の線形行列 T0)。種別分けは `_el_cod_step` の分岐そのまま:
  DRIFT (RFSW OFF の空洞を含む) → qtwiss の DRIFT の式、薄い多極・薄い四極 → `elkern.k_thin`、QUAD → `elkern._quad_full`、BEND → `elkern.k_bend`、
  それ以外 (MARK/MONI/SOL/空洞 RFSW ON…) → 転送行列 + 定数キック、MULT と F1 付き BEND は Python。
- `turn4` が要素 i0 から最後まで 1 周を回し、要素行列 mats[i]・入口 cod cods[i]・一周の積 Mtot・出口 cod を埋める。
  Python に戻すときは (i, why) を返す: why=1 は設計軌道 (delta=0, cod=0, ずれなし) の線形行列 T0 が未計算 (Python が `el.transfer_matrix(0)` で作って再開)、
  why=2 は要素 i を Python (`_el_cod_step_cached`) で処理してから i+1 で再開 (MULT、解析ヤコビアンが作れない状態など)。
- pack は (要素数, delta, RFSW) ごとに保持し、`_qcod_periodic` が Newton 反復の前に 1 回だけ検証する (`_get_turn_pack`):
  要素ごとに `__dict__` の写しと比べ、別のオブジェクトか属性が変わった要素だけ詰め直す (3563 要素の検証 ≈ 4 ms)。
  属性に配列や辞書を持つ要素は毎回詰め直す。`_EL_CACHE` (要素ごとのキャッシュ) は Python に戻った要素でだけ使う。
- 結果: 要素行列・入口 cod・出口 cod は Python 版と **同一** (同じカーネル・同じ式)。一周の積 Mtot だけ 6×6 の積の順序 (numpy の BLAS と numba の `_mm`) で相対 1e-13 の差。
  設計軌道の近道 (T0) と −0 の扱い (DRIFT の `transfer_matrix(0)` と RFSW OFF 空洞の `_drift_matrix_at_cod`) も Python 版と同じにした。
- 切り替え: 環境変数 `PYSAD_NO_TURNKERN=1` (または `lattice._USE_TURNKERN = False`) で要素ごとの Python ループに戻る。`PYSAD_NO_ELKERN=1` のときも使わない。

## 3. 検証

- `pysad/tests/test_elkern.py` に追加: 要素 11 種 + MARK/MONI + MULT + K1=0 QUAD の列 (Python に戻る要素を含む) で delta=0/1e-3、cod 4 個 (設計軌道 0 を含む) の
  要素行列・cod が同一、Mtot が 1e-12; 要素を書き換えたら pack が詰め直されること; 誤差入りの 3563 要素のリングの 1 周 (3 条件) と qcod (Newton 反復込み) の一致 (1e-9)。
- pytest: 246 passed (19 s; 新しい 3 件を含む) (test_bend_feeddown/test_coord/test_sext_radiation/test_track6n_ext は従来どおり除外)。`QUICK=1 run_regression.py`: FAIL 0 (50 ベンチ、116 s)。

### 速度 (FEL と同居・IDLE)

| 対象 | 性能 2 の後 | 今回 |
|---|---|---|
| 誤差入り 4D 1 周 (3563 要素、Newton 反復中) | 41 ms | 6 ms |
| qcod 1 回 (反復 7 回) | 0.29 s | 0.04 s |
| マイクロベンチ mb_emiterr (誤差なし emit + 誤差入り calc + emit) | 2.8 s | 0.9 s (calc 誤差入り 0.13 s、emit 誤差入り 0.20 s) |
| 誤差入りの cell calc emit 300 実現 | 2343 s | **221 s** (1 実現 0.74 s、SAD は 37 s)。300 実現のうち 292 が SAD と一致 (相対 4e-7 以下)、残り 8 実現は前回と同じ不安定なもの (記録のみ) |

## 4. 変更したファイル

- 新規: `pysad/turnkern.py`。
- `pysad/lattice.py`: `_el_has_offset` (旧 `_has_offset` のラムダを関数に)、`_TurnPack`/`_pack_el`/`_get_turn_pack`/`_design_matrix`、`_one_turn_cod_matrix` (pack 引数、カーネル + Python 再開のループ)、`_qcod_periodic` (pack を 1 回作る)。
- `pysad/tests/test_elkern.py`: `test_turn_kernel_matches_python_line`、`test_turn_kernel_matches_python_(誤差入りのリング)`。
- 作業フォルダ (リポジトリに含まれない): 差し替えの記録 (patch3.py・patch4_resume.py)、順に 1 本ずつ走らせる補助スクリプト、出力とログ。

## 2026-11-16 の確認・追記

- 2026-11-16 確認: 変化なし。根拠: この話題の記述は work_unimpl.md と docs/22〜35 に出てこない (grep で確認)、関係するモジュールに 11/06 以後この話題の変更の記録がない。この文書は 11/06 時点の記録のまま。
