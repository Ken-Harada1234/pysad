# 15_csr の検証の詳細 (gs との比較)

docs/15_csr.md の「検証」の詳細。旧 docs/40_csr_guide.md §5 と旧 docs/36_spin_csr.md §2 (CSR) の照合の記録を元の文脈つきで移した。

## 旧 40 §5 gs との一致の度合い


すべて SAD K64 (gs、2026-10-19 ビルド) の出力との比較 (試験は [`pysad/tests/test_csr.py`](../pysad/tests/test_csr.py))。

| 項目 | 一致 |
|----|----|
| `CSRTrack` (RF あり・なし、発散後 0、z・dp のその場の更新) | 相対 1e-11 |
| `CSRHaissin` | 相対 1e-11 |
| `CSRSetupOY` の J・ω・rhoj | 1e-9〜1e-7 (周回の離散計算のため) |
| `CSRMatrix` mode 22・33・23・32・110・220・330 | 全要素 |
| `CSRConvert` | 一致 |
| `CSRInit` (pipe 8 列と 5 列) | 相対 1e-10 |
| `CSRMatrix` mode 1000 ({2,2,2} と {5,5,5,5,5,5}) | 相対 1e-12〜1e-13 |
| `CSR.n`: `SkinDepth`、`EffectiveSkinDepth`、`EffectiveCond`、`PathDifference`、`StandardKList`、`VariableMesh`、`ZThre`、`Tweak` | gs と一致 |
| `CSRData` (8 列のメッシュ) の固有値 | 相対 1e-15 (固有ベクトルは縮退で並びが違うが物理量は同じ) |
| `Ez` | 相対 5e-10 |
| `MakeZL` (`{{10.,1.}}`、2000〜6000 /m の 3 点、PipeWidth 0.04、BunchLength 3e-4) | 相対 1e-10 以内 |

この文書の作成時に `py -3.12 -m pytest pysad/tests/test_csr.py` を `PYSAD_SAD_PACKAGES` 付きで実行して 6 通過・1 skip (skip は `PYSAD_SLOW=1` が要る MakeZL) を確かめた。

試験で使った量は小さい値の人工的な ZL・断面で、実機のパラメータ (エネルギー、実際のダクト、周長) での物理的な妥当性の検証ではない。ZL の `MakeZL` の照合は 3 周波数・1 区間だけである。


## 旧 36 §2 CSR の関数ごとの状態


`CSRTrack`、`CSRHaissin` (csrtrack.f)、`CSRSetupOY`、`CSROYMatrix` (csroy.f)、`CSRMatrix`、`CSRConvert` (csrtest.f の一部)。実FFT (`tcftr`・`trftr`・`tftrr`) と自然スプライン (`spline1`) も写した。

| 関数 | 状態 |
|----|----|
| `CSRTrack[z, dp, zl, {dz, ndiv, inten, r56, r65, frf, phirf, damp, sigmae}, nturn]` | gs と 11 桁一致 (乱数は ffs と共通の SadRandom、z・dp はその場で更新、発散後は 0、RF あり・なし) |
| `CSRHaissin[zl, {dz, inten, r65, sigz}]` | gs と 11 桁一致 (Newton 反復。tsolvg は pysad/tsvd.py) |
| `CSRSetupOY[rho, zl, {dz, nj, nphi, inten, r56, r65, zspan}]` | gs と J・ω 1e-9、rhoj 1e-9 (周回の離散計算のため)。dpjn の記憶域の重なり (B59) を再現 |
| `CSROYMatrix[...]` | ソースどおりに翻訳。手で作った入力 (nr=64、nj=3、nphi=2) の行列は gs と 4e-20 で一致。`CSRSetupOY` の出力を渡すと gs は wrongtype (原因は特定: B60) |
| `CSRMatrix[pipe, {dx,dy,f0..f3,tilt}, omega, rho, center, mode]` | mode 22・33・23・32・110・220・330 は gs と全要素一致。mode 1000 (らせん断面 `csrhelmat`) は実装済みで 2 断面が相対 1e-13 で一致 (§2.1) |
| `CSRConvert[v or m, pipe]` | gs と一致 (K64 の tfcsrvec と同じ) |
| `CSRInit[{{dx,dy}, center, rho, pipe}]` | 実装済み、gs と相対 1e-10 (§2.1) |


## 旧 36 §2.1 CSRInit・mode 1000・CSR.n の流れ (照合の記録)


- `CSRInit[{{dx,dy}, center, rho, pipe}]` は実装済み (csrtest.f の csrphimesh・csr2dpoisson。境界の電荷は `tsolva` = `tsvd(.false.)` の最小ノルム解で、pysad の `tsvd` を使う。
  `csr2dpoisson` の電場は 3 回の繰り返し)。gs と全要素の和・最大値・個別の値が相対 1e-10 で一致 (pipe 8 列と 5 列の 2 例)。実数 1 つの dx は K64 のバグ
  (dx を読まず 0 のまま) で常に wrongtype。試験は [`test_csr.py`](../pysad/tests/test_csr.py) の `test_csrinit_matches_gs`。
- `CSRMatrix` mode 1000 (らせん断面。`csrhelmat` + `csr2dpoisson1`) は実装済み。境界の点列を中点で 2 倍にして境界要素法で各セルの電場への応答 (2nm×2nm) を求める。
  `tsolvm(svd=.true.)` は最小ノルム解 (特異値 < 最大×1e-4 を捨てる) として numpy の SVD で解いた。K64 は `phix(0:nphimax)` を `csrphimesh` の `phix(*)` に渡すため、
  点の添字が 1 つずれる (これも写した)。gs と {2,2,2} と {5,5,5,5,5,5} の 2 つの断面で全要素の和・個別の値が相対 1e-13 で一致。{4,4,4,4} などは gs が SIGABRT/SIGSEGV で落ちる (参照なし)。
- `CSR.n` の流れが動く: `CSR[…]` の生成、`MakePipe[]`、`CSRData` (CSRMatrix → RealEig → CSRInit → CSRConvert)、`Ez`、`MakeZL` (周波数の表からインピーダンス)。
  8 列のメッシュの `CSRData` は固有値が gs と相対 1e-15、`Ez` が相対 5e-10 で一致 (固有ベクトルは縮退で並びが違うが物理量は同じ)。
  `CSR[NParallel->1, Nk->4, BunchLength->3e-4, BendingRadius->10, PipeWidth->0.04, Center->0.02]@MakeZL[{{10.,1.}}, 2000., 6000., 3]` の ZL の表 (3 周波数) は
  gs と相対 1e-10 以内で一致 (pysad で約 7 分、試験は `PYSAD_SLOW=1` で動く `test_csr_makezl_matches_gs`)。
  動かすための対応: `Do[…, {Infinity}]`・`{i, Infinity}` (Break するまで続ける)、`OpenRead["!コマンド"]` (コマンドの出力を読む。`top` のない環境 (Windows) では
  `FreeMem` 用に空き記憶領域を旧形式の top の行で作る)、開けないファイルは `$Failed`。
- gs の参照を取るには、WSL に `top` を入れるか、PATH に `top` の代用 (`Mem: …k total, …k used, …k free, …k buffers` を出すスクリプト) を置き、`OSTYPE=linux` を設定する
  (非対話の WSL では OSTYPE が空で、`FreeMem` の `Switch` が未評価になり OpenRead が失敗する)。
- `Packages/CSR.n` (1088 行) は `pkgload.py` の許可リストにあり、読み込める (描画の `ListPlot.n`・`Plots.n` なども許可リストにある)。
- **未照合 (gs と比べていない)**: `CSRSetupOY` の出力から作る `CSROYMatrix` (gs が wrongtype を返す: B60。手作りの入力は gs と一致)、`CSR.n` の `ScanTrack`・`OYStability` などの長い流れ、実用の大きさのメッシュの `MakeZL` (試験の 3 周波数・2000〜6000 の範囲以外)、
  mode 1000 の {4,4,4,4} など gs が落ちる断面 (参照なし)、実数 1 つの `CSRInit` (K64 のバグで常に wrongtype)。`CSR.n` の `SkinDepth`・`PathDifference`・`StandardKList`・`VariableMesh`・`ZThre` は gs と一致 (docs/history/work_unimpl.md の記録。試験は MakeZL だけ)。

- 2026-11-16 夜にソースで確認: `CSR.n` の `TrackAnimate` は `ListPlot` を呼んだあと `System["import -silent -window … -format -jpeg …"]` (ImageMagick) と `System["mkdir …"]` を呼ぶので、Windows の pysad (System は実行しない) では画像は作れない。
  `Wait[mem]` は `FreeMem[]` (top) と `Sleep` を使う。gs はこの環境で top が無く FreeMem が OpenRead エラーになる。pysad の `FreeMem` は Python (psutil) で代替して 1.0 を返す。並列を増やすと速くなるかは確かめていない (pysad は fork しない)。

K64 の癖は B59・B60。
