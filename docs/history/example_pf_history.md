# example/pf の履歴 (整理で md から取り除いた記述)

元の文書の文脈つきで、取り除いた文章をそのまま入れてある。

## 共通: 「gs と pysad の比較」の冒頭の v2 の注記 (22 文書)

各文書の「gs と pysad の比較」の冒頭にあった段落 (全文書で同文)。現在の md では「**v2 モードと k64 モード**」の見出しの段落に言い直した (日付・再実行・出力の作り直しの記述を除いた)。

```
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。保存してあった旧い出力は `PYSAD_COMPAT=k64` の出力と同じ (確認済み)。v2 は K64 のバグを直したので (`docs/47_v2_changes.md`。主に `Emittance[]`・EMIT・追跡の Jx の端の項 B75 と `MomentumCompaction` B77)、この例の出力のうち `Emittance[]` の表示 (εx が数 % 大きい、εz・減衰率・α・バケット高さなど) と、それに依る量が gs と違う。以下の「一致」の判定は k64 モードの結果。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。
```

対象の文書:

- `pf/injection/calc_traj_01_PF.md`
- `pf/injection/calc_traj_02_PF.md`
- `pf/injection/calc_traj_03_PF.md`
- `pf/injection/calc_traj_03b_PF.md`
- `pf/matching/FFS_FIT_chrom_PF.md`
- `pf/matching/FFS_FIT_injbump_PF.md`
- `pf/optics/FFS_cell_calc_emit_PF.md`
- `pf/optics/FFS_disp_type_PF.md`
- `pf/optics/FFS_draw_OpticsPlot_PF.md`
- `pf/optics/FFS_draw_topdrawer_PF.md`
- `pf/optics/calc_chromaticity_PF.md`
- `pf/optics/calc_dispersion_PF.md`
- `pf/optics/calc_respx_PF.md`
- `pf/optics/fn_PolynomialFit_PF.md`
- `pf/optics/fn_Twiss_LINE_PF.md`
- `pf/tracking/fn_DynamicApertureSurvey_DAsv_PF.md`
- `pf/tracking/fn_DynamicApertureSurvey_PF.md`
- `pf/tracking/fn_DynamicApertureSurvey_da_xy_PF.md`
- `pf/tracking/phase_PF.md`
- `pf/tune_radiation/FFS_amp_dep_ts_PF.md`
- `pf/tune_radiation/FFS_mom_dep_ts_PF.md`
- `pf/tune_radiation/calc_RadInt_PF.md`

## pf/injection/calc_traj_03_PF.md

```
**○ 一致 (2026-11-17 に pysad を直した)**: 計算の値 (`GUI-start`・`after Weaker`) も、窓の Canvas の項目数 (114、再描画後 115) も gs と pysad で一致する。`.pysad.txt` にエラー表示は無い。以前は窓の部分が 2 つの理由で動かなかった:
1. `Canvas$Widget = c1` のあと Tcl が `Tkinter::tclerror: 1501 (invalid command name "plotframe")` で失敗した。原因は 2 つ: (a) `e=Entry[w]`・`c1=Canvas[plotframe,…]` の値が `Widget[…]` にならず (`Tkinter.n` の入れ子の規則で `With` が呼び出し元の引数名を見ていた)、(b) `Canvas$Widget = c1` の UpValue (Canvas.n の `(Canvas$Widget=Widget[c_,_Canvas])^:=…`) が、遅延した定義を持つ記号 `The$CanvasObject := …` を照合のために評価して、既定の Canvas の窓を作っていた。いまは gs と同じく、描く先が `c1` になる (docs/22)。
2. `OpticsPlot` に自前のデータ列を渡すと `General::equalleng: Equal lengths of elements in Thread are expected` で失敗した。原因は、OpticsPlot.n の `OpticsList1$[pr_,{f_}?(Length[Dimensions[#]]>=3&),opt___]` を、`OpticsList1$[{pr1_,pr2_},f_,opt___]` より先に試すこと (gs は定義の順。pysad は特殊さの点数で後にしていた) と、`Legend -> False` を `MapThread` に渡したとき、リストでない要素を全部に繰り返すこと (gs の `MapThread[List,{{1,2},3}]` は `{{1,3},{2,3}}`)。試験 `test_opticsplot_with_own_data_lists_like_gs` (自前データの x は要素番号で、この例は `TrackParticles` が返す位置番号)。
```

## pf/injection/calc_traj_03b_PF.md

```
`calc_traj_03_PF` (教材の書き方) は、2026-11-17 に pysad でも窓が描けるようになった (そちらの `.md` 参照)。
```

```
`OpticsPlot` のラティスの絵は、以前の pysad ではデータ列と一緒に描けなかった。2026-11-17 の修正のあとは `calc_traj_03_PF` で描ける。この例ではラティスの絵を付けていない。
```

## pf/optics/FFS_draw_OpticsPlot_PF.md

```
pysad は v2.4.1 から `draw` を matplotlib で描く (保存先 `PYSAD_DRAW_FILE` か対話表示。docs/42 §7b。ここの `.pysad.txt` は v2.4.0 までの実行で、未対応のメッセージが残っている)。
```

## pf/lattice/PF_history.md

```
- 1981 の台本は、元のままでは pysad が二重の反転 `-RQUAD` (中に `-RQUADM`) を誤って展開し、リングが不安定になっていた (gs は元から正しい)。`pysad/parser.py` の `-LINE` の展開を直した (反転の文脈の中の `-LINE` は正順に戻す)。直したあと k64 で gs と一致する。K64 側の誤りではなく pysad 側の不具合。docs/21・docs/47 への記載は保留 (他の担当が編集中)。
```

