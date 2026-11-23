# 19_source_analysis から取り除いた履歴の記述

`docs/19_source_analysis.md` を現在のスナップショットだけの解説にするため取り除いた文章。元の文脈 (節の見出し) つきで、文章はそのまま。

## # 19. SAD ソース解析 総合リファレンスと pysad 対応表

取り除いた部分: 「(現状のスナップショット。K64 ソースの構成 (ファイル数・行数) は 11/06 に数えたまま。pysad から呼ばれる部分の変化 (COORD・TCAVI・UND・MULT ANGLE≠0・SOL・MAP・WAKE・BEAMBEAM・POL・CSR・RADLIGHT など) は 11/16 に直した)」

元の行:

> 最終更新: 2026-11-17 (現状のスナップショット。K64 ソースの構成 (ファイル数・行数) は 11/06 に数えたまま。pysad から呼ばれる部分の変化 (COORD・TCAVI・UND・MULT ANGLE≠0・SOL・MAP・WAKE・BEAMBEAM・POL・CSR・RADLIGHT など) は 11/16 に直した)

## ## 11. 未実装・未解決の解説

取り除いた部分: 「11/06 までのベンチ (Level 1〜31) の値で、11/16 に足した行の桁は各要素の文書 (docs/22〜37) の gs との比較による。」

元の行:

> 11. **この文書の §5・§7 の照合桁**は 11/06 までのベンチ (Level 1〜31) の値で、11/16 に足した行の桁は各要素の文書 (docs/22〜37) の gs との比較による。実用のリングや強い非線形での確認は別。

## ### 5.5 インタープリタ 🔵C / グラフィックス 🟢D / C ソース・sim・inc・Packages・extensions・contrib

取り除いた部分: 「11/06 は「pysad では不要 (描画は matplotlib: `plot_optics.py`)」としたが、**11/16 の現状は例外あり**: GUI」

元の行:

> - 🟢D: tfTkInter_.c、tfXlib_.c、sad_tcltk.c、sad_xlib.c、geodrw.f、Packages の Plots 系。11/06 は「pysad では不要 (描画は matplotlib: `plot_optics.py`)」としたが、**11/16 の現状は例外あり**: GUI (Tkinter.n・KBFrame.n) と描画 (ListPlot・Plot・OpticsPlot など) は SAD 自身の .n を pysad のインタプリタで解釈し、tkinter のキャンバスに描く (`tclbridge.py`・`tkcanvas.py`・`pkgload.py`。docs/22・35)。tfTkInter_.c のキャンバス関数は `tkcanvas.py` が tkinter の Tcl に置き換える。matplotlib の描画 (`plot_optics.py`) は別に残る。

## ### 5.5 インタープリタ 🔵C / グラフィックス 🟢D / C ソース・sim・inc・Packages・extensions・contrib

取り除いた部分: 「11/06 は init.n / emit.n / Fit.n / dapert.n / Optics.n の挙動だけを参照していたが、11/16 は `pkgload.py` が」

元の行:

> - Packages: init.n (NPARA=4 など既定値)、Optics.n (CALC6D, NormalCoordinates)、emit.n (Emittance[] の規則名、TouschekLifetime)、dapert.n、Fit.n、correction.n、optimize.n、extmap.n。11/06 は init.n / emit.n / Fit.n / dapert.n / Optics.n の挙動だけを参照していたが、11/16 は `pkgload.py` が correction.n・bump.n・corutility.n・optimize.n・D.n・Plots.n・ListPlot.n・CSR.n などを SADScript のまま読む (`PYSAD_SAD_PACKAGES` が要る。docs/34)。

## ### 7.2 tturne1 (tturne.f select case L470-675)

取り除いた部分: 「UND・BEAMBEAM は 11/16 から ✓」

元の行:

> | — | default | 無作用 (WIG, MARK)。UND・BEAMBEAM は 11/16 から ✓ (`und.py`・`beambeam.py`) | ✓ |

## ### 9.2 中 (将来)

取り除いた部分: 「✓ (11/16 に済み。INS は追加の実装が不要)」

元の行:

> | ttcav/ttcave, tcoord/tcoorde, temap, tins (追跡) | ✓ (11/16 に済み。INS は追加の実装が不要) | docs/25・26・30 |

## ### 9.2 中 (将来)

取り除いた部分: 「SOL・tmultae は ✓ (11/16)。」

元の行:

> | emit6d の SOL、tmultae (MULT ANGLE≠0)、MULT の SK/F1 | SOL・tmultae は ✓ (11/16)。MULT の SK0・F1・TAPER・VOLT は ✗ | docs/28・29 |

## ### 9.3 低・不要

取り除いた部分: 「自前実装。11/16 は SAD の .n パッケージもそのまま動かす)」

元の行:

> 🔵C インタープリタ (sadscript で自前実装。11/16 は SAD の .n パッケージもそのまま動かす)、🟢D グラフィックス (11/16 は SAD の .n を tkinter で動かす。docs/37)、⚫E ビルド、🟤G contrib、デッドコード 10 件。

## ### 9.3 低・不要

取り除いた部分: 「(11/16 は SAD の .n を tkinter で動かす。docs/37)」

元の行:

> 🔵C インタープリタ (sadscript で自前実装。SAD の .n パッケージもそのまま動かす)、🟢D グラフィックス (11/16 は SAD の .n を tkinter で動かす。docs/37)、⚫E ビルド、🟤G contrib、デッドコード 10 件。

## ## 10. pysad 現状

取り除いた部分: 「- 11/16 の追加 (要素・機能別の文書は docs/22〜37、大項目の解説は docs/38〜45): COORD」

元の行:

> - 11/16 の追加 (要素・機能別の文書は docs/22〜37、大項目の解説は docs/38〜45): COORD・TCAVI・UND・MULT ANGLE≠0・MAP・DELCOR・SOL の Twiss/Emittance[]・WAKE・BEAMBEAM・POL (スピン)・CSR・RADLIGHT・GaussianCoulomb・特殊関数・SAD パッケージ・GUI/描画。回帰は 50 ベンチ FAIL 0 (`QUICK=1`)。
