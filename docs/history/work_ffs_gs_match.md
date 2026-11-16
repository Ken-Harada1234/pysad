# FFS・光学・ラティス側の gs との差の修正 (2026-11-17)

`example/01〜07` の比較で見つかった差を、gs の実動作を基準に直した。試験は `pysad/tests/test_ffs_gs_match.py` (10 件)。

## 直した点 (gs で確認)
1. FFS 内の `MOMENTUM=…`: `USE` のたびに MAIN 側の値で上書きしていたのが原因。MAIN 側の値が変わったときだけ取り込む。追跡・Emittance[] (SIGE) が gs と一致 (1 GeV: x 2.237707E-4、SIGE 3.916263E-4)。
2. `Twiss["LENGTH"]` は実装済みだった。実際に落ちていたのは `Twiss["GAMMA",1]` (tftwiss.f:271 で GAMMA・GAMMABETA・S・SIZE は LINE と同じ処理)。未知の名前のメッセージも gs の `General::wrongval` に。
3. `USE [NO]EXPAND`、`CALC|GO [NO]EXPAND`、`NOCALEXP`: calexp は KEEPEXP で残る。USE の NOEXPAND は構築済みのラインを再利用 (成分ごとの値が残る)。
   副産物: gs の TrackParticles は成分ごとの値を揃えない (K1・K0・DX を `LINE[…,"QF.2"]=` で 1 台だけ変えて追跡すると効く)。pysad は追跡の前に揃えていたのをやめた。CALC・Emittance[]・SAVE・EXPAND は揃える (gs で確認。SAVE と EXPAND は tfsave で同じ)。
4. `RENUMBER` と `LINE["MULT"]` (tfinimult の直訳。同名 1 台だけは 0)。
5. TYPE の K5 (DODECA): 長さあたりの内部値のまま返していた。`_PER_LENGTH` と get/set に K5 を追加。
6. `DISP` の GA・Z・P・A・B (BEAM)・D と REFERENCE・DREFERENCE、`REFERENCE`、`Twiss["R…"/"D…"/"RALL"/"DALL"]`。見出しの文字の置換は tfdisp.f の trefheader/tdrefheader のとおり。
7. `DUMP`: gs の `word1` が `character*20` で、20 文字で切れる (例の D1.1・B.1 が出ない原因)。DL の drift は gs でも入る。DDK は K0 キーワードのある要素 (BEND・MULT・TCAVI) のみ (QUAD には入れない)。合う要素がなければ全要素。
8. `OrbitGeo[$$$]` の GX: Optics.n の式が `Twiss["DZ","$$$"]` (一周の経路長のずれ -5.338e-5) を進行方向に足す。pysad は CELL の 4D で最後の行の DZ を 0 にしていた。`lattice._twiss_table_impl` で周回出口の DZ を入れる。
9. `QF.3 K1 0.9;`: `名前.番号` で始まる文は式として評価されて何も設定しない (pysad は family 全体を変えていた)。出力の `Out[n]:=` 行は再現しない。

## 直さなかった点
- `PolynomialFit` (x の平均が 0 のとき最高次以外の係数が壊れる): pysad は SAD の `Packages/Fit.n` をそのまま評価しており (`sadscript._PRELUDE` が同一のコード、`If[xa,…]` の `xa=0` が偽になる癖もそのまま)、SAD と同じ振る舞い。直さない。
- `NOKEEPEXP`: tfvcomp が未翻訳。
- RENUMBER の INDX 属性による初期値、`Twiss["R…"]` の小数位置。
- USE 後も EMIT の効果 (追跡の結果) が gs では残る (EMIT 後に `USE` して TrackParticles を呼ぶと gs は EMIT 前と違う値)。pysad は USE で emit を捨てる。原因未調査 (MOMENTUM の件を切り分けているときに見つけた)。
- `Out[n]:=` の表示。

## 確認
`py -3.12 -m pytest pysad/tests -q --ignore=…` 558 通過・3 skip。回帰 `QUICK=1 run_regression.py` FAIL 0。例の `.pysad.txt` を再実行して `.md` の比較を更新 (`FFS_USE_EXPAND`、`FFS_RENUMBER`、`FFS_DISPLAY`、`FFS_DUMP`、`FFS_REFERENCE`、`fn_Twiss`、`fn_LINE`、`fn_OrbitGeo_GeoBase`、`elem_SEXT_OCT_DECA_DODECA`、`FFS_BEAMSIZE`、`flag_SPAC` と各 README)。
