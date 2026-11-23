# 16_fortran_status.md から取り除いた履歴の記述

16_fortran_status.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 16. SAD Fortran ファイル別 pysad 翻訳状況

整理前の 3 行目:

````text
最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット。内容は `history/work_unimpl.md` と `pysad/*.py` の docstring、v2 の作業の記録 `history/v2_changes_*.md` で裏を取った。§9 は [46_unimplemented_and_unverified.md](46_unimplemented_and_unverified.md) と同じ状態)
````

整理前の 5 行目:

````text
**互換モード**: 表の「状態」は、特記がなければ **k64 モード** (K64 (gs) と同じ動作。gs との照合はこのモード) のもの。v2 モード (v2.0.0 から既定) で変わった動作は「v2 モード:」と書いた (一覧は [47_v2_changes.md](47_v2_changes.md))。v2 の値は gs ではなく解析式・独立な数値計算で確かめたもので、gs と比較できないものがある (MEASURE・TRCOD、空間電荷の SPAC、領域内の SOL の要素など)。
````

## ## 0. 凡例と現在の pysad モジュール

整理前の 23 行目:

````text
v2.0.0 で追加: compat (互換モード) / phsrot / radtaper / spch (空間電荷)。`measure` は v2 で MEASURE・TRCOD ほかの実装を足した。
````

## ## 1. 粒子トラッキング

整理前の 60 行目:

````text
| spkick.f、tspch_.c・spch.c・psn.c、tintrb.f twspac | 空間電荷 SPAC・PSPAC・WSPAC | `spch.py` | **v2 モードで翻訳** (k64 モードは入れない)。K64 の SPAC は力が 1/γ 倍に弱い (B89) ほか B90〜B98 の誤りがあり、v2 は直した。PSPAC は gs と 3〜4 桁、WSPAC は gs と 1% (Δν)、SPAC は gs の誤りのため解析式 (包絡方程式・Laslett の式) と比較。Emittance[]・EMIT・CALC6D・CALC (RFSW ON) の WSPAC は v2.1.0 で実装 (Emittance[]・EMIT は k64 も gs と全桁。K64 の最後のスライスのレンズが 8 倍強い誤り: B127)。追跡の WSPAC のビーム行列は空間電荷なしのまま | test_v2_C_spch.py、test_v2_C2_spch.py、docs/08 §9.7 |
````

## ## 5. COD 補正

整理前の 114 行目:

````text
| DUMP コマンド | 要素の DK・DDK などの表 | ffscmd `_gwd` | 翻訳済み・照合済み (パターンは 20 文字まで、DDK は K0 のある要素だけ。誤差入りリングで値は gs と一致。DK・DDK 列の書式 (Fortran の 1p g12.5: 0 → `0.0000`、小さい値 → `4.22829E-04`) を gs に合わせた) | test_ffsextra.py `test_dump_format_matches_gs` |
````

## ## 8. 寿命・その他

整理前の 165 行目:

````text
| gamma.f、hg.f、tfunc*.f (特殊関数)、tffindroot.f、D.n | 特殊関数・Fit・FindRoot・NIntegrate | specfn.py、fitroot.py、miscfn.py | 458 式が gs と一致 (K64 が誤値を返す B43〜B47、B53、B54 は真値)。Fit・FindRoot は数値微分で 5e-10〜1e-9。NIntegrate は複素数の被積分関数を扱い (gs と虚部 1e-15)、`LogGamma1` もある。複素 FindRoot は v2 モードで実装 (gs が壊れた値: B61。k64 モードは TypeError)、複素引数の Bessel・超幾何・Airy (実数のパラメータ) と PolyLog の非整数次の大きい \|z\| (相対 1e-11、\|z\| ≤ 1e12) は両モードで実装。k64 と値が違う関数は InverseErf (B47) だけで、v2 は真値。複素引数の PolyLog・HypergeometricU・Zeta・HurwitzZeta は v2.1.0 で両モードで実装 (`specfn_cx.py`。gs が誤る点でも真値: B118〜B120) | docs/33 |
````

## ## 9. 未翻訳・一部未対応・未照合の一覧 (2026-11-18、v2.1.0 の後)

整理前の 168 行目:

````text
## 9. 未翻訳・一部未対応・未照合の一覧 (2026-11-18、v2.1.0 の後)
````

整理前の 171 行目:

````text
各行の「モード」は **k64** (K64 (gs) と同じ動作)、**v2** (v2.0.0 から既定の、K64 のバグを直し機能を足した動作)、**両** (モードによらない)。v2 の機能の根拠は解析式・数値積分で、gs と比較できないものが多い (gs が動かない・K64 が誤っている)。
````

整理前の 172 行目:

````text
実装・解消済みのもの (APERT 要素の追跡 (SOL 領域内を除く)、FFS 内の `MOMENTUM=`・`[NO]EXPAND`・RENUMBER・DISP の各モード・REFERENCE・DUMP・OrbitGeo、MakeBump・Couple の ElementValues・Rational・NIntegrate の複素数・LogGamma1・ParametricPlot・ListContourPlot・OptimizeOptics・Entry/ListBox/Scale のメソッド、LOSSMAP・RADTAPER・MULT のスキュー成分・ResponseMatrix の TRPT・NOKEEPEXP・RENUMBER の INDEX・`BeginPackage`/`Begin` の文脈・`Shared`。**v2.1.0 で**: Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE と ANGLE=0 の BEND・MULT の F1、光学・Emittance への空間電荷 (WSPAC)、光子 (tphrec) の BEND・QUAD・MULT・SEXT・OCT・ドリフト・BEAMBEAM + BEND、tdgeo の領域内の四極の DX・DY の解析微分、CALC6D・RADCOD つきの CALC・INTRA と SOL、傾いたソレノイドの POL の量、POL + TRPT + 加速 (領域の外)、`MASS`・`CHARGE`、`NPARA=n;`、複素の PolyLog・HypergeometricU・Zeta、`Parallelize` (POSIX)、成分ごとの FREE、`$Context` の代入・`WithContext`、VISIT/BYE、BEND の F1 + RAD の追跡・CAVI の DX・DY・EMIT のあとの USE の差) は、ここでは未対応に入れない (46 §6)。
````

## ### 9.1 未翻訳 (理由つき)

整理前の 178 行目:

````text
| wjfit、wfres、wfbin、wiord、wjdraw | 両 | MAIN の TRACK の Fourier 解析で、FFS から使えない (gs も SIGSEGV)。チューンの抽出は [`pysad/naff.py`](../pysad/naff.py) (v2.1.0) |
````

## ### 9.1a モードで状態が違うもの (k64 モードでは未対応・v2 モードで実装、または結果が違う)

整理前の 199 行目:

````text
| Emittance[] の POL の量 (傾いたソレノイド) | K64 の式 (B130) | v2.1.0 で修正 (追跡と 1e-13〜5e-7) |
````

整理前の 200 行目:

````text
| RADINT の K0 の極・四極のずれ・BMAG (B115・B126・B107)、ANGLE=0 の BEND の端の角 (B108)、tradke の符号 (B111)、WSPAC のレンズ (B127)、K1≠0 の BEND の追跡の端 (B117)、ずれた BEND の CALC (B114) | K64 のまま | 修正 (47 の v2.1.0 の (a)) |
````

## ### 9.3 未照合 (コードはあるが gs と数値比較していない、または gs が参照を出さない)

整理前の 237 行目:

````text
| v2 モードの機能全般 (MEASURE・TRCOD、SOL 領域内の要素、UND・PHSROT、空間電荷の SPAC、Polarization の有限値、RADLIGHT の 3 次元の軌道、v2.1.0 の修正の値) | v2 | 値の根拠は解析式・独立な数値計算。gs との比較は無い (gs が動かない・K64 が誤っている)。k64 モードの再現は gs と比べた |
````

## ### 9.4 優先度

整理前の 254 行目:

````text
| 中 | Polarization の追跡との差 (非線形の軌道)、v2 の Emittance[] と RADINT の矩形磁石の差、SynchroBeta の大きなリング、9.3 の未解決の差、v2 の要判断 (docs/47 v2.1.0 の (d)) | 原因の切り分け |
````

