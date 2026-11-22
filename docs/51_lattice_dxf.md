# 51. ラティスの CAD 図 (DXF) — 平面図と磁石の立体 (3D) の書き出し (v2.4.1)

最終更新: 2026-11-22 (v2.4.1 開発中)

ラティスを上から見た図 (平面図) を、CAD で開けるテキストの DXF (ASCII、AutoCAD R12 相当) に書き出す。オプション (`mode='3d'`・`Mode->"3D"`) で、磁石の立体の外形 (基準軌道に沿った直方体、偏向磁石は入口・出口の面が E1・E2 に傾いた楔形) も同じ DXF の形式 (3DFACE) で書ける (§2d)。K64 の SAD に同じ機能はなく (v2 の新機能。互換モードによらず同じ動作)、
SAD の `3GSR/esrf/drawmagvw3d.n` (K. Oide 作。`DrawMagnets3D[]` が Vectorworks の VectorScript を書く) の幾何の考え方を参考にした自前版である。
Vectorworks が一般的でないので、出力を DXF にした。既定は平面図 (2D) で、3D は選んだときだけ。`drawmagvw3d.n` は参照しただけで、コピー・改変していない。

| 入口 | 場所 |
|---|---|
| Python | `pysad.dxf_lattice.write_dxf(ffs, path, region=None, scale=1.0, ..., mode='plan', parts='body')`、`load_ffs(台本, 'RING')`、`read_dxf(path)` (検証用の最小パーサ。3DFACE も読む)、`check_solids(entities)` (3DFACE の立体が閉じているかの検査) |
| SADScript (gs でも pysad でも走る) | `DrawMagnetsDXF[file, options]` ([`example/dxf_lattice/DrawMagnetsDXF.n`](../example/dxf_lattice/DrawMagnetsDXF.n)。`Get["DrawMagnetsDXF.n"]` で読む) |
| 例 | [`example/dxf_lattice/`](../example/dxf_lattice/README.md) (FODO・SOL・ASTRID2・KEK-LS の 1 セルと全周。3D は FODO・KEK-LS の 1 セルと FODO の全周、`*_3d*`) |
| 試験 | `pysad/tests/test_dxf_lattice.py` (28 件のうち 3D は 11 件。ezdxf があれば動く 1 件を含む) |

## 1. 使い方

```python
from pysad.dxf_lattice import load_ffs, write_dxf
ffs = load_ffs('example/00_lattice/FODO_ring.sad', 'RING')       # USE RING まで (CALC は不要)
write_dxf(ffs, 'ring.dxf')                                       # 全周、単位 m
write_dxf(ffs, 'part_mm.dxf', region=('QF.1', 'QF.4'), scale=1000, text_height=0.1,
          spec={'B': {'Width': 0.5, 'PoleWidth': 0.15, 'Type': 'Ctype1'}, 'QUAD': {'Width': 0.36}})
write_dxf(ffs, 'cell.dxf', text_height=0.3, overlap='remove')    # 要素名の重なりを消す (台本が MagnetSpec[] を定義していれば寸法も読む)
```

```
FFS;
USE RING;
Get["DrawMagnetsDXF.n"];
DrawMagnetsDXF["ring.dxf"];
DrawMagnetsDXF["part_mm.dxf", Region->{"QF.1", "QF.4"}, Scale->1000, TextHeight->.1,
               Spec->{"B"->{Width->.5, PoleWidth->.15, Type->Ctype1}, "QUAD"->{Width->.36}}];
DrawMagnetsDXF["cell.dxf", TextHeight->.3, Overlap->"remove"];
DrawMagnetsDXF["cell_3d.dxf", Mode->"3D", Spec->{"BEND"->{Height->.3}}];                 (磁石の立体、§2d)
DrawMagnetsDXF["cell_3d_split.dxf", Mode->"3D", Parts->"split"];                         (ヨークとポールに分ける)
```
```python
write_dxf(ffs, 'cell_3d.dxf', mode='3d', spec={'BEND': {'Height': .3}})
write_dxf(ffs, 'cell_3d_split.dxf', mode='3d', parts='split')
```

| 引数 (Python / SADScript) | 既定 | 意味 |
|---|---|---|
| `region` / `Region->{始点, 終点}` | 全周 | 要素の番号か `NAME.n`。**終点の要素も含む**。複数に合う名前は始点が最初、終点が最後 |
| `scale` / `Scale` | 1 | m からの倍率。1 → `$INSUNITS`=6 (m)、1000 → 4 (mm)、100 → 5 (cm)、他は 0 |
| `text_height` / `TextHeight` | 0.15 m | 要素名の文字の高さ。`names=False` / `Names->False` で要素名なし |
| `mark_size` / `MarkSize` | 0.1 m | MARK の線の長さ。`mark_text` / `MarkText->True` で MARK の名前も出す |
| `arc_step` / `ArcStep` | 0.05 rad | ANGLE≠0 の要素の軌道の折れ線の最大の角度 |
| `text_side` / `TextSide` | -1 | 要素名を置く側 (-1: 進行方向の右 = 反時計回りのリングの外側、+1: 左、0: 磁石の上) |
| `spec` / `Spec` | なし | 磁石の寸法 `{要素名または型名: {Width, PoleWidth, Type, PhysicalLength, Height (3D)、Gap・PoleHeight・LegWidth・Aperture (3D の split)}}`。要素名が型名 (`"BEND" "QUAD" "SEXT" "OCT" "SOL" "CAVI" "MULT" "DECA" "DODECA" "TCAV"`) より優先。Type は `Htype` (左右対称)・`Ctype1`・`Ctype2` |
| `magnet_spec` / `UseMagnetSpec` | 台本の定義を読む | 台本の `MagnetSpec[MagnetType[要素名]]` を寸法に使う (§2b)。Python は `None` (既定。台本が定義していれば自動で読む)・`False` (読まない)・dict `{キー: {Width.. }}` (直接与える。キーは `magnet_type(要素名)`、なければ要素名)。SADScript は `UseMagnetSpec->True` (既定)・`False`。優先順は量ごとに **`Spec` の要素名 > `Spec` の型名 > `MagnetSpec` > 既定値** |
| `overlap` / `Overlap` | `"keep"` | 要素名の重なり (§2c)。`"keep"` は全部出す (従来の出力と全行同じ)、`"remove"` は外接矩形が重なる名前を偏向磁石 > 四極 > 6 極 > 8 極 > その他 > MARK の優先で間引く |
| `mode` / `Mode` | `"plan"` (SADScript は `"2D"`) | `"plan"` (`"2d"`) は平面図 (出力は 3D のオプションを足す前と 1 行も変わらない)。`"3d"` (SADScript は `"3D"` か `"3d"`) は磁石の立体 (§2d)。Python は大文字小文字を区別しない |
| `parts` / `Parts` | `"body"` | `mode="3d"` のときだけ。`"body"` は磁石 1 台 = 6 面体 1 個、`"split"` は偏向磁石 (ANGLE≠0 の MULT も) と四極をヨークとポールに分ける (§2d) |

## 2. 何を描くか

座標は、各磁石を入口の幾何 (`LINE["GEO", i]` = 位置 (gx, gy) と chi1) の局所座標 (u, v) = (進行方向, 左) で作り、実験室の (X, Y) に写す。ビームの x は −v なので ANGLE>0 は左に曲がる (chi1 が増える)。幾何は pysad の `ensure_geo` (= `LINE["GEO"]`) の値そのもので、軌道の円弧は `LINE["GEO", i+f]` (要素内の分数位置) の折れ線にする。

| 要素 (`LINE["TYPE"]`) | 描くもの | レイヤー (色) |
|---|---|---|
| 基準軌道 | 開いた POLYLINE 1 本。直線部は曲がり・COORD の前後の頂点だけ、ANGLE≠0 の要素は `ArcStep` 以下の折れ線 | TRAJECTORY (7) |
| BEND (2)、ANGLE≠0 の MULT (22) | 4 角形 (閉じた POLYLINE) と要素名 (TEXT) | BEND (1)、MULT (30) |
| QUAD (4)・SEXT (6)・OCT (8)・DECA (10)・DODECA (12)・MULT (22)・CAVI (31)・TCAV (32) | 長さ × 幅の矩形と要素名 | QUAD (5)・SEXT (3)・OCT (6)・MULT (30)・CAVI (2) |
| SOL (20) | 長さ 0 の境界 (BOUND=1)。BZ≠0 の SOL から次の SOL までを矩形 (名前は始点の SOL) | SOL (4) |
| MARK (41) | 軌道に直角な短い LINE (名前は `MarkText`) | MARK (8) |
| DRIFT ほか | 軌道の線だけ (長さ 0 の MONI・COORD なども) | — |
| 要素名 | 磁石の中心から `text_side` の側に Width/2 + 文字の高さだけ離す。向きは弦の方向 (読める向き)。重なりは `overlap` で消せる (§2c) | TEXT (7) |

偏向磁石の外形は `drawmagvw3d.n` の `DrawBend` の式に従う。連続する同じ名前 (と型) の要素を 1 台の磁石にまとめ (SPLIT した磁石や間に長さ 0 の補正磁石 `HD` が入る磁石)、入口 `is`・出口 `if+1` の GEO から弦 lc0 と偏向角 φ (chi1 の差) を求める。
物理長 la (`PhysicalLength`、なければ弧の長さ φ/(2 sin(φ/2)) lc0)、弦 lc = 2 sin(φ/2)/φ · la、半サジッタ hs = (la/φ) sin²(φ/4)。左右の半幅 wl・wr は Htype で Width/2 ずつ、Ctype1 は φ≥0 のとき wr = PoleWidth/2・wl = Width − wr (φ<0 では左右逆)、Ctype2 は φ<0 のとき wr = PoleWidth/2 (φ≥0 では左右逆)。
頂点は、入口の左 p1 = r1 (−sin e1, cos e1)・右 p2 = r2 (sin e1, −cos e1)、出口の p3・p4 (出口の向き φ に対して −e2 回した面)。r は辺が弦に平行になる距離 (左は (wl−hs)/cos(e−φ/2)、右は (wr+hs)/cos(e−φ/2))。e1 = E1·ANGLE (DIR=1。DIR=−1 なら E2)、e2 = E2·ANGLE (E1 と E2 は ANGLE に対する比。SAD の定義)。物理長が短いときは軌道の弦の中央に置く。

`drawmagvw3d.n` との違い (意図したもの): (1) 連続する要素を、型だけでなく名前も同じときだけまとめる (違う名前の磁石が隣り合っても別々に描く)。(2) 磁石が終わる要素がすぐ次の別の磁石のときも、その磁石を落とさない (原版は落とす)。(3) 出口の面の左右の半幅を入口と同じ (wl・wr) にして、辺が弦に平行になる (原版の r3・r4 は非対称の磁石で wl・wr が入れ替わっている)。(4) 物理長の磁石の位置を弦の方向の中央にそろえる。(5) 偏向角は chi1 の差 (原版は z 軸の 3 次元の角度の ArcCos)。(6) 立体は `mode='3d'` の押し出した直方体・楔形だけ (§2d。原版の `Xtrd`・`Set3Drot` による 3 次元の回転 (垂直面の chi2・ROTATE) は使わない)。コイル・端板・ベローズは持たない。`MagnetSpec[MagnetType[要素名]]` は読む (Width・PoleWidth・Type・PhysicalLength・Height と、3D の split の Gap など。§2b・§2d)。(7) 基準軌道・MARK・要素名・SOL の領域は原版にないもの。

## 2b. 磁石の寸法の定義 `MagnetSpec[]` を読む

`drawmagvw3d.n` と `draw_cad.sad` は、磁石の寸法を `spec = MagnetSpec[MagnetType[LINE["element", is]]]` で得る (定義は台本側のユーザー関数で、`modules_general.n` などに書く。参照した `drawmagvw3d.n` の版の同じフォルダの `modules_general.n` には定義がなく、書式は呼ぶ側の使い方から読んだ)。`MagnetType[要素名]` が磁石の種類 (キー。文字列でも Symbol でもよい) を返し、`MagnetSpec[キー]` が **規則のリスト** `{Width->.5, Height->.4, PoleWidth->.13, Type->Ctype1, PhysicalLength->.8, ...}` を返す。未定義のとき `MagnetSpec[...]` は未評価のまま残り、`drawmagvw3d.n` は `Head[spec]===MagnetSpec` で空リストに直す。`DrawMagnetsDXF[]`・`write_dxf` も同じ形を読む。

```
MagnetType[n_] := Switch[n, "QFA1H" | "QDA2", "QCELL", "BID5" | "BID4", "BLOW", _, n];
MagnetSpec["BLOW"] = {Width -> .50, Height -> .30, PoleWidth -> .12, Type -> Htype};
MagnetSpec["QCELL"] = {Width -> .32, Height -> .32};
```

| キー (規則の左辺) | 使い方 |
|---|---|
| `Width`・`PoleWidth`・`Type`・`PhysicalLength` | `Spec` と同じ意味 (§1)。`PhysicalLength->None` は指定なし |
| `Height` | 平面図では読んで `info['magnet_spec_used']` に残すだけ。`mode='3d'` では磁石の高さ (§2d) |
| `Gap`・`PoleHeight`・`LegWidth`・`Aperture` | `mode='3d'` の `parts='split'` だけが読む (§2d。`Gap` は原版の `Gap` と同じ名前・既定値、他の 3 つは本実装の名前) |
| `CoilGap`・`CoilWidth`・`CoilHeight`・`CurvedPole`・`CoilExtension`・`FlatCoil`・`DrawCoils`・`EndPlate`・`CoilBendRadius`・`CrossSection` | 読まない (コイル・端板・曲がった磁極・断面は未対応) |

- **キーの決め方**: `MagnetType[]` が定義されていればその戻り、なければ要素名 (`LINE["ELEMENT"]` の `.n` のない名前)。`MagnetSpec[キー]` がリストでなければ `MagnetSpec[Symbol[要素名]]` も試す。どれも無い要素は、規則が空で従来の既定値になる (磁石の型の既定: §5)。規則でない項 (リストの中の他の式) は無視する。
- **優先順** (1 つの量ごと。たとえば Width は `Spec` から、PoleWidth は `MagnetSpec` から、のように混ざる): **`Spec` の要素名 > `Spec` の型名 > `MagnetSpec` > 既定値**。呼び出しで明示した `Spec` を、台本のデータ (`MagnetSpec`) より強くした。注意: `Spec->{"QUAD"->{Width->.2}}` のように型名で幅を与えると、`MagnetSpec` が与える個々の四極の幅を**すべて上書きする** (型名の `Spec` も明示の指定のため)。`MagnetSpec` を効かせたい量には型名の `Spec` を与えないこと。定義を使いたくないときは `UseMagnetSpec->False` / `magnet_spec=False`。
- **Python**: `write_dxf(ffs, ..., magnet_spec=None)` (既定) は、`ffs` を作った台本が `MagnetSpec`・`MagnetType` を定義していれば、pysad のインタプリタで `MagnetSpec[MagnetType["要素名"]]` を評価して読む (要素名ごとに 1 回)。台本が `Exit;` で終わっていても読める。定義がなければ何も変わらない (従来と全行同じ出力)。`magnet_spec={キー: {規則}}` と `magnet_type=関数` で Python から直接与えることもできる。
- 値の型: 数はそのまま、Symbol (`Ctype1` など) は文字列、`None` は None。
- **書式の限界**: 別の書式 (`MagSpec[MagnetWidth->..., MagnetLength->...]` を返す定義ファイル) は読まない。キーの大文字小文字は SAD と同じく区別する。

## 2c. 要素名 (TEXT) の重なりを消す

`overlap="remove"` (SADScript は `Overlap->"remove"`) で、名前の TEXT の外接矩形を文字高と位置から見積もり、重なる名前を間引く。`"keep"` (既定) は全部出す。

- **外接矩形**: 中心は TEXT の位置 (中央寄せ 72=1・73=2 と同じ)、幅 = 文字数 × 文字高 × 0.6、高さ = 文字高、向きは TEXT の回転。0.6 は CAD の一般的な字体の目安で、**実際の字体の幅ではない** (フォントを測らず定数にしたのは、gs と pysad で同じ結果にするため。`OpticsPlot` の `RemoveOverlap` は Tk のフォントに依存して gs と pysad で残る名前が違う)。
- **判定**: 向きのついた 2 つの矩形が面積を持って重なるか (分離軸の判定。接するだけは重ならない)。近傍だけを格子 (一辺 = 最大の外接円の直径) で調べるので、全周でも速い。
- **優先順** (`OpticsPlot` の `RemoveOverlap` と同じ発想): 偏向磁石 (ANGLE≠0 の MULT も) > 四極 > 6 極 > 8 極 > その他の磁石 (DECA・DODECA・MULT・SOL・CAVI・TCAV) > MARK (`mark_text` の名前)。優先の高い名前から順に、すでに残した名前と重ならないものだけを残す。同じ優先なら描く順 (進行方向) が先のものを残す。貪欲法なので、残る名前の数が最大とは限らない。
- 磁石・基準軌道との重なりは調べない (名前どうしだけ)。磁石の図形・MARK の線・基準軌道は `keep` と同じ。消した名前は `write_dxf` の戻り値の `text_removed` (`ntext_candidates`・`ntext_removed`) にある (SADScript は画面に出さない)。
- **既定を `"keep"` にした理由**: 既存の出力 (v2.4.1 の最初の版の DXF 8 本。gs との全行一致の基準) を変えないため。また、消した名前は情報の欠落なので、使う側が選ぶ。
- 例 (KEK-LS の標準セル `ARCA`、`example/dxf_lattice/`): 文字高 0.12 m で 52 → 48 個、0.3 m で 52 → 24 個。全周 (`RING`、TextHeight 0.3 m) は 1048 → 481 個で、gs と pysad の台本の出力が**バイト単位で一致**した。図は `KEKLS_cell_v2_*.png` (matplotlib の描画は固定の文字の大きさなので、DXF の文字高とは見え方が違う)。

## 2d. 磁石の立体 (3D、`mode='3d'`・`Mode->"3D"`)

平面図の代わりに、磁石の立体の外形を `3DFACE` で書く (AutoCAD R12 にある図形。外部ライブラリなし)。**コイル・端板・ベローズ・曲がった磁極・ヨークの細かい輪郭は書かない** (未対応)。

- **外形 (`parts='body'`、既定)**: §2 の 4 角形 (BEND・ANGLE≠0 の MULT・QUAD・SEXT・OCT・DECA・DODECA・MULT・CAVI・TCAV・SOL の領域。2D で閉じた POLYLINE にしたもの) を、基準軌道の面 (z = 0) を中心に z = −Height/2 … +Height/2 に押し出した 6 面体にする。偏向磁石は入口・出口の面が E1・E2 に傾いた楔形 (側面・端面は鉛直、上面・下面は水平)。2D の閉じた POLYLINE の代わりに 3DFACE を書き、基準軌道 (POLYLINE)・MARK (LINE)・要素名 (TEXT) は平面図と同じ (z = 0)。レイヤー・色は 2D と同じ規則 (BEND 1・QUAD 5・SEXT 3・OCT 6・MULT 30・SOL 4・CAVI 2)。磁石の垂直方向の位置は基準軌道の面が中心で、ずれ・傾き (chi2・ROTATE) は §5 のとおり無視。
- **3DFACE の並び**: 1 つの立体 = 6 枚で、上面 (z = +Height/2。頂点は 2D の外形と同じ順: 左入口・右入口・右出口・左出口。上から見て反時計回り)、下面 (逆順)、側面 4 枚 (入口の面・右・出口の面・左の辺の順)。法線は外向き (右手系)。8 個の頂点は一度だけ数値にして使うので、共有する頂点の文字列は同じ。上面・下面の (X, Y) は 2D の外形の頂点と同じ値 (同じ式・同じ丸め)。R12 の 3DFACE には立体のグループがないので、**立体の区切りは「ファイルの順に 6 枚ずつ」** (`check_solids`・`example/dxf_lattice/plot_dxf3d.py` はこの約束で読む。`parts='split'` の部品も 1 部品 = 6 枚)。
- **Height (高さ)**: `Spec` の要素名 > `Spec` の型名 > `MagnetSpec` の `Height` > 既定 (量ごと。§2b)。既定は **BEND (ANGLE≠0 の MULT も) 0.406 m** (`drawmagvw3d.n` の `DrawBend` の既定)、**それ以外は幅 Width と同じ** (正方形の断面。原版の `DrawQuad` は Height の既定が 0.45 m で Width とは独立だが、本実装は Width を変えたとき高さも追随する)。Height ≤ 0 は Python が `ValueError` (SADScript は画面にメッセージを出して既定の高さ)。
- **ヨークとポールに分ける (`parts='split'`)**: 偏向磁石 (ANGLE≠0 の MULT も) と四極だけ。他の磁石 (SEXT など) は 6 面体 1 個のまま。ポールはレイヤー `POLE` (色 9。`split` のときだけレイヤー表に足す)、ヨークは型のレイヤー。**模式図であり、実機のヨーク・ポールの輪郭ではない**。全部品は軸に平行な直方体 (偏向磁石は楔形の帯) で、互いに重ならない。
  - 偏向磁石 (H 型・C 型。左右の半幅 wl・wr は §2 の `Type`・`PoleWidth`・`Width` から。軌道の面を z = 0、zi = Gap/2 + PoleHeight): ヨークの上の板 (横 −wr … wl、z = zi … Height/2) と下の板、側脚 (半幅が PoleWidth/2 より大きい側に、横の幅 LegWidth で z = −zi … zi。H 型は両側、C 型は広い側の 1 本。ANGLE の符号と `Ctype1`・`Ctype2` で側が決まる)、ポール 2 個 (軌道の左右 PoleWidth/2、上は z = Gap/2 … zi、下は −zi … −Gap/2)。部品の順は上板・下板・左脚・右脚・上ポール・下ポール (脚は無い側を除く)。端面は 2D の外形と同じ E1・E2 の傾き、側面は弦に平行 (曲がった偏向磁石でも、帯の 4 角形は 2D と同じ式で作る。全幅の帯が 2D の外形と 1e-12 以内で一致)。
  - 四極 (模式図。実機は 45 度の向きの 4 極): ヨークの 4 本の枠 (上・下: 横 ±Width/2、厚み LegWidth。左・右: 上下の枠の間) と、ボア (直径 Aperture) に向かう 4 つのポール (横幅 PoleWidth、ボアの縁 Aperture/2 から枠の内面まで)。部品の順は上枠・下枠・左枠・右枠・上ポール・下ポール・左ポール・右ポール。
  - 寸法 (`Spec`・`MagnetSpec` のキー。無ければ既定の割合。**既定の割合は本実装が決めた仮の値**):

| キー | 対象 | 既定 | 意味 |
|---|---|---|---|
| `Height` | 全磁石 | BEND 0.406、他は Width | 高さ (z の全幅) |
| `Gap` | 偏向磁石 | 0.04 m (`drawmagvw3d.n` の既定) | 上下のポールの間隔 |
| `PoleHeight` | 偏向磁石 | (Height − Gap)/4 | ポールの厚み (z)。ヨークの上下の板の厚みは (Height − Gap)/2 − PoleHeight |
| `PoleWidth` | 偏向磁石・四極 | 偏向磁石 0.13 m (§2 と同じ)、四極 Aperture | ポールの幅 |
| `LegWidth` | 偏向磁石・四極 | 偏向磁石 (Width − PoleWidth)/4、四極 (min(Width, Height) − Aperture)/4 | 側脚の幅・四極の枠の厚み |
| `Aperture` | 四極 | 0.4 × min(Width, Height) | ボアの直径 |

  - 寸法が入らないとき (Python は `ValueError`、SADScript は画面にメッセージを出してその磁石を `body` で描く): 偏向磁石は Gap ≤ 0、PoleHeight ≤ 0、Gap/2 + PoleHeight ≥ Height/2、PoleWidth が Width より大きい (左右の半幅が PoleWidth/2 より小さい)、LegWidth が大きくて脚がポールに当たる。四極は Aperture/2 ≥ min(Width, Height)/2 − LegWidth、PoleWidth > Aperture (上下のポールと左右のポールが重なる)。
- **戻り値** (`write_dxf`): `info['nsolids']` (立体の数)、`info['solids']` (`(要素名, 部品 'body'|'yoke'|'pole', レイヤー, 高さ)` の並び)、`info['counts']['3DFACE']`。
- **確認用の関数**: `read_dxf` は 3DFACE を `{'type': '3DFACE', 'layer', 'pts': [(x, y, z) × 4]}` で読む。`check_solids(entities)` は 6 枚ずつを 1 立体として、すべての有向辺が 1 回ずつ現れ逆向きの辺も 1 回ずつ現れる (= 各辺が 2 面で共有され向きがそろう)・体積が正 (外向きの法線)・bbox を返す。matplotlib の mplot3d の図は `example/dxf_lattice/plot_dxf3d.py`。
- 使い方は §1。SADScript と Python は同じ出力になる (§4)。2D の既定の出力は 1 行も変わらない (§4)。

## 3. DXF の形式

AutoCAD R12 (`$ACADVER` = AC1009) 相当の ASCII。外部のライブラリは使わない。HEADER (`$ACADVER`・`$INSUNITS`)、TABLES (LTYPE の CONTINUOUS、LAYER 11 本 (`0` を含む。3D の `parts='split'` だけ POLE を足した 12 本)、STYLE の STANDARD)、ENTITIES、EOF。
エンティティは `POLYLINE`+`VERTEX`+`SEQEND` (閉じた形は flag 70=1)、`LINE`、`TEXT` (中央寄せ 72=1・73=2)、3D の立体は `3DFACE` (頂点 4 個。コード 10〜13・20〜23・30〜33。不可視辺の 70 は付けない)。LWPOLYLINE は R13 以降なので使わない。色はレイヤーの色 (BYLAYER)。先頭に `999` のコメント 1 行 (3D は `3D magnet solids of the lattice ...`)。座標は scale 後の値を小数 6 桁に丸める (|値|<5e-7 は 0)。`$INSUNITS` は R2000 の変数だが、R12 の読み手は読み飛ばす。

## 4. 確認したこと (2026-11-22)

- **読み戻し** (自前の最小パーサ `read_dxf`。ezdxf は環境になく、あれば `test_dxf_lattice.py` の追加試験で使う): FODO (24 セル) で POLYLINE 121 (磁石 120 + 軌道 1)・TEXT 120・LINE 1、レイヤー 11 本と色番号、`$INSUNITS`、EOF。レイヤー別の個数は要素数と一致 (BEND 24・QUAD 48・SEXT 48 …)。
- **偏向磁石の外形** (SAD の幾何に依らない独立な確認): 入口・出口の GEO の点が端面の上にある、E1・E2 = 0.5 の矩形で端面が弦に直角、E1・E2 を変えると端面が `E·ANGLE` だけ回る、辺が弦に平行で間隔 = Width、軌道が箱の中心線から入口・出口で hs 左・中央で hs 右、ANGLE<0 で逆、Ctype の左右の半幅 (1e-5 以内)。
- **閉じたリング**: 基準軌道の始点と終点の差は FODO 1.7e-14 m (折れ線の長さは周長の 0.9995 倍以上)、KEK-LS (RING、570.72 m) 8e-8 m。ASTRID2 は 1.9e-5 m で、これは台本の `ANGLE = 0.523599` (12 台の和が 2π より 3e-6 rad 小さい) の丸めで、描画の誤差ではない。
- **gs との一致**: 同じ台本 `DrawMagnetsDXF[]` を gs (SAD 1.4.4.2k64、WSL) と pysad で走らせ、8 本の DXF (FODO 3、SOL 1、ASTRID2 2、KEK-LS 2) が**全行 (文字列として) 一致**した。pysad の Python 関数 `write_dxf` の出力も同じ。再現は `py -3.12 example/dxf_lattice/run_dxf_examples.py`。台本と Python の一致は pysad だけで走る試験 (`test_script_matches_python`、3 通りのオプション) にもある。
- **目視**: DXF を読み戻して matplotlib で描いた図 (`example/dxf_lattice/*.png`)。
- **3D (立体、§2d。2026-11-22)**:
  - 3DFACE の立体の検査 (自前のパーサ `read_dxf`・`check_solids`): FODO の全周 (120 台、720 面)・KEK-LS の標準セル (52 台・312 面)・FODO の 1 セル (5 台・30 面) の `body` で、すべての立体が閉じている (各有向辺が 1 回、逆向きも 1 回 = 各辺が 2 面で共有され向きがそろう)、符号つきの体積が正で **2D の外形の面積 × Height に相対 7e-13 以内**、上面・下面の (X, Y) が **同じ引数の平面図の外形の頂点と完全に一致**。`split` (FODO 1 セル 23 部品、KEK-LS 1 セル 280 部品) も全部品が閉じて体積が正、1 台ごとの部品の体積の和が `body` 以下、bbox が `body` の中。
  - 試験 (`test_dxf_lattice.py` の 3D の 11 件): 直線の偏向磁石・四極の `split` のヨーク・ポールの体積が解析式 (直方体の和) と 1e-6 で一致・部品どうしが重ならない (AABB)、C 型の側脚が 1 本、曲がった偏向磁石の部品が 2D の外形の中、楔形の端面の向きが E1·ANGLE、側面が鉛直、全幅の帯が 2D の外形と 1e-12 以内、Height の優先順 (`Spec` の要素名 > 型名 > `MagnetSpec` > 既定)、入らない寸法が `ValueError`、`check_solids` が面の欠けと向きの反転を見つける。
  - **gs との一致**: 同じ台本 (`DXF_FODO_3d.sad`・`DXF_KEKLS_3d.sad`) を gs (SAD 1.4.4.2k64、WSL。`local/wake/gsrun.sh` の呼び出し) と pysad で走らせ、5 本の 3D の DXF (FODO の 1 セル `body`・`split`・全周、KEK-LS の 1 セル `body`・`split`) が**全行一致** (gs は行末が LF、pysad は CRLF の違いだけ)。pysad の台本と Python 関数 `write_dxf(..., mode='3d')` も全行一致 (`run_dxf_examples_3d.py run`)。試験にも pysad の台本と Python 関数の一致 (`body`・`split`・`MagnetSpec`・`Scale`・`Overlap`) がある。
  - **2D の既定の出力が変わらないこと**: `mode` を足した後の `DrawMagnetsDXF.n` と `write_dxf` で、既存の DXF 8 本と v2 の 5 本 (台本のとおりに作り直したもの) が gs・pysad の台本とも**全行一致**。
  - 目視: mplot3d の図 `example/dxf_lattice/*_3d*.png` (DXF を読み戻して描いた。軸の縮尺は x・y・z 同じ)。
- **`MagnetSpec[]` の読み込みと重なり消し** (2026-11-22、試験は `test_dxf_lattice.py` の追加分):
  - `MagnetSpec`・`MagnetType` を定義した台本で、幅・物理長・PoleWidth・Ctype1 が図形に出る (Ctype1 の左右の半幅は端面の傾きを含む式と 3e-5 以内)。優先順 (`Spec` の要素名 > 型名 > `MagnetSpec` > 既定)、`PhysicalLength->None`、`magnet_spec=False`、dict で与えた場合 (台本から読んだのと同じ DXF)、`MagnetType` なしの文字列・Symbol のキー、定義のない台本の出力が従来と同じ、を確かめた。
  - 重なり消し: 残った名前どうしが重ならない・消した名前は優先が同じか高い残った名前と重なる (独立な多角形の判定で、24 セルの FODO 192 個の名前)。優先順 (偏向磁石 > 四極 > 6 極、同じ優先は描く順、MARK は最後)。矩形の判定は独立な実装とランダム 400 組で一致。
  - **gs との一致**: KEK-LS の標準セルの 5 本 (`MagnetSpec`・`Overlap`・`UseMagnetSpec` を組み合わせた `DXF_KEKLS_cell_v2_*.dxf`) と全周 1 本 (`Overlap->"remove"`) が、gs・pysad の台本・Python 関数 (台本から読む場合と dict で与える場合。全周は Python 関数と pysad の台本) で全行一致。pysad の台本と Python 関数の一致は 4 通りのオプションの試験にもある。再現は `py -3.12 example/dxf_lattice/run_dxf_examples_v2.py` (既存の DXF・PNG は上書きしない)。
  - 既定のまま (`overlap="keep"`、`MagnetSpec` なし) の出力は、既存の `DXF_KEKLS_cell.dxf` と全行一致 (`test_example_v2_files`)。

## 5. 限界・未対応・未確認

- **既定は平面図**。立体 (3D、§2d) は `mode='3d'` を選んだときの磁石の外形だけで、**コイル・端板・ベローズ・曲がった磁極・ヨークとポールの細かい輪郭・真空ダクト・架台は未対応** (`split` の分け方も模式図。既定の割合は仮の値)。垂直面の幾何 (chi2、ROTATE のある偏向磁石や垂直ステアリング) は無視し、平面の chi1 だけを使う (ROTATE=90 度の偏向磁石は偏向角 0 の矩形に見える)。
- **磁石の寸法**: `MagnetSpec[MagnetType[要素名]]` (§2b) か `Spec`/`spec` で与える (3D の Height の既定は BEND 0.406 m・他は Width で、これも実機の寸法ではない)。どちらにも無い量の既定は、BEND の Width 0.63 m・PoleWidth 0.13 m・Htype (`drawmagvw3d.n` と同じ)、QUAD 0.45 m (同じ)、SEXT・OCT・DECA・DODECA・MULT 0.30 m、SOL・CAVI・TCAV 0.40 m (この 3 組は本実装が決めた仮の値で、実機の寸法ではない)。**実機の図にするときは `MagnetSpec` か `Spec` を与えること**。`example/dxf_lattice/MagnetSpec_KEKLS.n` の値は説明用の仮の値。`Height` は平面図では読むだけ (3D では高さ)。`MagnetSpec` の書式は呼ぶ側 (`drawmagvw3d.n`・`draw_cad.sad`) の使い方から読んだもので、`MagnetSpec` を定義している実際の台本では確かめていない (`MagSpec[...]` の別の書式は読まない)。
- **要素名の重なり消し**: 外接矩形は文字数 × 文字高 × 0.6 の見積もりで、CAD の実際の字体では多少違う (消しすぎ・消し残しがありうる)。名前どうしだけを調べ、磁石や軌道との重なりは調べない。貪欲法 (優先順、同じ優先は描く順) で、残る数が最大とは限らない。実行時間: KEK-LS の全周 (3224 要素、名前 1048 個) の台本で、`Overlap->"remove"` は pysad 約 31 秒・gs 約 9 秒 (磁石の走査を 2 回通すため。`"keep"` の pysad は約 12 秒。どちらもラティスの読み込みを含む)、Python 関数は 0.1 秒以下 (GEO を計算した後)。
- **BEND の再現の範囲**: 4 角形 (辺は直線で弦に平行、端面は E1・E2) だけ。コイル (CoilExtension・CoilWidth・CoilHeight・DrawCoils)・端板 (EndPlate)・曲がった磁極 (CurvedPole)・端面の曲がりは描かない (磁石の高さ Height は 3D の `mode='3d'` だけ)。E1・E2 の符号は SAD の定義に従うが、CAD で実機の図と重ねて確かめてはいない。SPLIT した偏向磁石の E1 は分割した要素ごとの値と ANGLE の積で読む (原版と同じ)。
- 長さ 0 の磁石 (薄いレンズの MULT など) は描かない (MARK のような線も出さない)。ANGLE≠0 の MULT は E1=E2=0 のセクター形、他の型の ANGLE は無視。
- 要素名は `LINE["ELEMENT"]` (`.n` の番号を付けない名前)。逆向きの要素 (`-ARCAH`) は GEO と DIR に従う。
- CAD ソフト (AutoCAD・LibreCAD・Vectorworks・DraftSight など) で実際に開いてはいない (R12 の形式の記述と自前のパーサ、matplotlib の描き直しだけ。ezdxf も環境になく未実施)。TEXT の中央寄せ・`$INSUNITS` の扱いは CAD によって違うかもしれない。**3D は 3DFACE を 6 枚ずつ並べただけで、CAD 側では 1 つの立体 (3D ソリッド) にならない** (面の集まり。立体のグループ・隠線処理・陰影・法線の向きの扱いは CAD で未確認)。`3DFACE` は R12 にあるが、CAD が不可視辺・面の向きをどう扱うかは確かめていない。
- KEK-LS の全周の DXF は 0.7 MB (要素名つき) なので追跡対象にしない (台本で作る)。他は 9〜90 KB。3D は FODO の全周で 0.19 MB (120 台、720 面)、KEK-LS の標準セルの `split` で 0.37 MB (280 部品、1680 面) と、2D より大きい (3DFACE は 1 枚が 28 行・約 0.25 KB。8 頂点を共有できる polyface mesh にすれば小さくなるが、本実装は 3DFACE だけ)。
- 3D の `split` の部品の数と割合 (`Gap` 0.04 m など) は仮で、`MagnetSpec` に `Gap`・`PoleHeight`・`LegWidth`・`Aperture` を書いた実際の定義ファイルでは確かめていない (キー名は `Gap` が原版と同じで、他は本実装の名前)。偏向磁石の `split` はポールを軌道の左右対称 (±PoleWidth/2) に置く。SADScript は入らない寸法のとき画面にメッセージを出してその磁石を `body` で描くが、Python 関数は `ValueError` (出力の一致の試験は入る寸法だけ)。
- 実行時間: KEK-LS (RING、3224 要素) の台本の `DrawMagnetsDXF[]` は、セルと全周の 2 本で pysad が約 10 秒、gs は 1 秒未満 (Python 関数 `write_dxf` は数秒)。
