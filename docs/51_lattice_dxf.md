# 51. ラティスの CAD 図 (DXF) — 平面図の書き出し (v2.4.1)

最終更新: 2026-11-22 (v2.4.1 開発中)

ラティスを上から見た図 (平面図) を、CAD で開けるテキストの DXF (ASCII、AutoCAD R12 相当) に書き出す。K64 の SAD に同じ機能はなく (v2 の新機能。互換モードによらず同じ動作)、
SAD の `3GSR/esrf/drawmagvw3d.n` (K. Oide 作。`DrawMagnets3D[]` が Vectorworks の VectorScript を書く) の幾何の考え方を参考にした自前版である。
Vectorworks が一般的でないので、出力を DXF にし、立体 (3D) は作らない (平面図だけ)。`drawmagvw3d.n` は参照しただけで、コピー・改変していない。

| 入口 | 場所 |
|---|---|
| Python | `pysad.dxf_lattice.write_dxf(ffs, path, region=None, scale=1.0, ...)`、`load_ffs(台本, 'RING')`、`read_dxf(path)` (検証用の最小パーサ) |
| SADScript (gs でも pysad でも走る) | `DrawMagnetsDXF[file, options]` ([`example/dxf_lattice/DrawMagnetsDXF.n`](../example/dxf_lattice/DrawMagnetsDXF.n)。`Get["DrawMagnetsDXF.n"]` で読む) |
| 例 | [`example/dxf_lattice/`](../example/dxf_lattice/README.md) (FODO・SOL・ASTRID2・KEK-LS の 1 セルと全周) |
| 試験 | `pysad/tests/test_dxf_lattice.py` (9 件) |

## 1. 使い方

```python
from pysad.dxf_lattice import load_ffs, write_dxf
ffs = load_ffs('example/00_lattice/FODO_ring.sad', 'RING')       # USE RING まで (CALC は不要)
write_dxf(ffs, 'ring.dxf')                                       # 全周、単位 m
write_dxf(ffs, 'part_mm.dxf', region=('QF.1', 'QF.4'), scale=1000, text_height=0.1,
          spec={'B': {'Width': 0.5, 'PoleWidth': 0.15, 'Type': 'Ctype1'}, 'QUAD': {'Width': 0.36}})
```

```
FFS;
USE RING;
Get["DrawMagnetsDXF.n"];
DrawMagnetsDXF["ring.dxf"];
DrawMagnetsDXF["part_mm.dxf", Region->{"QF.1", "QF.4"}, Scale->1000, TextHeight->.1,
               Spec->{"B"->{Width->.5, PoleWidth->.15, Type->Ctype1}, "QUAD"->{Width->.36}}];
```

| 引数 (Python / SADScript) | 既定 | 意味 |
|---|---|---|
| `region` / `Region->{始点, 終点}` | 全周 | 要素の番号か `NAME.n`。**終点の要素も含む**。複数に合う名前は始点が最初、終点が最後 |
| `scale` / `Scale` | 1 | m からの倍率。1 → `$INSUNITS`=6 (m)、1000 → 4 (mm)、100 → 5 (cm)、他は 0 |
| `text_height` / `TextHeight` | 0.15 m | 要素名の文字の高さ。`names=False` / `Names->False` で要素名なし |
| `mark_size` / `MarkSize` | 0.1 m | MARK の線の長さ。`mark_text` / `MarkText->True` で MARK の名前も出す |
| `arc_step` / `ArcStep` | 0.05 rad | ANGLE≠0 の要素の軌道の折れ線の最大の角度 |
| `text_side` / `TextSide` | -1 | 要素名を置く側 (-1: 進行方向の右 = 反時計回りのリングの外側、+1: 左、0: 磁石の上) |
| `spec` / `Spec` | なし | 磁石の寸法 `{要素名または型名: {Width, PoleWidth, Type, PhysicalLength}}`。要素名が型名 (`"BEND" "QUAD" "SEXT" "OCT" "SOL" "CAVI" "MULT" "DECA" "DODECA" "TCAV"`) より優先。Type は `Htype` (左右対称)・`Ctype1`・`Ctype2` |

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
| 要素名 | 磁石の中心から `text_side` の側に Width/2 + 文字の高さだけ離す。向きは弦の方向 (読める向き) | TEXT (7) |

偏向磁石の外形は `drawmagvw3d.n` の `DrawBend` の式に従う。連続する同じ名前 (と型) の要素を 1 台の磁石にまとめ (SPLIT した磁石や間に長さ 0 の補正磁石 `HD` が入る磁石)、入口 `is`・出口 `if+1` の GEO から弦 lc0 と偏向角 φ (chi1 の差) を求める。
物理長 la (`PhysicalLength`、なければ弧の長さ φ/(2 sin(φ/2)) lc0)、弦 lc = 2 sin(φ/2)/φ · la、半サジッタ hs = (la/φ) sin²(φ/4)。左右の半幅 wl・wr は Htype で Width/2 ずつ、Ctype1 は φ≥0 のとき wr = PoleWidth/2・wl = Width − wr (φ<0 では左右逆)、Ctype2 は φ<0 のとき wr = PoleWidth/2 (φ≥0 では左右逆)。
頂点は、入口の左 p1 = r1 (−sin e1, cos e1)・右 p2 = r2 (sin e1, −cos e1)、出口の p3・p4 (出口の向き φ に対して −e2 回した面)。r は辺が弦に平行になる距離 (左は (wl−hs)/cos(e−φ/2)、右は (wr+hs)/cos(e−φ/2))。e1 = E1·ANGLE (DIR=1。DIR=−1 なら E2)、e2 = E2·ANGLE (E1 と E2 は ANGLE に対する比。SAD の定義)。物理長が短いときは軌道の弦の中央に置く。

`drawmagvw3d.n` との違い (意図したもの): (1) 連続する要素を、型だけでなく名前も同じときだけまとめる (違う名前の磁石が隣り合っても別々に描く)。(2) 磁石が終わる要素がすぐ次の別の磁石のときも、その磁石を落とさない (原版は落とす)。(3) 出口の面の左右の半幅を入口と同じ (wl・wr) にして、辺が弦に平行になる (原版の r3・r4 は非対称の磁石で wl・wr が入れ替わっている)。(4) 物理長の磁石の位置を弦の方向の中央にそろえる。(5) 偏向角は chi1 の差 (原版は z 軸の 3 次元の角度の ArcCos)。(6) 立体・コイル・端板・`MagnetSpec[]` の呼び出し (modules_general.n のユーザー関数) は持たない。(7) 基準軌道・MARK・要素名・SOL の領域は原版にないもの。

## 3. DXF の形式

AutoCAD R12 (`$ACADVER` = AC1009) 相当の ASCII。外部のライブラリは使わない。HEADER (`$ACADVER`・`$INSUNITS`)、TABLES (LTYPE の CONTINUOUS、LAYER 11 本 (`0` を含む)、STYLE の STANDARD)、ENTITIES、EOF。
エンティティは `POLYLINE`+`VERTEX`+`SEQEND` (閉じた形は flag 70=1)、`LINE`、`TEXT` (中央寄せ 72=1・73=2)。LWPOLYLINE は R13 以降なので使わない。色はレイヤーの色 (BYLAYER)。先頭に `999` のコメント 1 行。座標は scale 後の値を小数 6 桁に丸める (|値|<5e-7 は 0)。`$INSUNITS` は R2000 の変数だが、R12 の読み手は読み飛ばす。

## 4. 確認したこと (2026-11-22)

- **読み戻し** (自前の最小パーサ `read_dxf`。ezdxf は環境になく、あれば `test_dxf_lattice.py` の追加試験で使う): FODO (24 セル) で POLYLINE 121 (磁石 120 + 軌道 1)・TEXT 120・LINE 1、レイヤー 11 本と色番号、`$INSUNITS`、EOF。レイヤー別の個数は要素数と一致 (BEND 24・QUAD 48・SEXT 48 …)。
- **偏向磁石の外形** (SAD の幾何に依らない独立な確認): 入口・出口の GEO の点が端面の上にある、E1・E2 = 0.5 の矩形で端面が弦に直角、E1・E2 を変えると端面が `E·ANGLE` だけ回る、辺が弦に平行で間隔 = Width、軌道が箱の中心線から入口・出口で hs 左・中央で hs 右、ANGLE<0 で逆、Ctype の左右の半幅 (1e-5 以内)。
- **閉じたリング**: 基準軌道の始点と終点の差は FODO 1.7e-14 m (折れ線の長さは周長の 0.9995 倍以上)、KEK-LS (RING、570.72 m) 8e-8 m。ASTRID2 は 1.9e-5 m で、これは台本の `ANGLE = 0.523599` (12 台の和が 2π より 3e-6 rad 小さい) の丸めで、描画の誤差ではない。
- **gs との一致**: 同じ台本 `DrawMagnetsDXF[]` を gs (SAD 1.4.4.2k64、WSL) と pysad で走らせ、8 本の DXF (FODO 3、SOL 1、ASTRID2 2、KEK-LS 2) が**全行 (文字列として) 一致**した。pysad の Python 関数 `write_dxf` の出力も同じ。再現は `py -3.12 example/dxf_lattice/run_dxf_examples.py`。台本と Python の一致は pysad だけで走る試験 (`test_script_matches_python`、3 通りのオプション) にもある。
- **目視**: DXF を読み戻して matplotlib で描いた図 (`example/dxf_lattice/*.png`)。

## 5. 限界・未対応・未確認

- **平面図だけ**。立体 (高さ・Height・コイルの上下) は作らない。垂直面の幾何 (chi2、ROTATE のある偏向磁石や垂直ステアリング) は無視し、平面の chi1 だけを使う (ROTATE=90 度の偏向磁石は偏向角 0 の矩形に見える)。
- **磁石の寸法**: 原版の `MagnetSpec[]` (modules_general.n のユーザー関数) は読まない。`Spec`/`spec` で与える。無いときの既定は、BEND の Width 0.63 m・PoleWidth 0.13 m・Htype (`drawmagvw3d.n` と同じ)、QUAD 0.45 m (同じ)、SEXT・OCT・DECA・DODECA・MULT 0.30 m、SOL・CAVI・TCAV 0.40 m (この 3 組は本実装が決めた仮の値で、実機の寸法ではない)。**実機の図にするときは `Spec` を与えること**。
- **BEND の再現の範囲**: 4 角形 (辺は直線で弦に平行、端面は E1・E2) だけ。コイル (CoilExtension・CoilWidth・CoilHeight・DrawCoils)・端板 (EndPlate)・曲がった磁極 (CurvedPole)・磁石の高さ (Height)・端面の曲がりは描かない。E1・E2 の符号は SAD の定義に従うが、CAD で実機の図と重ねて確かめてはいない。SPLIT した偏向磁石の E1 は分割した要素ごとの値と ANGLE の積で読む (原版と同じ)。
- 長さ 0 の磁石 (薄いレンズの MULT など) は描かない (MARK のような線も出さない)。ANGLE≠0 の MULT は E1=E2=0 のセクター形、他の型の ANGLE は無視。
- 要素名は `LINE["ELEMENT"]` (`.n` の番号を付けない名前)。逆向きの要素 (`-ARCAH`) は GEO と DIR に従う。
- CAD ソフト (AutoCAD・LibreCAD・Vectorworks・DraftSight など) で実際に開いてはいない (R12 の形式の記述と自前のパーサ、matplotlib の描き直しだけ)。TEXT の中央寄せ・`$INSUNITS` の扱いは CAD によって違うかもしれない。
- KEK-LS の全周の DXF は 0.7 MB (要素名つき) なので追跡対象にしない (台本で作る)。他は 9〜90 KB。
- 実行時間: KEK-LS (RING、3224 要素) の台本の `DrawMagnetsDXF[]` は、セルと全周の 2 本で pysad が約 10 秒、gs は 1 秒未満 (Python 関数 `write_dxf` は数秒)。
