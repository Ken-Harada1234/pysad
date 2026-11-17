# 44. ソレノイド (SOL) の解説

最終更新: 2026-11-17

SOL 領域の仕組み、置ける要素、Twiss・Emittance・追跡・幾何それぞれの対応、K64 の癖、動く SAD スクリプト、gs との一致の度合い、限界を、初めて読む人向けにまとめる。
詳細な検証表は `docs/29_sol_optics.md` (Twiss・Emittance[]・幾何)、追跡の式は `docs/08_tracking.md` §3.7。

## 1. 何ができるか

| 機能 | 状況 | 実装 |
| --- | --- | --- |
| TrackParticles (RAD OFF/ON、FLUC) | 対応 (DRIFT・QUAD のみ領域内に置ける) | `track6.py`、`track6n_ext.py` |
| CALC / INS の Twiss (CELL も) | 対応 | `solmap.py` |
| Emittance[] / EMIT (放射込み) | 対応 | `solemit.py` |
| 幾何 (S・GX・GZ・GCHI2・GCHI3) | 対応 (領域内に DRIFT・QUAD・SOL) | `sgeo.py` (tsgeo) |
| POL のスピン (Emittance[]・EMIT・追跡) | 対応 (領域内は DRIFT・QUAD) | `docs/36_spin_csr.md` |
| CALC6D、RADCOD つきの CALC | 未対応 (SADError) | `ffscmd._check_sol_optics` |
| Touschek (INTRA) との組み合わせ | 未対応 | |
| 領域内の BEND・MULT・CAVI | 未対応 | §5 |
| BOUND=0 で BZ≠0 の SOL (領域の外のソレノイド) の Twiss | 未対応 (NotImplementedError) | `solmap.annotate` |
| ずれ・傾き (DX DY DPX DPY DZ CHI) のある SOL の追跡・Twiss・Emittance[] | 未対応 (幾何だけ対応) | |

## 2. SOL 領域の仕組み

SAD の SOL は「ソレノイドの境界」を表す長さ 0 の要素である。

- `BOUND=1` の SOL を **対** で置く。最初の SOL から次の BOUND≠0 の SOL までが 1 つの領域。
- 領域の中は、BZ を持つソレノイド内の **正準運動量** (力学的運動量は `px + bz/2·y`、`py − bz/2·x`) で動き、横方向が結合する。
- BZ の単位は T。規格化した `bz = BZ / (p0/c)` が使われる (K64 の `tfbzs`)。BZ は領域の始端の SOL の値が領域に入り、領域内の SOL (BOUND=0、`SM = (L=0 BZ=0.8)` のような) は途中で BZ を変える。終端の SOL は BZ を 0 に戻す (pysad は終端の BZ を 0 として扱う)。
- 領域内の SOL の `DISFRIN=1` はフリンジなし。
- **GEO=1 は対の片方だけ** に付ける。両方に付けると K64 は ` Duplicated GEO of Solenoid` を出して領域を扱わない (pysad も同じ警告)。BOUND の対がないと ` Missing BOUND of Solenoid`。
- L は 0 にする (SOL 自体は長さを持たない)。領域の長さは中の DRIFT・QUAD の長さの和。

### 動く SAD スクリプト (最小例)

`docs/29_sol_optics.md` と `pysad/tests/test_sol_optics.py` の s1 と同じ (gs で参照を取り、チューン 1e-15 で一致)。

```
MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.3) QD = (L=0.3 K1=-0.25);
BEND B = (L=1 ANGLE=0.2);
SOL SA = (L=0 BZ=0.3 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF DR B DR2 QD DR SA DR2 QF DR2 QD DR2 SB DR B DR);
FFS;
USE BL;
CELL; NORAD; NOFLUC;
CALC;
Print["TUN ", Twiss["NX", "$$$"]/(2 Pi), " ", Twiss["NY", "$$$"]/(2 Pi)];
Exit;
```

BZ=0.3 T、1.7 GeV の領域である。ただしこの例は GEO=1 を付けていない (§4 の K64 の癖。Twiss は GEO に影響されないが、S 座標を使う量は影響される)。実用では `SA = (L=0 BZ=0.3 BOUND=1 GEO=1)` とする。

### 追跡と放射 + Emittance の例

`bench_level9.sad` (リング `RINGS`、2.5 GeV、`SB = (BZ=2 BOUND=1 GEO=1)`、`SE = (BZ=2 BOUND=1)`、領域内に `D1 QD`):

```
SOL   SB = (BZ=2 BOUND=1 GEO=1)  SE = (BZ=2 BOUND=1);
LINE CELLS0 = (ORG CA QF DQS SF DS B SB D1 QD SE DQS SD DS B D1);
...
USE RINGS; RFSW; CELL; CALC; emit;
NORAD;  b100 = TrackParticles[beam0, 1, 1, 100];
RAD; NOFLUC; b100 = TrackParticles[beam0, 1, 1, 100];
```

領域内の追跡は `pysad/tests/test_track6n_ext.py` (`test_track6n_vs_track6`: SOL を含む格子で numba 版と Python 版が NORAD・RAD・FLUC とも座標 1e-13 以内、乱数の最終状態まで一致、`test_thin_mult_and_sol_not_ignored`、`test_apert_with_solenoid`)。
Emittance[] は `pysad/tests/test_sol_emit.py` (e1: 12 セルのリング、BZ=2、回転四極 QS、領域内の SOL SM、QFE の DX DY で軌道が領域を斜めに通る。RFSW)。幾何は `test_sol_geo.py`。

## 3. 領域内に置ける要素

| 要素 | 追跡・Twiss・Emittance | 備考 |
| --- | --- | --- |
| DRIFT | 対応 | bz≠0 では `tdrift_solenoid` (回転しながら進む) |
| QUAD (K1≠0) | 対応 | ずれ (DX DY)・F1・FRINGE なしの場合。ROTATE は可 |
| L=0 の要素 (MARK など) | 対応 (恒等) | |
| 領域内の SOL (BOUND=0) | 対応 | BZ の段差 (フリンジ `tsfrin`)。ずれ・回転・F1 なし |
| BEND | **未対応** | 追跡では NotImplementedError、幾何では警告を出して直進扱い |
| MULT | **未対応** | NotImplementedError (SADHelp は置けるとしているが pysad は未対応) |
| CAVI | **未対応** | NotImplementedError |
| BEAMBEAM | 未対応 | K64 の tsol が SIGSEGV で落ちるので比較できない |

領域の外 (BZ=0) の要素は通常どおり。`bench_level9.sad` の冒頭コメントは SADHelp どおり「DRIFT・直線 BEND・QUAD・MULT」と書いているが、実際に動くのは上の表の範囲。

## 4. 各機能の対応 (pysad の構成と K64 の対応)

### 追跡 (`track6.py`、`track6n_ext.py`)

K64 の `tsol.f` の直訳。`track6.py` の `_tsol_region` (Python)、`track6n_ext.py` の `compile_sol_region` と `K_SOL_FRIN`・`K_SOL_DRIFT`・`K_SOL_QUAD` (numba)。
RAD ON の四極は `tsolqur` (スライスを外側、粒子を内側、各スライスの後に全粒子 `tradk`)。スライス数は `ndiv = 1 + int(|L·hypot(K1·L/L, bz)|/eps)`、eps=0.2 (EPS 指定あり: 0.2×EPS)、RAD のときは `ndivrad` の大きい方 (上限 1000)。
numba 版で作れない要素 (領域内の BEAMBEAM など) は `NotImplementedError` で Python 版にフォールバックするが、損失の処理は Python 版にない。

### Twiss (CALC / INS) (`solmap.py`、K64 の `tsole1`)

追跡と同じ numba カーネルを 1 粒子に使い、COD での **数値微分 (中心差分 1e-7)** を行列にする。K64 の `tsolque` の解析行列は写していない。
`CALC` のたびに `solmap.annotate` が領域内の要素に `_solmap` を付ける。領域内の Twiss の座標は領域内の正準運動量での値 (K64 の `tsole1` と同じ)。幾何 (GX など) は従来どおり。

### Emittance[] / EMIT (`solemit.py`)

K64 の `tsole1` の放射の経路。写像は追跡カーネルの数値微分、放射は bz 付きの `emit6d._State.tradke(al, 0, bzh)`。要素ごとの対応は `docs/29_sol_optics.md` の表。
`tradke` の bz 付き (正準運動量 → 力学的運動量、`bzhr0`、`dbeam` の `b=bzh/2` の項) は `perfkern.tradke_db` の `bzh`・`bzhr0` 引数。

### CALC6D

未対応。SOL (BZ≠0 または BOUND≠0) を含む CALC6D と RADCOD つきの CALC は `SADError` で拒否する (黙って SOL を無視して結合 0 の値を出さないため)。

### 幾何 (`sgeo.py`、K64 の `tsgeo.f`)

GEO≠0 の SOL (始端なら前向き、終端なら後ろ向き) の DX DY (位置) と DPX DPY (向き) から軌道を作り、領域内の DRIFT・QUAD・SOL を通して `geo`・`pos` を作る。
軌道が軸上 (DX DY DPX DPY が 0) のときは直進と同じで従来の経路。`LINE["S"|"GX"|"GZ"|"GCHI2"|"GCHI3", k]` が gs と一致する (§6)。
**fit の `GX..` の解析微分 (K64 の `tdgeo`) の領域内の幾何は未翻訳** (領域幾何つきの fit は中心差分に落ちる扱い。ここは未確認)。

## 5. K64 の癖

1. **SOL (BOUND≠0) の長さを周長に入れない**。周回周波数 (omega0)・空洞の位相・MomentumCompaction・偏極時間に効く。
   pysad は `emit6d._emit_circ` で同じにしている (追跡の `circ` もこれを使う)。PLANAR + SOL で周長 45.6 のときと偏極時間が 15 桁一致した (`docs/history/work_unimpl.md`)。
2. **GEO が無いと S 座標が途切れる**。K64 の `tsgeo` は ` Missing GEO of Solenoid` を出し、領域の S 座標がつながらない。MomentumCompaction が約 100 倍、偏極時間が約半分など、周長を使う量が狂う。
   pysad は同じ警告を出して直進扱い (つなぐ) にするので、**この狂いは再現しない**。参照の gs はすべて GEO=1 付きで作った。
3. **`tsolqur` の最後の半ステップ**: RAD ON の bz≠0 の四極では、スライスごとの `bzp = bz/(1+δ_i)` を粒子ごとに設定するが、最後の半ステップは再計算せず、**最後に `tzsetparams` を呼んだ粒子 (最後のスライスの最後の粒子) の bzp を全粒子に使う** (K64 `tsolqu.f:519-531`)。
   結果が粒子の並びに依存し、差は δ_i − δ_last に比例する。pysad の `tsolqu_bz_all` が再現している (`track6n_ext.py`)。RAD OFF の `tsolqu` は粒子ごとの値。bench_level9 の追跡で確認。
4. `dbeam(8)` が `beam(2)−beam(10)+b·beam(4)` と書かれている (力学的運動量の共分散なら `beam(9)` のはず): gs に合わせてそのまま写した (`docs/21_sad_bug_report.md` B39。bzh=0 では効かない。試験の ε は直しても全桁同じで、どちらが正しいかは決まらない)。
5. 領域内の軌道が 0 でも K64 は `itgetqraddiv` の下限 (pmin) で刻み数が増える。pysad も写している。
6. `tsolqum` (MULT 用) の z の dv 項は `tsolqu` と式が違う (SAD の不整合を再現。MULT の追跡の話で領域内ではない)。

## 6. gs との一致の度合い

数値は `docs/29_sol_optics.md` と各試験から (いずれも gs = K64 2026-10-19 ビルドの出力を参照に取ったもの。特記なければ NORAD または RAD ON、領域内の要素は DRIFT・QUAD・SOL)。

| 項目 | 試験 | 条件 | 一致 |
| --- | --- | --- | --- |
| Twiss (チューン、β、分散) | `test_sol_optics.py` 4 本 | BZ=0.3 (1.7 GeV)、BZ=1.2 (3 GeV、回転四極・領域内 SOL)、INS・RAD ON・DP≠0 | チューン 1e-15、β 3e-14、領域内の β α 分散 結合 1e-11 |
| Emittance[] e1 | `test_sol_emit.py` | BZ=2、3 GeV、12 セルのリング、RFSW | ε_x 2e-12、ε_y 4e-10、ε_z 1e-7、チューン 1e-13、軌道 1e-11 (相対) |
| Emittance[] e2 | 同上 | BZ=5、K1=0.3 L=0.5、軌道 2.7 mm | ε_x 4e-9、ε_y 2e-8、ε_z 3e-6、チューン 1e-11 |
| 幾何 | `test_sol_geo.py` (sg2 sg3) | GEO が始端 (前向き) と終端 (後ろ向き)、DX DY DPX DPY あり、3 GeV | 18 要素の S GX GZ GCHI2 GCHI3 が相対 1e-9 (絶対 1e-13) |
| 追跡 NORAD | `bench_level9` (Level 9) | BZ=2、2.5 GeV、領域に D1 + QD | 5e-11 以下 |
| 追跡 RAD (NOFLUC) | 同上 | 同上 | 1.5e-10 以下 (**ただし保存済みの参照に対して。下の限界 1 を参照**) |
| スピン追跡 | `test_spin_track.py` | SOL 領域内の DRIFT (BZ=0.05) | 4e-12 以下 (四極は照合できていない) |

(追跡の NORAD・RAD の 5e-11・1.5e-10 は `docs/08_tracking.md` の Level 9 の記録)

## 7. 限界・未確認

1. **現在の gs と bench_level9 の保存値が SOL + 四極 + RAD で食い違う (どちらが正しいかは断定できないが、現在の gs が誤っている疑いが強い。ソースのどこかは未特定)**。
   `bench_level9.result` (gs、2026-10-14 の実行) の `S R T1 P1` の x は -0.0027867898 で pysad と一致する。現在の gs (2026-10-19 ビルド) は同じ入力で -0.0014539780 を返す。NORAD は両方で一致する。
   pysad は保存済みの参照 (と `tsolqur` の癖の読み) に合わせている。回帰の bench_level9 は保存値との比較で、現在の gs では一致しない。したがって **SOL (BZ≠0) 内の四極 + RAD の追跡は、現在の gs に対しては未確認**。SOL 領域内の四極のスピン追跡も同じ理由で未確認。
   **2026-11-16 の切り分け** (単一粒子の 1 周追跡、2.5 GeV、RAD NOFLUC、emit のあと、SOL の BZ を 0 から 0.01 T に変える):
   (a) 領域内がドリフトだけ (`SB D1 SE`): gs = pysad (BZ=0 も 0.01 も 9 桁以上)。
   (b) 領域内が四極 (`SB QD SE`): BZ=0 は gs = pysad。BZ=0.01 では gs が x=6.33e-4、y=−5.6e-16 (横の結合なし)、pysad は x=4.04e-4、y=1.1e-8 (連続)。
   (c) `SB D1 QD SE`、BZ=0.01: gs の結果が四極なしの `SB D1 SE` と 8 桁で一致する (x=1.194038e-4、px=−5.363477e-4)。**gs は RAD ON かつ BZ≠0 の領域内で四極の収束・発散力を入れていないように見える**。BZ=0.01 と 0.2 で gs の x はほぼ同じ (2.95e-5。BZ=0 は −2.23e-4) で、BZ に連続でない。
   (d) 現在の gs の 2 つのビルド (10-19 の Linux6 と、11-14 に -O0 で作り直した `~/sad_dbg`) は同じ値 (x=2.954816161E-5、10 桁一致) なので最適化の違いではなく、現在の K64 ソースから作った gs の振る舞い。10-15 の古い実行ファイルは symbol lookup error で動かず、保存済みの参照 (10-14) を出した版かは確認できない。
   結論: pysad が壊れているのではなく、現在の gs が SOL 領域内の四極 + RAD を物理的に誤って扱う疑いが強い (pysad は物理と保存済みの参照に合う)。ただし確証ではない。条件は 2.5 GeV、RAD NOFLUC、BZ は 0〜0.2 T の少数点のみ。
2. **領域内の BEND は未対応**。K64 の `tsol.f` は双極キック付きの `tdrift(al, bzs, phiy, phix, …)` を呼ぶので直訳を試したが (`track6._tdrift_kick` を bz≠0 に一般化)、現在の gs は偏向をほぼ無視する (BZ=0.5 T の領域で ANGLE=0.1 でも位置・角度が 1e-5 しか変わらない。pysad の直訳は ANGLE × 距離だけ曲がる)。原因は未確認 (K0・ANGLE が領域内の BEND に入らない K64 側の癖かもしれない)。合わないので取り下げた。MULT・CAVI の領域内は未着手。
3. Emittance[] は SOL では `emit_extra` の規則 (Orbit・Matrix・BeamMatrices) を使わず `run_emit` の値に落ちる。EMT・TUN・ALP・COD は合うが **Stable は 0、HarmonicNumber は NaN、OneTurnTransferMatrix などの規則は出ない**。
   理由は `emit_extra`・`calc6d` が Python 版の `_tquade` を使い、軌道のある QUAD で gs と 1e-5 (ε_y 1e-4) ずれるため (SOL と無関係の既存の差。`emit6d.tturne` は gs と 1e-12)。
4. 行列が数値微分なので、解析行列の gs と 1e-10 程度の差が残る。
5. RADCOD の閉軌道探索中は K64 の `tsolque` が内部の分割の間でパラメータを固定するが、`solemit` は各スライスの δ で計算し直す (差は 1 スライスあたり 1e-9 程度で 2 次。未確認)。
6. QUAD の FRMD の向き (逆向きの LINE の `direlc<0` での入れ替え) は写していない。
7. SOL の DX DY DZ CHI1..3 と F1 は、追跡 (numba・Python)・Twiss・Emittance[] で未対応 (tsconv の CHI を使う `tsrote`)。幾何だけ対応。
8. 領域内の L=0 / K1=0 の四極の幾何は未対応。
9. 領域内の MULT・BEND・CAVI の `tmulte1` と、`tsgeo`/`tdgeo` の領域幾何 (上記のうち tdgeo)、CALC6D、Touschek は未翻訳。
10. 上の一致の数値は、GEO=1 を付けた SOL (GEO 無しは S 座標が途切れる K64 の癖があり、pysad は再現しない。§5)、エネルギー 1.7〜3 GeV (追跡は 2.5 GeV)、BZ は 0.3〜5 T の範囲だけで取った。これ以外 (低エネルギー、BZ が極端に大きい場合) は未確認。

## 8. 関連文書

- `docs/29_sol_optics.md`: Twiss・Emittance[]・幾何の検証の詳細。
- `docs/08_tracking.md` §3.7: `tsol` の追跡の式と実装の対応。
- `docs/36_spin_csr.md`: POL のスピンと SOL。
- `docs/21_sad_bug_report.md`: B39 (`dbeam(8)`)。
- `docs/history/work_unimpl.md`: 未解決の差の記録。
