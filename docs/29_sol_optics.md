# SOL 領域の Twiss (CALC / INS) と Emittance[] / EMIT

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `tsole.f` (`tsole1`、`qsol`)。`pysad/solmap.py`。追跡 (`TrackParticles`) は以前から対応 (docs/08。スピン追跡は docs/36)。Twiss の写像・行列は `solmap.py`、放射込みの 1 周 (Emittance[]・EMIT) は `pysad/solemit.py`。

```
SOL SA = (L=0 BZ=0.3 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
LINE BL = (ST QF DR B DR2 QD DR SA DR2 QF DR2 QD DR2 SB DR B DR);
```

BOUND≠0 の SOL から次の BOUND≠0 の SOL までが領域。領域の中は BZ のあるソレノイド内の正準運動量で動き、横方向が結合する。

## 作り方

領域内の要素の写像は、追跡 (K64 の `tsol.f`) と同じ座標なので、追跡の numba カーネル (`track6n_ext` の `K_SOL_FRIN` `K_SOL_DRIFT` `K_SOL_QUAD`) を
1 粒子に使い、COD での数値微分 (中心差分、1e-7) を行列にする。`tsolque` の解析行列は写していない。
`CALC` のたびに `solmap.annotate` が領域内の要素に `_solmap` を付け、`_el_cod_step` と `transfer_matrix` (オフモーメンタム・線形の経路) が使う。

## Emittance[] / EMIT (`solemit.py`、K64 の `tsole1` の放射の経路)

領域内の各要素を K64 の順に、写像は上と同じ追跡カーネルの数値微分、放射は bz 付きの `emit6d._State.tradke(al, 0, bzh)` で入れる。

| 要素 | K64 | solemit |
| --- | --- | --- |
| DRIFT (bz≠0) | `tsetr0(bzs/2)` → `tdrife0` (`tsoldz`) → `tradke(al, 0, bzs/2)` | `tdrift_solenoid1` の数値微分 → `tradke` |
| QUAD | `tsolrote` → `tsetr0` → `tqfrie` → `tsetr0` → `tsolque` を `itgetqraddiv` 回 → `tqfrie` → `tsolrote` | 同じ順。`tsolque` は内部分割 (EPS×0.1) ごとに [半/全ドリフト + 本体] → `tradke` (最後以外)、最後に半ドリフト → `tradke` |
| SOL の端・領域内の SOL | `tsfrie` (F1=0 なので放射なし) | `tsfrin1` の数値微分 |

- `tradke` の bz 付き: 正準運動量 → 力学的運動量 (`px+bzh·y`、`py−bzh·x`)、前回の基準点の `bzhr0` (`_State.bzhr0`)、`tr2` と `tr1` の補正、
  `dbeam` の `b=bzh/2` の項 (`perfkern.tradke_db` の `bzh`, `bzhr0` 引数。bz=0 のときは従来と同じ)。
- K64 の `dbeam(8)` は `beam(2)-beam(10)+b*beam(4)` と書かれている (`B12 - B44 + b B13`)。力学的運動量の共分散なら `B34` のはずで、誤記に見えるが、gs と合わせるためそのまま写した (bzh=0 では効かない。`beam(9)` に直しても試験の ε は全桁同じで、どちらが正しいかは試験では決まらない。docs/21 B39)。
- 領域内の軌道が 0 でも K64 は `itgetqraddiv` の下限 (pmin) で刻み数が増える。写しているので同じ刻み数になる。

## 検証 (gs との比較)

### Twiss (`pysad/tests/test_sol_optics.py` 4 本)

- 単純な領域 (BZ=0.3): チューン 1e-15、β 3e-14、分散 1e-14。
- 強い場 (BZ=1.2、領域内に回転四極 QS、領域内の SOL SM) の CELL: チューン 4e-16。領域内 8 要素の β、α、分散、結合 (R2 R3) が 1e-11 以内。
- INS・RAD ON・DP≠0: チューン 1e-15、領域内の Twiss が 1e-11 以内。

### Emittance[] (`pysad/tests/test_sol_emit.py` 3 本。12 セルのリングの 1 セルに領域、RFSW)

- e1 (BZ=2、回転四極、領域内の SOL、QFE の DX DY で軌道が領域を斜めに通る): ε_x 2e-12、ε_y 4e-10、ε_z 1e-7、チューン 1e-13、軌道 1e-11 (相対)。
- e2 (BZ=5、長い四極 K1=0.3 L=0.5、大きい軌道 2.7 mm): ε_x 4e-9、ε_y 2e-8、ε_z 3e-6、チューン 1e-11。
- e1z (BZ=0 の領域): 普通の要素になり gs と 1e-10。
- 要素ごとの切り分け (SA SB のみ、SA D3 SB、SA QS SB、回転なし、DISFRIN=1) もすべて一致した。

## 限界・未確認

- **SOL には GEO=1 が要る** (K64 `tsgeo` が "Missing GEO of Solenoid" を出し、なければ S 座標がそこで途切れて周長 (MomentumCompaction など) が狂う)。pysad の `tsgeo` は翻訳済み (`sgeo.py`、末尾の「幾何」)。
  ただし Emittance[] の周長は GEO の有無で変わらない (pysad は S 座標を途切れさせない。gs は GEO なしだと Emittance[] の αp が 100 倍以上になる)。参照の gs は GEO=1 で作った。`tdgeo` (fit の GX.. の領域内の幾何微分) は領域の部分を翻訳していない。
- **K64 は SOL (BOUND≠0) の長さを周長に入れない**。omega0・MomentumCompaction・偏極時間・追跡の周回周波数 (空洞の位相) に効く。gs で確かめた: PLANAR のリングに GEO のある SOL (L=0.3、0.6、BZ=0) を足しても、
  周長 45.6 のときと偏極時間が 15 桁一致する。pysad は `emit6d._emit_circ` で同じにした (Emittance 系と、`fn_track` の `circ`)。確認は偏極時間 (`test_spin.py` の `test_solenoid_ring_matches_gs`) による。Twiss (CALC) だけの量 (チューンなど) への影響は未確認。
- **CALC6D・RADCOD つきの CALC は未対応** (SADError)。Touschek (INTRA) と SOL の組み合わせも未対応。
- Emittance[] は `emit_extra.emittance_rules` (Orbit/Matrix の規則、BeamMatrices) を SOL では使わず、`run_emit` (emittance6d) の値に落ちる
  (`ffscmd.emittance` の except)。EMT・TUN・ALP・COD は上のとおり合うが、**Stable は 0 (不安定扱い)、HarmonicNumber は NaN、OneTurnTransferMatrix などの規則は出ない**。
  理由: 同じ経路の `emit_extra` / `calc6d` の要素ループは Python 版の `_tquade` を使い、軌道のある QUAD で gs と 1e-5 (ε_y は 1e-4) ずれる
  (SOL と無関係。SOL なしの QFE にずれのあるリングで確認)。`emittance6d.tturne` (elkern) は gs と 1e-12 で合う。そのためフックを入れなかった。
- RADCOD の閉軌道探索中は K64 の tsolque が内部の分割の間でパラメータを固定するが、solemit は各スライスの δ で計算し直す (δ の変化は 1 スライスあたり 1e-9 程度で影響は 2 次。確認はしていない)。
- QUAD の FRMD の向き (逆向きの LINE の `direlc<0` で FRMD を入れ替える所) は写していない。
- 領域内に置ける要素は追跡と同じ: DRIFT、QUAD (K1≠0、ずれ・F1・FRINGE なし)、L=0 の要素、SOL (ずれ・回転・F1 なし)。MULT・BEND・CAVI は NotImplementedError。
- **領域内の BEND** は K64 の `tsol.f` どおりに `tdrift(al, bzs, phiy, phix, …)` で直訳 (`track6._tdrift_kick` を bz≠0 に一般化) して試したが、現在の gs は偏向をほぼ無視する
  (BZ=0.5 T の領域で ANGLE=0.1 でも位置・角度が 1e-5 しか変わらない。直訳は ANGLE×距離だけ曲がる) ので合わず、実装を取り下げた。原因は未確認 (tpara で K0・ANGLE が領域内の BEND に入らないなど K64 側の癖かもしれない)。
  領域内の MULT・CAVI は未着手 (docs/history/work_unimpl.md)。
- **追跡の食い違い (SOL BZ≠0 の領域内の四極 + RAD)**: 保存した gs の出力 (`bench_level9.sad`、2026-10-14) は pysad と一致するが、現在の gs (2026-10-19 ビルド) は同じ入力で別の値を返す
  (S R T1 P1 の x が -.0027867898 から -.0014539780)。NORAD は両方で一致する。保存した出力と現在の gs のどちらが正しいかは、2026-11-16 夜の単一粒子 (粒子の並びの癖が効かない) の 1 周追跡 (2.5 GeV、RAD NOFLUC、emit のあと、SOL の BZ を 0 → 0.01 T) で絞り込んだ:
  (1) 領域内がドリフトだけ (`SB D1 SE`): gs = pysad (BZ=0・0.01 とも 9 桁以上で一致)。
  (2) 領域内が四極 (`SB QD SE`): BZ=0 は gs = pysad。BZ=0.01 で gs は x=6.33e-4、y=−5.6e-16 (横の結合なし)、pysad は x=4.04e-4、y=1.1e-8 (BZ=0 から連続)。
  (3) `SB D1 QD SE` BZ=0.01 の gs の結果は、四極を除いた `SB D1 SE` と 8 桁で一致する (x=1.194038e-4、px=−5.363477e-4)。BZ=0.01 と 0.2 でも gs の x はほぼ同じ (2.95e-5。BZ=0 は −2.23e-4) で、BZ に連続でない。
  以上から、**gs は RAD ON かつ BZ≠0 の領域内で四極の収束・発散力を入れていない疑いが強い** (pysad は物理と保存済みの参照に合う)。ソースのどこかは未特定で、確定ではない。
  (4) 現在の gs の 2 つのビルド (10-19 の Linux6 と、11-14 に -O0 で作り直した `~/sad_dbg`) は同じ値 (x=2.954816161E-5、10 桁一致) → 最適化の違いではなく、現在の K64 ソースから作った gs の振る舞い。10-15 の古い実行ファイル (`obj/AMD64-Linux5/sad1.exe`) は symbol lookup error で今は動かず、保存済みの参照 (10-14) を出した版かどうかは確認できない。`bench_level9` は回帰 (`run_regression.py`) のベンチの 1 つで、保存した参照を使う。この食い違いは追跡 (TrackParticles) の話で、上の Emittance[] の試験とは別の経路 (照合した gs のビルドの日付は未確認)。
- BOUND=0 で BZ≠0 の SOL (領域の外のソレノイド) の Twiss は NotImplementedError。
- 行列が数値微分なので、解析行列の gs と 1e-10 程度の差が残る。
- 領域内の要素の Twiss の座標は領域内の正準運動量での値 (tsole1 と同じ)。幾何 (GX など) は従来どおり。

## 幾何 (tsgeo、`pysad/sgeo.py`)

SOL にずれ・傾き (DX DY DPX DPY) があるとき、K64 の幾何 (LINE["S"|"GX"|…]) は領域内の軌道に沿って決まる (`tsgeo.f`)。
GEO≠0 の SOL (始端なら前向き、終端なら後ろ向き) の DX DY (位置) と DPX DPY (向き) で軌道を作り、領域内のドリフト・四極・領域内の SOL を通して `geo`・`pos` を作る
(フリンジは `tsfrin1`、四極は `k_sol_quad` を 1 粒子に使う)。軌道が軸上 (DX DY DPX DPY が 0) のときは直進と同じなので従来の経路。

- 検証 (`pysad/tests/test_sol_geo.py`): GEO が始端 (前向き、回転四極 K1>0、領域内の SOL) と終端 (後ろ向き、K1<0) の 2 通りで、S・GX・GZ・GCHI2・GCHI3 の 18 要素が gs と 1e-9 (絶対 1e-13) で一致。
- `Solenoid` の DPX DPY DX DY DISFRIN DBZ を parser が読むようにした (以前は BZ BOUND GEO だけ)。
- 未対応: 領域内の BEND・MULT・CAVI (NotImplementedError → 警告を出して直進扱い。BEND は上の「限界」のとおり gs と合わず取り下げた)、領域内の L=0 / K1=0 の四極。
- **ずれ・傾きのある SOL の追跡・Twiss・Emittance[] は未対応のまま** (tsconv の CHI1..3 を使う `tsrote`)。幾何だけが gs と合う。
- `GEO` がない SOL 領域は K64 と同じ文面の警告 (` Missing GEO of Solenoid …`) を出して直進扱い (K64 は S 座標が途切れる: 上の限界を参照)。

解説: docs/44_sol_guide.md (SOL)
