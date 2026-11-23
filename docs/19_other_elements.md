# その他の要素 (COORD・TCAVI・UND・曲げつき MULT・MAP)

この文書で分かること: 専用の文書を持たない 5 種類の要素 (座標変換 COORD、横偏向空洞 TCAVI、アンジュレータ UND、曲げつき多極 MULT (ANGLE≠0)、外部写像 MAP) が、pysad のどの計算 (追跡・光学・エミッタンス・幾何) でどう扱われ、何ができないか。
前提: SAD の要素定義 (`要素名 = (キーワード=値 …);`、[02_lattice_main.md](02_lattice_main.md)) と、追跡・光学の基本 ([07_calc_optics.md](07_calc_optics.md)、[10_tracking.md](10_tracking.md)、[09_emittance.md](09_emittance.md))。
読む順: 下の対応表で要素を探し、該当の章を読む。各章は「定義 → 使われる所 → 検証 → 限界」の順。

用語: **K64** は SAD のソース (Fortran、K. Oide ら)。**gs** は K64 から作った SAD の実行体。**v2** は pysad の既定動作、**k64** は `PYSAD_COMPAT=k64` で K64 の動作に合わせるモード ([28_v2_vs_k64.md](28_v2_vs_k64.md))。**Emittance[]** は放射減衰・量子励起を含む 6 次元のエミッタンス計算 ([09_emittance.md](09_emittance.md))。**RFSW** は RF の ON/OFF。**COD** は閉軌道。

## 0. 対応表

次の表は、各要素の K64 での元のルーチンと、pysad の実装ファイルを示す。

| 要素 | 内容 | K64 のルーチン | pysad の実装 |
| --- | --- | --- | --- |
| COORD | 座標の平行移動と回転 | `tcoord.f` (追跡)、`tcoorde.f`・`qcoord.f` (6 次元の写像と行列) | [`pysad/coord.py`](../pysad/coord.py) |
| TCAVI | 横偏向空洞 (クラブ空洞) | `ttcav.f` (追跡)、`ttcave.f` (6 次元の写像と行列) | [`pysad/tcav.py`](../pysad/tcav.py) |
| UND | アンジュレータ | `undulator.f` (K. Ohmi と E. Forest、KEK Report 92-14) | [`pysad/und.py`](../pysad/und.py) |
| MULT (ANGLE≠0) | 曲げの上に多極をのせる | `tmulta.f` (追跡 `tmulta`、Twiss・Emittance[] の `tmultae`) | [`pysad/multa.py`](../pysad/multa.py) (追跡)、[`pysad/multae.py`](../pysad/multae.py) (COD 写像と行列) |
| MAP | ユーザーが SADScript で与える写像 | `temap.f` (追跡 `temap`、Emittance[] `temape`、Twiss `qemap`、幾何 `tgmap`) | [`pysad/elements.py`](../pysad/elements.py) の `Map`、`ffscmd._install_map_hooks` |

## 1. COORD (座標変換)

K64 の `tcoord.f` (追跡)、`tcoorde.f` / `qcoord.f` (6 次元の写像と行列)。幾何 (`LINE["GX"]` など) は `geo.py` (`tfgeo1` の icCOORD)。

```
COORD CO = (DX=1e-3 DY=-2e-3 DZ=5e-3 CHI1=0.01 CHI2=-0.02 CHI3=0.03 DIR=0);
```

- DIR=0 (既定): 入口で DX DY DZ だけずらして (位置の変換は 1 を引く向き) から CHI1・CHI2・CHI3 で回す。DIR≠0: 逆変換 (回して、出口で戻す)。
- DX DY DZ は誤差の属性ではなく変換のパラメータ。pysad は `Coord.params` に置き、ずれ・回転の誤差のコードに混ざらないようにしている (BEAMBEAM も同じ)。

### 使われる所

次の表は、COORD がどの計算のどの実装で扱われるかを示す。

| 計算 | 実装 |
| --- | --- |
| 追跡 (TrackParticles) | `track6n` の `K_COORD` (`coord.k_coord`)。`tturn` の icCOORD は `exit` なので、tapert までの要素数 (la) を減らさない |
| Twiss (CALC、INS、CELL) | `lattice._el_cod_step` の Coord (qcoord → tcoorde の COD 写像と解析行列)。COD が設計軌道上でも動く |
| Emittance[]、CALC6D、Touschek | `emit6d._tcoord` (`piece`、z 行のシンプレクティック組み直しなし) |
| 幾何 | `geo._coord` |

### 検証

K64 (gs) と一致する (追跡は相対 0、Twiss と幾何は 1e-11 以下、Emittance[] は 5e-10 以下)。[`pysad/tests/script_tests/test_coord_elem.py`](../pysad/tests/script_tests/test_coord_elem.py) 12 本。詳細は [history/19_other_elements_validation.md](history/19_other_elements_validation.md)。

### 限界

- track6 (低速追跡) は COORD を扱わない (NotImplementedError)。
- COORD を含むリングでの FIT (解析微分の経路) は gs と全桁一致を確認した (2026-11-16: COORD の対を含むリングで NX・NY=0.2、QF K1・QD K1 を FREE)。COORD のパラメータそのものを FIT 変数にする場合は未確認 (数値微分にもならない可能性がある)。
- TrackParticles の `dv` は COORD で更新されない (`tcoord.f` のまま)。

## 2. TCAVI (横偏向空洞・クラブ空洞)

K64 の `ttcav.f` (追跡)、`ttcave.f` (6 次元の写像と行列)。

```
TCAVI TC = (L=0.4 K0=1e-4 FREQ=4.5e8 PHI=0.1 DX=2e-4 DY=-1e-4 ROTATE=0.05);
```

次の表は TCAVI の主なキーワード。

| キーワード | 意味 |
| --- | --- |
| K0 | 偏向の強さ (ak、kick = ak sin(位相)) |
| FREQ / HARM | RF の周波数。HARM≠0 は ω0·HARM/c。HARM の小数部は周回ごとに位相をずらす (追跡) |
| PHI | 位相 (電荷の符号を掛ける) |
| DX DY ROTATE | ずれと回転 (入口と出口の TENT/TEXIT) |
| RAD | 0 (既定) で放射する (追跡は RAD と FLUC の両方が要る) |
| RANKICK RANPHASE | 乱数の揺らぎ (tgauss)。**追跡は未対応** (NotImplementedError) |

RFSW OFF (CALC の既定) では長さ L のドリフト (tturn、ttcave とも)。

### 使われる所

| 計算 | 実装 |
| --- | --- |
| 追跡 | `track6n` の `K_TCAV` (`tcav.apply_tcav`)。TENT/TEXIT の回転は K64 の式の順 (更新した x を使う) のまま |
| Emittance[]、CALC6D、Touschek | `emit6d._ttcave` (tchge → [tsetr0] → 半ドリフト → キック (解析行列) → 半ドリフト [tradke] → tchge) |
| Twiss (CALC) | RFSW OFF はドリフト。**RFSW ON の CALC は未対応** (NotImplementedError。Twiss が momentum を持たないため) |

### 検証

追跡・Emittance[] とも K64 (gs) と一致する (追跡は相対 0。Emittance[] は L=0 でエミッタンス 1e-11)。追跡と 6 次元の写像が同じ写像であることも確認済み (1 粒子で 1e-19)。[`pysad/tests/script_tests/test_tcav.py`](../pysad/tests/script_tests/test_tcav.py) 8 本。詳細は [history/19_other_elements_validation.md](history/19_other_elements_validation.md)。

### K64 の癖・不具合

- **L≠0 の TCAVI を含むリングの Emittance[] は gs が壊れる** ([30_sad_bug_report.md](30_sad_bug_report.md) B36)。K0=0 でも EnergyLossU0 が NaN、縦のチューンが 0、エミッタンスが桁違い (K0=1e-6 で 1.3e-13 など) になる。ドリフトに置き換えると正常。`ttcave` の 2 番目のドリフトが `tdrife0` に放射の旗 `krad` を渡し、長さ 0 の放射 (alr=0) で `tradke` を呼ぶことが疑わしいが、RAD=1 でも NaN が出るので、原因は特定できていない。pysad は gs の NaN を再現せず、物理として自然な値を出す (ドリフトと同じ結果から K0 の効果が加わる)。L=0 は gs と一致する。
- 追跡は RAD と FLUC の両方 (`rad .and. RAD_TCAV==0`) のときだけ放射する。Emittance[] は `RAD_TCAV==0` だけを見る (NORAD でも TCAVI は放射の旗が立つ)。pysad は Emittance[] の `enarad` に従う。

### 限界

- RANKICK・RANPHASE (追跡)。
- RFSW ON の Twiss (CALC、INS)。
- wake との併用、SOL 領域内の TCAVI (track6n の領域外の要素のみ対応) は確認していない。

## 3. UND (アンジュレータ)

K64 の `undulator.f` (K. Ohmi と E. Forest、KEK Report 92-14)。

```
UND U1 = (L=1 KX=0.1 QY=0.2 POLE=10 SLICE=4 PHI=30);              ! 磁場なし (k64 と v2 の両方で書ける)
UND U2 = (L=1 FBX=1.8 POLE=20 SLICE=2000);                         ! 磁場つき (v2 だけ)
```

### 3.1 モードによる違い

**K64 の UND は磁場を設定できない** ([30_sad_bug_report.md](30_sad_bug_report.md) B37): `initb1.f` は `kytbl(kwFBX,icUND)=ky_FBX_UND` と書くが、`kwFBX` は未定義の変数 (=0) で、`FBX` `FBY` `BX` `BY` はどれも "Undefined keyword" になる。そのため F0 = G0 = 0 で、undulator は磁場なしの (正準運動量の近似の) ドリフトを刻みごとに行う。
さらに `undulator.f` は陰的な刻みの右辺に dpds を二重に掛けている ([30_sad_bug_report.md](30_sad_bug_report.md) B78)。磁場があっても収束力が刻み幅に比例して消える (200 刻みと 2000 刻みで py が 10 倍違う)。

次の表は k64 モードと v2 モードの UND の違い (v2 だけが磁場を受け、光学・放射にも効果が入る)。

| | k64 モード (K64 と同じ) | v2 モード (既定) |
| --- | --- | --- |
| `FBX`・`FBY` | NotImplementedError | 受ける (磁場つきの UND) |
| 追跡 | 磁場なしのドリフト | 磁場つきのカーネル (`und.und_seg1`。dpds は 1 回だけ掛ける) |
| Twiss・Emittance[]・CALC6D | 恒等 (長さは幾何にだけ入る) | 追跡の写像の数値微分が入る。**磁場がなくても長さ L のドリフトとして光学に入る** (追跡と揃えた。`SLICE=0` は恒等) |
| 放射 (RAD・FLUC) | なし | 1 周期を 64 個以上に分けて `tradke`。RAD で放射する |

磁場なしの UND の追跡は両モードで同じ (gs と一致)。光学が違うのは「磁場なしで v2 はドリフト、k64 は恒等」の 1 点で、v2 のこの扱いは未決のまま (K64 の「光学では恒等」を仕様と見るなら、磁場つきのときだけにする選択もある)。

| キーワード | 意味 |
| --- | --- |
| L | 長さ |
| FBX | v2: 垂直磁場のピーク [T] (x 方向に蛇行)。軸上の By = (FBX·c/p) cos(kz s) |
| FBY | v2: 水平磁場のピーク [T]。Bx = (FBY·c/p) cos(kz s + PHI) |
| KX、QY | 横方向の波数 (kx = KX、ky = √(kz²+KX²)、qy = QY、qx = √(kz²+QY²)) |
| PHI | 位相 [度] |
| SLICE | 刻み数 (ds = L/SLICE) |
| POLE | 極の数 (kz = 2π POLE / L)。0 だと kz=0 で NaN になる (K64 も同じ) |

### 3.2 使われる所

| 計算 | k64 モード | v2 モード |
| --- | --- | --- |
| 追跡 (`TrackParticles`、MAIN の TRACK、DA。`track6n`) | `K_UND` (`und.k_und`)。入口・出口の TENT/TEXIT は dx=dy=θ=0 で何もしない | 磁場つきのカーネル |
| Twiss・Emittance[]・EMIT・CALC6D・Touschek | 恒等 (K64 の `qtwiss`・`tturne` は `case default` で何もしない) | 追跡のカーネルの数値微分。放射は 1 周期を 64 個以上に分けた `tradke` |
| 幾何 | `geo` が `UND` の長さだけ直進 | 同じ |
| スピン追跡 (track6 の経路。POL + RAD) | NotImplementedError | 刻みごとに [写像 → 放射とスピンの歳差]。縦磁場 Bs (軸外) と、周期の途中で終わる端のポテンシャルの飛び (薄い蹴り) の歳差も入る |

### 3.3 検証

- k64 モードは gs と一致する (追跡は相対 0、INS の Twiss と幾何は 2e-10 以下)。[`pysad/tests/script_tests/test_und.py`](../pysad/tests/script_tests/test_und.py) 4 本。
- v2 モード (gs は磁場を設定できないので比べていない) は、平面アンジュレータの解析式 (軌道、行路長、垂直の収束、放射損失、リングでのチューン変化・U0・減衰率・量子励起) と 1e-3〜1e-7 で一致する ([`test_v2_C_und.py`](../pysad/tests/script_tests/test_v2_C_und.py) 9 本)。スピンは Thomas-BMT の積分と 2e-6 で一致 ([`test_v2_C2_spin.py`](../pysad/tests/script_tests/test_v2_C2_spin.py))。個々の許容誤差は [history/19_other_elements_validation.md](history/19_other_elements_validation.md)。
- 例: [`example/01_lattice_elements/elem_UND.sad`](../example/01_lattice_elements/elem_UND.sad)。

### 3.4 使える範囲と条件 (v2)

- `SLICE` は `POLE` の 20 倍以上 (刻みは 1 次の陰的な蛙跳び。刻みが周期の整数倍だと sin(kz s) = 0 で磁場が見えない)。試験は 3 GeV (スピンは 1 GeV)、1 T・1.8 T、λ = 5 cm、20 周期。
- 実機の大きなリング・挿入光源のあるリングでは確かめていない。

### 3.5 限界

- k64 モード: 磁場がないので、アンジュレータの放射・軌道への効果は計算されない (K64 と同じ)。
- v2 モード: `BX`・`BY` という名前は K64 にも無い (ValueError)。端の磁極の整形は無い (軌道は x' = −(F0/kz) sin(kz s) で、中心が −F0/kz² ずれる)。PHI≠0 の入口では正準運動量をそのままつなぐので、ベクトルポテンシャルの飛びが蹴りになる (K64 の式のまま)。
- track6 (低速追跡) は、k64 モードでは UND を扱わない。v2 では POL + RAD のスピン追跡の経路で扱う。

## 4. MULT (ANGLE≠0、曲げつき多極)

K64 の `tmulta.f` (追跡 `tmulta`、Twiss・Emittance[] の `tmultae`)。

```
MULT M1 = (L=1 ANGLE=0.05 K1=0.1 K2=0.5 E1=0.5 E2=0.5);
```

曲げ磁石 (ANGLE) の上に、K0 (双極の過不足)、K1、K2 … と SK0 … (スキュー) の多極をのせる要素。`ndiv` 個のスライスに割り、
「半スライスの曲げ → 曲率つきの多極キック → 全スライスの曲げ → キック → … → 半スライスの曲げ」の順に進む。
曲げは BEND と同じ `tbend` (K1=0 は tbendcore、K1≠0 は tbendi) をスライスの引数で呼ぶ。キックは `gkn` の表 (`gkninit`) を使う解析式。

### 使われる所

| 計算 | 扱い |
| --- | --- |
| 追跡 (track6n) | `K_MULTA` (`multa.k_multa`)。放射 (RAD)、揺らぎ (FLUC)、TRPT の ndivrad ×10 も BEND と同じ |
| Twiss (CALC) | `lattice._el_cod_step` → `multae.lattice_step`。COD 写像は解析、行列は各段の解析行列の積 |
| Emittance[]・CALC6D・Touschek | `emit6d._tmulte` → `multae.emit` (誤差なしの BEND の `_tbende_core` を段ごとに呼ぶ) |
| 幾何 | BEND と同じ (ANGLE を使う) |

### 検証

追跡 (K1・K0・スキュー・E1 E2・フリンジ・ずれ・放射・大角度など 11 通り) は gs と相対 0、Twiss はチューン 1e-15・β α 1e-13、Emittance[] (12 セルのリング) は ε 1e-11・チューン 1e-13 で一致する。[`pysad/tests/script_tests/test_mult_angle.py`](../pysad/tests/script_tests/test_mult_angle.py) 21 本。詳細は [history/19_other_elements_validation.md](history/19_other_elements_validation.md)。

### 限界・未確認

- 未対応 (NotImplementedError): BZ≠0 (K64 自体が abort)、VOLT (RF)、F1、TAPER、スキュー双極 SK0 (K64 の ndiv 計算が `fact(-1)` を読む。結果が定まらない)、L=0。
- 追跡の光子 (PHOTONS) は未対応 (ffscmd が SADError で拒否)。偏極 (CALPOL・POL) のスピンは、ANGLE=0 の MULT は [16_spin_pol.md](16_spin_pol.md) のとおり gs と一致し、ANGLE≠0 の MULT (`tmultae`) も gs と一致した (2026-11-16: 偏極時間 14 桁、スピンチューン 1e-14。平面リングなので偏極ベクトルは 0,1,0 で、縦成分のある偏極ベクトルは確かめていない。TrackParticles のスピン追跡での ANGLE≠0 の MULT は未確認)。
- Twiss・Emittance[] で FB1/FB2 (tblfre) と CHI2 は未対応 (NotImplementedError)。追跡は対応。
- DROTATE≠0 の Emittance[] は ε で 5e-6 ずれる (チューン 1e-9)。追跡は一致。tchge の DROTATE の行列が数値微分のため (BEND と同じ)。
- K0 + DX DY + ROTATE を同時に与えたリングの Emittance[] は gs と一致しない (δ の閉軌道が 4e-3 と 5e-4 に分かれる)。水平チューンが 1e-3 で整数に近く、閉軌道の探索が不安定なためとみて追っていない。2 つずつの組み合わせはすべて一致する。
- tmultae の tradke は m ループの中にあり 22 回/スライス呼ばれる。結果に影響しない (1 回でも同じ) が、そのまま写してある (`multae.TRADKE_CALLS`)。
- 追跡と tmultae は K64 の書き方どおり違う: tmultae は ndiv に放射の項がなく、キックを (1+δ) で割らない。

## 5. MAP (ExternalMap)

K64 の `temap.f` (追跡 `temap`、Emittance[] `temape`、Twiss `qemap`、幾何 `tgmap`)。

```
MAP P1 = (L=0.5);
ExternalMap["TRACK", LINE["POSITION","P1"], nt_, x_] := (...);
ExternalMap["EMIT",  LINE["POSITION","P1"], cod_]    := {cod1, trans [, dtrans, dbeam]};
ExternalMap["OPTICS", LINE["POSITION","P1"], cod_]   := {cod1, trans};
```

MAP はユーザーが SADScript の関数 `ExternalMap` で写像を与える要素。**定義がなければ恒等で、長さも 0** (`L` は K64 では使われない。`tgmap` の既定は `pos(l+1)=pos(l)`)。
MAP を長さ L のドリフトにすると周長が L だけ増えるので、pysad は長さ 0 にしている。

### 使われる所

| 計算 | 呼び出し | pysad |
| --- | --- | --- |
| TrackParticles | `ExternalMap["TRACK", l, nt, x]`。x は (7, 粒子数): x, px(1+δ), y, py(1+δ), z, δ, 生存フラグ (番号順) | `track6n.K_MAP` → `ffscmd._map_track` |
| CALC・INS (Twiss、閉軌道) | `ExternalMap["OPTICS", l, cod]` → `{cod1, trans}` | `lattice._map_step` (`_el_cod_step`、`Map.transfer_matrix`) |
| Emittance[]・EMIT | `ExternalMap["EMIT", l, cod]` → `{cod1, trans}` か `{cod1, trans, dtrans, dbeam}` | `emit6d._tmap` |

- TRACK: 同じ形のリストが返れば座標を置き換える。返らなければ何もしない。フラグ 0 の粒子は落とす。**落ちている粒子もフラグ≠0 なら生き返る** (K64 の `doinject`。ヘルプには「復活できない」とあるが、コードは生き返らせる。gs で確認)。
- EMIT の 4 つ組: `T ← trans·T`、`D ← trans·D`、`B ← (trans+dtrans) B (trans+dtrans)ᵀ + dbeam (上三角)`、`D += dtrans (T+D)` (temape.f の順)。2 つ組: `T, D, B` に `trans`。
- **軌道 `cod1` はユーザーが返す。`trans` から作らない。** `cod1` と `trans` が食い違うと (例: `cod+kick` に結合のある `trans`) K64 の閉軌道は物理的でない値で止まる。pysad は固定点を求めるので別の値になる。試験は `cod1 = trans.cod + kick` で作った。

### 検証

CALC (OPTICS)、Emittance[] (EMIT、4 つ組も)、TrackParticles (TRPT、座標の書き換え・フラグ 0 での落下・復活を含む) とも gs と一致する (チューン 1e-14、エミッタンス 1e-11、追跡 10 桁)。既定 (ExternalMap なし) は恒等で、S は MAP の L を数えない (gs と同じ)。[`pysad/tests/script_tests/test_map.py`](../pysad/tests/script_tests/test_map.py) 5 本。詳細は [history/19_other_elements_validation.md](history/19_other_elements_validation.md)。

### 限界・未確認

- `ExternalMap["GEO", l, geo, pos]` (幾何) は未対応 (MAP は直進、長さ 0)。`CompiledMap` (Fortran の動的コンパイル) も未対応。
- TRACK は TrackParticles (`TrackRun`) だけ。MAIN の `TRACK` 文・`track6.py` へのフォールバック・DAPERT・Touschek の経路は MAP を恒等として通す (ExternalMap を呼ばない)。
- 最初からフラグ 0 で渡した粒子の復活は未対応 (SADError)。
- ExternalMap が返す形が違うときのメッセージは近い文面にしたが、gs の文面とは同一でない。
- CALC4D では trans の 4×5 だけが効く (ヘルプ) が、pysad は 6×6 をそのまま使う (z の行は 4D の結果に影響しないので差は出ないはず。未確認)。
- Emittance[] の Orbit・Matrix 規則 (`emit_extra`) は `Map` を通す (`_tmap` を入れた) が、gs との照合は基本の規則だけ。
