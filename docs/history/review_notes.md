# docs レビューメモ (2026-10-15)

対象: 01_overview, 02_lattice_main, 03_ffs_commands, 04_sadscript, 05_calc_optics, 06_matching, 07_emittance, 08_tracking, 09_dynamic_aperture, 11_sad_quirks (10/12 は執筆中のため未読)。
照合先: `<SAD>/src`, `Packages/`, `pysad/*.py` (コードは変更していない)。

## 1. 修正した箇所

| ファイル | 箇所 | 内容 | 根拠 | 修正 |
|---|---|---|---|---|
| 11_sad_quirks.md | A7 | 「parser.py:349-365 は推測規則のまま、要修正」は古い。現在の parser.py:353-354 は `E1 = E1*angle + AE1`, `E2 = E2*angle + AE2` (tbende.f:448 `pe1=psi1*phi0+apsi1`、tffs.f:1261 と同じ)。逆向きは parser.py:385-386 で E1↔E2 | parser.py, tbende.f, tffs.f | 済 |
| 11_sad_quirks.md | E 節 | 上記 A7 の残課題行を削除 | 同上 | 済 |
| 01_overview.md | §6 未対応表「その他」 | 「parser の E1/E2 規則 (→ 11 A7)」を削除 | 同上 | 済 |
| 01_overview.md | §2 FFS 層の解釈順 | 「tflag → tgetfv → terror → コマンド → tfgetv → SADScript 式」は順序が違う。tffsa.f は tfprint (:245, 式に見える語を SADScript として評価) → UNTIL/REPEAT → tfif → tflag (:293) → tgetfv (:297) → terror (:306) → tffile (:315) → tfgetlineps (:330, 要素名なら 7000 tfgetv) → コマンド群 → 7000 tfgetv (:1049) の順。03 §1 の表と整合させた | tffsa.f | 済 |
| 04_sadscript.md | §1 | 「FFS コマンドに当てはまらなければ SADScript として評価」は SAD の順序と逆 (SAD は tfprint が先)。pysad が逆向きの判定をしている旨に書き換え | tffsa.f:245 | 済 |
| 01_overview.md | §6 到達点 | 「線形光学 ≤1e-10」「追跡 ≤5e-11」は 05 §5 (Level 2/3/5 は ≤1e-7) と 08 §7 (Level 7-9 の一部 ≤6e-10) の表と食い違うので、括弧で注記を追加 | 05 §5, 08 §7 | 済 |
| 03_ffs_commands.md | §2 フラグ表 RFSW | 既定 F (defflg) だけだと 11 C12 / 01 §2 の「FFS に入ると RFSW=ON」と矛盾して見える。track.f:146-148 (`trpt=.false.; rfsw=.true.`)、:104 (`calc6d=.false.`) で強制される旨を注記 | track.f | 済 |
| 03_ffs_commands.md | §5 USE | 「CONVERGENCE=1e-9 も tffsa で設定」は誤り。tffsa.f:141-177 で設定するのは scale/NP=1000/ExponentOfResidual=2 で、CONVERGENCE=1E-9 は track.f:35 の `tfevals` | tffsa.f:160-177, track.f:35 | 済 |
| 05_calc_optics.md, 08_tracking.md | Level 12 の要素数 | 「3564 要素」と 01/07 の「3563 要素」が不一致。bench_level12.result は NLAT=3564 ($$$ 込み)、要素は 3563。「3563 要素、NLAT=3564」に統一 | bench_level12*.result | 済 |
| 08_tracking.md | §5 SeedRandom[] の戻り | 「4 語の dump」と 11 B3 の「5 要素」が食い違って見える。dump は uint32 4 語 (seed, filled, next の下位/上位) + 名前 "SAD" で 5 要素 (tfRandom_.c:177-195, random_plugin_sad.c:50-55)。注記を追加 | tfRandom_.c, random_plugin_sad.c | 済 |
| 11_sad_quirks.md | C5, C13, C14, C15 の ffscmd.py 行番号 | 642-647 → 700- (`fn_ffs`)、998 → 1060 (MCOD/TCOD/MON)、667 → 729 (EnergyLossU0)、1501-1539 → 1563- (`type_cmd`)、`_shortest` は :83。約 60 行ずれていた | ffscmd.py | 済 |
| 05_calc_optics.md | §1 INS | `ffscmd.py:1134` (BXI) → 1196 | ffscmd.py | 済 |

## 2. ソースで確認して正しかった主な記述 (修正なし)

- 乱数: `seed ← 48828125·seed`、Gauss の φ=seed·π/2³¹, ρ=√(−2log(seed/2³²+0.5))、SeedRandom[n] の `|n/2|·2+1`、初期種 17 (random_plugin_sad.c:50-68, 107-140)。
- tffsmatch の定数 (flim, aimp, badc, amedc, alit=0.75, wlmin=0.009, rtol1=1.02, tstol=1e-6, amtol=1e-9, eps=1e-5, eps1=1e-8)、`rl = CONVERGENCE·max(nfcol,1)` (tffsmatch.f:3-6, 98)、res2r (tffs.f:1655)、OffMomentumWeight の `sqrt(max(1,|mfitp|))` (tffscalc.f:339)、tweigh の重み表 (tffsmatch.f:835-)。
- qcod: itmax=10、`conv=min(1e-10, 1e-18(1+(δ/1e-3)²+nlat))` (qtwiss.f:820-837)、tsolvg の SVD 閾値 1e-8 (tsolvg.f:8)。qcell の cosμ>1→INTRES / <−1→HALFRES (qcell.f:214-230)、xyth=0.375 (tffs.f:210)。
- iprior 表 (tfstk.f:5170-5187): 04 §2.1 の値はすべて一致 (単項 − = 10 も)。
- nlist 64 個 (tffs.f:1622-)、nflag=53 (tffs.f:333)、MAXITERATION 既定 40 (tffsa.f:107)、NP=1000 (tffsa.f:168)、DRIFT L の vlim(1)=0.1 (tfinit.f:28)。
- 既定値: LOSSAMPL=1, LOSSDZ=100, DAPWIDTH=7, GCUT=1e300, NPARA=1 (initb1.f:1716-1815) を init.n:219 の NPARA=4 が上書き。ky_PHI_CAVI=4。defflg: RAD/RFSW/RADCOD/TRPT/GAUSS OFF、COD/FLUC/RADPOL/EMIT ON。
- 追跡: npparamin=9, npnlatmin=3000 (tftrack.f:4)、粒子リストの列数 7-9 / 9-11 (tftrack.f:106-116)、xlimit=3·LOSSAMPL (tturn.f:225)、tsolqu の epsdef=0.2 / tsolque の 0.1 (tsolqu.f:14, tsolque.f:381)、tbendi の `ndiv=1+int(|phi0/eps|)` with eps=1e-3·EPS (tbendi.f:40-41)、tbend の ndivmax=1024、tsetfringep の `−|akk·f1²|/24`, `|akk|·f2` (tffs.f:960-981)、ndivrad (tffs.f:641-651)。
- 放射: cave=8/(15√3), cuu=11/27, gmin=−0.999999 (photons.f:381)、tradkn の `dg=−de/(1+2de)`, `ddpx=−½dpx·dg`, `x+=ddpx·al1/3` (photons.f:1191-1196)、crad/cuc/anrad は tphyzp (tfsetparam.f:105-122)。
- tcod: codw0=(1e-6,1e-7,1e-7,1e-8,1e-5,1e-6), conv0=1e-4 (tcod.f:3-5)。tcave の ndiv 式 (tcave.f:47)。qcav の `rfsw=rfsw.and.trpt` (qcav.f:19)。
- DA: dapert.n の n1=Min[201,…], NP=Min[4000, n1·Min[DAPWIDTH,51]·2] (:34-35)、tfdapert の n1=min(200), np0=max(2000, n1, n1·min(51,DAPWIDTH)) (:56, 94)、trackd の n1p0=256, n2p=51, nw=16, minnp=16、emx/emz (trackd.f:3-4, 161-162)。
- emit.n: SaveEMIT->False で EMITX 等を復元 (:80)、EnergyLossU0/RfVoltageVc (:47)、TransferMatrices (:76)。Fit.n の `xl=If[xa, Table[…], 1]` (:49-50)。
- 02_lattice_main: 型コード (tfstk.f:277-286: SOL=20, ST=21, MULT=22, CAVI=31, TCAV=32, MARK=41, MONI=42, APRT=43, LINE=CELL=100)、MOMENTUM 既定 30e9 / PBUNCH 1e10 / CHARGE 0 (initb1.f:1697-1760)、単位表 (initb1.f:1618-)、tffsinitparam の DP→0.01, BX/BY/BZ→1 (tffs.f:2904-)、doelem の同型再定義は rdkwdl 再実行 (部分上書き)・別型はエラー (doelem.f:57-66)、MONI のキーワードに L なし (initb1.f:1576-1582)、E1/E2 の比→rad と ori による入替 (tffs.f:1258-1262)。parser.py 側: `_SAD_ELEMENT_TYPES` に DECA/DODECA なし、`_SAD_UNITS` は DEG/MRAD/MM/CM/…のみ、ANGLE=0 BEND→Drift (`_sad_type`)、`('REPEAT',n,id)`/`('REVERSE',id)`、`define_element` が dict を置換 — いずれも 02 §3 の記述どおり。02 §3.2 の E1/E2 変換は 11 A7 (修正後) と整合。
- pysad: `_num_jacobian` eps=1e-7 (lattice.py:1258)、`NJ = dict(cache=True, error_model='numpy', nogil=True)` と NAN_MODE=1 (track6n.py:36-44)、PageWidth 既定 131 (sadscript.py:939)、Cavity.PHI=0.0 (elements.py:445)、`_raw` (ffscmd.py:357)、`_fuzz`/`_log10` (match6.py:112-119)、`Matched.`/`Residual =` (match6.py:781-788)、INTRES 判定 (match6.py:405-432)、分散の qcell1 直訳 (lattice.py:1805-)、tests/test_pysad.py 49 件。

## 3. 確認できず修正しなかった点・要確認

| ファイル | 箇所 | 内容 |
|---|---|---|
| 03_ffs_commands.md §6 / 08_tracking.md §1 | TrackParticles の POL ON 列順 | 03 は `[,{sz},{psi}]`、08 は「phis, sz が加わる」。tftrack.f 自体が :344 のメッセージでは `{sz},{psis}`、:502 のコメントでは `(…, delta, phis, sz, alive)` と食い違う。POL は未対応なので保留 |
| pysad (コード、文書外) | `ffscmd._FLAGS_DEFAULT['RFSW'] = False` | SAD は FFS 進入時に RFSW=ON に強制する (track.f:148、11 C12)。pysad の FFS 層は既定 OFF のままなので、スクリプトが RFSW/NORFSW を明示しない場合に Emittance[]/TrackParticles の挙動が SAD と異なる可能性がある。ベンチはすべて明示しているため未検出。コードは変更していない |
| pysad (コード、文書外) | `ffscmd._TYPECODE` の MULT=20 / SOL=22 | 02 §3.3 項目 8 の指摘どおり逆 (tfstk.f:280-281 は SOL=20, MULT=22)。`LINE["TYPE",…]` の戻り値が SOL/MULT で SAD と食い違う。コードは変更していない |
| 02_lattice_main.md §1.3 | 「MAIN の `ON RAD;` は FFS フラグにならない (未検証)」 | 文書自身が未検証と明記。toplvl/doflag と tflag の関係は未照合のまま |
| 11_sad_quirks.md D5 | 非コンパイル 122 本 | mk/sad.obj.mk との突き合わせは未実施 |
| 04_sadscript.md §2.2 | tfetok の `=.`/`/.` + 数字の規則 | tfetok.f の該当分岐は未照合 (記述の真偽は不明) |
| 05_calc_optics.md §3 端数位置 | qfraccomp の E1/E2・FRMD の扱い | 未照合 |
| 06_matching.md §1 | tdfun の関数別変換表、nstab の数え方 (INS で 2) | qdcell.f/tffscalc.f の該当箇所は未照合 |
| 07_emittance.md §2-3 | tecalc/tnorm の並べ替え・alphap・bh の式 | 未照合 (式の形は 16b/16e で全桁一致しているので実害は小さい) |
| 09_dynamic_aperture.md §6 | 16a の `bench_level16a_run.result` が純 Python 版のまま | 文書自身が「再実行して更新するとよい」としている。未実行 |

## 4. 用語の揺れ (修正せず、参考)

- 「閉軌道」と「COD」、「オフモメンタム」と「off-momentum」(11 A13)、「放射あり」と「RAD ON」が混在するが、いずれも文脈で判別できる。
- 追跡側の `tsolrot`/`tbshift`/`tbrot` と CALC/emit 側の `tsolrote`/`tchge`/`tbrote` は別ルーチンなので、表記の違いは誤りではない (08 §4 の注記どおり)。
- Level 12 の「3563 要素」(要素数) と「NLAT=3564」($$$ 込み) は上記のとおり統一した。00_RESUME.md:203 にも「3564 要素」があるが、対象外のため未修正。
