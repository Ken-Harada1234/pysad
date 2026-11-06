# 作業ログ: MAIN レベルの TRACK コマンド (C-14)

担当ファイル: `pysad/parser.py` (MAIN 命令の読み取り)、`pysad/ffscmd.py` (命令の実行)、`pysad/maintrack.py` (新規、TRACK 本体)、
`pysad/track6n.py` (損失報告の記録だけ追加)。日付: 2026-10-16。作業用ファイルは `track_work/` (捨ててよい)。
SAD 側の確認は WSL の `bin/gs` (K64 版)。

## 0. 結論 (先に)

1. **SAD K64 の MAIN の `EMIT` と `TRACK` は動かない。** `EMIT USE=RING;` は `Use EMIT_TANCE or Emittance[] within FFS.` と言うだけで何も計算せず、
   続く `TRACK USE=RING NPART=.. TURNS=..;` は ヘッダ (フラグ一覧) の表示のあと **SIGSEGV** で落ちる。
   落ちる場所は `tracka_` → `__trackdlib_MOD_tinip` (初期分布の `tmap` が emit の正規モード行列 `temw` の `r` を使うが、MAIN の EMIT が何も計算しないので未設定)。
   `NX=(1,2)` のリストは K64 の MAIN では `list is not allowed` で拒否される (`NX=1` と書けば通るが、そのあと同じ場所で落ちる)。
   FFS 内の `EMIT` は `r` を作るが、FFS から MAIN に戻って TRACK を実行する経路は見つからなかった (`END;` のあとに入力が読まれない、`GetMAIN["file"]` は SADScript として読む)。
   → **MAIN の TRACK の出力そのものを SAD で取ることはできない** (§3)。
2. ユーザースクリプトの TRACK は古典 SAD の書き方 (`ON/OFF` → `EMIT USE=RING;` → `PHSPLOTS=n;` → `TRACK USE=RING NPART= TURNS= NX=(..) NY=(..) NZ=(..);` → `STOP;`)。
   ソース (`ActTra.f` → `track.f` → `tracka.f` と `tinip` `tsmear` `wjfit` `wtunem` `tcftr` `tapert`) を読んで直訳し、`pysad/maintrack.py` に実装した。
   描画 (`PHSPLOTS` の TopDrawer) 以外の数値出力は、フラグ一覧・損失の報告・smear の表・tune (FFT) の表・NX/NY/NZ の確認表示・`TunesByTracking` を含めて全部出る。
3. 照合: K64 で動かないので、**同じ手順 (ActTra/tracka) を SADScript で書き、追跡だけ `TrackParticles` (= tturn) に任せた版** (`bench_level27a.sad`, `bench_level27b.sad`) を SAD で走らせ、
   pysad の MAIN TRACK の全桁データと比べた。RFSW ON/RAD OFF、RAD ON/FLUC OFF、RAD ON/FLUC ON (乱数列も) の 3 ケース (27a) と、粒子が落ちる 2 ケース (27b) で
   smear (全行)・最終座標・TunesByTracking の 8 本×3 面の周波数と振幅・損失した粒子と周が **SAD と一致 (座標 ≤5e-11、tune ≤5e-11、smear/振幅の相対 ≤7e-9 = SAD の印字桁)**。
   (初期分布の規格化振幅の割り当て (NX→x, NY→y, NZ→dp) は SAD ソースの読みで、SAD の MAIN からは確認できない。DA の `trackd`/`tinip1` (照合済み) と同じ式。)

## 1. ユーザースクリプトでの使われ方 (`ユーザーの旧スクリプトのディレクトリ`, `sad`)

`^ *TRACK +USE` の行を持つファイル **725 本** (VSX_Ring 492、PF_before_STRUG 145、DBA 26、Super_SOR 16、sad/PF 14、example 6、PLS 6、sad/3GSR 6 …)、TRACK 文 1727 個。
FFS (`FFS;`) も含むものは 9 本だけで、残りは MAIN だけのスクリプト (ラティスファイルと連結して `gs lattice script`)。

TRACK 文のキーワード (出現数): `USE` 1725、`NPART` 1726、`TURNS` 1726、`NX` 1740、`NY` 1573、`NZ` 1575。**これ以外 (XEMI YEMI SIGS SIGE ELEC POSI PROTON RADI SPAN CENTER AFTER AXIS SYNC) は 1 回も使われていない。**
値はすべて数値定数 (`NPART=21 TURNS=1000` が 662 回、`NPART=20/150/41/25/10/11/200 TURNS=1000` など)。リストは `NX=(1,10,20,…)` (コンマ) と `NX=(1 5 10 …)` (空白) の両方。

TRACK の前後の MAIN 命令の並び (文の先頭語の列、多い順):
`ON OFF ON EMIT TRACK STOP` 182、`ON ON EMIT PHSPLOTS TRACK STOP` 90、`ON ON EMIT TRACK STOP` 75、`ON ON EMIT PHSPLOTS TRACK PHSPLOTS TRACK PHSPLOTS TRACK STOP` 61、
`ON ON EMIT TRACK TRACK STOP` 42、`ON ON EMIT TRACK EMIT TRACK EMIT TRACK STOP` 22 … → 使われる MAIN 命令は **ON/OFF、EMIT、`PHSPLOTS=n;`、TRACK、STOP だけ**。
`PHSPLOTS=` は 801 回 (描画ファイルの装置番号)、`OFF SMEAR` / `ON FOURIE` は example の 3 本だけ。

## 2. SAD の経路とコマンド表 (src/)

```
toplvl.f (icACT) → doACT → ActTra.f   TRACK 文の引数 (initbl.f: USE TURNS NPART XEMI YEMI SIGS SIGE ELEC POSI PROTON RADI SPAN CENTER AFTER AXIS SYNC NX NY NZ)
                            ↓ 粒子 j (0 始まり) の規格化振幅 x=NX[j mod n], y=NY[j mod n], g(dp)=NZ[j mod n]、px=py=z=0。リストが無ければ x=y=1, g=0
                   track.f  フラグ一覧 (RFSW RADCOD RAD FLUC INTRA POL COD DAPER EMIOU K64 FOURI SMEAR) を表示 → tracka
                   tracka.f 初期分布 (tinip, trackd.f) → 周ごとの記録 → tturn で 1 周 → 損失 (tapert) → 終了後 wjfit/tphplt
                   track.f  `Tracking end:  CPU time = …` → ActTra.f が NX= NY= NZ= を表示
```

| 項目 | SAD の動作 (ソース) | pysad (`maintrack.py`) | 状況 |
|---|---|---|---|
| `ON/OFF フラグ…` | doflag.f。MAIN の ON/OFF で RFSW RAD FLUC SMEAR … を切り替える | 書かれた順に FFS のフラグ辞書へ (`FFS.main_exec`) | 実装 |
| `EMIT USE=名前;` | 古典 SAD: emit の結果 (正規モード行列 r、閉軌道 codin、EMITX/Y/Z、SIGE) を作る。**K64: 何もしない** | `USE` + `run_emit()` (報告も出す)。r, codin, EMITX.. はこれを使う | 実装 (古典の意味) |
| `PHSPLOTS=n;` | 描画ファイルの装置番号 (大域変数) | 値を大域変数に入れるだけ | 描画は無視 |
| `TRACK USE=` | ラインを選ぶ。emit と別のラインなら r は古いまま | 名前が違えば `use()` (emit は消える → 次の TRACK で黙って EMIT) | 実装 |
| `NPART=` `TURNS=` | 粒子数・周回数。TURNS=0 は `Use EMIT_TANCE or Emittance[] within FFS.` を出して終わる (K64) | 同じ | 実装 |
| `NX= NY= NZ=` | σ 単位の規格化振幅 (上記)。リストは粒子数より短ければ繰り返し | `initial_distribution` | 実装 |
| 初期分布の大きさ | `emx=√\|EMITX+EMITY\|` (x,px,y,py に共通)、`emz=√\|EMITZ\|` (RFSW ON、0 なら \|SIGE\|)、RFSW OFF は `\|SIGE\|`。物理座標 = r·(振幅×(emx,emx,emx,emx,emz,emz)) + codin + (TDXI,TDYI,TDZI) | 同じ | 実装 |
| 追跡 | tturn (要素外側・粒子内側、15 要素ごと・周末に tapert、RAD/FLUC は大域の乱数列) | `track6n.turn` (TrackParticles と同じ。Level 6〜11 で照合済み) | 実装 |
| 損失 | `P. id lost in n-th turn at l(name), amplitudes: x y z` (`tapert_report_dropped`、`1x,'P. ',i5,' lost in ',i5,a,' turn',' at ',i5,'(',a,'), amplitudes:',3(1x,g13.7)`) | `track6n.k_apert` が (ID, l, x, y, z, 周) を記録して同じ書式で表示 | 実装 |
| smear (SMEAR ON、既定) | n が nsmear=`max(100,TURNS/50)/100*100` の倍数か n=TURNS+1 のとき、`Smear at n th turn:` と粒子ごとの 2Jx 2Jy 2Jz (周平均) とそのゆらぎ、`Avrg` (tsmear.f) | `smear_report` | 実装 |
| tune (SMEAR ON) | wjfit (FOURIE OFF → wtunem): 粒子ごとに ` P.  i:` と x,y,z 面の FFT パワー (重み .5,1,.5 で平滑) の極大 8 本 (放物線補間の周波数と振幅)、落ちた粒子は `,lost:`。大域 `TunesByTracking` にも入る | `wtunem`, `tune_table` | 実装 |
| `NX=` `NY=` `NZ=` の確認表示 | `autofg(x,'S7.4')` を 7 桁ずつ、70 桁を超えたら改行 (ActTra.f) | `echo_lists` | 実装 |
| `Tracking end:  CPU time = …` | `1X,2A,F10.3,A` (タイトルは 20 桁) | 同じ書式 (時間は実測) | 実装 |
| `PHSPLOTS` の TopDrawer 出力 (tphplt, tdplt$) | 描画ファイル | — | 無視 (警告なし) |
| `ON FOURIE` (wjfit の maxord>0: 多項式フィットの Fourier 解析) | 追加の長い出力 | — | **未対応** (使用 1 本、無視) |
| `ResultsOfTracking[n_,x_]` フック (周ごとに 7 本のリスト) | tracka.f の `RESULTSOFTRACKING` | — | 未対応 (example の注釈だけ) |
| `TURNS<0` (trpt)、`NPARA>1`、`QUICK`、`SPECT/PLOT` 系の旧コマンド、ELEC/POSI/PROTON の電荷・質量の細部、wake・spin | — | ELEC (電荷 −1)・PROTON (質量) は渡すだけ。他は未対応 | 使用なし |

### 2.1 SAD のソースから読み取った癖 (実装に入れたもの)

- **tsmear の正規化は常に nsmear**: TURNS=10 のように短いときも n=TURNS+1 で `nsmear`=100 で割る (11 周ぶんの和を 100 で割った小さい値が出る)。n=TURNS+1 が `nsmear` の倍数+1 なら (nsm=1) 表は出ない。表を出すと和は 0 に戻る。
- **周の頭に記録する**: n=1 は追跡前の初期座標、n=TURNS+1 は追跡後。zbuf の長さ nts は TURNS+1 以上の最小の 2 のべき (残りは 0 のまま FFT に入る)。
- **FFT の符号**: `tcftr(.., conj=.false.)` は Σ z_n exp(+2πi nk/N)。z = x + i px は正の tune で exp(−iφ) なので peak は +ν に出る (テスト `test_wtunem_pure_rotation`)。
- **損失の報告の位置 l**: `l_track=l` は tapert の呼び出しの**後**で設定されるので、要素ループの頭の判定で報告される位置は 1 つ前の要素 (周の先頭では前の周の最後 = nlat。名前は `$$$` とした)。周末の判定は nlat。
  1 周目の最初の判定の l_track は未初期化 (0 のはず) で、そこで落ちる粒子は印が付かない癖があるが、再現していない (初期振幅が LOSSAMPL=1 m 超の場合だけ)。
- **表示される amplitudes は limittrack で切った値**: |x|,|y| は LOSSAMPL×3 (= 3.0)、z は 1e10 に切られた値 (`3.000000  3.000000  1.0000000E+10`)。
- **落ちた粒子の px, py の出力値**: TrackParticles の出力は tconvm を生きている粒子 (npa 個) にだけ掛けるので、落ちた粒子の px, py は pysad の値と違う (SAD の方が大きい)。
  TRACK が表示する x,y,z と tune・smear には影響しない (比較は x,y,z,dp だけ)。
- **smear の表示順は粒子 ID 昇順** (落ちた粒子の入れ替えがあってもスロット順にならない)、`Avrg` は生きている粒子の数で割る。

## 3. K64 で TRACK が動かないことの記録

```
ON RFSW RADCOD COD;
EMIT USE=RING;                      → Use EMIT_TANCE or Emittance[] within FFS.      ($ECH end: CPU time = 0.002 sec)
TRACK USE=RING NPART=5 TURNS=10;    →  RFSW RADCOD RAD  FLUC  INTRA POL   COD  DAPER EMIOU K64   FOURI SMEAR
                                          T     T     F     T     F     F     T     F     F     T     F     T
                                       Program received signal SIGSEGV: Segmentation fault - invalid memory reference.
TRACK USE=RING NPART=5 TURNS=10 NX=(1,2) NY=(0) NZ=(1);   →  *** $TRACK *** list is not allowed (NX の ( ) で、NY NZ も同様) の後に同じ SIGSEGV
TRACK USE=RING NPART=1 TURNS=-1;    →  同じ SIGSEGV (trpt でも)
```

- SIGSEGV のバックトレースのアドレスを `nm` (obj/AMD64-Linux5/sad1.exe) で引くと `toplvl_ → master.0.toplvl_ → doact_ → acttra_ → track_ → tracka_ → __trackdlib_MOD_tinip` (+0x1db)。
  `tinip` の `tmap(xa,-1)` が `use temw, only: rn=>r` の行列を使うので、EMIT が r を作っていないと落ちる。
- 実スクリプト (`VSX_Ring/memory17/v18_67` + `phase`) も K64 では同じ結果 (リストの拒否 → ヘッダ → 落ちる)。
- FFS 内で `USE RING; CELL; CALC; EMIT;` すると r ができるが、MAIN の TRACK を呼ぶ手段がない: `END;` で FFS を出ると `In[2]:=` (SADScript の対話) になり次の行が読まれず、
  `GetMAIN["file"]` (FFS 内) はファイルを SADScript として読んで `LALR(1) parser: syntax error`。
- 結論: **MAIN の TRACK の参照出力は K64 では取れない。** 古典の SAD (KEK 本流の MAIN) では動くはずだが手元にない。
  そこで §4 の SADScript 版で「同じ手順を TrackParticles で」照合した。

## 4. 照合 (bench_level27a / bench_level27b)

- `bench_level27a.sad`: ラティスは Level 6 と同じ (CELL0 に空洞、六極あり)。FFS で `emit` → `Emittance[Matrix->True]` から `NormalCoordinates` (ri) と `OrbitAtEntrance` (codin)、グローバル EMITX/EMITY/EMITZ/SIGE。
  ActTra/tracka の手順を SADScript で書く: 振幅 = (NX,0,NY,0,0,NZ)×(emx,emx,emx,emx,emz,emz)、初期座標 = Inverse[ri]·振幅 + codin、周ごとに `TrackParticles[b,1,1,1]` (1 周)、
  周の頭に ri·(x−codin) から zbuf と 2J の和、nsmear の倍数で smear、最後に `Fourier[]` (= tcftr/√N) から wtunem の直訳 (SADScript) でピーク 8 本。
  6 粒子 `NX=(1,3,5) NY=(1,2) NZ=(1,2,0)`、250 周 (nts=256、nsmear=100 → n=100, 200, 251 に表)。(1) RFSW ON/RAD OFF (N1) (2) RAD ON/FLUC OFF (R1) (3) RAD ON/FLUC ON, SeedRandom[17] (F1)。
- `bench_level27a_main.sad`: 同じラティス + MAIN の命令 (`ON RFSW; OFF RAD; EMIT USE=RING; TRACK USE=RING NPART=6 TURNS=250 NX=(1,3,5) …;` を 3 回、2 回目は空白区切りのリスト)。
  `bench_level27a_pysad.py` が `run_sad_string` で実行し、`ffs.track_results` の全桁データ (emx, emz, codin, smear 全行, 最終座標, TunesByTracking, 損失) を SAD の `.result` と比べる。
- `bench_level27b.sad` / `bench_level27b_main.sad`: `NX=(1,60,100,150,200,300) NY=(1) NZ=(1)` で粒子 4,5,6 が 1 周目、3 が 2 周目 (RAD OFF)、3 周目 (RAD+FLUC) に落ちる。
  SAD 側は `TrackParticles` の生存フラグで落ちた周を取り、落ちた粒子を smear・FFT から外す。

結果 (`python bench_level27a_pysad.py a|b`):

| ケース | emx, emz | codin | smear (全行) | 最終座標 | tune (8×3×6) | 振幅 | 損失 (粒子, 周) |
|---|---|---|---|---|---|---|---|
| 27a N1 (NORAD) | 5.7e-9 (SAD の印字 8 桁) | 2.8e-16 | 18 行 6.6e-9 | 4.4e-11 | 4.9e-11 | 4.4e-9 | なし |
| 27a R1 (RAD, NOFLUC) | 同上 | 同上 | 18 行 5.2e-9 | 5.0e-11 | 5.0e-11 | 3.2e-9 | なし |
| 27a F1 (RAD+FLUC) | 同上 | 同上 | 18 行 4.6e-9 | 4.7e-11 | 5.0e-11 | 3.5e-9 | なし |
| 27b N1 | 同上 | 同上 | 6 行 5.5e-9 | 4.1e-10 (落ちた粒子は x,y,z,dp) | 4.9e-11 | 1.9e-9 | (3,2) (4,1) (5,1) (6,1) = SAD と一致 |
| 27b F1 | 同上 | 同上 | 6 行 4.0e-9 | 4.3e-11 | 4.7e-11 | 2.0e-9 | (3,3) (4,1) (5,1) (6,1) = SAD と一致 |

(smear・振幅は相対差、座標・tune・codin は絶対差。差は SAD の `$FORM=15.10` の印字桁。)
`run_regression.py` に 27a, 27b を追加 (各 約 10 秒)。

照合できていないもの: 損失報告の l (要素番号と名前、上記の l_track の読み)、ヘッダ・書式の行 (gfortran の `1p g12.4` などは gfortran 11 と照合した `fmt_g` で、
ActTra の `autofg` は `fmt_real`)、`Tracking end` の時間、MAIN の EMIT の報告。これらは SAD K64 の MAIN が動かないため、ソースの読みだけ。

## 5. 実装の要点

- `parser.parse_main_commands(src)`: MAIN 部 (FFS; より前) の文を `;` で切って (`()[]{}` と文字列の中は切らない、`!` と `(* *)` は注釈)、
  ON/OFF・EMIT・TRACK・`NAME=値` (PHSPLOTS EMITX … TDXI などの大域変数)・STOP 類を順に返す。`track[...]:=…` (SADScript の関数定義) や `emit = Emittance[]` (代入) は命令としない。
  `parse_track_statement`: `KEY=値` を拾い、リストは `( )` `{ }` の中をコンマ・空白区切りで数値に (簡単な式 `2*3` や `1.5GEV` も可)。同じキーワードは後のものが勝つ (ActTra の引数解析と同じ)。
- `ffscmd.run_sad_string`: MAIN に TRACK 文があるスクリプトだけ、MAIN 命令を書かれた順に `FFS.main_exec` で実行する (ON/OFF が EMIT/TRACK の前後で効くため)。それ以外は従来どおり (ON/OFF を先に反映)。
- `maintrack.run_track`: 上の表のとおり。追跡は `track6n.Prog` + `track6n.turn` を 1 周ずつ呼ぶ (DA の `dapert.py` と同じ作り)。乱数は大域の `track6.RNG` を続けて使う (SAD は乱数列が 1 本)。
  `PYSAD_PROGRESS=100` を環境変数に付けると 100 周ごとに生存数と経過時間を表示する (長い追跡で止まったのか遅いのかを区別するため)。
- `track6n.k_apert` に損失報告の記録 (`log`, `nlog`, `lcur`) を追加 (行数 0 の `log` なら何もしない = DA・TrackParticles は従来どおり)。`turn(apert=...)` は 5 要素 (従来) か 7 要素 (+log, nlog) のタプル。
- 速度: ESRF リング (960 要素) の 21 粒子 1000 周 RAD+FLUC で約 21 秒 (EMIT 別)。

## 6. 回帰とテスト

- 作業前に 13a/13b/14a/16a の `run_sad.py` 出力を保存 (`base_track/`)、作業後 (`after_track/`) と `cmp` で全部同一。
- `QUICK=1 run_regression.py`: 結果は §7 に追記。
- pytest: `pysad/tests/test_main_track.py` 8 本 (命令の読み取り、書式、wtunem の符号と補間、確認表示の 70 桁改行、smear の書式、初期分布、小リングの端から端まで)。

## 7. 結果の最終確認

- `QUICK=1 PYTHONIOENCODING=utf-8 python run_regression.py`: **35 ベンチ FAIL 0** (27a, 27b を含む。初回は 27a の比較スクリプトが SAD 側の結果ファイルの旧形式 (生存フラグ列なし) を読めず FAIL → SAD 側の 27a をフラグ付きの版で取り直して解消)。
- pytest (`--ignore` 3 本): **146 passed** (作業前 138 + test_main_track.py 8)。
- 13a/13b/14a/16a の `run_sad.py` 出力: 作業前 (`base_track/`) と作業後 (`after_track/`) で `cmp` 同一。
- tar.gz: `pysad_20261016.tar.gz` を更新 (bench_level27a/27b の sad・main・result・pysad.py を含む)。
- 既存の動作への影響: MAIN に TRACK 文を持つスクリプトだけ MAIN の命令を順に実行する新経路に入る。それ以外 (ラティス + FFS の bench) は従来どおり ON/OFF を先に反映する経路のまま。
