# 作業ログ: TrackParticles の z/dp と戻りの先頭要素 (2026-11-05)

最終更新: 2026-11-16

TrackParticles の z/dp (経路長・空洞) が入っていなかった件 (00_RESUME.md 11/04 のまとめ 1 位)。
pysad の追跡 (track6n の k_cav・経路長) 自体は 10/16 から正しく、原因は **MAIN 部の `OFF RFSW` の扱い**と **TrackParticles の引数 (終点) の解釈**だった。

## 原因と修正

| # | 症状 | 原因 (SAD の仕様、WSL の gs で実測) | 修正 |
|---|---|---|---|
| 1 | MAIN に `OFF RAD FLUC COD RFSW RADCOD;` を書いた 6 次元追跡: SAD は 1 周で z=−1.395e-3・dp=.0099960 (経路長と空洞の作用)、pysad は z=0・dp=初期値 | SAD K64 では **MAIN 部の `ON/OFF RFSW` は FFS の RFSW に届かない**。MAIN の ON/OFF は整数フラグ `$RFSW$` 等を立て、FFS 入口 (track.f:80-107) で `$RAD$ $FLUC$ $COD$ $RADCOD$ $INTRA$ …` を FFS フラグに写すが、RFSW だけは `isynch=$RFSW$` を読んだ後 **`rfsw=.true.` に固定して tffs を呼ぶ (track.f:149)**。FFS 入口の表 (`RFSW RADCOD RAD …`) の RFSW 列は `$RFSW$` の値 (既定 F) で、FFS の `FLAG[{"RFSW"}]` は 1。pysad は MAIN の `OFF RFSW` を FFS の RFSW に反映していたので、`tturn` の規則 (RFSW OFF なら 1 周ごとに z=0、空洞は tlim=0 で位相差 0) が働いて z=0 のままだった | `run_sad_string` と `main_exec` で MAIN の RFSW は `ffs.main_flags['RFSW']` (= `$RFSW$`、MAIN の TRACK 用) に入れ、`ffs.flags['RFSW']` は触らない。`maintrack.run_track` は `main_flags['RFSW']` (TURNS=−1 なら真、track.f:160) を使う |
| 2 | `TrackParticles[before, 1]` の戻りの先頭要素が SAD 179 (=nlat)、pysad 1 | tftrack.f: 第 2 引数は**終点 ld** (既定 nlat)、第 3・4 引数が周番号 nt..nend (mt=nend−nt+1)。ls, ld は nlat−1 を法に正規化 (ls∈[1,nlat]、ld∈[2,nlat]; ld=1 は nlat)。ld<=ls なら mt+1。mt>1 なら ls→nlat、(mt−2) 周、1→ld; mt==1 なら ls→ld。戻りの先頭は `merge(ld, ls, normal)` = 正規化した ld。pysad は第 2 引数を周回数と読み、先頭に始点をそのまま返していた | `fn_track` を書き直し: 始点・終点を `_resolve_loc` で解決して正規化、区間 (`els[a-1:b-1]`) ごとに `track_particles_fast` を呼ぶ (区間でも空洞の周波数が全周長で決まるよう `circ=` を追加: track6n/track6)。RFSW OFF の z=0 は区間 (tturn) ごと |

## SAD の実測 (WSL の gs)

- (MAIN の ON/OFF → FFS 入口で `FLAG[{"RFSW","RAD","FLUC","COD","RADCOD","INTRA",…}]`):
  何も書かない → {1,0,1,1,0,0}; `OFF RAD FLUC COD RFSW RADCOD;` → {1,0,0,0,0,0}; `OFF RFSW;` → {1,…}; `ON RAD; ON RFSW; OFF COD;` → {1,1,1,0,0,0};
  `OFF FLUC;` → {1,0,0,1,0,0}; `ON RADCOD; ON INTRA;` → {1,0,1,1,1,1}。**RFSW 以外は MAIN の ON/OFF が届き、RFSW は常に 1。** (既定: FLUC と COD は ON)
- 3564 要素のリング、dp=0.01 の 1 粒子: FFS 内の `NORFSW;` で z=0・dp=.01 (pysad の旧出力と同じ値)、`RFSW;` で z=−1.3952845e-3・dp=.0099959616。
  VOLT=0 なら z=−1.3955e-3 のまま dp=.01 (dp の変化は空洞の位相 ∝ z による)。区間 `{179,..},179` / `{179,..}` / `{179,..},100`、2 周 `,1,1,2` も取った。
- 注意: FFS 内の `OFF RFSW;` は RFSW を**立てる** (10/16 D 項: FFS 内の ON/OFF は無効で後ろの語が普通に読まれる)。最初の確認はこれで RAD/RADCOD まで立ってしまい、放射損失 (dp0 の粒子が 1 周で −1.05e-4) が混ざった。

回帰: `QUICK=1 py -3.12 run_regression.py` FAIL 0 (50 ベンチ)、pytest 156 passed。
MAIN の TRACK (bench_level27a/b_main.sad) は `ON RFSW;` を書いているので影響なし。MAIN に `ON RFSW` を書かない MAIN TRACK は SAD (track.f) と同じく RFSW OFF になる (従来は ON だった)。

## 変更したファイル

- `pysad/ffscmd.py`: `fn_track` (書き直し)、`main_flags` (`__init__`)、`main_exec` と `run_sad_string` の RFSW の扱い。
- `pysad/track6n.py` `track_particles_fast(..., circ=None)`、`pysad/track6.py` `track_particles(..., circ=None)`。
- `pysad/maintrack.py`: RFSW を `main_flags` から。
- docs: 02_lattice_main.md §1.3 (MAIN の ON/OFF が FFS に届く範囲)、03_ffs_commands.md (TrackParticles)、19_sad_bug_report.md B30、00_RESUME.md。


## 2026-11-16 の確認・追記

- 2026-11-16 確認: 変化なし。根拠: z/dp と TrackParticles の終点の解釈に関する変更は無い。TrackParticles は 11/14〜11/16 に RADLIGHT の軌道と POL のスピン追跡の戻り (7・8 列) が加わったが (work_unimpl.md)、ここの記述 (z/dp、終点) には影響しない。
