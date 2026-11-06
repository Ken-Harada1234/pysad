# 作業ログ: ユーザー実スクリプト (esrf / pf2) の再点検 (2026-10-16)

10/15〜10/16 に入った多数の修正 (解析微分、INS、K0、ROTATE、欠落関数、CHROMATICITY、REPEAT/IF、CalculateOptics、Element 系、tsvd、不安定続行、RADINT、DISP 範囲、FSHIFT、EMIT 報告など。docs/history/00_RESUME.md 参照) の後で、
ユーザーの実スクリプトの対応状況を再点検した。**コードは変更していない** (不具合は記録のみ)。`ユーザーのディレクトリ` (ユーザーのディレクトリ) では何も実行していない。
結果の表は `docs/12_user_scripts_status.md` の §3・§4 (現状に差し替え、旧表は 3.A・4.A として 10/15 時点のまま残した)。ここには実行ログ・判定基準・追加の観察を書く。

## 1. 対象と成果物の場所

- 対象 86 本: esrf 52 本 (`<work-esrf>/` の disp 系、chrom、matching 系、tune_survey 系、calc_*、dx_sv/dy_sv、correct_cod*、amp_dep_ts*、print_*、touscheklife、tH/tH2/tw_frac、go_1。lattice 込みの `amp_dep_ts_q.sad` を含む。`matching1`、`chrom1` は `matching`、`chrom` と同一、`matching1_dbg*` はデバッグ用なので除外)
  + pf2 34 本 (`<work-pf2>/pysad_jobs.txt` の 32 本 = `_r` 31 本 + chrom_r、および disp、matching)。
- 一覧は `<work-recheck>/manifest.txt` (group、name、lattice、script、SAD 参照出力)。
  esrf のラティスは `pysad/tests/kek-ls_v03_12.lat` (= `<work-esrf>/v03_12`)、correct_cod 系だけ `<work-esrf>/v05_61` (v03_12 には PM/HD/VD が無い)。pf2 は `pysad_jobs.txt` のラティス (v30_10、v31_24 ほか)。
- 出力はすべて `<work-recheck>/`:
  - `<group>.<name>.pysad.out` (pysad)、`run_log.txt` (rc・秒・行数)
  - `<group>.<name>.sad.out` と `<group>.<name>_full.sad` (今回新しく取った SAD)、`sad_log.txt`
  - `nodraw/` (draw を無効化して取り直した SAD。disp_b、disp_arc_optics)
  - `run_one.sh`、`sad_missing.sh` (実行に使ったスクリプト)

## 2. 実行ログ

### 2.1 SAD 参照出力 (手順 1)

- 既存の SAD 出力があるもの 59 本はそのまま使用 (esrf: `disp.sad.out`、`chrom1.sad.out`、`calc_*.sad.sad.out`、`calc_chromaticity.sad_r.out`、`amp_dep_ts.sad_r.out`、`amp_dep_ts_q.sad_r.out`、`correct_cod_v05.sad_r.out`、`dx_sv.sad_r.out` など。pf2: `<name>.sad.out`、matching は `matching.sad_r.out`)。
- 無かった esrf の 27 本は `cat <lattice> <script>` を `NPARA = 1` に書き換えて、`<work-esrf>` の中で `wsl -d Ubuntu -- gs` で取得 (300 s 制限、1 本ずつ)。
  `sad_log.txt` の通り 24 本は 1〜39 s で終了、`disp2` は 296 s、`correct_cod_nomon_dq` は 247 s、`disp_da`・`disp_test`・`print_cod_magnitude` は 300 s の制限で打ち切り (SAD も終わらない)。
  gs はパスに空白があると読めないので、連結ファイルを `<work-esrf>/rc1016_<name>_full.sad` に一時コピーして相対パスで渡し、実行後に消した。
  最初の試行は `wsl` が while ループの標準入力を食って 1 本で止まった (`< /dev/null` で解決)。
- FITP>1 の GO で SAD が落ちる (Run time Environment から再起動) 件は、`chrom1.sad.out` に 61 回、`chrom_b.sad.out` に 3 回、`<work-pf2>/chrom_r.sad.out` に 47 回出ている。最終解が取れないので比較は「再起動の途中で出た値との照合」になる。
- `draw` (WriteGraphs 内) は WSL の SAD が core dump する (`disp_b`、`disp_arc_optics`、`plot_inj_r` の SAD 出力が `dumped core`)。disp_b と disp_arc_optics は `modules_general.n` の `FFS["draw …"]` を `Print["nodraw …"]` に置換して取り直した (`<work-recheck>/nodraw/`)。
- 上記以外に、`go_1.sad` の SAD は `System["go …"]` が `sh: go: not found` になる (pysad は実行しない旨を表示)。

### 2.2 pysad の実行 (手順 2)

- コマンド: `PYTHONIOENCODING=utf-8 python -u run_sad.py <lattice> <script> > <work-recheck>/<group>.<name>.pysad.out 2>&1` を `timeout 600` 付きで、`xargs -P 3` で 3 並列。
- 86 本全部が終わるまで約 80 分 (12:20 ごろに全件終了)。打ち切り (rc=124) は esrf 10 本、例外終了 (rc=1) は esrf tune_survey_momentum と pf2 plot_inj_r の 2 本。
- 最長は touscheklife (pf2、542 s)、matching3_01 (559 s)、amp_dep_ts_x_octa (462 s)、chrom (465 s)、amp_dep_ts_x (441 s)。10/15 に 900 s で打ち切られた amp_dep_ts (今回 312 s で完走)、pf2 chrom (今回 164 s で完走) は大きく改善。

### 2.3 比較 (手順 3)

- 既存の `cmp_tags.py`/`cmp_lines.py`/`cmp_blocks.py` はタグ付き (T001、L001 など) のベンチ用なので、ユーザーのスクリプトの出力 (タグ無し) には合わない。
  代わりに、10/15 の `work_pf2_survey.md` と同じ方法の小さな比較スクリプト (scratchpad の `cmpgen.py`。リポジトリには入れていない) を使った:
  SAD の出力から入力エコー (連結した lattice・script・modules_*.n の行と一致する行)、バナー、`Out[n]:=` を除き、数値を `#` に置換した骨格を difflib で整列して、対応行の数値を比較する。絶対差 1e-9 未満は無視。
  `Momentum compact. alpha` の行は表示桁の違い (pysad は小さい値で 6 桁) があるので最後の数だけ除いた。
- 数値行の一致度は目安 (整列は骨格一致ベース)。判断は各 `.out` を直接見て行った (disp_b・chrom_b・matching・correct_cod・disp_R56 などは差分を読んだ)。
- 分類は docs/12 §3 の冒頭に定義がある。例外は 16 本で、うち 4 本 (disp_L、disp_eta、disp_optics、calc_H_2) は USE 先 (ARCAL、TEST2、ARCINJ など) がラティスに無く、SAD も無意味な出力を出すもの (スクリプトがこのラティス用ではない)。

## 3. 結果の要約 (表は docs/12 §3.2・§3.3)

| 対象 | 本数 | 数値一致 | 一部一致 | 動くが不一致 | 例外 | 打ち切り |
|---|---|---|---|---|---|---|
| esrf | 52 | 17 | 13 | 1 | 11 | 10 |
| pf2 | 34 | 14 | 13 | 2 | 5 | 0 |
| 合計 | 86 | 31 | 26 | 3 | 16 | 10 |

### 3.1 10/15 との前後比較

- pf2 の `_r` 31 本 (work_pf2_survey.md): 例外なしで完走 13 → 26、例外 16 → 5、打ち切り 2 → 0、数値まで一致 4 → 12 (一部一致 12、不一致 2)。
  数値一致に加わったもの: calc_chromaticity、amp_dep_ts_dp (10/15 は不安定で SADError 停止)、amp_dep_ts_x (νx 2.6e-3 のずれが解消)、bump2、disp_response、disp_kick_effect (K0 の効き)、matching02 (GO 反復列が一致)、matching_BMNCELL (打ち切り → 6 s)。
  例外 → 完走 (一部一致): da_xy (Repeat/Until)、matching_v30_12・tune_sv_v30_24_1 (2 段目の発散が解消、Residual は SAD と一致)、touscheck・touscheklife (ON INTRA)、kick (K0 構文)、correct_cod_9_16 (長さ 0 四極)。完走・不一致 → 一部一致: calc_bump (打ち切りからも)、bump3、sad2ele (E1 の返値)。
  動くが不一致に残る: error_ring (TypeError は解消したが誤差表が違う)、print_all_error (CalculateOptics は解消したが値が違う)。
  例外のまま: disp_q_K1、beamsize、disp_bt、matching_geo_ring、plot_inj (下の §5)。
- esrf の一括実行 40 本 (`<work-esrf>/batch_log.txt`、10/15 朝、240 s 制限。今回は 600 s 制限なので打ち切りの比較は厳密ではない): rc=0 17 → 29、rc=1 12 → 1、打ち切り 11 → 10 (うち 3 本 chrom・matching3_01・dy_sv は完走になり、print_all_error・print_cod_magnitude は例外から打ち切りになった)。
  数値一致 13、一部一致 8、動くが不一致 1 (disp_R56)、例外 8、打ち切り 10。10/15 に rc=1 だった calc_H、calc_RadInt は数値一致になった。
- docs/12 の代表 11 本 (旧表): 数値一致 4 → 6 (matching、pf2 matching が加わる)、一部一致 4 (chrom、amp_dep_ts、dx_sv、pf2 chrom)、例外 1 (correct_cod)。
  chrom は 15 分 → 465 s、amp_dep_ts は 900 s 打ち切り → 312 s で完走 (i≤181 の tune が一致。INS 初期オフセットの値は 47.5800386898 で SAD と一致)、pf2 chrom は 900 s 打ち切り → 164 s。
  esrf matching は最終 K1 が 14 桁一致 (10/15 は QFA1H 0.478 と 0.400 の別解)。

## 4. 追加で分かったこと (不具合の記録。コードは未変更)

1. **CodCorr の Dot** (esrf correct_cod、correct_cod_nomon、correct_cod_nomon_dq、print_emittance の 4 本): pysad は `Monitor["PM*"]`・`Steer[…]`・`CalculateOptics[…]`・`CalculateCODResponse[1]` (160×200 の応答行列) まで通り、
   `CodCorr[matx,befx,IEX]` の中で `???General::equalleng: Equal lengths of arguments are expected in Dot` で止まる。`Twiss["DX","PM*.*"]` は 200 要素、`LINE["NAME","HD*.*"]` は 160 要素と長さは正しいので、原因は CodCorr 内の別の Dot (Eigensystem の戻りの形か `ivect.laminvmat.ivec`)。
   小さな行列で再現しようとしたが、FFS モード内での定義の都合で未達 (次に調べる点)。SAD 側は `correct_cod` (v05 版) では `MON 0 0 0 0 PM*` が "Missing value for PM*" で停止し、`nomon` 版は 270 s で完走する。ただし SAD は誤差 (DELX など) が COD に効かない (1e-17) ので、pysad (誤差が効く) とは比較できない。
2. **FFS コマンドのエラーの扱い**: SAD は `FFS["go;",6]` が "No variable" などで失敗すると `FFS::error` で外側の `Do` ごと中断し、`go` の後の出力を出さない。pysad は `?No variable.` を表示して Do を続行する (chrom_b、disp_b で 10 回分余分に出る)。invset (`LINE["ANGLE","BID3"]={}` など、matching21、dx_sv、dy_sv) でも同様で、SAD は中断、pysad は継続または別のエラー。
3. **TrackParticles** (disp_R56): RFSW・RAD を OFF にしても SAD は 1 周で z (R56) と dp (R66) が変わる (z=-1.39e-3、dp=.0099660)。pysad は z=0、dp=初期値 (経路長・空洞を入れていない)。2005 行中 2000 行が不一致になる。R16〜R46 は 1e-5 級の差。matching_test では戻りの先頭要素 (SAD 179、pysad 1) だけが違う。
4. **誤差の実現**: print_all_error (esrf・pf2)、print_cod_magnitude、error_ring は `SeedRandom[187]` と `Random[1000]` から作るシードは SAD と一致 (125188771、60490425、145796967)。しかし誤差を入れた後の COD/ETA の値が違う (pf2: .0101577 と .0062326)。
   さらに `Twiss["NX","***"]` が pysad では 1 要素のリスト (出力が `{{ .0062326272}`) になる。SAD では数値。
5. **不安定ラティス**: SAD は 1.59e199 や 3.18e200 を出して進む (amp_dep_ts の i≥201、dx_sv)。pysad は値 (40.78/25.72) を返すか不安定ラティス例外で止まる。安定域は全て一致。
6. **K0 の fit 変数**: `free K1 K0 K2 K0` (bump3)、`free KK1..4 K0` (calc_bump)。K0 を値として与えた COD は一致するようになった (bump2、disp_response、disp_kick_effect)。fit で更新した K0 が bump3 では反映されず (要素表示に `ANGLE =0` と `K0` が出ない)、calc_bump では別解。
7. **ON INTRA** (touscheck、touscheklife): FFS 内の `ON INTRA` は通り、Touschek 表は SAD と全点一致。INTRA ON の後の Emittance[] 報告ブロック (MINCOUP を反映した Emittance Y=8.66e-11 や Beam size eta=9.7e-6) が違う/出ない。touscheck は 298 s、touscheklife は 542 s かかる。
8. **fit の別解**: pysad の 1 回目の解は SAD と全桁一致する場合がある (pf2 chrom_r の K2 SD -26.645487186032284 / SFH 23.91648140839749)。esrf chrom、matching3_01、pf2 matching_v30_12 の 2 段目では別解 (多解)。SAD の FITP>1 の GO は再起動を繰り返すので、最終解が取れない。
9. **UI の細部**: `Momentum compact.` が小さい値で 6 桁 (.00321、SAD は .0032133、2e-3〜1e-3 の相対差)。BEND の要素表示のキーワード順が SAD は `E1 E2 K1`、pysad は `K1 E1 E2`。未定義名のエラー文言が違う (SAD: `Missing value for MKSFA*`、`Incomplete expression`、pysad: `?Undefined command or element`)。
10. **FileNotFoundError** (tune_survey_momentum): `OpenRead` の失敗が Python 例外のまま漏れる (SAD は `Read` が wrongtype エラーを出して Do を中断)。

## 5. 残る未対応機能 (頻度順)

順位表は docs/12 §4 を参照 (本数は 86 本のうち主因または副因のスクリプト数)。上位は次の通り。

1. 性能 (打ち切り 10 本 + 300 s 超 5 本): tune_survey 系、DAsv を含む disp2/disp_da/disp_test/disp_bunchlength、print_all_error、print_cod_magnitude、chromaticity_plot、touscheck 系。
2. fit の最終解が SAD と違う (多解、SAD の FITP>1 の再起動): chrom、chrom_b、chrom_r、matching3_01、matching_v30_12_r (5 本)。
3. CodCorr の Dot 長さ不一致 (4 本): esrf の COD 補正の本体を塞ぐ。
4. 不安定ラティスの扱い (8 本)、FFS エラー時の Do 中断 (5 本)、K0 の fit 変数 (3 本)、誤差の実現 (6 本)。
5. TrackParticles の z・dp (2 本)、GX/GY・BEAMSIZE・裸の要素名・USE LinacPFBT・OpenRead (各 1〜2 本)、小さな数値差 (5 本)。

## 6. 次にやると良いこと

- CodCorr の Dot の原因を特定する (4 本が通れば esrf の correct_cod 系が SAD の 270 s の run と比較できる。ただし SAD は誤差が効かないので、pysad と SAD の比較用に `FFS["MON …"]` を使わず `DELX` が効く版を探す必要がある)。
- 性能: cell calc と DAsv の高速化が 15 本を左右する。
- fit の多解は、SAD の restart の最終解 (SAD は FITP>1 で落ちる) が取れないので、`FITP 1` に落とした版で解の一致を確認する方が実用的。
- `Twiss["NX","***"]` のスカラー化、TrackParticles の z・dp、FFS エラー時の Do 中断は小さな修正で効果が見える。
