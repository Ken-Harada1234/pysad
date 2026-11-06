# 作業ログ: INS (非周期) で先頭 MARK に初期軌道オフセットを与えたときの Twiss/チューン (2026-10-15)

担当範囲: `pysad/lattice.py`, `pysad/ffscmd.py` (`match6.py` は変更なし)。
課題: esrf の `amp_dep_ts` (`<work-esrf>/amp_dep_ts*.sad`, `modules_amp_dep_ts.n`) の 1 点目が pysad 47.5924241431 / SAD 47.5800386899。
**結果: 47.5800386898 (SAD と 1e-10)、NY 17.6200386132 (SAD と一致)。2 点目も 47.5784879349 vs 47.5784879354。**

## 1. SAD の INS の経路 (ソースで確認)

| 段階 | ソース | 内容 |
|---|---|---|
| (a) スクリプトの設定 | `modules_amp_dep_ts.n` `CalcAmpDepTS` | `ininame=LINE["NAME",1]` (= `MKS`) に対し `LINE["DX",ininame]=…; LINE["DPX",…]; LINE["DY",…]; LINE["DPY",…]` で先頭 MARK の値を書き、`FFS["ins calc"]`。`SetInitialParameters` は `FFS["cell calc emit fix"]` の後 `Twiss["AX",1]` 等を読む。BXI コマンドや Element は使っていない |
| (b) CALC の入口 | `tffscalc.f:93` | `twmov(1,twiss,nlat,ndim,.true.)`: 先頭 MARK の値 (cmp%value(1:ntfun)) を `twiss(1,0,:)` に読む。ntfun は `orbitcal`(既定 .true., track.f:128) なら 28 全部。`NX=NY=0`、`DDP += dp(0)` (dp(0)=0) |
| (c) CELL | `qcell.f` `qcell1` | `cod=twiss(lb,idp,DX:DDP)` を **qcod の初期値**にして閉軌道を探索、周期解の Twiss を `twiss(lb,0,:)` に作る |
| (d) CELL の後 | `tffscalc.f:119` | 安定 (`stabx,staby,stabz`) なら `twmov(1,twiss,nlat,ndim,.false.)` で周期解 (AX,BX,AY,BY,EX,EPX,EY,EPY,R1–R4,DETR,DX,DPX,DY,DPY,DZ,DDP) を **先頭 MARK の値に書き戻す**。これが「CELL → INS に切り替えた直後の初期値」であり、`LINE["BX","MKS"]` で読める値になる |
| (e) INS | `qcell.f:281-289`、`qtwiss.f` `qtwiss1` | qcod は呼ばない。`cod=twiss(ipa,DX:DDP)` (MARK の値) から各要素を通し、**要素ごとに軌道を更新して軌道上の 4×5 行列 `trans`** で `qmat2twiss` を呼ぶ。R1–R4 は MARK の値、`normal = DETR<1`、`|値|<1e-30` は 0 に丸める |
| (f) 軌道オフセットの行き先 | 同上 | DX/DPX/DY/DPY は**初期軌道 `cod(1:4)`** に入る (分散 EX には入らない)。DZ は cod(5)、DDP は cod(6) = DP0 (`DP0` は `LINE["DDP",1]` の別名: track.f:39) |
| (g) 行列 | `qtwiss1` `case` 各要素 | 要素ごとの `q*` ルーチンが **解析式の行列**を返す (BEND は `tbende` の `trans1`)。結合 (`coup`) があれば `qmat2twiss` の結合の式、なければ R だけ `t1..t4` で更新する式 |

## 2. 原因 (pysad の何が違っていたか)

1. **INS の入口が CELL の周期解でなかった**: `ffscmd._ins_rows` は先頭 MARK の raw 値 (ラティスファイルの `BX=3.89`…) を読むだけで、CELL CALC は MARK に何も書き戻していなかった (`lattice._twmov_to_marker` は属性を書くだけで raw 値に反映されず、未使用)。amp_dep_ts の `Print` 前の値は、周期解 (βx=4.643) ではなく MARK の 3.89 から伝搬していたため、チューンが 0.012 ずれた。
2. **CELL の qcod が MARK の軌道ではなく直前の `_cod_prev` から始まっていた**: SAD は (c) の通り MARK の DX..DDP から探索を始める。収束判定が緩い (conv≈1e-14 → cod の誤差 ~1e-7) ので、**初期値が違うと Twiss が 1e-7 程度変わる**。amp_dep_ts は INS 用に MARK に軌道を書いた後すぐ `cell calc emit fix` するので、この「ノイズ」が SetInitialParameters → 次の振幅に伝わる。
3. **Twiss の伝搬が結合なしの式だった**: `_transport_twiss` は R1–R4 を一定のまま、x と y を独立に伝搬していた。軌道が y 方向にずれると六極が斜め四極として効いて結合し (DY=5e-4 で NY が 0.009 ずれる)、R と分散 (EY, EPY) も要素ごとに変わる。
4. **BEND の行列が数値微分で、SAD の解析式とずれていた**: SAD の `tbedge`/`tbendebody0` の行列は写像の厳密な微分ではなく (px, py≠0 で 2e-8 ずれる)、軌道が大きい INS では EX/BX が 1e-7〜1e-6 違った。K1≠0 の BEND も、`emit6d._tbende` が入口・出口エッジを数値微分にしていたので 4e-9 ずれた。
5. 小さな抜け: INS の初期軌道の DZ、MARK の DETR (`normal` 判定)、`DP0` ↔ `LINE["DDP",1]` の別名。

## 3. 修正

`pysad/ffscmd.py`
- `FFSState._twmov_to_mark(rows)` (新規): CELL CALC (4D と CALC6D 両方) の後に周期解の入口 Twiss と閉軌道を先頭 MARK の raw 値に書く (tffscalc.f:119 の直訳。PSIX=PSIY=0)。`_TWMOV_ROW` に対応表。
- `calc()` の CELL: 閉軌道探索の初期値を先頭 MARK の DX/DPX/DY/DPY にする (tffscalc.f:93)。
- `_ins_rows`: 入口は MARK の Twiss + R + **DETR**、初期軌道は DX/DPX/DY/DPY/**DZ**。
- `get_value`/`set_value`: 先頭 MARK の `DDP` と `DP0` を同じ値にする (track.f:39 の別名)。
`pysad/lattice.py`
- `_transport_twiss`: **SAD `qmat2twiss` の直訳に置き換え** (結合 coup/非結合、normal/入れ替えの 4 通り、R1–R4・DETR・EX..EPY・位相を更新)。要素行列に x-y 結合ブロックがあれば結合の式。CELL の `twiss()` もこれを使う。
- `ins_twiss`: 入口で qtwiss1 と同じ `|値|<1e-30` の丸めと DETR の作り直し (normal=DETR<1)。
- BEND の解析行列 (新規): `_tbedgemaxwell_mat`, `_tbedgedrift_mat`, `_tbedgebody_mat`, `_tbendebody0_mat` (tbfrie.f / tbende.f の `trans`/`trans1` の直訳), `_bend_edge_chain`, `_bend_cod_matrix` (K1=0, F1=0, 据え付け誤差なし)。K1≠0 は `_bend_k1_cod_matrix` で入口・出口エッジだけ解析行列に差し替え (本体は `emit6d` の `_tbendecorr`/`_tbendebody`)。誤差・F1 フリンジありは従来どおり数値微分。z 行は 4D の CALC では使われないので、`tbendebody0` の `h0/h1emit**3·aln` 項は省略。
- `_twmov_to_marker`: Marker オブジェクトの属性は書かない (Lattice.twiss() はオフモメンタムや微分でも呼ばれる。MARK への書き戻しは ffscmd 側の本計算のみ)。`cell_initial` は参考として残す。

`emit6d.py` と `match6.py` は変更していない (`_tbende` の K1≠0 は lattice 側で組み直した)。

## 4. 照合 (SAD は同じスクリプトを WSL の `gs` で実行)

ベンチは `bench_level22*.sad` (先頭にラティス `pysad/tests/kek-ls_v03_12.lat` と同じ内容を連結)。`.result` が SAD、`_pysad.result` が pysad。差は `max(|SAD|,1)` で割った値の最大 (0 近傍の値は絶対差)。

| ベンチ | 内容 | 最大差 (pysad vs SAD) |
|---|---|---|
| 22a | RAD OFF。CELL → INS (オフセットなし、DX=1e-3、DY=5e-4+DPX=1e-4、DX=2.6e-5)、BXI/AXI 変更、CELL に戻る (MARK の軌道から探索) | 2.4e-10 |
| 22b | CELL を一度も実行しない INS: ラティスの MARK、AXI/BXI、EX/EPX/EY/EPY、R1–R4 (結合)、軌道を順に与える | 8.7e-11 |
| 22c | RAD ON COD RFSW RADCOD。CELL → INS に切り替えただけ (DP0=0 と 0.003) | 1.5e-10 |
| 22d | RAD ON、小さなオフセット (DX=4e-5, DY=8e-5) | 2.9e-11 |
| 22e | DP0=0.004、四極オフセット (閉軌道≠0)、2 回目の CELL、DP0=0 に戻して INS、DZ=1e-3 | INS 側 5e-10。CELL (DP0≠0, 閉軌道≠0) の βy 等が 2e-9 (SAD の qcod 収束ノイズの範囲) |
| 22f | INS (DX=1e-3) の全要素 k=1..60 の NX BX AX EX DX NY BY | 3.4e-10 (k≥4。k=1–3 は AX≈0 の相対) |
| 22g | 単独要素 (BID5/BID4/BID3H/BPI6 K1≠0/QFA1/SDA/DRIFT) に軌道・Twiss を与えた INS | 2.0e-12 |
| amp_dep_ts 1 点目 (<work-esrf>/amp_dep_ts_q.sad: `{i,1,22,20}` の 2 点) | NX | pysad 47.5800386898 / SAD 47.5800386899 (修正前 47.5924241431)。NY 17.6200386132 (修正前 17.6831884572、SAD と一致)。2 点目 47.5784879349 vs 47.5784879354、NY 17.6192760580 vs 17.6192760583 |

修正の寄与 (切り分けは作業中の途中結果から): 初期値 (MARK 書き戻し) だけ直した段階の 22f (DX=1e-3) は k=60 で NX 4e-9、BX 5e-8、EX 4e-7 の差が残り、BEND を解析行列にして ≤3e-10 になった。y 軌道が 0 でないケース (22a INS3: DY=5e-4, DPX=1e-4) は、初期値だけの段階で NY が SAD と 0.009 違い、結合の式 (`qmat2twiss`) と BEND の解析行列を入れて 1.4e-10 になった。22a の「INS5 の後の CELL に戻る」は MARK の軌道 (DX=2.6e-5) から qcod を始める修正で 2.7e-9 → 9e-13 (SAD も MARK の値から始めるため、初期値で Twiss が 1e-7 違った)。

副産物: **bench_level14a (esrf calc_chromaticity、DP0 スキャン 21 点) が SAD に一致するようになった**。修正前はチューンが SAD と 1e-7 (最大 2.5e-7) 違い、3 次当てはめの係数が相対 2e-5 違ったが、修正後は 21 点とも 1e-10 以内、当てはめの 1 次係数は -854.629590701 (SAD -854.629595042、相対 5e-9)。
切り分け: 修正前のコードに「qcod を MARK の軌道から始める」変更だけ入れても NX の差は 2.5e-7 のまま (別に作った `before2` で確認)。つまり主因は **BEND の解析行列** で、DP0≠0 では BEND の行列が軌道上 (δ≠0) の値になり、数値微分との 2e-8 のずれが 1e-7 のチューンの差になっていた。以前「オフモメンタムの Twiss は SAD の qcod 収束ノイズで 1e-6 揺らぐ」(docs/05) としていたものの大部分はこの行列の差だった。

## 5. 回帰

- 修正前後の出力を比較 (修正前は `scratchpad/ins_work/before/pysad` に退避した修正前のコード):
  `bench_level13a/13b/16a/18a/19a` は絶対値 1e-9 以上の数値がすべて同じ (AX≈1e-14 などの 0 近傍の丸めのみ違う)。`14a` は上のとおり SAD に近づいた。
- `python -m pytest -q pysad/tests`: **57 passed** (修正前も同じ)。ただし `test_bend_feeddown.py` (収集時に `sys.exit`)、`test_coord.py`、`test_sext_radiation.py` (`track_tthin()` の引数) は**修正前から収集に失敗する**ので `--ignore` で外した (今回の変更とは無関係。`test_bend_feeddown.py` は `Bend.constant_kick` の DX フィードダウンの期待値が現行の elements.py と合わない)。
- `QUICK=1 python run_regression.py`: 21 ベンチ **FAIL 0** (KNOWN は 3c の結合時の分散 EX/EPX [差 3.48e-7 → 2.84e-7 に改善]、14a の EMIT dz 表示桁のみ。いずれも従来から)。
- 遅いもの `python -u run_regression.py 7 12t 13a 16a`: **FAIL 0** (L7 1.7e-10、12t 3.5e-11、13a は TWB 4.7e-8 の KNOWN (13A-2 解なし、従来から)、16a は DA の損失周数 1 個 ±1 の KNOWN (従来から))。QUICK と合わせて 25 ベンチ全部 FAIL 0。

## 6. 残ること・注意

- **CELL の収束ノイズ**: 閉軌道が 0 でない CELL (22e) は SAD と 1e-9〜1e-8 違う。qcod の停止点が SAD の解析行列と pysad の数値微分 (四極・六極・MULT は中心差分) の違いで変わるため。四極/薄肉多極の解析行列を直訳すれば詰まるはず (今回は不要)。
- **z (DZ) の伝搬**は SAD と違う (22d の `Twiss["DZ",100]` が 1.6%、BEND の dvemit 項などの z 行を省いたため)。4D の光学には影響しない。INS の初期 DZ は MARK から読むだけ。
- **MARK の書き戻しは raw (名前単位)**: 先頭が `MKS`、同名の `MKS.2` が別にあるとき、pysad は同じ raw を共有するので `LINE["BX","MKS.2"]` も周期解になる (SAD は先頭の MARK だけ)。INS/CELL の計算には影響しない。
- `Matcher._calc_optics`/`_matcher_optics` (GO の評価中) は MARK を書き換えない。GO 後の最終 CALC (`go()` の末尾) で書かれる (SAD は評価ごとに書く)。matching の結果への影響は確認していない (回帰で変化なし)。
- INS + CALC6D フラグ (`qcell6d` の非周期版) は未実装のまま。
- 作業中の一時ファイルは `scratchpad/ins_work/` (cmpnum.py, diffnum.py, runcmp.sh, before/)。ベンチは `bench_level22a..g.sad` と `.result`/`_pysad.result` を残した。


## 7. 2026-10-16 追記: §6 の「CELL の収束ノイズ」の解決 (docs/history/work_emit_report.md)

22e の CELL (閉軌道≠0、DP0=0.004) の SAD との差 5.3e-10 (表示 12 桁) は、四極・薄肉多極の数値微分の雑音 (相対 1e-10) が原因だった。
tquade (tsolrote、tqfrie、tqlfre、tsolque) と tthine の解析行列を直訳して lattice の CALC 経路に入れたところ **1.4e-11** になった。
INS の `Twiss["DX"..]` 等の行も ≤1.4e-11。「SAD の qcod 収束ノイズの範囲」とした説明は半分正しく (最後に評価した軌道で行列を作る)、半分は行列そのものの雑音だった。
