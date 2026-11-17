# DELCOR (相関のあるずれ)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `tcorr.f` と `tvcorr.f` (FFT は `tcftr.f`)。`ffscmd._delcor`、`_tvcorr`。

```
DELCOR [ADD] corlength delx dely element [element1 ...];
```

リングに沿った相関のあるランダムなずれを、指定した要素の DX・DY に入れる。`corlength` は相関長 [m]、`delx` `dely` は大きさ。ADD が付けば足し、なければ置き換える。

## k64 モードと v2 モードの違い

以下の 3 点が違う (それ以外は両モードで同じ)。k64 モードは K64 (gs) と同じ動作、v2 モード (既定) は修正した動作。

| | k64 モード | v2 モード |
|---|---|---|
| MARK・CAVI の書き込み先 (B40) | MARK は DETR・DPX、CAVI は HARM・PHI (符号も逆) | MARK・CAVI も DX・DY に、ほかの磁石と同じ符号 (+delx) で入れる |
| corlength≠0 の波数 −1 の成分 (B62) | `cv(n)` が未設定 (2 回目以降は古い記憶域の値)。負の波数の重みが 1 つずれる | 正負の波数で対称に置く (乱数を引く順は同じ) |
| corlength≠0 の規格化 (B99) | `s = 1/sqrt(4Σw²+2)`。複素数 c の絶対値の 2 乗平均が 1 なので、DX の rms は delx/√2 | `s = 1/sqrt(2Σw²+1)`。実部・虚部それぞれの分散が 1。DX の rms = delx、DY の rms = dely |

corlength = 0 の枝 (ランダムウォーク、1 m あたり分散 1) は両モードでビット単位で同じ。

## 作り方

- 周長 `al = pos(nlat)` (INS は 2 倍)。点数 `nc = 128` から `nc/2 < al` になるまで 2 倍 (最大 65536)。
- `tvcorr`: `corlength=0` は 2 つのガウス乱数を足していくランダムウォークで、両端が合うように直線を引く。≠0 は周波数空間に `exp(-(k·2π cl/(2 al))²/4)` の重みのガウス乱数 (実部・虚部) を置き、
  `tcftr` (exp(+2πi jk/n) の和 = `n·ifft`) で戻して規格化の定数を掛ける (k64 は `1/sqrt(4Σw²+2)`、v2 は `1/sqrt(2Σw²+1)`)。乱数は `tgauss` と同じ列 (`SadRandom`、GCUT 適用)。呼び出しの最初に必要な分をまとめて引く。
- 要素ごとに、中心位置 s の周期的な線形補間の値 (実部 × delx、虚部 × dely) を入れる。同じ要素の連なり (master) には同じ値、ワイルドカード以外は最初の 1 つだけ。最後に `tfadjst` (隙間なく並んだ磁石は先頭と同じずれ)。
- 対象は BEND QUAD SEXT OCT (と MARK、CAVI)。DECA・DODECA・MULT・MONI・ドリフトは `Undefined element` 扱い (K64 の型の判定どおり)。

## 検証 (gs との比較、`pysad/tests/test_delcor.py` 3 本。2 本は CELL と INS の組み合わせ、1 本は CAVI)

- CELL の 12 セル: 相関長 2 m (QF QD)、相関長 0 (B)、ADD とワイルドカード (Q*) の 3 回の DELCOR のあとのチューン 1e-14、閉軌道 1e-10。
- INS (4 セル、SEXT・BEND を含む、`*` で全要素): チューンと軌道が一致 (乱数の列が gs と同じ)。

## v2 モードの根拠と使える範囲 (`test_v2_C_delcor.py`、`test_v2_C2_delcor.py`。gs との比較ではなく統計・解析式)

- 規格化: 場の rms (500 m、512 点、400 標本) は x0 = 3、15、5000 m で実部・虚部とも 1 ± (2、4、8)% 以内 (k64 は 0.708、0.703)。DELCOR コマンド (12 セルのリング 62.4 m、四極 24 台、corlength 5 m、150 回) で rms(DX)/delx = 1.021、rms(DY)/dely = 1.043 (統計誤差 約 2.6%。k64 は 0.713、0.746)。
- 相関関数 ⟨c(s) c*(s+Δ)⟩/⟨|c|²⟩ = exp(−Δ²/(2x0²)) が 5% 以内 (虚部 2%)。相関長 1000 km ではリング全体が同じ値 (差 < 1e-9) で、rms 1.0 ± 5% のガウス分布。
- B40: 相関長 ≫ 周長なら全要素のずれが同じ値・同じ符号になる (QF・RFC・ST が 5% 以内。相関長 5 km では全桁同じ)。
- 使える範囲: v2 モード、corlength > 0。rms は多数回の平均としての値で、GCUT を掛ければその分小さい。相関長が格子間隔 (周長/nc、nc ≥ 128、INS は周長の 2 倍) より短いと線形補間で rms が下がる (間隔 = 相関長/4 で約 1%)。corlength = 0 の枝は別の量 (分散 = delx²·s(1−s/周長) [m]) で、delx の意味は「1 m あたりの rms」。
- 未確認・判断: 規格化 (delx を DX の rms と読む) は作者の文書による裏づけが無い推定 (SADHelp に DELCOR の項は無い)。K64 の値が必要なら k64 モード、または delx に 1/√2 を掛ける。CAVI の符号 (+delx) は決定済みだが、K64 が SOL と同じ形の式 ((設計の DX) − delx) を写した意図があるなら逆。
- 例: `example/v2_C/C2w_delcor_rms.sad`。

## 見つけたこと・限界 (k64 モード。gs との比較)

- **MARK には DX・DY ではなく DETR と DPX が入る (k64 モード。v2 は DX・DY)。** K64 は `latt(j)+15` と `+17` に書くが、値の添字では `ky_DETR_Mark=15`、`ky_DPX_Mark=17` (`ky_DX_Mark=16`、`ky_DY_Mark=18`)。
  一つずれた書き方に見える (docs/21 B40)。gs と合わせるためそのまま写した。INS の始点 (MARK の DETR・DPX) が変わる。
- SOL は master が 0 なので K64 でも対象にならない (何も起きない)。**CAVI** は (k64 モードでは) K64 の添字のずれで DX・DY ではなく **HARM と PHI** に書く (tcorr.f の `le=latt+3`、`le1=le+1`)。add でなければ設計の DX・DY (`idvalc+13`, `+14`) から delx・dely を引いた値、add なら -delx・-dely を足す。gs で確かめた (`CAVI RFC=(VOLT=1.5e6 HARM=64 DX=1e-5 DY=2e-5)`、`DELCOR 3 1e-4 1e-4 RFC` で HARM=-.00015145432812411818、PHI=.00011480572629164664。pysad は同じ値)。HARM は実数のまま持つ。
- CAVI の書き込みは `test_delcor_cavity_writes_harm_and_phi_like_gs` が gs の値と比べる。`docs/history/work_unimpl.md` の一覧に「CAVI は未対応」の古い記述が残るが、現在は実装済みで、この文書が現状。
- 各要素の DX・DY の添字 (QUAD は 5、BEND は 9 など) は gs の結果で合うことを確認 (QUAD・SEXT・BEND)。OCT は対象に含めた DELCOR で、チューンと QF の DX が gs と 12 桁で一致した (2026-11-16 夜)。ただし OC 自体の DX は gs も pysad も 0 を返すので、OCT への効果そのものを見る強い比較ではない。
- 'DELCOR' 以外の Rule 系の誤差コマンド (ERROR の DELX など) は従来どおり (docs/11)。
- **2 回目以降の `DELCOR` で corlength≠0 の場合、gs と値が一致しない (乱数の列は同じ)。k64 モードのこと。v2 は cv(n) を対称に置くので、この不一致は起きない (gs とは違う値になる)。** K64 の `tvcorr` の corlength≠0 の枝は `cv(n)` (最後の要素) を設定しない
  (`cv(2..n/2)` と `cv(n/2+1..n-1)` と `cv(1)` だけ) ので、`tcftr` は `ktaloc` で得た記憶の古い内容を読む。1 回目は 0 で pysad (0) と一致するが、2 回目以降は前の確保の中身が入って
  数 % 違う (corlength=0 の枝と ADD の回は一致。確かめた 4 通りの組み合わせ: 0→3 と 3→2 と 3→3 は不一致、3→0 と 0→0 は一致)。再現できないので pysad は 0 のままにし、docs/21 に B62 として残した。
