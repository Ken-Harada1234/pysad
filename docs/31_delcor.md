# DELCOR (相関のあるずれ)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `tcorr.f` と `tvcorr.f` (FFT は `tcftr.f`)。`ffscmd._delcor`、`_tvcorr`。

```
DELCOR [ADD] corlength delx dely element [element1 ...];
```

リングに沿った相関のあるランダムなずれを、指定した要素の DX・DY に入れる。`corlength` は相関長 [m]、`delx` `dely` は大きさ。ADD が付けば足し、なければ置き換える。

## 作り方

- 周長 `al = pos(nlat)` (INS は 2 倍)。点数 `nc = 128` から `nc/2 < al` になるまで 2 倍 (最大 65536)。
- `tvcorr`: `corlength=0` は 2 つのガウス乱数を足していくランダムウォークで、両端が合うように直線を引く。≠0 は周波数空間に `exp(-(k·2π cl/(2 al))²/4)` の重みのガウス乱数 (実部・虚部) を置き、
  `tcftr` (exp(+2πi jk/n) の和 = `n·ifft`) で戻して `1/sqrt(4Σw²+2)` を掛ける。乱数は `tgauss` と同じ列 (`SadRandom`、GCUT 適用)。呼び出しの最初に必要な分をまとめて引く。
- 要素ごとに、中心位置 s の周期的な線形補間の値 (実部 × delx、虚部 × dely) を入れる。同じ要素の連なり (master) には同じ値、ワイルドカード以外は最初の 1 つだけ。最後に `tfadjst` (隙間なく並んだ磁石は先頭と同じずれ)。
- 対象は BEND QUAD SEXT OCT (と MARK)。DECA・DODECA・MULT・MONI・ドリフトは `Undefined element` 扱い (K64 の型の判定どおり)。

## 検証 (gs との比較、`pysad/tests/test_delcor.py` 3 本。2 本は CELL と INS の組み合わせ、1 本は CAVI)

- CELL の 12 セル: 相関長 2 m (QF QD)、相関長 0 (B)、ADD とワイルドカード (Q*) の 3 回の DELCOR のあとのチューン 1e-14、閉軌道 1e-10。
- INS (4 セル、SEXT・BEND を含む、`*` で全要素): チューンと軌道が一致 (乱数の列が gs と同じ)。

## 見つけたこと・限界

- **MARK には DX・DY ではなく DETR と DPX が入る。** K64 は `latt(j)+15` と `+17` に書くが、値の添字では `ky_DETR_Mark=15`、`ky_DPX_Mark=17` (`ky_DX_Mark=16`、`ky_DY_Mark=18`)。
  一つずれた書き方に見える (docs/21 B40)。gs と合わせるためそのまま写した。INS の始点 (MARK の DETR・DPX) が変わる。
- SOL は master が 0 なので K64 でも対象にならない (何も起きない)。**CAVI** は K64 の添字のずれで DX・DY ではなく **HARM と PHI** に書く (tcorr.f の `le=latt+3`、`le1=le+1`)。add でなければ設計の DX・DY (`idvalc+13`, `+14`) から delx・dely を引いた値、add なら -delx・-dely を足す。gs で確かめた (`CAVI RFC=(VOLT=1.5e6 HARM=64 DX=1e-5 DY=2e-5)`、`DELCOR 3 1e-4 1e-4 RFC` で HARM=-.00015145432812411818、PHI=.00011480572629164664。pysad は同じ値)。HARM は実数のまま持つ。
- CAVI の書き込みは `test_delcor_cavity_writes_harm_and_phi_like_gs` が gs の値と比べる。`docs/history/work_unimpl.md` の一覧に「CAVI は未対応」の古い記述が残るが、現在は実装済みで、この文書が現状。
- 各要素の DX・DY の添字 (QUAD は 5、BEND は 9 など) は gs の結果で合うことを確認 (QUAD・SEXT・BEND)。OCT は対象に含めた DELCOR で、チューンと QF の DX が gs と 12 桁で一致した (2026-11-16 夜)。ただし OC 自体の DX は gs も pysad も 0 を返すので、OCT への効果そのものを見る強い比較ではない。
- 'DELCOR' 以外の Rule 系の誤差コマンド (ERROR の DELX など) は従来どおり (docs/11)。
- **2 回目以降の `DELCOR` で corlength≠0 の場合、gs と値が一致しない (乱数の列は同じ)。** K64 の `tvcorr` の corlength≠0 の枝は `cv(n)` (最後の要素) を設定しない
  (`cv(2..n/2)` と `cv(n/2+1..n-1)` と `cv(1)` だけ) ので、`tcftr` は `ktaloc` で得た記憶の古い内容を読む。1 回目は 0 で pysad (0) と一致するが、2 回目以降は前の確保の中身が入って
  数 % 違う (corlength=0 の枝と ADD の回は一致。確かめた 4 通りの組み合わせ: 0→3 と 3→2 と 3→3 は不一致、3→0 と 0→0 は一致)。再現できないので pysad は 0 のままにし、docs/21 に B62 として残した。
