# FFS_DUMP (DUMP と誤差コマンド)
種別: FFS コマンド。フォルダ: 02_optics_cod。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。

## 何を見る例か
誤差コマンド `DELX DELY DTHETA DL DK DDK 量 パターン` で成分ごとの誤差 (ずれ・回転・長さ・強さ) を入れ、`DUMP パターン` で入っている誤差を表で見る。`UNIFORM` を付けると乱数を使わず指定した量をそのまま入れる (`RANDOM` は 量×ガウス乱数)。`SeedRandom[n]` で乱数の種を固定できる。

| コマンド | 入る先 | 対象 |
|---|---|---|
| DELX / DELY | DX / DY | BEND QUAD SEXT … |
| DTHETA | ROTATE (BEND は DROTATE) | QUAD … |
| DL | L | DRIFT |
| DK | 既定キーワード (相対誤差 1+量) | 磁石類 |
| DDK | K0 (BEND) | BEND SOL |

## スクリプト
```
DELX UNIFORM 1E-4 QF.1|QF.2;    DELY UNIFORM -2E-4 QF.1|QD.3;    DTHETA UNIFORM 1E-3 QD.2;
DL UNIFORM 1E-3 D1.1;           DK UNIFORM 1E-3 QF.3;            DDK UNIFORM 1E-5 B.1;
DUMP QF.1|QF.2|QF.3|QD.2|QD.3|D1.1|B.1|Q*.4;
SeedRandom[7];    DELX 1E-4 SF.*;    DUMP SF.1|SF.2|SF.3;
```
(成分は `名前.番号`。複数の指定は `|` でつなぐ。)

## 結果の読み方
`DUMP` の列: `element DL(mm) DK DDK DTHETA(mrad) DELX(mm) DELY(mm)`。誤差の入っている成分の行だけが出る。
- gs: `QF.1` に DELX=0.1 mm、DELY=-0.2 mm、`QD.2` に DTHETA=1 mrad (1e-3 rad)、`QF.3` に DK=1e-3。乱数の `SF.2` の DELX=0.10616 mm、`SF.3` の DELX=-0.02600 mm。

## 物理・解析式
なし (入力どおりの値が出る)。乱数 (`SeedRandom[7]` のガウス乱数 × 1e-4 m) が gs と pysad で全桁一致するのが確認点。

## gs と pysad の比較
**○ 一致 (2026-11-17 に対応)。** 差の原因は gs の `DUMP` が語を 20 文字 (`character*20`) で切ること: 長い `QF.1|QF.2|QF.3|QD.2|QD.3|D1.1|B.1|Q*.4` は先頭の `QF.1|QF.2|QF.3|QD.2|` だけが使われるので D1.1 (DL) と B.1 (DDK) は出なかった (短い `DUMP QF.1|D1.1|B.1;` は gs も出す)。pysad も 20 文字で切る。どの要素にも合わない語のときは全要素 (誤差のあるもの) を出す。
もう 1 点: `DDK` は K0 のキーワードのある要素 (BEND・MULT・TCAVI) だけに入る (gs で QUAD に入らないのを確認。pysad は QUAD にも入れていた)。`QF.1|QF.2` の 2 つ目に誤差が入らないのは、誤差コマンドが `*%{}` を含まない語では最初の一致だけに作用するため (gs・pysad とも)。

## 未確認・限界
- `RANDOM` と `COHERENT`/`INCOHERENT`、`ADD`/`PUT` のオプション、`DBEAM` は未確認 (`DBEAM` は pysad が読み飛ばす、docs/03)。誤差がある状態での閉軌道は flag_COD (errors_ring) で扱う。
