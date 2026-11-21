# FFS cell / calc / emit (PF 1981, 460 nm rad): 履歴ラティスの周期解と平衡パラメータ
種別: FFS コマンド。加速器: PF 2.5 GeV リング、建設時のオプティクス (`lattice/pf_1981_460nmrad.sad`、周長 187.0741 m、チューン 5.2987 / 4.1638)。フラグ: `ON LOG; ON RAD FLUC COD RFSW RADCOD;`。説明は [`lattice/PF_history.md`](../lattice/PF_history.md)。

## 何を見る例か
`FFS_cell_calc_emit_PF` と同じ台本の形 (ラティスを先頭に貼って `USE RING; cell calc emit;`)。PF の履歴ラティスの 1 本を、周長・チューン・エミッタンス・要素の数・直線部と標準セルの Twiss・クロマティシティで調べる。4 本 (1981・1986・1997・2005) は同じ台本で、ラティスと、Twiss を出す MARK の名前だけが違う。

## スクリプト (`FFS;` 以降の要点)
```
USE RING;  cell calc emit;
Print["circumference [m] = ", LINE["S","$$$"]];  ...            ! 周長・要素数・チューン・EMITX・SIGE・SIGZ
CountEl[];                                                      ! 偏向磁石・四極・6 極・八極・空洞の数
ShowAt["NAME"];                                                 ! その場所の BX BY EX AX と NX NY/2pi (直線部の MARK、標準セルの偏向磁石)
ChromFD[0.001];                                                 ! 色収差 = 運動量 ±0.1% のチューンの差分 (設計の6 極あり)
SF* 0; SD* 0; FFS["cell calc"]; ChromFD[0.001];                 ! 6 極を 0 にした自然値
```
- 台本の関数名は `ShowAt[]` (`At[]` と名づけると pysad は何も出さない。詳しくは `lattice/PF_history.md` §g)。

## 結果の読み方
| 量 | 値 | 出所 |
|---|---|---|
| 周長 | 187.0741 m | `LINE["S","$$$"]` |
| チューン | 5.2987 / 4.1638 | `Twiss["NX","$$$"]/2/Pi` |
| εx (k64 = gs) | 4.0337e-7 m rad | `EMITX` |
| εx (v2) | 4.6018e-7 m rad | `.pysad.txt` (v2 の既定) |
| 名前のエミッタンス | 460 nm rad | ファイル名 |
標準セル・直線部・クロマティシティ・放射積分の値は `.gs.txt` と、`lattice/PF_history.md` の表。
**ヘッダ・名前との対応**: ヘッダの NUX = 8.454、NUY = 3.300 とは合わない (要確認)。名前の 460 nm rad は v2 の値 (460.2) に近い、k64 (403.4) は −12%。

## 物理・解析式
放射積分からの εx の照合は `PF_history.md` §c (εx ≈ Cq γ² I5 / (Jx I2) が v2 の Jx で 0.3% 以内)。

## gs と pysad の比較
`.pysad.txt` は既定の v2 モードの出力 (リポジトリの他の例と同じ)。`PYSAD_COMPAT=k64` で走らせると、`.gs.txt` と、Print の出力 (周長・要素数・チューン・EMITX・SIGE・SIGZ・要素の数・各点の Twiss と位相・クロマティシティ) と emit の表の主な行が、表示桁 (10 桁) で一致する (差は 2005 の SIGE の最後の桁 1e-10 だけ。`alpha` の表示桁数は gs が 5 桁、pysad が 3 桁。εy とビームの傾きは 0 の雑音)。`tools/pf_history/cmp_gs_pysad.py` で確かめられる。v2 は `Emittance[]`・EMIT の Jx の端の項 (B75) などで、emit の表 (εx・σE・減衰・α など) が k64 と違う。光学 (周長・チューン・Twiss・クロマティシティ) は同じ。

## 未確認・限界
- 設計値のラティス・誤差なし。垂直エミッタンスは結合なしで 0。
- `.gs.txt` の `LINE["LENGTH"]` の数 (リングの要素数 + 終端) は pysad と同じ。
- ヘッダ・名前との食い違いの理由は調べていない (要確認)。
