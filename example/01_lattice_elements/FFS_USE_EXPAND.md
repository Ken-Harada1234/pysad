# FFS_USE_EXPAND (USE と EXPAND / NOEXPAND)
種別: FFS コマンド。フォルダ: 01_lattice_elements。ラティス: `00_lattice/FODO_ring.sad` (24 セル、周長 105.6 m、3 GeV)。フラグ: CELL、`KEEPEXP`/`CALEXP` (既定 CALEXP)。

## 何を見る例か
**`USE ラインの名前`** で FFS が使うラインを決める (このとき光学の条件などは消える)。`EXPAND` は「要素 (名前) の値を、同じ名前の全成分に配る」操作で、**成分ごとの個別のずれ (誤差) は消える**。
`LINE["K1","QF.5"] = 0.6` で 1 台の成分だけ値を変え (誤差を想定)、`CALC`・`CALC NOEXPAND`・`USE EXPAND`・`USE NOEXPAND` のどれでそのずれが残るかを見る。

## スクリプト
```
LINE["K1","QF.5"] = 0.6;                          ! 5 番目の QF 成分だけ
CALC;                       ! 既定は EXPAND: ずれは消える → K1=0.7 に戻る
LINE["K1","QF.5"] = 0.6;   CALC NOEXPAND;         ! ずれを残して計算
EXPAND;                                             ! 全部消す
LINE["K1","QF.5"] = 0.6;   USE EXPAND RING; CELL; CALC;     ! 消える
LINE["K1","QF.5"] = 0.6;   USE NOEXPAND RING; CELL; CALC;   ! 残る
```

## 結果の読み方 (gs)
| 操作 | QF.5 の K1 | NX |
|---|---|---|
| 設計値 | 0.7 | 6.696329 |
| `CALC` (既定の EXPAND) | 0.7 (ずれが消える) | 6.696329 |
| `CALC NOEXPAND` | 0.6 (残る) | 6.637402 |
| `EXPAND` のあと | 0.7 | |
| `USE EXPAND` | 0.7 | 6.696329 |
| `USE NOEXPAND` | 0.6 | 6.637402 |

NOEXPAND にすると `???-FFS-Info-Element values are not expanded.` が出て、以後の `CALC` の表示が `NOCALEXP` になる (この設定は `KEEPEXP` で次の CALC にも引き継がれる)。

## 物理・解析式
QF 1 台の K1 を 0.7→0.6 に弱めたときのチューンの変化は、線形近似で $\Delta\nu_x\approx\frac1{4\pi}\beta_x\,\Delta K_1$ ($\beta_x$=QF.5 の中心、約 6.8 m、$\Delta K_1=-0.1$) $=-0.054$ (実測 -0.05893 は同符号・同程度。β の平均値を使った見積もりなので桁だけ確認)。

## gs と pysad の比較
**○ 一致 (2026-11-17 に対応)。** `USE EXPAND/NOEXPAND RING`、`CALC NOEXPAND`・`CALC EXPAND`・`GO NOEXPAND` が gs と同じ動作になった: NOEXPAND で成分ごとのずれ (QF.5 の K1=0.6、NX=6.637402) が残り、`???-FFS-Info-Element values are not expanded.` と `NOCALEXP` の表示が出る。`NOEXPAND` は KEEPEXP なので以後の `CALC` にも残る (`CALC EXPAND` か `USE` で戻る)。`SAVE`・`EXPAND` は代表 (最初の 1 台) の値を全成分に揃える。TrackParticles は成分ごとの値をそのまま使う (CALC・Emittance[] で揃う)。

## 未確認・限界
- `QF.3 K1 0.9;` のように `名前.番号` を先頭にした値設定は、gs では式として評価されて何も設定されない。pysad も同じ (何もしない) にした。成分だけの設定は `LINE["K1","QF.3"]=0.9`。
- `NOKEEPEXP` (CALC のたびに tfvcomp で EXPAND を決める): 最初の版では未対応だった。v2 で gs に合わせて対応した (`docs/47_v2_changes.md`。gs の実測と一致と記録されている)。ただし pysad の `tfvcomp` は常に偽 (成分ごとの FREE が無い)。この例では使っておらず、再実行していない。
- `KEEPEXP`/`CALEXP` フラグ、`USE BeamLine[...]` は未確認。
