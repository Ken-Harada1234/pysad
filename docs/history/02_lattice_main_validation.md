# 02_lattice_main の検証の記録 (K64 との比較の詳細)

docs/02_lattice_main.md §4 (旧 02 §4「実例: 3 GeV 20 セルのリング (Level 12)」) から移した、K64 との要素ごとの比較の元の記述。

## 4. 実例: 3 GeV 20 セルのリング (Level 12)

- ラティス: `<試験用ラティス>.lat` (209 行)。`MOMENTUM = 3GEV;` の後に DRIFT/BEND/QUAD/SEXT/OCT/CAVI/MARK/MONI を定義 (式 `L=.413/2`, `K3=-6000*0.0125`, `VD=(ROTATE=-90 DEG)`, MKS/MKE に Twiss・EMITX・DP)、`LINE ARCAH … ARCAH2 … ARCA=(ARCAH -ARCAH2) RING=(MKS 20×ARCA + RFC MKE);`、その後 BEND/QUAD/SEXT (+DRCT) を 3 回再定義 (最後が有効。再定義は全キーワードを書いている)。
- SAD 側: `bench_level12.sad` = 上記 + FFS 部 (`USE RING; LINE["NAME"|"TYPE"|"L"|"K1"|"ANGLE"|"K2"|"K3"|"E1"|"ROTATE","*"]` を出力し `CELL; CALC;` → Twiss、`RFSW; emit;`)。結果 `bench_level12.result`。
- pysad 側: `bench_level12_pysad.py` が `parser.load_sad_file(<試験用ラティス>.lat)` → `build_lattice('RING')` し、SAD の `LINE[]` 出力と要素ごとに比較する (pysad の K1·L, K2·L, K3·L を SAD の積分値と、Bend.E1 を SAD の E1×ANGLE と比較)。結果 `bench_level12_pysad.result`:

| 項目 | pysad | SAD | 差 |
|---|---|---|---|
| 要素数 / 周長 | 3563 / 603.120000 m | NLAT−1 = 3563 / 603.120000 m | 一致 |
| 要素名 (順序込み) | 不一致 0 | | |
| L, ANGLE, K1, K2, K3, E1, ROTATE | 最大差 K1 4.9e-11, ANGLE 4.0e-11, K2 4.5e-11 | (`$FORM="15.10"` の 10 桁出力) | 表示桁の丸め |
| Twiss 末尾 `$$$` (BX, AX, BY, AY, EX, EPX, NX, NY) | | | ≤ 4e-11 |
| Twiss `MKC` | | | ≤ 4e-11 |
| emit (εx, σE, σz, α, U0, Qx/Qy/Qz, 減衰) | | | 相対差 ≤ 2e-5 (SAD 5 桁出力) |

ANGLE=0 の `VD`/`HD` の扱いは §3.2 (Bend のまま保持) のとおりで、要素数・名前・周長は SAD と一致する。反転ライン `-ARCAH2` の中の要素も (`REV_EL` 経由で) 同じ順序になることをこの比較が保証している (ARCAH2 の BEND に E1/E2 は無いので入替の効果は Level 12 では見えない)。


(元の §4 の末尾の文)
ANGLE=0 の `VD`/`HD` の扱いは §3.2 (Bend のまま保持) のとおりで、要素数・名前・周長は SAD と一致する。反転ライン `-ARCAH2` の中の要素も (`REV_EL` 経由で) 同じ順序になることをこの比較が保証している (ARCAH2 の BEND に E1/E2 は無いので入替の効果は Level 12 では見えない)。
