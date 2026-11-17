# COORD (座標変換)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `tcoord.f` (追跡)、`tcoorde.f` / `qcoord.f` (6 次元の写像と行列)。`pysad/coord.py`。幾何 (`LINE["GX"]` など) は既存の `geo.py` (`tfgeo1` の icCOORD)。

```
COORD CO = (DX=1e-3 DY=-2e-3 DZ=5e-3 CHI1=0.01 CHI2=-0.02 CHI3=0.03 DIR=0);
```

- DIR=0 (既定): 入口で DX DY DZ だけずらして (位置の変換は 1 を引く向き) から CHI1・CHI2・CHI3 で回す。DIR≠0: 逆変換 (回して、出口で戻す)。
- DX DY DZ は誤差の属性ではなく変換のパラメータ。pysad は `Coord.params` に置き、ずれ・回転の誤差のコードに混ざらないようにしている (BEAMBEAM も同じ)。

## 使われる所

| 計算 | 実装 |
| --- | --- |
| 追跡 (TrackParticles) | `track6n` の `K_COORD` (`coord.k_coord`)。`tturn` の icCOORD は `exit` なので、tapert までの要素数 (la) を減らさない |
| Twiss (CALC、INS、CELL) | `lattice._el_cod_step` の Coord (qcoord → tcoorde の COD 写像と解析行列)。COD が設計軌道上でも動く |
| Emittance[]、CALC6D、Touschek | `emit6d._tcoord` (`piece`、z 行のシンプレクティック組み直しなし) |
| 幾何 | `geo._coord` |

## 検証 (gs との比較、`pysad/tests/test_coord_elem.py` 12 本)

- 追跡: DX/DY/DZ のみ、CHI のみ、すべて、DIR=1、CHI3 のみ (200 粒子、3 周): 相対 0。
- Twiss と幾何 (INS、BEND を含む並びで 12 要素の BX BY NX R1 R3 EX DX DY と GX GY GZ GCHI1-3): 1e-11 以下。
- Emittance[] (リングの中に COORD を 2 個: 順変換と逆変換、または非対称): エミッタンス・チューン・COD・運動量収縮因子が 5e-10 以下。

## 限界

- track6 (低速追跡) は COORD を扱わない (NotImplementedError)。
- COORD を含むリングでの FIT (解析微分の経路) は gs と全桁一致を確認した (2026-11-16 夜: COORD の対 C1・C2 (DX DY DZ CHI3 付き) を含むリングで NX・NY=0.2、QF K1・QD K1 を FREE。K1・チューンとも全桁一致)。
  COORD のパラメータそのものを FIT 変数にする場合は未確認のまま (数値微分にもならない可能性がある)。
- TrackParticles の `dv` は COORD で更新されない (tcoord.f のまま)。
