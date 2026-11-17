# 04_emittance_radiation 一覧
ラティスは `00_lattice/FODO_ring_V4MV.sad` (FODO_ring の RF 電圧だけ 4 MV に変更)。`RINGQ` を使う (RADINT の NaN を避けるため)。式の検算は `formula_check.py` → `formula_check.txt`。

| ファイル | 説明 | gs と pysad |
|---|---|---|
| FFS_EMITTANCE | EMIT の報告と放射積分の式による検算 (U0, τ, εx, σδ, νs, σz) | ○ |
| FFS_RADINT | RADINT の I1〜I5 と自前の積分 | ○ |
| fn_Emittance | Emittance[] の規則 | ○ |
| flag_RFSW_RAD_RADCOD_FLUC | 4 つのフラグの効き方 | △ (Emittance[] は一致、FLUC の追跡は乱数列が違う) |
| flag_RADTAPER | RADTAPER の前後 | × pysad 未対応 (値が変わらない) |
| flag_CALC6D | CALC6D の縦の光学 | ○ |
| fn_SynchroBetaEmittance | シンクロ・ベータ結合の emitx (平面リング) | ○ (結合のあるラティスは gs が落ちる) |

## このフォルダで分かったこと
- RADINT の $J_x$ (0.9415) と EMIT の $J_x$ (0.9651) は 2.4% 違う。$I_4$ の端の項を 2 倍にすると EMIT に近づく (FFS_EMITTANCE.md)。
- `FODO_ring` の RF 1.5 MV は $U_0$ より小さく、バケットがない。
