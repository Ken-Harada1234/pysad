# 07_spin_solenoid_special (スピン・ソレノイド・特殊)

使ったラティスは `00_lattice/` の `FODO_ring_V4MV.sad` (POL)、新規の `sol_ring.sad`、`multa_ring.sad`、`coord_ring.sad`、`radlight_line.sad`。
各例は 4 ファイル。一致の記号: ○ 一致、△ 差あり、× gs 失敗または pysad 未対応。

| 名前 | 1 行説明 | gs と pysad |
|---|---|---|
| [flag_POL_RADPOL](flag_POL_RADPOL.md) | 偏極時間・スピンチューン・Sokolov-Ternov の立ち上がり | ○ (k64 モード。v2 の既定ではスピン追跡の `spin angle` が gs と違う: 原因は未確認) |
| [elem_SOL](elem_SOL.md) | ソレノイド領域の 4×4 伝達行列 (結合) | ○ |
| [elem_MULT_ANGLE](elem_MULT_ANGLE.md) | ANGLE 付き MULT (BEND と同じ・複合機能型) | ○ |
| [elem_COORD](elem_COORD.md) | 座標変換 (ずらし・回転)、挟んだリングと打ち消し | ○ |
| [flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum](flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum.md) | 軌道の記録とスペクトル (ガウスパルスと照合) | △ (軌道・スペクトルは k64 モードで ○。v2 の既定は BEND の中も記録するので点数が 13 で gs の 11 と違う: B65。**RadiationField は gs が SIGSEGV で失敗**、pysad は動く) |
| [fn_SymplecticJ](fn_SymplecticJ.md) | シンプレクティック行列 J と M^T J M = J | ○ |

図: `flag_RADLIGHT_spectrum.png`。

## このフォルダで分かったこと
- スピンチューン $a\gamma=6.808148$、Sokolov-Ternov の偏極時間 (143.3652 s、式 143.3651 s)、`sy` の立ち上がりは式と一致し、gs と pysad も同じ。`NORADPOL` で立ち上がりが止まる。
- ソレノイドの結合 (`SC`、`S²/K` の項) は式の 4×4 行列と 8 桁一致。SOL 内の四極と RAD の組み合わせは gs が信用できない (癖) ので避けた。
- COORD の逆変換 (`DIR=1`) は、`DX` を同じ値にし、回転角の符号を反転するだけで元に戻る (`tcoord.f` の `DIR≠0` の分岐では `DY`、`DZ` の符号が反転する)。
- `RadiationField` は gs (この K64 ビルド) が SIGSEGV で落ちるため、値は比べられない。
