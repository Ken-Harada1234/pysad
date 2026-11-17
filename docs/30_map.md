# MAP (ExternalMap)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `temap.f` (追跡 `temap`、Emittance[] `temape`、Twiss `qemap`、幾何 `tgmap`)。`pysad/elements.py` の `Map`、`ffscmd._install_map_hooks`。

```
MAP P1 = (L=0.5);
ExternalMap["TRACK", LINE["POSITION","P1"], nt_, x_] := (...);
ExternalMap["EMIT",  LINE["POSITION","P1"], cod_]    := {cod1, trans [, dtrans, dbeam]};
ExternalMap["OPTICS", LINE["POSITION","P1"], cod_]   := {cod1, trans};
```

MAP はユーザーが SADScript の関数 `ExternalMap` で写像を与える要素。**定義がなければ恒等で、長さも 0** (`L` は K64 では使われない。`tgmap` の既定は `pos(l+1)=pos(l)`)。
MAP を長さ L のドリフトにすると周長が L だけ増えるので、pysad は長さ 0 にしている。

## 使われる所

| 計算 | 呼び出し | pysad |
| --- | --- | --- |
| TrackParticles | `ExternalMap["TRACK", l, nt, x]`。x は (7, 粒子数): x, px(1+δ), y, py(1+δ), z, δ, 生存フラグ (番号順) | `track6n.K_MAP` → `ffscmd._map_track` |
| CALC・INS (Twiss、閉軌道) | `ExternalMap["OPTICS", l, cod]` → `{cod1, trans}` | `lattice._map_step` (`_el_cod_step`、`Map.transfer_matrix`) |
| Emittance[]・EMIT | `ExternalMap["EMIT", l, cod]` → `{cod1, trans}` か `{cod1, trans, dtrans, dbeam}` | `emit6d._tmap` |

- TRACK: 同じ形のリストが返れば座標を置き換える。返らなければ何もしない。フラグ 0 の粒子は落とす。**落ちている粒子もフラグ≠0 なら生き返る** (K64 の `doinject`。ヘルプには「復活できない」とあるが、コードは生き返らせる。gs で確認)。
- EMIT の 4 つ組: `T ← trans·T`、`D ← trans·D`、`B ← (trans+dtrans) B (trans+dtrans)ᵀ + dbeam (上三角)`、`D += dtrans (T+D)` (temape.f の順)。2 つ組: `T, D, B` に `trans`。
- **軌道 `cod1` はユーザーが返す。`trans` から作らない。** `cod1` と `trans` が食い違うと (例: `cod+kick` に結合のある `trans`) K64 の閉軌道は物理的でない値で止まる。
  pysad は固定点を求めるので別の値になる。試験は `cod1 = trans.cod + kick` で作った。

## 検証 (gs との比較、`pysad/tests/test_map.py` 5 本)

- CALC (OPTICS): チューン 1e-14、軌道 1e-11。
- Emittance[] (EMIT、`mm.cod+kick` の整合した写像、x-y 結合の回転 + せん断): ε_x ε_y ε_z、チューン、減衰率、軌道、σδ、σz が 1e-11 以内。4 つ組 (dtrans、dbeam つき) も同じ。
- TrackParticles (TRPT、20 粒子、4 周、MAP を 2 か所): 座標を書き換える写像、フラグ 0 で落とす写像、落ちた粒子の復活を含めて 10 桁一致。ExternalMap が受け取る nt、粒子数、座標の和も一致。
- 既定 (ExternalMap なし): 恒等、S は MAP の L を数えない (gs と同じ 4.1)。

## 限界・未確認

- `ExternalMap["GEO", l, geo, pos]` (幾何) は未対応 (MAP は直進、長さ 0)。`CompiledMap` (Fortran の動的コンパイル) も未対応。
- TRACK は TrackParticles (`TrackRun`) だけ。MAIN の `TRACK` 文・`track6.py` へのフォールバック・DAPERT・Touschek の経路は MAP を恒等として通す (ExternalMap を呼ばない)。
- 最初からフラグ 0 で渡した粒子の復活は未対応 (SADError)。
- ExternalMap が返す形が違うときのメッセージは近い文面にしたが、gs の文面とは同一でない。
- CALC4D では trans の 4×5 だけが効く (ヘルプ) が、pysad は 6×6 をそのまま使う (z の行は 4D の結果に影響しないので差は出ないはず。未確認)。
- Emittance[] の Orbit・Matrix 規則 (`emit_extra`) は `Map` を通す (`_tmap` を入れた) が、gs との照合は基本の規則だけ。
