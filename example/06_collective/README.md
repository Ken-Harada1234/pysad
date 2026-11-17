# 06_collective (集団効果)

使ったラティスは `00_lattice/` の `FODO_ring.sad`、`FODO_ring_V4MV.sad`、新規の `bb_ring.sad` (BEAMBEAM)、`wake_ring.sad` (MARK WK)、`map_ring.sad` (MAP MP)、`FODO_ring_V4MV_0p5GeV.sad` (0.5 GeV)。
各例は 4 ファイル。一致の記号: ○ 一致、△ 差あり、× gs 失敗または pysad 未対応。

| 名前 | 1 行説明 | gs と pysad |
|---|---|---|
| [elem_BEAMBEAM](elem_BEAMBEAM.md) | 弱-強ビーム・ビームのキックとチューンシフト | ○ |
| [fn_WakeFunction](fn_WakeFunction.md) | ウェイク関数の表と、ウェイクによるキックの最大値 | ○ |
| [flag_LWAKE_TWAKE_WAKEOPT](flag_LWAKE_TWAKE_WAKEOPT.md) | LWAKE / TWAKE / WAKEOPT の組み合わせ | ○ |
| [fn_CSRInit_CSRTrack](fn_CSRInit_CSRTrack.md) | CSRInit と CSRTrack (縦の線形写像) | ○ |
| [fn_TouschekLifetime](fn_TouschekLifetime.md) | Touschek 寿命 (INTRA + EMIT) | ○ (k64 モード。v2 の既定では寿命が約 2% 違う: Jx の修正 B75) |
| [fn_GaussianCoulomb_BBBrem](fn_GaussianCoulomb_BBBrem.md) | ガウス分布の力、BBBrem1 の事象 | ○ |
| [flag_INTRA](flag_INTRA.md) | バンチ内散乱 (IBS) による平衡値の増加 | ○ |
| [fn_ExternalMap](fn_ExternalMap.md) | MAP 要素の ExternalMap (TRACK/OPTICS) | ○ |

図: `elem_BEAMBEAM_tuneshift.png`、`fn_GaussianCoulomb_BBBrem_force.png`。

## このフォルダで分かったこと
- 8 例すべてで gs と pysad が一致した。式 (キックの傾き、ウェイクの最大値、ガウスの力、薄いレンズのチューンなど) で照合できた量は、式とも合う。
- ウェイクの表は、追跡では最後の 1 行が使われない (K64 の癖)。表の使い方は fn_WakeFunction.md。
- 光学の WAKEOPT の縦のキック (DDP) は軌道の数が奇数だと K64 が 1 つずらす癖があり、式と合わない (flag_LWAKE_TWAKE_WAKEOPT.md)。
- CSRTrack の分割は、対称分割の式と 1e-4〜1e-3 で一致するが、K64 の厳密な分割は未確認。
