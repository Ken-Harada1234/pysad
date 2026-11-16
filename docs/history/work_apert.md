# APERT 要素の追跡 (2026-11-17)

担当: pysad/track6n.py、track6.py、parser.py、ffscmd.py、track6n_ext.py、tests/test_apert.py、docs/08_tracking.md §2.4

## 問題
pysad は追跡で APERT を L=0 のドリフトとして扱い、粒子を落とさなかった。`k_apert` は LOSSAMPL/LOSSDZ の判定だけだった。gs は APERT の楕円・矩形・ずれ・回転で落とす (`example/05_tracking/elem_APERT_LOSSAMPL`: gs {1,1,1,0,0,1,0}、旧 pysad は全部 1)。

## 直したこと (K64: tapert.f 全体、tturn.f:417-423、SADHelp の APERT)
- `track6n.k_apert_gen`: tapert の直訳 (楕円 `AX AY`、矩形 `min/max(DX1,DX2)` ほか、`DX DY`、`ROTATE`、`eli` で矩形の幅が 0 なら ±1e100、`|z|≤LOSSDZ`、JDPX/JDPY/DP の再入射、落ちた粒子を末尾の生きている粒子と入れ替え、報告の記録)。`k_apert` (LOSSAMPL) は同じ関数の薄い包み。`k_apert_el` が APERT 要素用。
- `Prog._compile` の `K_APERT`、`turn` で APERT の位置で判定して `la = 15` (la は減らさない)。TRPT では報告を APERT のときだけにする (`_rep`)。
- parser は APERT を `Drift(L=0)` + `_apert=True` + キーワードの属性で作る (光学は従来どおり)。
- track6 (Python): `_apert_pass`、`track_particles(alive=, zlost=)`。ffscmd の fn_track は track6 に落ちたときも alive を持ち回る (戻りの flag)。
- SOL 領域内の APERT は numba・Python とも NotImplementedError (tsol.f:789 の tapert1 は未翻訳)。

## 照合 (gs。3 GeV の FODO リング 24 セル、RFSW/RAD なし、16 粒子、NPARA=1)
楕円、矩形、楕円+矩形、DX DY、ROTATE、全キーワード 0、APERT 2 個、3 周: flag と落ちた粒子の座標が numba・Python 版とも一致 (座標差 0)。JDPX JDPY DP: 粒子 ≤11 個で一致。
`tests/test_apert.py` に参照値。例 `example/05_tracking/elem_APERT_LOSSAMPL` は gs と一致に更新。

## 未確認
- JDPX/JDPY: 粒子 12 個以上で gs の値が合わない (gs が別の乱数列を使う。pysad は粒子 ID 順に 2 個ずつ引く。原因は未調査)。
- RAD・FLUC つき、TRPT、MAIN の TRACK の損失報告、DAPERT での APERT は gs と比べていない。
- track6 は LOSSAMPL の判定 (最初の要素・15 要素ごと) を持たない (従来どおり)。
