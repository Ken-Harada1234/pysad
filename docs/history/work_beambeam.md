# BEAMBEAM 実装の経緯 (2026-11-14)

最終更新: 2026-11-16

- 根拠: ビルドに使われるのは bbstrhl.f (mk/sad.obj.mk)。bb.f は古い版で未使用。SADHelp に記述なし。
- 実装: `pysad/beambeam.py` (setup = bbinit、bb_kick = beambeam、cerrf、gauinv)、`elements.BeamBeam` (Marker の子)、parser (BEAMBEAM → BeamBeam)、track6n (K_BB)。track6 (低速) は NotImplementedError。
- それまで parser は BEAMBEAM を黙ってドリフトにしていた (L=0 なので実害は追跡で相互作用が無視されること)。
- 検証: local/wake/genbb.py (b1..b12) を gs と pysad で実行して比較。10 ケースは完全一致、傾き (b6, b9) は 7e-10。試験は pysad/tests/test_beambeam.py (14 本)。全体 347 passed、6 skipped。
- 発見: 傾き回転の cost/sint の取り違え (docs/21 B34)。直訳で再現。
- 続き: BSTRL (ビームストラールング) を実装 (乱数は合わせなくてよい、という指示)。pysad の tran (SAD と同じ LCG) を使うと、SeedRandom[7]・330 粒子以下では gs と乱数列まで一致した。330 粒子を超えると gs の列が変わる (360 粒子で最初の粒子から違う。gs は再現する) ので、20000 粒子で分布を比べた (放射の割合、KS 検定、平均)。試験は bs1・bs10 (乱数列まで一致)。
- 訂正: 最初の比較の `TrackParticles[..., 5, 1]` の 5 は終点の要素番号で、周回数ではなかった (BB を 1 回通るだけ、「5 周」は誤り)。`10, 1, 5` (全周を 5 周) で取り直し、全ケースが一致 (傾き 3e-18、BSTRL は乱数列まで)。
- 続き (未対応の対応): LuminosityMonitor (gs と一致)、BEAMBEAM の SOL 領域内 (gs が SIGSEGV: B35、pysad は NotImplementedError)、wake との併用 (一致)、光学 (gs も恒等、一致)、PHOTONS (BEAMBEAM の光子の PhotonList と stdout の write、乱数は 1 光子ごとに 2 回)。
- 発見: RAD FLUC のあとの `Random[]` が gs と違った → gs は TRPT で ndivrad が 10 倍 (tffs.f)、pysad は落としていた (四極の分割 gs 4、pysad 1)。直して一致。TRPT の放射ありの結果は gs に近づく。
- 未対応: 他の要素の光子 (tphrec、四極・BEND など)。そのため RAD FLUC で他の要素が光子を出すと PhotonList の中身と、あとの乱数列が gs とずれる。

## 2026-11-16 の確認・追記

- 2026-11-16 確認: 本文は変化なし (`pysad/beambeam.py` の更新日時は 11/14 17:02 で、この文書の最後の記述より後の整理のみと見られるが、差分の内容は未確認)。11/16 に BEAMBEAM の追加の作業は無い。未対応のまま: 他の要素の光子 (tphrec)、BEAMBEAM の SOL 領域内 (gs が SIGSEGV: B35)。現状の解説は docs/24 と、11/16 に新設の docs/39。根拠: work_unimpl.md に BEAMBEAM の新しい項目が無い、docs/24 §5。
