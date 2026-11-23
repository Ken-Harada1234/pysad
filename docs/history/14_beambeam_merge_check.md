# 14_beambeam の統合の照合

docs/14_beambeam.md を、旧 docs/39_beambeam_guide.md (土台)、旧 docs/24_beambeam.md、旧 docs/32_gaussiancoulomb_bbbrem.md (BBBrem1 の部分) から統合したときの節の対応。

## 旧 39_beambeam_guide.md

| 旧の節 | 新 14 での場所 |
|---|---|
| §0 先に知っておくこと | 1 |
| §1 パラメータ表 | 2.1 |
| §2 bbstrhl のキック | 3.1 |
| §3 BSTRL | 3.2 |
| §4 PHOTONS | 4.1 (台本は 2.3) |
| §5 LuminosityMonitor | 4.2 (意味)、2.4 (台本) |
| §6 光学での扱い | 4.3 |
| §7 動く台本 (7.1〜7.4) | 2.2・2.3・2.5 |
| §8 BBBrem | 4.4 |
| §9 tlum | 4.5 |
| §10 gs との一致の表 | history/14_beambeam_validation.md。本文は 6 の要約 |
| §11 K64 の癖 | 5 |
| §12 限界・未確認 | 7 |
| §13 関連文書 | 9 と 8 (コード) |

## 旧 24_beambeam.md

| 旧の節 | 新 14 での場所 |
|---|---|
| §1 使い方・キーワード表 | 2.1 (39 §1 と同内容) |
| §2 計算の流れ | 3.1 (39 §2 と同内容) |
| §2b ビームストラールング | 3.2 (39 §3 と同内容。PHOTONS・LuminosityMonitor の記述は 4.1・4.2) |
| §3 検証 | history/14_beambeam_validation.md。内容は 6 に要約 |
| §4 K64 の癖 (ndivrad の式を含む) | 5 (ndivrad の式と、光子の放射角を 39 に無かった事実として追加) |
| §5 未対応・限界 | 7 (STURN・COUPLE、光子の write を追加。ほかは 39 §12 と同じ) |

## 旧 32_gaussiancoulomb_bbbrem.md

- BBBrem1 の節: 4.4 (戻り値、乱数 tran nran=3 を追加)。検証 (相対 1e-15) は 6。
- GaussianCoulomb 系の節: 本文書の担当外 (docs/06_sadscript_library.md が受ける)。4.4 から docs/06 へリンクを置いた。

## 意図して削った・移したもの

- 落としたものは無い。gs との一致の表と旧 24 §3 の検証の一覧は history/14_beambeam_validation.md に移した。

## 追記: BBBrem1 の検証と限界

旧 32 の BBBrem1 の検証 (SeedRandom[3] のあと 6 事象と、√s = 10〜40 GeV の 4 事象が相対 1e-15 で一致) と限界 (棄却時の診断 print を写していない、`Parallelize` での動作は未確認) を、14 の 4.4 (本文)、6 (一致の 1 行)、7 (限界) に入れた。一致の桁の詳細は history/14_beambeam_validation.md の末尾。
