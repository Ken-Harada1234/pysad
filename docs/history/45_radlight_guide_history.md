# 45_radlight_guide の履歴 (整理で本文から取り除いた記述)

本文 (`docs/45_radlight_guide.md`) を版のないスナップショットにするとき、取り除いた版・経過の記述を、元の文脈つきで残したもの。元の行をそのまま引用する (引用した行のうち、履歴にあたるのは版・経過を述べた部分)。

## 元 3 行 (節: 45. RADLIGHT と放射場 (RadiationField・RadiationSpectrum) の解説)

> 最終更新: 2026-11-18 (v2.1.0 の後)

## 元 21 行 (節: 1. 何ができるか)

> | PHOTONS と `PhotonList` | BEAMBEAM の光子に加え、v2.1.0 で **BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND の光子** に対応 (§6) | `beambeam.py`、`track6.py`、`ffscmd` |

## 元 22 行 (節: 1. 何ができるか)

> | SOL 領域内の要素・COORD・TCAVI・MAP のあるラティスの光子 | 未対応 (SOL 領域内は、v2.2.0 から v2 が DRIFT・QUAD・MULT だけ記録する。BEND・SEXT・OCT と k64 は記録せず乱数だけ引く) | |

## 元 95 行 (節: 3. K64 の癖 (tlinit / tlstore / tlresult。k64 モードの動作。v2 モードは §3a))

> 5. 同じ位置の点が 2 つ並ぶことがある (区間の始めの点と、直前の要素の後の点。試験の S2 で 7.8 が 2 回出る)。**複数区間 (複数周) の軌道では区間の始めの点が重なり (t が同じ点が 2 つ)、以前の pysad の `RadiationField` は `x must be strictly increasing sequence` で失敗した** (v2.1.0 で直した: 両モード。完全に同じ点を除く。時刻が戻る・同じ時刻で位置が違う点は SADError)。gs は `RadiationField` が SIGSEGV するので比較できない。

## 元 100 行 (節: 3a. v2 モードの軌道 (B65))

> K64 の §3 の癖 1〜3 (BEND の中を記録しない、横方向の座標が入らない、BEND で向きが変わらない) を、v2 モードでは直した。試験 [`test_v2_C_radlight.py`](../pysad/tests/test_v2_C_radlight.py) (3 件)、例 [`example/v2_C/C5_radlight.sad`](../example/v2_C/C5_radlight.sad)。

## 元 106 行 (節: 3a. v2 モードの軌道 (B65))

> - **使えない・未確認**: SOL 領域を含む区間 (SADError。両モード)、COORD による幾何の変化、曲げのある MULT の中の点 (出口だけ)。複数周の軌道の `RadiationField` は v2.1.0 で動く (重なった点を除く。k64・v2。`TrackParticles[d, 9, 1, 2]` の軌道 22 点・32814 点で場を返す)。高エネルギー・多数粒子の規模は未確認。

## 元 132 行 (節: 5. gs との一致の度合い)

> 注意: `RadiationField` は **gs と比べていない** (gs は SIGSEGV)。上の 6e-3 は物理の解析解との比較。**v2.1.0 で Liénard-Wiechert の厳密な数値計算と比べた** (C。独立な計算 [`pysad/tests/_c3_indep.py`](../pysad/tests/_c3_indep.py): 円弧・直線の軌道 (解析的な位置・速度・加速度) で、観測時刻 τ ごとに遅延時間 t' を `brentq` で解き、E = q/(4πε0)[(n−β)/(γ²(1−n·β)³R²) + n×((n−β)×β̇)/(c(1−n·β)³R)] を評価)。pysad の `RadiationField` は K64 の Feynman-Heaviside の式 (軌道の点列と自然スプラインの微分)。試験 [`test_v2_C3_radfield.py`](../pysad/tests/test_v2_C3_radfield.py) (8 件)、例 [`example/v2_next/C3_radfield_bend.sad`](../example/v2_next/C3_radfield_bend.sad):

## 元 151 行 (節: 6. PHOTONS と光子)

> - **対応する要素 (v2.1.0。両モード)**: BEND (K1=0: `tbend` の刻みごと、K1≠0: `tbendi` の刻みごと + 最後。ROTATE・DX・DY・DROTATE)、QUAD (`tquad`、`tsolqur` の刻み)、MULT (K0・K1・K2・DX・DY・ROTATE)、SEXT・OCT (DX・ROTATE)、磁場つきドリフト、BEAMBEAM (ビームストラールング)、BEAMBEAM + BEND (TRPT つきでも)。`tdusrnpl` (`tdusrn` と同じ乱数の引き方で光子ごとのエネルギーと位置も返す)・`radangle1` (放射角と偏光。`tgauss` を 2 回)・`tphrec` (放射点 xr = x − r(px − dpx r/2)al、`tphchge`・`tforbitgeo`・`tfgeofrac`) の直訳で、K64 と同じ乱数の順序・`pcvt%fr0` の進め方 (tbend.f:367、tbendi.f:141,178、tquad.f:71,108、tsolqu.f:453)。PhotonList の行 = [E (eV)、gx gy gz、ux uy uz、xi1 xi2 xi3、粒子番号、要素番号]。BEND・QUAD・MULT・SEXT・OCT があるラティスは Python の追跡 (`track6`) になり遅い (24 セルで 1 粒子 1 周 1 ms 程度)。

## 元 155 行 (節: 6. PHOTONS と光子)

> - **使えない・未確認**: SOL 領域内の要素・磁場なしのドリフトの光子は記録しない (乱数だけ K64 と同じに引く。位置の情報 `tsetpcvt` を領域内の要素に足し、領域の座標系の `geo_frac` の使い方を確かめる必要)。COORD・TCAVI・MAP のあるラティスは従来の高速追跡で、BEND があっても光子を記録しない。POL の追跡との併用は未確認。QUAD の光子の位置は gs と未照合 (BEND は全桁)。強い6 極・8 極で不安定になる条件。track6 の `_ndivrad` が TRPT の 10 倍を見ていなかったのを直した (POL・RADPOL・光子の追跡の TRPT の結果が gs に近づく)。

## 元 160 行 (節: 7. 限界・未確認)

> - **RadiationField の gs 参照がない** (gs は SIGSEGV で落ちる。B49。引数を変えても再確認できず)。等速運動の解析解と、v2.1.0 の Liénard-Wiechert の厳密な数値 (§5) との比較。K64 の式との照合は読みによる。

## 元 161 行 (節: 7. 限界・未確認)

> - **k64 モードでは BEND の中を記録しない**ので、そのままの BEND では放射場の主要な部分 (偏向の円弧) が入らない (K64 の仕様どおり。v2 モードは §3a)。Oide の例の流れ (B を ANGLE=0、K0=0.1) は 1 周の軌道なら pysad で通る (§2)。複数区間・複数周の軌道の `RadiationField` は v2.1.0 で動く (§3 の 5)。gs は `RadiationField` で落ちるため通しの照合は無い。

## 元 166 行 (節: 7. 限界・未確認)

> - **SOL 領域内の光子 (v2.2.0)**: v2 は領域内の DRIFT・QUAD・MULT (ANGLE=0) の光子を軌道上の位置に記録する (docs/44 の v2.2.0 の追加。3 GeV、bz = 0.1 1/m、傾き・ずれのない領域で確認)。領域内の BEND・SEXT・OCT・傾いた領域は記録しない。k64 は K64 のまま (領域内の QUAD は記録しない: B150)。

