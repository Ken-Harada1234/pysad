# flag_SPAC
フラグ・05_tracking・ラティス `00_lattice/spac_ring.sad` (新規: FODO_ring の DRIFT に `RADIUS=0.03`、運動量 20 MeV)・`NORFSW NORAD NOFLUC`

## 何を見る例か
`SPAC` を ON にすると追跡で空間電荷が入る (実粒子数 `PBUNCH`、マクロ粒子数 `NP`)。ダクトは円筒対称とし、半径を `DRIFT` と `MULT` の `RADIUS` で与える: RADIUS>0 なら空間電荷を入れ、かつその半径で損失の絞りを置く。RADIUS<0 なら空間電荷なし、半径 |RADIUS| の絞りだけ。RADIUS=0 は空間電荷なし (SADHelp)。

## スクリプト
```
NP = 200; PBUNCH = 1e13;
NOSPAC;  r = TrackParticles[{1, beam}, 1, 1, nt];    ! nt = 1, 20, 100 (A)
SPAC;    r = TrackParticles[{1, beam}, 1, 1, nt];
Element["RADIUS","D1"] = -4e-4; ...                   ! (B) 負の RADIUS = 絞りだけ
```
- 低エネルギー (20 MeV、$\gamma=40$) にするため、ラティス側の `MOMENTUM = 20MEV` を使う。`FFS` に入った後で `MOMENTUM=20e6` と書くと gs では効くが pysad は追跡に反映しない (下)。
- 粒子は $x\sim3\times10^{-4}$ sin 型の決定的な 200 粒子 (乱数を使わない)。

## 結果の読み方
| | gs | pysad |
|---|---|---|
| (A) SPAC と NOSPAC の違い (PBUNCH=1e13、1・20・100 周) | **違いなし** (rms x、y が同じ) | 違いなし |
| (A) 初期ラティスの `Element["RADIUS","D1"]` | 0.03 | 0 (pysad は RADIUS を持たない) |
| (B) RADIUS=−0.4 mm、SPAC、2 周: 生存数 | 84 / 200 | 200 / 200 |
| (B) NOSPAC で同じ半径: 生存数 | 200 / 200 | 200 / 200 |

(A) の gs の結果: rms x [m] = 2.59080e-4 (1 周)、2.77989e-4 (20 周)、1.62915e-4 (100 周) (SPAC あり・なしとも同じ)。

## 物理・解析式
空間電荷の収束力の強さ (丸いガウスビーム、$\beta\approx1$) は $K_{sc}=\dfrac{r_e\lambda}{\gamma^3\sigma^2}$ ($\lambda=N/(\sqrt{2\pi}\sigma_z)$ は線密度、$r_e=2.82\times10^{-15}$ m)。$N=10^{13}$、$\sigma_z\approx7$ mm、$\sigma=1.5\times10^{-4}$ m、$\gamma=40$ を入れると $K_{sc}\sim10^3$ m$^{-2}$ で、チューンシフト $\Delta\nu\simeq\frac{1}{4\pi}\oint\beta K_{sc}\,ds\gg1$ (桁外れに大きい)。空間電荷が働いているなら追跡結果が大きく変わるはずだが、gs の (A) には全く出ていない (係数は未確認で、桁の見積もりだけ)。

## gs と pysad の比較
× **pysad は未対応**: `SPAC` は値を保持するだけで追跡に効かず、RADIUS も読まない (docs/08 §7 の未対応一覧)。(B) の絞りも効かない (生存数 200)。
さらに gs 側も、(A) の `RADIUS=+0.03` で空間電荷の効果は何も出なかった。RADIUS>0 で小さい値 (0.4 mm など、粒子のビームが半径内に収まらない設定) では gs は不正終了 (SIGSEGV 相当のバックトレース) した。したがって gs でも「空間電荷が働く」ことは確認できていない。働くのは負の RADIUS の絞りだけ。

別件: `MOMENTUM = 20e6;` を FFS の中 (`USE` の前) で再設定した最初の版では、gs は新しい運動量を使い (NOSPAC でも rms x が 3 GeV のときと 0.07〜14% 違った)、当時の pysad は 3 GeV のままの結果を返した (`Emittance[]` も同様に変わらなかった: SIGE=1.1914e-3 のまま)。この例ではラティスの先頭 (MAIN レベル) に運動量を置いて避けた。この差は 2026-11-17 に直した (FFS 内の MOMENTUM は USE で MAIN 側の値に戻らない)。

## 未確認・限界
- 空間電荷が gs で働く条件 (PBUNCH、RADIUS、NP の組み合わせ) を見つけられなかった。`SPAC` の追跡は、K64 では `tspac.f` 系で行われるが、この例では確認していない。
- 20 MeV、200 マクロ粒子 (決定的な分布) のみ。
