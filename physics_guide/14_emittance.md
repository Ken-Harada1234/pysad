# 14. 平衡エミッタンス: $\Sigma = M\Sigma M^\top + \Delta$

**この章で分かること**: SAD の `EMITTANCE`・`Emittance[]` (K64 の `temit`、`tecalc`) が、放射の 1 周の写像 $M$ (減衰を含む) と励起 $\Delta$ から、平衡のビーム行列 $\Sigma$ (6×6 の 2 次モーメント) を求める手順。(1) $\Sigma = M\Sigma M^\top + \Delta$ の導出と、$\Delta$ を各要素の励起の輸送として組む方法 (`trans(6,12)`、`beam(42)`)。(2) `tecalc`: 固有値・正規座標・減衰の分離、チューンと減衰率、分配数 $J_x, J_y, J_z$、各モードのエミッタンス。(3) 同期運動と RFSW。(4) MINCOUP。放射積分の標準の式との比較で数値を確かめる。

**前提**: [12_linear_optics.md](12_linear_optics.md) (1 周の行列、結合)、[13_radiation.md](13_radiation.md) (励起 $\Delta$ の導出)。

**読む順**: §1 (結論) → §3 → §4 → §6 (確認)。

## 1. 結論

- **方程式** (ヘルプ equilibrium-beam-envelope): $\langle X_iX_j\rangle = M\langle X_iX_j\rangle M^\top + \Delta_{ij}$。$X = (x, p_x, y, p_y, z, \delta)$ は閉軌道からのずれ、$M$ は減衰を含む 1 周の行列、$\Delta$ は 1 周の放射の励起。1 周を線形な写像 $X_{n+1} = MX_n + \xi_n$ ($\xi_n$ は励起の確率的な擾乱、$\langle\xi\xi^\top\rangle = \Delta$) と見ると、平衡の $\Sigma = \langle X_nX_n^\top\rangle$ がこの式を満たす (§3)。
- **$\Delta$ の組み立て** (`tradke`): 要素ごとの励起 $\Delta_i$ (13 章 (13.7)) を、残りの要素の行列 $T_{i\to\mathrm{end}}$ で 1 周の終わりまで運ぶ: $\Delta = \sum_iT_{i\to\mathrm{end}}\Delta_iT_{i\to\mathrm{end}}^\top$。K64 は `beam(42)` に走る和を持ち、各要素の行列 `tr1` について `tmulbs(beam, tr1)` ($\mathrm{beam} \leftarrow T\,\mathrm{beam}\,T^\top$) をしてから、その要素の励起 `dbeam` を足す。
- **`tecalc`**: $M$ (`trans(:,1:6)`) の固有値・固有ベクトルから正規座標の行列 $R$ を作り (`teigenc`、`tnorm`、`tsymp`)、チューン ($\mathrm{Im}\log\lambda/2\pi$)、減衰率 ($\mathrm{Re}$、`rd = ri*trans(:,7:12)*r` から)、分配数 $J_k = 4\,a_k/\sum a$ を出す。正規座標の $\Sigma_n = R^{-1}\Sigma R^{-\top}$ は、各モードの 2×2 のブロックの対角が $\varepsilon_k$ になる。モード $k$ のブロックの方程式 $\Sigma_k = e^{-2a_k}P(\mu_k)\Sigma_kP(\mu_k)^\top + \Delta_k$ のトレースから、**$\varepsilon_k = \bar\Delta_k/(1 - e^{-2a_k})$** ($\bar\Delta_k = (\Delta_{n,XX} + \Delta_{n,PP})/2$)。
- **確認** ([emit_01_equilibrium.out.txt](outputs/emit_01_equilibrium.out.txt)): 8 セルの FODO リング (3 GeV、周長 41.6 m、CAVI 3 MV) で、$\varepsilon_x = 5.2996\times10^{-5}$ m rad は、放射積分の式 $C_q\gamma^2I_5/(J_xI_2)$ と 2.2% 以内、$\sigma_\delta$ は 0.5%、$J_x$ ($0.1340$) は $1 - I_4/I_2 = 0.1364$ (1.8%)、$U_0 = 2.3449$ MeV は完全に一致、$\alpha_p$ は 0.2%、$\nu_s$ は 0.4%。$\Delta_{\delta\delta}$ は 16 台の $N_\gamma\langle u^2\rangle$ の和と 0.2% 以内。各モードの $\varepsilon_k$ の式は $0.03\%$ 以内。
- **注意**: $J_x$ が小さい ($0.13$、分散が大きくて $I_4/I_2 \approx 0.86$) 例では、$\varepsilon_x \propto 1/J_x$ なので、2% の差も $J_x$ の 2% の差から来る。$\varepsilon_y = 0$ は、垂直の励起も結合もないため (実機では MINCOUP のような下限や誤差が要る、§5)。

## 2. 記号と単位

| 量 | 意味 | K64 の Fortran の変数 | pysad (`calc6d` の `emit`) |
|---|---|---|---|
| $M$ | 1 周の行列 (減衰を含む) | `trans(6,12)` の 1:6 列、`rx` | `one_turn_matrix` |
| $\Delta$ | 励起の行列 | `beam(42)` (21 成分の対称行列 × 2: 前半が励起の和、後半 `beam(22:42)` が `calint` 用) | `beam_matrix` から $\Sigma - M\Sigma M^\top$ |
| $\Sigma$ | 平衡のビーム行列 | `beam`、`beamn` (正規座標)、`beamp` | `beam_matrix` |
| $R$、$R^{-1}$ | 正規座標への変換 | `r`、`ri`、`rx` | `normal_r`、`normal_ri` |
| 固有値 | $\lambda_k$ | `ceig` | `cd` |
| チューン | $\nu_k = \mathrm{Im}\log\lambda_k/2\pi$ | `params(ipnx:ipnz) = imag(cd(4:6))/m_2pi` | `tunes` |
| 減衰率 | 1 周あたりの減衰 (振幅) $a_k$ | `params(ipdampx:ipdampz) = dble(cd(1:3))` | `damping_rate` |
| 分配数 | $J_k = 4a_k/\sum a$ | `params(ipjx:ipjz)` | |
| エミッタンス | $\varepsilon_k$ | `emitn(21)` (正規座標の 2 次モーメント) | `emittance_x`・`_y`・`_z` |
| $\alpha_p$ | 運動量圧縮率 | `alphap` | `momentum_compaction` |
| $\omega_z$ | シンクロトロン角振動数 | `omegaz` | `omegaz` |
| $U_0$ | 1 周のエネルギー損失 | `u0` ($p_0$ 単位)、`u0*pgev` | `energy_loss_eV` |

21 成分の添字 `iaidx(i,j)` は、6×6 の対称行列の下三角 ($i \ge j$) を並べたもの ($k = i(i-1)/2 + j$)。

## 3. $\Sigma = M\Sigma M^\top + \Delta$ の導出

ビームの 1 粒子の 1 周の写像を、閉軌道のまわりで線形化すると、放射 (光子の放出) が確率的な擾乱 $\xi$ ($\langle\xi\rangle = 0$、共分散 $\Delta_1$) を足す:

$$
X_{n+1} = M X_n + \xi_n,\qquad \langle\xi_n\xi_m^\top\rangle = \Delta\,\delta_{nm},\qquad \langle\xi_nX_n^\top\rangle = 0
\tag{14.1}
$$

($M$ の平均的な減衰 (平均の放射損失の線形化) は $M$ に入っている。擾乱は各周で独立で、粒子の位置と相関しない。) 両辺の外積の期待値をとって、$\Sigma_n = \langle X_nX_n^\top\rangle$ とすると

$$
\Sigma_{n+1} = M\Sigma_nM^\top + \Delta,\qquad\text{平衡 }\Sigma_{n+1} = \Sigma_n \equiv \Sigma\ \Rightarrow\ \Sigma = M\Sigma M^\top + \Delta
\tag{14.2}
$$

(クロス項 $\langle MX_n\xi_n^\top\rangle = 0$ を使った。) $M$ の固有値の絶対値が 1 未満 (減衰) なので (14.2) は一意な解を持つ。

### 3.1 $\Delta$ の組み立て

励起は 1 周の途中の位置 $s_i$ で生じる ($\Delta_i$、13 章 (13.7))。位置 $s_i$ での擾乱は、リングの始点まで残りの行列 $T_{i\to\mathrm{end}}$ で運ばれて、始点の座標への擾乱 $T_{i\to\mathrm{end}}\xi_i$ になる。1 周ぶんの擾乱の共分散は

$$
\Delta = \sum_iT_{i\to\mathrm{end}}\,\Delta_i\,T_{i\to\mathrm{end}}^\top
\tag{14.3}
$$

(各位置の擾乱は独立なので共分散が足し合わされる。) 漸化式で書くと、始点から終点に向けて要素を 1 つずつ通るとき、$B_0 = 0$、$B_{i} = T_iB_{i-1}T_i^\top + \Delta_i$ ($T_i$ は要素 $i$ の行列) で、最後の $B_N = \Delta$ である。K64 の `tmulbs(beam, trans1, calint)` が $B \leftarrow TBT^\top$ を、`beam(1:21) = beam(1:21) + dbeam` が $+\Delta_i$ を実行する。例えば `tradke` は、要素の行列を作って `tmulbs` してから `dbeam` を足す。`trans(6,12)` の 7:12 列も同じ考え方 (7:12 列に放射の微分による減衰) で、減衰は $M$ に取り込まれる (§4.1)。

**確認** ([emit_01_equilibrium.out.txt](outputs/emit_01_equilibrium.out.txt) (c)): RF をラインの先頭に置いた例 (励起は RF のあとの放射だけ) で、$\Delta = \Sigma - M\Sigma M^\top$ の $\delta\delta$ 成分は $6.742\times10^{-9}$、16 台の偏向磁石の $N_\gamma\langle u^2\rangle$ の和 (13 章。$\gamma$、$\rho$ は 13 章の例と同じ) は $6.757\times10^{-9}$ (比 $0.9978$)。$\Delta$ は対称で、固有値は $\ge -3\times10^{-25}$ (半正定値)。RF を末尾に置くと、RF の $z$–$\delta$ 結合を通って $\Delta_{\delta\delta}$ が変わる (別の計算、比 $0.90$)。

## 4. `tecalc`: 正規座標と各モード

### 4.1 固有値と減衰

$M$ (6×6、実、シンプレクティックでない: 減衰を含む) の固有値は 3 組の複素共役ペア $\lambda_k = e^{-a_k}e^{\pm i\mu_k}$ ($a_k > 0$: 減衰、$\mu_k = 2\pi\nu_k$)。K64 は `teigenc` で固有値・固有ベクトルを求め、`tnorm` で固有ベクトルを (複素共役ペアを実の 2 本に直して) 正規化し、`tsymp` でシンプレクティックな成分 (減衰がないとしたときの形) に直す。この $R$ について、$R^{-1}MR$ は 3 つの 2×2 ブロック $e^{-a_k}P(\mu_k)$ (回転 × 減衰) になる。

減衰率は、減衰を除いたシンプレクティックな $R$ を使って、`rd = ri*trans(:,7:12)*r` の対角ブロック (`trans(:,7:12)` は放射の微分による減衰の部分: 13 章 §6 の `tr1`) から

$$
\mathtt{cd}(k) = \frac{(\mathtt{rd}_{ii} + \mathtt{rd}_{i+1,i+1})/2 + \mathrm{i}\,(\mathtt{rd}_{i,i+1} - \mathtt{rd}_{i+1,i})/2}{\mathtt{cc}(i)},\quad a_k = -\mathrm{Re}\,\mathtt{cd}(k)
\tag{14.4}
$$

(`cc` は $e^{i\mu_k}$。符号は SAD の約束で、`params(ipdampx) = dble(cd(1))` は負。) 減衰時間 $\tau_k = 2\pi/(\omega_0a_k)$ (`taurdx` など、$\omega_0$ は周回の角周波数)。

**分配数と Robinson の定理**: $J_k = 4a_k/(a_x + a_y + a_z)$ (`sr = 4/sum`、`params(ipjx:ipjz)`)。定義により $J_x + J_y + J_z = 4$ ((e)。Robinson の定理は、放射積分の式 $J_x = 1 - I_4/I_2$、$J_y = 1$、$J_z = 2 + I_4/I_2$ の和が 4 になる、$U_0$ の運動量依存性からの結果)。

**確認** ((a)(b)): 減衰率は $(5.2327\times10^{-5},\ 3.9074\times10^{-4},\ 1.1189\times10^{-3})$ で、$J = (0.1340, 1.0006, 2.8654)$。放射積分の式では $(0.1364, 1, 2.8636)$。差は 1〜2%。$J_y = 1.0006$ は 1 に近い (放射積分の式では厳密に 1)。差の原因 (四極の軸ずれの放射など) は未確認 (§7)。

### 4.2 正規座標の $\Sigma$ とエミッタンス

(14.2) を正規座標 $\Sigma_n = R^{-1}\Sigma R^{-\top}$、$\Delta_n = R^{-1}\Delta R^{-\top}$ で書くと、$\Sigma_n = \Lambda\Sigma_n\Lambda^\top + \Delta_n$、$\Lambda = \mathrm{diag}(e^{-a_k}P(\mu_k))$。モード $k$ の 2×2 のブロックは $\Sigma_k = e^{-2a_k}P\Sigma_kP^\top + \Delta_k$ ($\Delta_k$ は同じブロック)。回転 $P$ は $\mathrm{tr}$ を変えないので、$\mathrm{tr}\,\Sigma_k\,(1 - e^{-2a_k}) = \mathrm{tr}\,\Delta_k$。$\Sigma_k$ が円 ($\Sigma_k = \varepsilon_kI$) なら

$$
\varepsilon_k = \frac{\bar\Delta_k}{1 - e^{-2a_k}}\simeq\frac{\bar\Delta_k}{2a_k},\qquad \bar\Delta_k\equiv\frac{\Delta_{n,XX} + \Delta_{n,PP}}{2}
\tag{14.5}
$$

(モード間の相関 $\Sigma_{n,kl}$ ($k \ne l$) は、異なる周波数の回転で平均されて 0 に近い。共鳴 ($\mu_k \pm \mu_l \approx 2\pi\cdot$整数) の近くでは大きくなる。) K64 の `btr(21,21)` は、(14.2) を 21 成分の連立方程式 $(1 - \mathrm{btr})\,\Sigma_n = \Delta_n$ に書いたもので、モード間の相関も含めて解く。上の (14.5) は、その対角ブロックのトレースだけを取り出した式。

**確認** ((d)(e)): $\Sigma_n$ の対角は $(5.2997, 5.2995, 0, 0, 7.8982, 7.8982)\times10^{-5}$ ($X$ と $P_X$ が等しい: 円)、モード間の相関は最大 $2.5\times10^{-8}$ (対角の $3\times10^{-4}$)。(14.5) の予想は $\varepsilon_x = 5.2981\times10^{-5}$、$\varepsilon_z = 7.8984\times10^{-5}$ で、$\Sigma_n$ の $5.2996\times10^{-5}$、$7.8982\times10^{-5}$ と $0.03\%$ 以内 (近似 $\bar\Delta/2a$ は $5.2978\times10^{-5}$、$7.890\times10^{-5}$)。

## 5. 同期運動と RFSW、MINCOUP

**RFSW ON** (RF あり): 縦モード ($z$、$\delta$) はシンクロトロン振動で回転し、$\mu_z \ne 0$。$M$ の (5,6)・(6,5) 成分に RF の収束が入る。`tecalc` (`synchm`): チューンが $\nu_s$ のとき

$$
\alpha_p = -\frac{\mathrm{Im}\,\mathtt{cd}(6)\,|\mathrm{Im}\,\mathtt{cd}(6)|}{(c\,2\pi/\omega_0)\,(dv_{cacc}/p_{gev})},\qquad \omega_z = |\mathrm{Im}\,\mathtt{cd}(6)|\,\omega_0/2\pi
\tag{14.6}
$$

(`dvcacc` は RF の電圧の $z$ 微分。$\mathrm{Im}\,\mathtt{cd}(6) = \mu_z$。シンクロトロン運動の $\nu_s^2 = -\alpha_p\,(\text{RF の傾き})$ から $\alpha_p$ を逆算する式。) 標準の式 $\nu_s = \sqrt{h\alpha_pV\cos\varphi_s/(2\pi E)}$ ($\sin\varphi_s = U_0/V$) との比較: $I_1/C$ から $\alpha_p = 0.39857$ ($\mathtt{calc6d}$ は $0.39780$、0.2%)、$\nu_s = 0.0503$ (calc6d $0.0505$)。バケット高さ `bh` $= \sqrt{|V/(\pi|\alpha_p|h\,p_{gev})\cdot(2\cos\varphi_s - (\pi - 2\varphi_s)\,U_0p_{gev}/V)|}$ (標準のバケットの式、未確認)。

**RFSW OFF**: 縦運動は閉じない。`rx(6,1:5) = 0`、`rx(6,6) = 1` として $\delta$ の変化を止める。`synchm` が偽のとき、$\alpha_p = -R_{56}\cdot\omega_0\,h_0/(2\pi c\,p_0)$ (`rx(5,6)` から)、`omegaz` は RF の電圧から (`sqrt(abs(alphap*...))`)。縦のエミッタンス `beam(iaidx(5,1:6)) = 0` として扱わない (同期運動がないので $\varepsilon_z$ は定義されない)。

**MINCOUP**: 垂直のエミッタンスは、完全に結合のない理想リングでは 0 (本例: $\varepsilon_y = 0$) になり、IBS や空間電荷に使うと非現実的になる。ヘルプ (equilibrium-beam-envelope): $\varepsilon_{x,y} = \max[\varepsilon^0_{x,y},\ \mathtt{MINCOUP}\times(\varepsilon^0_x + \varepsilon^0_y)]$ ($\varepsilon^0$ はラティスだけが与えるエミッタンス)。この下限は、IBS (`INTRA`) と WSPAC の反復 (`tintraconv`) の中だけで使われる (読み。数値では確認していない)。

## 6. 確認スクリプト

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `emit_01_equilibrium.py` | (a) `calc6d` のチューン・減衰・エミッタンス、(b) 放射積分の式との比較 ($U_0$、$J_x$、$\varepsilon_x$、$\sigma_\delta$、$\alpha_p$、$\nu_s$)、(c) $\Delta = \Sigma - M\Sigma M^\top$ の $\delta\delta$ 成分と 13 章の $N_\gamma\langle u^2\rangle$、半正定値性、(d) 正規座標の $\Sigma_n$ の対角、(e) (14.5) と Robinson | [emit_01_equilibrium.out.txt](outputs/emit_01_equilibrium.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/emit_01_equilibrium.py` (十数秒)。

## 7. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV、8 セルの FODO リング (QF = 0.25、QD = −0.30 [1/m]、各 0.3 m、偏向 16 台 × 1.2 m の扇形、周長 41.6 m、$\nu = (1.6165, 1.1209)$、$\eta_x = 3.08$ m)、CAVI $V$ = 3 MV、harmonic 64、$\varphi_s = 0.897$ rad、$\nu_s = 0.0505$。偏向の端面なし、6 極なし、結合なし、IBS・空間電荷なし。

**未確認**:
- `tnorm`、`tsymp`、`teigenc` の導出 (固有ベクトルの正規化の符号、Krein 符号の扱い、縮退の処理)。
- `btr(21,21)` の係数の導出と、モード間の相関 (共鳴の近く) の確認。(14.5) の対角ブロックのトレースの式のみ確認。
- 結合のあるリング (垂直のエミッタンス)、MINCOUP の下限の数値の効果、IBS (`tintraconv`)、空間電荷 (`WSPAC`)。
- ソレノイドのある場合の `tmulbs` のうち、$bzh$ の項。
- 分配数 $J_x$ の $1 - I_4/I_2$ との差 (1.8%): K64 の $\Delta$ は四極の軸ずれの放射を含み ($J_y = 1.0006$)、放射積分の式はそれを含まないため、と読めるが、寄与を分けて確認していない。
- RF を末尾に置いたときの $\Delta_{\delta\delta}$ の違い (§3.1) の定量的な説明。
- SynchroBeta (縦方向の分散 $\mathtt{ZX}$ など、シンクロトロン・ベータトロン結合) の導出。ヘルプの拡張 Twiss ($V = PBR_6H$) との対応は [12_linear_optics.md](12_linear_optics.md) §8 のとおり未確認。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| equilibrium-beam-envelope | **原文で確認** (TeX 版・平文版とも): $\langle X_iX_j\rangle = M\langle X_iX_j\rangle M^T + \Delta_{ij}$、MINCOUP の式。「$M$ the one-turn transfer matrix including the damping」。(14.2) と同じ。 |
| ヘルプに無いもの | (1) (14.1)〜(14.3) の導出、(2) 正規座標での解法 (14.5)、(3) 減衰率と分配数の定義 (14.4)、(4) 同期運動の式 (14.6)。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

特になし (この章で確認した範囲)。分配数の差 (§7) は放射積分の式との定義の違いで、疑いではない。

## 付録 C. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `trans(6,12)` | 1:6 列が $M$ (または要素の行列)、7:12 列が放射の微分 | `temit.f`、`tradke` | |
| `beam(42)` | 21 成分の対称行列 ($\Delta$ の走る和) × 2 | `temit.f` | |
| `iaidx(i,j)` | 対称行列の下三角の添字 | `temw` | |
| `r`、`ri`、`rx` | 正規座標の行列、逆行列、`trans(:,1:6)` | `tecalc` | `normal_r`、`normal_ri` |
| `ceig`・`dceig` | 固有値 | `tecalc` | `cd` |
| `cd` | 減衰 + $\mathrm{i}\,\Delta\nu$ (1:3)、$\log\lambda$ (4:6) | `tecalc` | |
| `cc` | $e^{\mathrm{i}\mu}$ | `tecalc` | |
| `rd` | `ri*trans(:,7:12)*r` (減衰の行列) | `tecalc` | |
| `rdb`・`btr` | 正規座標の放射の行列、(14.2) の 21×21 の係数 | `tecalc` | |
| `emitn`・`beamn` | 正規座標の $\Sigma_n$ | `tecalc` | |
| `sr` | $4/(a_x + a_y + a_z)$ | `tecalc` | |
| `params(ipdampx:...)`・`ipjx` | 減衰率・分配数 | `tecalc` | |
| `taurdx` など | 減衰時間 | `tecalc` | |
| `synchm` | 縦運動が有効 (`rfsw` かつ $\mu_z \ne 0$) | `tecalc` | |
| `alphap`・`omegaz`・`bh` | $\alpha_p$、$\omega_z$、バケット高さ | `tecalc` | `momentum_compaction`、`omegaz` |
| `heff`・`vceff`・`wrfeff`・`dvcacc` | 実効 harmonic、実効電圧、RF の角周波数、電圧の $z$ 微分 | `tmacro` | `heff` |
