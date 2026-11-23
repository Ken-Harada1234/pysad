# 14. 平衡エミッタンス: $\Sigma = M\Sigma M^\top + \Delta$

**この章で分かること**: SAD の `EMITTANCE`・`Emittance[]` (K64 の `temit`、`tecalc`) が、放射の 1 周の写像 $M$ (減衰を含む) と励起 $\Delta$ から、平衡のビーム行列 $\Sigma$ (6×6 の 2 次モーメント) を求める手順。(1) $\Sigma = M\Sigma M^\top + \Delta$ の導出と、$\Delta$ を各要素の励起の輸送として組む方法 (`trans(6,12)`、`beam(42)`)。(2) `tecalc`: 固有ベクトル・正規座標 (`teigenc`、`tnorm`、`tsymp`)、チューンと減衰率、分配数 $J_x, J_y, J_z$、各モードのエミッタンス、21 成分の連立方程式 `btr`。(3) 同期運動と RFSW。(4) 有限のシンクロトロンチューンによる補正と、結合のあるリングの垂直エミッタンス。放射積分の標準の式との比較で数値を確かめる。

**前提**: [12_linear_optics.md](12_linear_optics.md) (1 周の行列、結合)、[13_radiation.md](13_radiation.md) (励起 $\Delta$ の導出)。

**読む順**: §1 (結論) → §3 → §4 → §6 (有限の $\nu_s$) → §8 (確認)。

## 1. 結論

- **方程式** (ヘルプ equilibrium-beam-envelope): $\langle X_iX_j\rangle = M\langle X_iX_j\rangle M^\top + \Delta_{ij}$。$X = (x, p_x, y, p_y, z, \delta)$ は閉軌道からのずれ、$M$ は減衰を含む 1 周の行列、$\Delta$ は 1 周の放射の励起。K64 は $M = T + D$ と分ける ($T$: 放射を除いた 1 周のヤコビアン = `trans(:,1:6)`、$D$: 放射による微分 = `trans(:,7:12)`)。1 周を線形な写像 $X_{n+1} = MX_n + \xi_n$ ($\xi_n$ は励起の確率的な擾乱、$\langle\xi\xi^\top\rangle = \Delta$) と見ると、平衡の $\Sigma = \langle X_nX_n^\top\rangle$ がこの式を満たす (§3)。
- **$\Delta$ の組み立て** (`tradke`): 要素ごとの励起 $\Delta_i$ (13 章) を、残りの要素の行列 $T_{i\to\mathrm{end}}$ で 1 周の終わりまで運ぶ: $\Delta = \sum_iT_{i\to\mathrm{end}}\Delta_iT_{i\to\mathrm{end}}^\top$。K64 は `beam(42)` に走る和を持ち、各要素の行列 `tr1` について `tmulbs(beam, tr1)` ($\mathrm{beam} \leftarrow T\,\mathrm{beam}\,T^\top$) をしてから、その要素の励起 `dbeam` を足す。減衰側 $D$ も同じ形の和 $D = \sum_iT_{i\to\mathrm{end}}D_iT_{0\to i}$ で、右から $T_{0\to i}$ (RF の kick を含む) が掛かる (§3.2)。
- **`tecalc`**: $T$ (`trans(:,1:6)`、シンプレクティック) の固有値・固有ベクトルから正規座標の行列 $R$ を作り (`teigenc`、`tnorm`、`tsymp`: §4.1)、チューン ($\mathrm{Im}\log\lambda/2\pi$) を出す。減衰率は、$D$ を 1 次の摂動として固有値の変化から求める ($\mathtt{rd}=R^{-1}DR$ の対角ブロック: §4.2)。分配数は $J_k = 4a_k/\sum a$。正規座標の $\Sigma_n = R^{-1}\Sigma R^{-\top}$ は、モード $k$ の 2×2 のブロックの対角が $\varepsilon_k$ になる。モード $k$ のブロックの方程式のトレースから、**$\varepsilon_k = \bar\Delta_k/(1 - e^{-2a_k})$** ($\bar\Delta_k = (\Delta_{n,XX} + \Delta_{n,PP})/2$: §4.3)。K64 は 21 成分の連立方程式 `btr` を解く (§4.4)。
- **確認** ([emit_01_equilibrium.out.txt](outputs/emit_01_equilibrium.out.txt)): 8 セルの FODO リング (3 GeV、周長 41.6 m、CAVI 3 MV) で、$\varepsilon_x = 5.2996\times10^{-5}$ m rad は、放射積分の式 $C_q\gamma^2I_5/(J_xI_2)$ と 2.2% 以内、$\sigma_\delta$ は 0.5%、$J_x$ ($0.1340$) は $1 - I_4/I_2 = 0.1364$ (1.8%)、$U_0 = 2.3449$ MeV は完全に一致、$\alpha_p$ は 0.2%、$\nu_s$ は 0.4%。各モードの $\varepsilon_k$ の式は $0.03\%$ 以内。
- **標準の式との 1〜2% の差の原因は、有限のシンクロトロンチューン** (§6、[emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt))。放射積分の式は $\nu_s\to0$ の極限の式で、$\nu_s = 0.05$ では水平のシンクロ・ベータ結合の補正が入る。3 GeV で加速電圧だけを変えると、$\varepsilon_x$ の比 $-1$ は $\nu_s^2$ にほぼ比例する ($(\text{比}-1)/\nu_s^2 = 8.7$)。RF の kick を 0 に取る極限では、$J_x$・$\varepsilon_x$ とも標準の式に $2\times10^{-4}$ で一致する。このリングには軸ずれも誤差も無く、四極の軸ずれによる放射は原因ではない。
- **注意**: $J_x$ が小さい ($0.13$、分散が大きくて $I_4/I_2 \approx 0.86$) 例では、$\varepsilon_x \propto 1/J_x$ なので、$J_x$ への補正が $\varepsilon_x$ に大きく効く。$\varepsilon_y = 0$ は、垂直の励起も結合もないため (結合があると §7、MINCOUP は §5)。

## 2. 記号と単位

| 量 | 意味 | K64 の Fortran の変数 | pysad (`calc6d` の `emit`) |
|---|---|---|---|
| $M$ | 1 周の行列 (減衰を含む) = $T + D$ | `trans(6,12)` の 1:6 列 + 7:12 列 | `one_turn_matrix` |
| $T$ | 放射を除いた 1 周のヤコビアン (シンプレクティック) | `trans(:,1:6)`、`rx` | `rx` |
| $D$ | 放射の微分による減衰の行列 | `trans(:,7:12)` | (`st.D`) |
| $\Delta$ | 励起の行列 | `beam(42)` (21 成分の対称行列 × 2: 前半が励起の和、後半 `beam(22:42)` が `calint` 用) | (`st.B`) |
| $\Sigma$ | 平衡のビーム行列 | `beam`、`beamn` (正規座標)、`beamp` | `beam_matrix` |
| $R$、$R^{-1}$ | 正規座標への変換 | `r`、`ri` | `normal_r`、`normal_ri` |
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

### 3.1 方程式の導出と解の存在

ビームの 1 粒子の 1 周の写像を、閉軌道のまわりで線形化すると、放射 (光子の放出) が確率的な擾乱 $\xi$ ($\langle\xi\rangle = 0$、共分散 $\Delta$) を足す:

$$
X_{n+1} = M X_n + \xi_n,\qquad \langle\xi_n\xi_m^\top\rangle = \Delta\,\delta_{nm},\qquad \langle\xi_nX_n^\top\rangle = 0
\tag{14.1}
$$

($M$ の平均的な減衰 (平均の放射損失の線形化) は $M$ に入っている。擾乱は各周で独立で、粒子の位置と相関しない。) $\Sigma_n = \langle X_nX_n^\top\rangle$ として、(14.1) の両辺の外積の期待値をとる。右辺を展開する:

$$
\langle X_{n+1}X_{n+1}^\top\rangle = M\langle X_nX_n^\top\rangle M^\top + M\langle X_n\xi_n^\top\rangle + \langle\xi_nX_n^\top\rangle M^\top + \langle\xi_n\xi_n^\top\rangle
$$

$\xi_n$ は $X_n$ と独立で $\langle\xi_n\rangle = 0$ なので、$\langle X_n\xi_n^\top\rangle = \langle X_n\rangle\langle\xi_n\rangle^\top = 0$ (第 2・3 項は消える)。第 4 項は $\Delta$。したがって

$$
\Sigma_{n+1} = M\Sigma_nM^\top + \Delta,\qquad\text{平衡 }\Sigma_{n+1} = \Sigma_n \equiv \Sigma\ \Rightarrow\ \Sigma = M\Sigma M^\top + \Delta
\tag{14.2}
$$

**解の存在と一意性。** 対称な $\Sigma$ の 21 成分に対する線形方程式である。行列の積を縦ベクトル化 ($\mathrm{vec}$、列を縦に並べる) すると $\mathrm{vec}(M\Sigma M^\top) = (M\otimes M)\,\mathrm{vec}\,\Sigma$ なので

$$
(1 - M\otimes M)\,\mathrm{vec}\,\Sigma = \mathrm{vec}\,\Delta
\tag{14.2a}
$$

$M\otimes M$ の固有値は $\lambda_i\lambda_j$ (全部の組)。したがって、$\lambda_i\lambda_j = 1$ となる組が無ければ解は一意。$M$ の固有値は $\lambda_k = e^{-a_k}e^{\pm i\mu_k}$ ($a_k > 0$) なので、$|\lambda_i\lambda_j| < 1$ で一意である。解は、$n$ 周前に入った擾乱が $M^n$ で減衰して残ったものの和でもある:

$$
\Sigma = \sum_{n=0}^{\infty}M^n\Delta\,(M^\top)^n
\tag{14.2b}
$$

($M\Sigma M^\top$ を (14.2b) に代入すると、$\Sigma - \Delta$ に一致することが、項をずらして確かめられる。級数は $M\otimes M$ のスペクトル半径 $\max|\lambda|^2 < 1$ で収束する。) 減衰が弱い ($a_k\sim10^{-4}$) ので、解は $\Delta/(2a)$ のように大きな因子を持ち、$1 - M\otimes M$ は 1 に近い量の差になって桁落ちしやすい。K64 が `btr` の対角を別の式で作る理由である (§4.4)。

### 3.2 $\Delta$ と $D$ の組み立て

励起は 1 周の途中の位置 $s_i$ で生じる ($\Delta_i$、13 章)。位置 $s_i$ での擾乱 $\xi_i$ は、リングの始点まで残りの行列 $T_{i\to\mathrm{end}}$ (要素 $i+1$ から最後までの行列の積) で運ばれて、始点の座標への擾乱 $T_{i\to\mathrm{end}}\xi_i$ になる。各位置の擾乱は独立なので、共分散は足し合わされる:

$$
\Delta = \sum_iT_{i\to\mathrm{end}}\,\Delta_i\,T_{i\to\mathrm{end}}^\top
\tag{14.3}
$$

漸化式で書くと、始点から終点に向けて要素を 1 つずつ通るとき、$B_0 = 0$、$B_{i} = T_iB_{i-1}T_i^\top + \Delta_i$ ($T_i$ は要素 $i$ の行列)。**(14.3) との一致** (帰納法): $B_1 = \Delta_1$、$B_2 = T_2\Delta_1T_2^\top + \Delta_2$、$B_3 = T_3T_2\Delta_1T_2^\top T_3^\top + T_3\Delta_2T_3^\top + \Delta_3$ と続き、$B_N$ の第 $i$ 項は $T_N\cdots T_{i+1}\Delta_iT_{i+1}^\top\cdots T_N^\top$ で (14.3) の第 $i$ 項に等しい。K64 の `tmulbs(beam, trans1, calint)` が $B \leftarrow TBT^\top$ を、`beam(1:21) = beam(1:21) + dbeam` が $+\Delta_i$ を実行する。

減衰側も同じ考え方で、要素 $i$ の放射による追加の写像の微分 $D_i$ (13 章 §6 の `tr1` の 7:12 列) が、始点から要素 $i$ までの行列 $T_{0\to i}$ の後ろに掛かり、残りの $T_{i\to\mathrm{end}}$ で運ばれる:

$$
D = \sum_iT_{i\to\mathrm{end}}\,D_i\,T_{0\to i}
\tag{14.3a}
$$

(pysad の `apply`: 要素の行列 $J$ を $T\leftarrow JT$、$D\leftarrow JD$、$B\leftarrow JBJ^\top$ と更新する。$D_i$ はその要素の中で足される。) **$D$ は RF の傾きに依存する。** 8 セルのリングで RF を先頭に置くと、すべての $T_{0\to i}$ が RF の行列 $K = 1 + \kappa\,e_\delta e_z^\top$ ($\delta \to \delta + \kappa z$、$e_\delta$、$e_z$ は $\delta$、$z$ 方向の単位ベクトル) を右から含み、$D = D'K$ ($D'$ は RF を含まない部分) となる。数値で、$D_{\delta z}/D_{\delta\delta} = 6.0293\times10^{-3}$ は $\kappa = 6.0293\times10^{-3}$ m$^{-1}$ に等しい ([emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (C))。一方 $B$ には $T_{i\to\mathrm{end}}$ だけが掛かり、RF (先頭) は含まれない。

**観測点の取り方と $\Delta$。** $\Delta$ は、ライン (1 周) の端の観測点での値である。RF をラインの先頭に置くと、観測点は RF の前 (リングの終わり) で、$\Delta_{\delta\delta}$ は 16 台の偏向磁石の励起の和に一致する ($\delta$ は偏向磁石、ドリフト、四極で変わらないので、$T_{\delta\delta}=1$、$T_{\delta j}=0$)。RF を末尾に置くと、観測点は RF のあとで、$B' = KBK^\top$ ($K$: RF の行列) となる。$\delta\delta$ 成分は

$$
B'_{\delta\delta} = B_{\delta\delta} + 2\kappa B_{z\delta} + \kappa^2B_{zz}
\tag{14.3b}
$$

数値 ([emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (D)): RF 先頭 $6.7418\times10^{-9}$、RF 末尾 $6.1109\times10^{-9}$ ((14.3b) の右辺も $6.1109\times10^{-9}$ で一致、比 $0.9044$)。これは計算の違いではなく、観測点の違い (正準変換 $\Sigma' = K\Sigma K^\top$、$M' = KMK^{-1}$、$\Delta' = K\Delta K^\top$。$K$ はシンプレクティックで、$M'\Sigma'M'^\top = KM\Sigma M^\top K^\top$ なので (14.2) が保たれる)。エミッタンスは観測点によらず、$\varepsilon_x = 5.2996\times10^{-5}$ (両方)、$\sigma_\delta$ は $1.23475\times10^{-3}$ と $1.23479\times10^{-3}$ ($\delta$ の定義点が RF の前後で $\kappa z$ だけ違う分)。

**確認** ([emit_01_equilibrium.out.txt](outputs/emit_01_equilibrium.out.txt) (c)): RF をラインの先頭に置いた例で、$\Delta = \Sigma - M\Sigma M^\top$ の $\delta\delta$ 成分は $6.742\times10^{-9}$、16 台の偏向磁石の $N_\gamma\langle u^2\rangle$ の和 (13 章。$\gamma$、$\rho$ は 13 章の例と同じ) は $6.757\times10^{-9}$ (比 $0.9978$)。$\Delta$ は対称で、固有値は $\ge -3\times10^{-25}$ (半正定値)。$\Delta_{\delta\delta}$ の `st.B` を直接見ても同じ値 ($6.7417\times10^{-9}$) なので、比 0.9978 は $\Sigma - M\Sigma M^\top$ の計算の誤差ではない。この $-2.2\times10^{-3}$ の差は、RFSW OFF・RADCOD OFF (閉軌道の $\delta$ が 0) でも同じ ($-2.213\times10^{-3}$) で、$E = 1, 2, 3$ GeV で $-8.1\times10^{-5}$、$-6.5\times10^{-4}$、$-2.2\times10^{-3}$ と $U_0/E_0$ ($2.9\times10^{-5}$、$2.3\times10^{-4}$、$7.8\times10^{-4}$) の約 $2.8$ 倍に比例する。閉軌道の $\delta$ ののこぎり歯の効果ではない。原因は未確認 (13 章の式は一定の $\gamma$ の古典的な式で、$U_0/E_0$ の 1 次の補正は含まない。`tradke` のスライス内の扱いに由来する可能性がある)。

## 4. `tecalc`: 正規座標と各モード

K64 の `tecalc` は、**$T$ (放射を除いた 1 周の行列) の固有ベクトルから正規座標 $R$ を作り、$D$ はその座標で 1 次の摂動として扱う** (`rx = trans(:,1:6)` を `teigenc` に渡し、`rd = ri*trans(:,7:12)*r`)。以下、この順に導出する。

### 4.1 正規座標: `teigenc`・`tnorm`・`tsymp`

**固有値と固有ベクトルの構造。** $T$ は実でシンプレクティック ($T^\top ST = S$、$S$ は各 (座標, 運動量) の 2×2 に $\begin{pmatrix}0&1\\-1&0\end{pmatrix}$ を置いた 6×6)。安定なら固有値は単位円上の 3 組の共役ペア $e^{\pm i\mu_k}$ で、固有ベクトル $v_k$ と $v_k^*$ を持つ。シンプレクティック性から、固有ベクトルの間の直交性が出る: $v_k^\dagger T^\top STv_l = \lambda_k^*\lambda_lv_k^\dagger Sv_l$ で、左辺は $v_k^\dagger Sv_l$。よって $(\lambda_k^*\lambda_l - 1)\,v_k^\dagger Sv_l = 0$、異なるモード ($\lambda_k^*\lambda_l\ne1$) では $v_k^\dagger Sv_l = 0$ (同様に $v_k^\top Sv_l = 0$、$\lambda_k\lambda_l\ne1$)。同じモードでは、$v_k^\dagger Sv_k$ ($S$ が実反対称なので純虚数) が残る。

**`teigenc`**: 実行列の固有値・固有ベクトルを求めるルーチン (`teigen.f`: 釣り合い `tbal` → Hessenberg 化 `thess` (高速 Givens 変換) → QR 法 `tqr`)。数値線形代数の標準の手順で、ここでは導出しない。結果の検査: $R^{-1}TR$ が 2×2 の回転ブロックになるか ([emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (E1): $x$ ブロックは $\begin{pmatrix}-0.74717&-0.66463\\0.66463&-0.74717\end{pmatrix}$、$\cos\mu_x = -0.74717$)。縮退した固有値 (RFSW OFF の $z$ の Jordan 形) では固有ベクトルが正しく出ない、と `teigen.f` の冒頭のコメントに書いてある。

**`tnorm`**: 固有ベクトルを、実の 2 本の列 $(R_X, R_P)$ で表し、規格化と向きを決める。モード $k$ の複素固有ベクトルを $r = (R_X - iR_P)/2$ とおくと ($x = 2\mathrm{Re}[r\,\zeta]$、$\zeta = X + iP$ とすると $x = R_XX + R_PP$)、実の列の対で $R^\top SR = S$ が成り立つ条件は、モード内の条件 $R_X^\top SR_P = 1$ (`tnorm` の `s` $= \sum_j(r_{j,i}r_{j+1,i+1} - r_{j,i+1}r_{j+1,i}) = R_X^\top SR_P$) と、モード間の直交性 (上の導出) である。$r^\dagger Sr$ を計算すると

$$
r^\dagger Sr = \tfrac14(R_X + iR_P)^\top S(R_X - iR_P) = \tfrac14\left(-iR_X^\top SR_P + iR_P^\top SR_X\right) = -\tfrac{i}{2}\quad (R_P^\top SR_X = -R_X^\top SR_P)
$$

したがって規格化は $i\,r^\dagger Sr = 1/2$ (> 0)。`tnorm` の手順:

1. $s = R_X^\top SR_P$ を計算する。$s>0$: $\sqrt s$ で両方の列を割る ($s\to1$)。$s<0$: 固有値を共役に取り替え (`ceig(i) = conjg(ceig(i))`)、$R_X$ を $\sqrt{|s|}$ で、$R_P$ を $-\sqrt{|s|}$ で割る (符号つきの割り算で $s\to+1$。この符号は $i\,v^\dagger Sv$ の符号、いわゆる Krein 符号が負の固有値を、共役の側に取り替えて正にすることに相当する)。$s = 0$ は不安定 (警告)。
2. モードの並び替え: $z$ モードは、成分 $(z, \delta)$ の「面積」$(R_{51}R_{62} - R_{52}R_{61})^2$ が最大のモード、残りで $(x,p_x)$ の面積 $(R_{11}R_{22} - R_{12}R_{21})^2$ が大きい方を $x$ モードとする。
3. 位相の自由度を固定する。固有ベクトルは位相 $e^{i\varphi}$ の不定性を持つ。実の列では $(R_X, R_P)$ の回転 $R_X' = R_X\cos\theta + R_P\sin\theta$、$R_P' = -R_X\sin\theta + R_P\cos\theta$ に対応し、$R^\top SR = S$ も $R^{-1}TR$ (回転 $P(\mu)$ と可換) も保つ。`tnorm` は $\tan\theta = R_{i,i+1}/R_{i,i}$ ($i$ はそのモードの座標の行) として、$R'_{i,i+1} = -R_{ii}\sin\theta + R_{i,i+1}\cos\theta = 0$ とする。実際に、$R$ の $x$ モードで $R_{1,2} = 5.1\times10^{-14}$ (ゼロ)。このとき $x = R_{11}X$、すなわち $X = x/\sqrt{\beta_x}$ ($R_{11} = \sqrt{\beta_x}$) の規格化座標になる。
4. 実固有値の対 (RFSW OFF の $z$、不安定面) には別の枝があり、列の比を $\sqrt[4]{s_a/s_b}$ ($s_a$、$s_b$ は 2 列の長さの 2 乗) で揃えて $R_{i,i+1}$ を消す。この枝の導出は未確認 (§9)。

**`tsymp`**: 固有ベクトルの数値誤差のため、$R$ は厳密にシンプレクティックではない。$t^s = -Jt^\top J$ (シンプレクティックな逆行列。$t$ がシンプレクティックなら $t^{-1}$ に等しい) として、$E = t^st - 1$ (シンプレクティックなら 0) とおき、

$$
t' = \tfrac32t - \tfrac12t\,t^st = t\left(1 - \tfrac12E\right)
\tag{14.4}
$$

とする (`tsymp`: `rx = t * tinv6(t)`、`-0.5*(rx*t) + 1.5*t`)。誤差 $E' = t'^st' - 1$ を求める。$(AB)^s = B^sA^s$ と、$E^s = (t^st)^s - 1 = t^st - 1 = E$ ($t^s t$ は自己随伴) を使うと ($(1 - E/2)^s = 1 - E/2$)

$$
t'^st' = \left(1 - \tfrac E2\right)t^st\left(1 - \tfrac E2\right) = \left(1 - \tfrac E2\right)^2(1 + E) = 1 - \tfrac34E^2 + \tfrac14E^3,\qquad E' = -\tfrac34E^2 + \tfrac14E^3
$$

($E$ 同士は可換。sympy で $(1 - E/2)^2(1 + E)$ の展開、$(t^st)^s = t^st$ (4×4 の記号行列)、$(AB)^s = B^sA^s$ を確認した: [emit_03_derivations.out.txt](outputs/emit_03_derivations.out.txt) (a))。**2 次収束** (Newton 法の形): ランダムな摂動を加えた $R$ で、$|E| = 9.6\times10^{-2}\to1.70\times10^{-3}$ (予想 $\tfrac34|E^2| = 1.69\times10^{-3}$)、$1.04\times10^{-2}\to6.3\times10^{-5}$ (予想 $6.3\times10^{-5}$)。`tecalc` の $R$ 自体は $|R^\top SR - S| = 3.3\times10^{-16}$。

### 4.2 減衰率: 固有値の 1 次の摂動

$M = T + D$ ($D$ は小さい) の固有値 $\lambda + \delta\lambda$ を、$T$ の固有値 $\lambda$、右固有ベクトル $r$ ($Tr = \lambda r$)、左固有ベクトル $\ell$ ($\ell T = \lambda\ell$、$\ell r = 1$) から 1 次で求める。$(T + D)(r + \delta r) = (\lambda + \delta\lambda)(r + \delta r)$ の 1 次の項は $T\delta r + Dr = \lambda\delta r + \delta\lambda\,r$。左から $\ell$ を掛け、$\ell(T - \lambda) = 0$ を使うと

$$
\delta\lambda = \ell Dr
\tag{14.5}
$$

$r = (R_X - iR_P)/2$ に対して $\ell = \ell_X + i\ell_P$ ($\ell_X$、$\ell_P$ は $R^{-1}$ の $X$、$P$ の行。$\ell r = \tfrac12(\ell_XR_X + \ell_PR_P + i(\ell_PR_X - \ell_XR_P)) = 1$) である。$\mathtt{rd} = R^{-1}DR$ の 2×2 ブロック (成分 $\mathtt{rd}_{XX} = \ell_XDR_X$ など) で書くと

$$
\ell Dr = \tfrac12\left[(\mathtt{rd}_{XX} + \mathtt{rd}_{PP}) + i(\mathtt{rd}_{PX} - \mathtt{rd}_{XP})\right]
$$

$\zeta = X + iP$ の 1 周の係数は $\lambda = e^{-i\mu}$ ($R^{-1}TR$ が $\begin{pmatrix}\cos\mu&\sin\mu\\-\sin\mu&\cos\mu\end{pmatrix}$ のとき、$\zeta' = (\cos\mu - i\sin\mu)\zeta$)。減衰率は $\delta\lambda/\lambda = -a + i\,\delta\mu$ の実部 (固有値の変化の比、絶対値の変化) である。K64 の式は、複素共役側 ($\lambda^* = e^{i\mu} = \mathtt{cc}$) で書いたもの:

$$
\mathtt{cd}(k) = \frac{(\mathtt{rd}_{ii} + \mathtt{rd}_{i+1,i+1})/2 + \mathrm{i}\,(\mathtt{rd}_{i,i+1} - \mathtt{rd}_{i+1,i})/2}{\mathtt{cc}(i)} = \overline{\left(\frac{\delta\lambda}{\lambda}\right)},\quad a_k = -\mathrm{Re}\,\mathtt{cd}(k)
\tag{14.6}
$$

(`cc` は $e^{i\mu_k}$。符号は SAD の約束で、`params(ipdampx) = dble(cd(1))` は負。虚部 `params(ipdnux)` は摂動によるチューンのずれ $\delta\mu$ の符号反転。) 減衰時間 $\tau_k = 2\pi/(\omega_0a_k)$ (`taurdx` など、$\omega_0$ は周回の角周波数)。**確認** ([emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (B)): $T$ の複素固有ベクトル ($i\,r^\dagger Sr = 1/2$) を直接作って $a_k = -\mathrm{Re}(\ell Dr/\lambda)$ を計算すると、$(5.232652\times10^{-5},\ 3.907412\times10^{-4},\ 1.118899\times10^{-3})$ で `tecalc` と 7 桁一致 (チューンも)。

**分配数と和の関係。** $J_k = 4a_k/(a_x + a_y + a_z)$ (`sr = 4/sum`、`params(ipjx:ipjz)`)。$\sum J_k = 4$ は定義により自明なので、検証にはならない。非自明な関係は、$\sum a_k = 2U_0/E_0$ (放射積分の式 $J_x = 1 - I_4/I_2$、$J_y = 1$、$J_z = 2 + I_4/I_2$ の和が 4 になることの元) である。6 個の固有値の対数の和は $\sum_{6}\log\lambda = \log\det(T + D)$ で、右辺は $\det T = 1$ から $\log\det(1 + T^{-1}D)\approx\mathrm{Tr}(T^{-1}D)$ (1 次)。左辺は共役ペアで虚部が消えて $-2\sum_ka_k$ なので

$$
\sum_ka_k = -\tfrac12\mathrm{Tr}(T^{-1}D) = -\tfrac12\sum_i\mathrm{Tr}\,D_i
\tag{14.6a}
$$

(最後の等号は (14.3a) と $\mathrm{Tr}(T_{0\to i}^{-1}T_{i\to\mathrm{end}}^{-1}\cdots)$ の巡回性による。$T^{-1}D = \sum_iT_{0\to i}^{-1}D_iT_{0\to i}$ の各項の跡が $\mathrm{Tr}\,D_i$ になる。) 各スライスの $D_i$ の跡は、$p_x$ と $p_y$ の縮み ($-u_i$ ずつ) と、損失の $\delta$ 依存 ($u\propto(1+\delta)^2$ による $-2u_i$) で $-4u_i$、したがって $\sum_ka_k = 2U_0/E_0$ が 1 次の理論の値。**数値**: $-\tfrac12\mathrm{Tr}(T^{-1}D) = 1.561966\times10^{-3}$ は `tecalc` の $\sum a_k$ と一致 (比 1.0000)、$2U_0/E_0 = 1.563266\times10^{-3}$ との比は $0.99917$。この差は $U_0/E_0$ に比例し (3 GeV で $-8.3\times10^{-4}$、1 GeV で $-3\times10^{-5}$、0.5 GeV で $-1\times10^{-5}$: 係数は約 $-1.06\,U_0/E_0$)、**$\nu_s$ には依らない** (3 GeV で加速電圧を 3〜24 MV に変えても 0.99917 で一定)。$U_0/E_0$ の 1 次の補正 (放射の 2 次の項) と見られるが、その導出は未確認。

**$J_y$ のずれ。** 垂直方向は結合がなければ $y$ モードの固有ベクトルが $(y, p_y)$ だけで、$D$ のうち効くのは $D_{p_yp_y} = -U_0/E_0$ だけ。したがって $a_y = U_0/(2E_0)$ で、実際 $a_y/(U_0/2E_0) = 0.99981$。$J_y = 4a_y/\sum a_k = 0.99981/0.99917 = 1.00064$ で、`tecalc` の $J_y = 1.00064$ に一致する。**$J_y - 1$ の原因は、上の $\sum a_k$ の $U_0/E_0$ の補正である。**

### 4.3 各モードのエミッタンス

固有ベクトルの規格化 ($i\,r^\dagger Sr = 1/2$) のもとで、$\zeta_k = X_k + iP_k = \ell_kX$ ($\ell_k$ は左固有ベクトル) の 1 周の更新は $\zeta_k' = \lambda_k\zeta_k + \ell_k\xi$ (他のモードとの混合を無視する)。両辺の絶対値の 2 乗の期待値をとると、$\langle\zeta_k\,\ell_k\xi\rangle = 0$ なので

$$
\langle|\zeta_k|^2\rangle = e^{-2a_k}\langle|\zeta_k|^2\rangle + \ell_k\Delta\ell_k^\dagger
$$

$|\zeta_k|^2 = X_k^2 + P_k^2$、$\varepsilon_k = \langle|\zeta_k|^2\rangle/2$ (円なら $X_k^2 = P_k^2 = \varepsilon_k$)、$\bar\Delta_k = \ell_k\Delta\ell_k^\dagger/2 = (\Delta_{n,XX} + \Delta_{n,PP})/2$ なので

$$
\varepsilon_k = \frac{\bar\Delta_k}{1 - e^{-2a_k}}\simeq\frac{\bar\Delta_k}{2a_k},\qquad \bar\Delta_k\equiv\frac{\Delta_{n,XX} + \Delta_{n,PP}}{2} = \tfrac12\ell_k\Delta\ell_k^\dagger
\tag{14.7}
$$

(正規座標の $\Delta_n = R^{-1}\Delta R^{-\top}$。ブロックの 2 次モーメントの方程式 $\Sigma_k = e^{-2a_k}P\Sigma_kP^\top + \Delta_k$ のトレースを取った式と同じ: 回転 $P$ は $\mathrm{tr}$ を変えない。) 固有ベクトルが $T$ のもので、$D$ によるモード間の混合を無視している点が近似である。

**確認** ((d)(e)、[emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (B)): $\Sigma_n$ の対角は $(5.2997, 5.2995, 0, 0, 7.8982, 7.8982)\times10^{-5}$ ($X$ と $P_X$ が等しい: 円)。(14.7) の予想は $\varepsilon_x = 5.2981\times10^{-5}$、$\varepsilon_z = 7.8984\times10^{-5}$ で、$\Sigma_n$ の $5.2996\times10^{-5}$、$7.8982\times10^{-5}$ と $0.03\%$ 以内 (近似 $\bar\Delta/2a$ は $5.2978\times10^{-5}$、$7.890\times10^{-5}$)。$x$ の $0.03\%$ の差の原因は、$D$ のモード間の混合 (次節) で、§4.4 の数値で確かめた。

### 4.4 21 成分の連立方程式 `btr`

**$\Sigma_n$ の方程式。** $\Lambda' \equiv R^{-1}(T + D)R = P + \mathtt{rd}$ ($P = \mathrm{diag}(P(\mu_1), P(\mu_2), P(\mu_3))$、各 $P(\mu) = \begin{pmatrix}\cos\mu&\sin\mu\\-\sin\mu&\cos\mu\end{pmatrix}$、$\mathtt{rd} = R^{-1}DR$。$R^{-1}TR = P$ による) と書くと、正規座標の (14.2) は

$$
\Sigma_n = \Lambda'\Sigma_n\Lambda'^\top + \Delta_n,\qquad \mathrm{vec}\,\Sigma_n = (1 - \Lambda'\otimes\Lambda')^{-1}\mathrm{vec}\,\Delta_n
\tag{14.8}
$$

K64 の `btr(k,l)` は、行 $k = (i,j)$ ($i\ge j$)、列 $l = (m,n)$ について `btr(k,l) = -(rdb(i,m)+rdb(i,m+6))*(rdb(j,n)+rdb(j,n+6))` ($\mathtt{rdb}(\cdot,1\!:\!6) = \mathtt{rd}$、$\mathtt{rdb}(\cdot,7\!:\!12) = P$ なので、$\mathtt{rdb}(i,m) + \mathtt{rdb}(i,m+6) = \Lambda'_{im}$) を作り、対角に 1 を足す。つまり `btr` は (14.8) の $(1 - \Lambda'\otimes\Lambda')$ を、対称な 21 成分に縮約したものである ($\Sigma_{mn} = \Sigma_{nm}$ と $\Lambda'_{im}\Lambda'_{jn}$ の足し合わせの重複は、`iaidx` の添字で畳み込まれる)。**検証**: $\Lambda' = R^{-1}(T + D)R$ から作った 36×36 の (14.8) を直接解くと、`tecalc` の結果と $\varepsilon_x$、$\varepsilon_z$ とも相対 $10^{-5}$ 以内で一致する ([emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (E2): 比 1.00000、1.00000)。

**モード $k$ のブロックの対角の置き換え。** $x$ の 2×2 ブロックの 3 成分 ($\Sigma_{XX}$、$\Sigma_{PP}$、$\Sigma_{XP}$) を、和 $u = (\Sigma_{XX} + \Sigma_{PP})/\sqrt2$、差 $v = (\Sigma_{PP} - \Sigma_{XX})/\sqrt2$、交差 $w = \sqrt2\Sigma_{XP}$ に取り替える (`tecalc` の $\sqrt{1/2}$ を掛ける操作)。ブロックの変換が $\Lambda_k = (1 - a)P(\mu)$ のとき、$\Sigma\to\Lambda_k\Sigma\Lambda_k^\top$ は

$$
\begin{pmatrix}u'\\v'\\w'\end{pmatrix} = (1-a)^2\begin{pmatrix}1&0&0\\0&\cos2\mu&-\sin2\mu\\0&\sin2\mu&\cos2\mu\end{pmatrix}\begin{pmatrix}u\\v\\w\end{pmatrix}
\tag{14.8a}
$$

(sympy で確認: [emit_03_derivations.out.txt](outputs/emit_03_derivations.out.txt) (d)。固有値は $(1-a)^2$ と $(1-a)^2e^{\pm2i\mu}$。) **和 (トレース) は他の 2 成分と分離し、減衰率 $(1-a)^2$ だけを持つ。差と交差は $2\mu$ で回転して平均が 0 になる** ($\mu\ne0,\pi$ のとき、2 次モーメントの楕円が回転で円に均される)。したがって (14.8) を解くと、和の成分に対する対角の係数が $1 - (1-a)^2 = 2a - a^2$ で (14.7) の $1 - e^{-2a}$ に一致する (差は $a^2$ のオーダー: $a_x = 5.2\times10^{-5}$ で相対 $2.6\times10^{-5}$)。

$\mathtt{rd}$ に、$-aP$ ではない対称な部分 $S = \begin{pmatrix}s_1&s_2\\s_2&-s_1\end{pmatrix}$ (トレースなし) があるとき ($\mathtt{rd} = -aP + S$)、`tecalc` が `btr(k,k)` に代入する式 (Fortran: $-(\mathtt{rd}_{ii}^2 + \mathtt{rd}_{i,i+1}^2 + \mathtt{rd}_{i+1,i}^2 + \mathtt{rd}_{i+1,i+1}^2)/2 - \cos\mu\,(\mathtt{rd}_{ii} + \mathtt{rd}_{i+1,i+1}) - \sin\mu\,(\mathtt{rd}_{i,i+1} - \mathtt{rd}_{i+1,i})$) は、sympy で展開すると

$$
\mathtt{btr}(k,k) = 2a - a^2 - s_1^2 - s_2^2
\tag{14.8b}
$$

で、$\Lambda' = P + \mathtt{rd}$ の (和,和) の係数を $1$ から引いたものに**厳密に一致する** (差 0: 同じ出力 (b))。この置き換えの理由: 一般の対角 `1 - (Λ'⊗Λ')` は、$1$ から $1 - 10^{-4}$ を引く形の桁落ち (相対 $10^{-12}$ まで) を持つので、$\mathtt{rd}$ の成分から直接 $2a$ 程度の量を作る式に替える。`ab` による行・列のスケーリング ($\sqrt{|\mathtt{btr}(k_1,k_1)/\mathtt{btr}(k_2,k_2)|}$) は、和と差の行の大きさを揃えて連立方程式を解く (`tsolva`、特異値分解で $10^{-8}$ 以下の特異値を捨てる) ときの条件を良くするためと読める (数値の工夫で、物理には影響しない。未確認)。

**モード間の相関の見積り。** 異なるモード $k\ne l$ の相関 $\Sigma_{n,kl}$ は、(14.8) を複素座標 $\zeta_k = X_k + iP_k$ で書くと、励起 $\Delta_{kl}$ の項だけでなく、$\mathtt{rd}_{kl}$ (放射によるモード間の混合) の項からも生じる: $\langle\zeta_k\zeta_l^*\rangle = (\Delta_{kl} + \text{混合項})/(1 - \lambda_k\lambda_l^*)$。分母は共鳴 $\mu_k\mp\mu_l\approx2\pi\cdot$整数 の近くで小さくなる。**数値** ((E2)): $x$–$z$ の $|1 - \lambda_x\lambda_z^*| = 1.73$、$|1 - \lambda_x\lambda_z| = 1.96$ (共鳴から遠い)。励起だけの見積り $|\Delta_{xz}|/|1 - \lambda\lambda^*| = 5.6\times10^{-9}$。$\mathtt{rd}$ のモード間ブロック ($|\mathtt{rd}_{zx}|\le6.6\times10^{-4}$、$|\mathtt{rd}_{xz}|\le7.1\times10^{-5}$、対角ブロックの減衰 $a_x = 5.2\times10^{-5}$ より大きい) を 0 にして (14.8) を解くと、$x$–$z$ の相関の最大は $5.9\times10^{-9}$ (見積りに一致) で、ブロックを残した全体の解は $2.45\times10^{-8}$ (`tecalc` の $2.45\times10^{-8}$)。**主な源は $D$ のモード間の混合** である。このブロックを 0 にしたときの $\varepsilon_x = 5.2977\times10^{-5}$ は (14.7) の $5.2981\times10^{-5}$ とほぼ一致し、全体の解 $5.2996\times10^{-5}$ との $0.03\%$ の差は混合による。

## 5. 同期運動と RFSW、MINCOUP

**RFSW ON** (RF あり): 縦モード ($z$、$\delta$) はシンクロトロン振動で回転し、$\mu_z \ne 0$。**$\alpha_p$ の式の導出。** 縦の 2×2 で、RF が $\delta\to\delta + \kappa z$ ($\kappa = \mathtt{dvcacc}/p_{gev}$、RF の電圧の $z$ 微分を $p_0c$ で割ったもの) の薄いレンズ、リングの残りが $z\to z + R_{56}^{\rm eff}\delta$ ($R_{56}^{\rm eff} = R_{56} + \sum_jR_{5j}\eta_j$、分散で経路が変わる分を含む) のとき、1 周の行列は $\begin{pmatrix}1&R_{56}^{\rm eff}\\\kappa&1 + \kappa R_{56}^{\rm eff}\end{pmatrix}$。行列式は 1、跡は $2 + \kappa R_{56}^{\rm eff} = 2\cos\mu_s$ なので

$$
\kappa R_{56}^{\rm eff} = 2\cos\mu_s - 2 = -4\sin^2(\mu_s/2) \simeq -\mu_s^2\ (\text{級数: } -\mu_s^2 + \mu_s^4/12 - \cdots)
\tag{14.9}
$$

(sympy で確認: [emit_03_derivations.out.txt](outputs/emit_03_derivations.out.txt) (c)。) $R_{56}^{\rm eff} = -\alpha_pC$ ($C = c\cdot2\pi/\omega_0$ は周長。$z$ の向きにより符号は変わる) とおくと、`tecalc` (`synchm`) の式:

$$
\alpha_p = -\frac{\mathrm{Im}\,\mathtt{cd}(6)\,|\mathrm{Im}\,\mathtt{cd}(6)|}{(c\,2\pi/\omega_0)\,(dv_{cacc}/p_{gev})},\qquad \omega_z = |\mathrm{Im}\,\mathtt{cd}(6)|\,\omega_0/2\pi
\tag{14.10}
$$

(`dvcacc` は RF の電圧の $z$ 微分。$\mathrm{Im}\,\mathtt{cd}(6) = \mu_z$。$\mu|\mu|$ は $\mu^2$ (符号つき) で、(14.9) の $4\sin^2(\mu/2)$ の小さい $\mu$ の近似に当たる。v2 は $\mu|\mu|$ を $4\sin^2(\mu/2)$ (符号つき) に替えて (14.9) と一致させる。K64 の式は $\mu_z^2/12$ (相対) の分だけ $\alpha_p$ が大きく出る。) 標準の式 $\nu_s = \sqrt{h\alpha_pV\cos\varphi_s/(2\pi E)}$ ($\sin\varphi_s = U_0/V$) は (14.9) の $\mu_s^2$ に $\kappa = hV\cos\varphi_s\cdot2\pi/(CE)$ ($= 2\pi h/C\cdot V\cos\varphi_s/E$) を入れたもの。比較: $I_1/C$ から $\alpha_p = 0.39857$ ($\mathtt{calc6d}$ は $0.39780$、0.2% 低い)、$\nu_s = 0.0503$ (calc6d $0.0505$)。0.2% の差は、$x$–$z$ の結合 (§6) が縦の固有値を動かすことによる。縦の 2×2 から決める $\kappa R_{56}^{\rm eff} = -0.099969$ と、実際の (結合を含む) $\mu_z$ から決まる $-4\sin^2(\mu_s/2) = -0.099776$ の比は $1.00194$ で、$\alpha_p$ の比 $0.99807$ はその逆数 ([emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (D))。

**バケット高さ `bh` の導出。** RF の 1 周の kick と位相の進みは、$\varphi$ (同期位相からのずれ) と $\delta$ で

$$
\frac{d\delta}{dn} = \frac VE\left[\sin(\varphi_s + \varphi) - \sin\varphi_s\right],\qquad\frac{d\varphi}{dn} = 2\pi h\alpha_p\delta
$$

($n$ は周回数、$E = p_{gev}$)。ハミルトニアン $H(\varphi,\delta) = \pi h\alpha_p\delta^2 + \frac VE[\cos(\varphi_s + \varphi) + \varphi\sin\varphi_s]$ ($\partial H/\partial\delta = 2\pi h\alpha_p\delta$、$-\partial H/\partial\varphi = \frac VE[\sin(\varphi_s + \varphi) - \sin\varphi_s]$) が保存する。不安定固定点は $\varphi_u = \pi - 2\varphi_s$ ($\sin(\varphi_s + \varphi_u) = \sin(\pi - \varphi_s) = \sin\varphi_s$)、$H_u = \frac VE[-\cos\varphi_s + (\pi - 2\varphi_s)\sin\varphi_s]$ で、これを通るセパラトリクスの $\varphi = 0$ (安定固定点) での $\delta$ が最大の $\delta_{\max}$:

$$
\pi h\alpha_p\delta_{\max}^2 + \frac VE\cos\varphi_s = H_u\ \Rightarrow\ \delta_{\max}^2 = \frac{V}{\pi h\alpha_pE}\left[-2\cos\varphi_s + (\pi - 2\varphi_s)\sin\varphi_s\right]
$$

右辺は ($\alpha_p$ の符号を逆にする向きでは) $\dfrac{V}{\pi|\alpha_p|h\,p_{gev}}\left[2\cos\varphi_s - (\pi - 2\varphi_s)\dfrac{U_0p_{gev}}{V}\right]$ ($\sin\varphi_s = U_0p_{gev}/V$、`tecalc` の `bh`)。数値: $\varphi_s = 0.8973$ rad、$\alpha_p = I_1/C$ で $1.5583\times10^{-3}$、calc6d の $\alpha_p$ の式で $1.5598\times10^{-3}$、`tecalc` の `bh` は $1.5598\times10^{-3}$ (一致)。

**RFSW OFF**: 縦運動は閉じない。`rx(6,1:5) = 0`、`rx(6,6) = 1` として $\delta$ の変化を止める。`synchm` が偽のとき、$\alpha_p = -R_{56}\cdot\omega_0\,h_0/(2\pi c\,p_0)$ (`rx(5,6)` から。$h_0/p_0 = 1/\beta_0$ が、$z$ と経路長の差の換算)、`omegaz` は RF の電圧から (`sqrt(abs(alphap*...))`)。縦のエミッタンスは `beam(iaidx(5,1:6)) = 0` として扱わない (同期運動がないので $\varepsilon_z$ は定義されない)。

**MINCOUP**: 垂直のエミッタンスは、完全に結合のない理想リングでは 0 (本例: $\varepsilon_y = 0$) になり、IBS や空間電荷に使うと非現実的になる。ヘルプ (equilibrium-beam-envelope): $\varepsilon_{x,y} = \max[\varepsilon^0_{x,y},\ \mathtt{MINCOUP}\times(\varepsilon^0_x + \varepsilon^0_y)]$ ($\varepsilon^0$ はラティスだけが与えるエミッタンス)。この下限は、IBS (`INTRA`) と WSPAC の反復 (`tintraconv`) の中だけで使われる (読み。数値では確認していない)。IBS の拡散・摩擦の項は $\Sigma$ に依存するので、(14.2) は $\Sigma = M\Sigma M^\top + \Delta_{\rm rad} + \Delta_{\rm IBS}(\Sigma)$ という非線形の方程式になる。`tintraconv` は、前の反復の $\Sigma$ から $\Delta_{\rm IBS}$ を作り、(14.2) を線形の連立方程式として解き直す不動点の反復と読める (ソースの読み。収束条件は未確認)。IBS の式自体の導出は [17_collective.md](17_collective.md) §9。

## 6. 有限のシンクロトロンチューンの補正 (シンクロ・ベータ結合)

放射積分の標準の式 ($J_x = 1 - I_4/I_2$、$\varepsilon_x = C_q\gamma^2I_5/(J_xI_2)$) は、縦運動 (シンクロトロン振動) が水平方向と結合しない極限 $\nu_s\to0$ の式である。`tecalc` は 6×6 の行列を直接扱うので、$\nu_s$ が有限のとき、水平モードの固有ベクトルが $z$、$\delta$ の成分を持ち、減衰も励起も変わる。

### 6.1 現象

3 GeV で加速電圧だけを変える ($\nu_s\propto\sqrt{V\cos\varphi_s}$。[emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (A)):

| $V$ [MV] | $\nu_s$ | $\varepsilon_x$ (calc6d) / 標準の式 | $(\text{比}-1)/\nu_s^2$ | $J_x$ (calc6d) | $J_x$ (式) | $J_y$ | $\sum a/(2U_0/E_0)$ |
|---|---|---|---|---|---|---|---|
| 3 | 0.0505 | 1.02215 | 8.69 | 0.13400 | 0.13643 | 1.00064 | 0.99917 |
| 4 | 0.0666 | 1.03863 | 8.71 | 0.13224 | 0.13643 | 1.00064 | 0.99917 |
| 6 | 0.0873 | 1.06632 | 8.71 | 0.12959 | 0.13643 | 1.00064 | 0.99917 |
| 12 | 0.1288 | 1.13774 | 8.30 | 0.12499 | 0.13643 | 1.00064 | 0.99917 |
| 24 | 0.1871 | 1.17739 | 5.07 | 0.13931 | 0.13643 | 1.00064 | 0.99917 |

見方: $\varepsilon_x$ の標準の式からの比 $-1$ は、$\nu_s\le0.09$ では $\nu_s^2$ に比例する (列 4 が 8.7 で一定)。$\nu_s$ がさらに大きいと高次の項で比例からはずれる ($\nu_s = 0.19$ では $J_x$ の向きが逆転)。$J_y$ と $\sum a/(2U_0/E_0)$ は $V$ に依らない (これは $\nu_s$ ではなく $U_0/E_0$ の効果。§4.2)。

エネルギーと電圧を同時に変えた場合 (2 GeV 以下は $U_0/E_0$ が小さい): 1 GeV・0.3 MV ($\nu_s = 0.0349$) で比 $1.0106$ ((比$-1$)$/\nu_s^2 = 8.76$)、1 GeV・0.1 MV ($\nu_s = 0.0197$) で $1.0035$ (8.95)、0.5 GeV・0.05 MV ($\nu_s = 0.0202$) で $1.0040$ (9.77)。$J_y - 1$ と $\sum a_k$ の不足は $U_0/E_0 = 2.9\times10^{-5}$ に比例して小さくなる ($J_y = 1.00003$、$\sum a/(2U_0/E_0) = 0.99997$)。

### 6.2 理論: 1 周の行列と RF の kick

RF をラインの先頭に置くと、1 周の行列は $T = AK$。$A$ は残りの全要素の行列 (RF を含まない)、$K = 1 + \kappa\,e_\delta e_z^\top$ は RF の薄い kick ($\delta\to\delta + \kappa z$)。$A$ は $\kappa$ に依らず (差は $9\times10^{-10}$)、$\kappa$ は電圧に比例する。$\kappa = \kappa_0 = 6.0293\times10^{-3}$ m$^{-1}$ が基準 (3 MV)。次の 3 点が使える:

1. **縦運動との関係**: $\kappa$ と $\nu_s$ は (14.9) $\kappa R_{56}^{\rm eff} = -4\sin^2(\mu_s/2)$ で結ばれる。$R_{56}^{\rm eff} = -16.581$ m ($= -I_1$、$\alpha_pC$)、$\kappa_0R_{56}^{\rm eff} = -0.09997$、$-4\sin^2(\mu_s/2) = -0.09978$ (比 1.0019 = 結合による $\mu_s$ の違い)。したがって $\kappa\simeq\mu_s^2/(\alpha_pC) = (2\pi\nu_s)^2/I_1$。
2. **$D$ の依存**: $D = D'K$ (§3.2)。$B$ は $\kappa$ に依らない (RF が先頭で、観測点は RF の前)。
3. **$\kappa\to0$**: RF なし (縦運動なし) の行列は $T_0 = A$。この極限で標準の式に戻るか?

**$\kappa = 0$ の極限が標準の式を再現する。** $T = A$ ($\delta$ は保存し、$z$ は $R_{56}$ で進む) で、水平モードの $a_x$、$\bar\Delta_x$、$\varepsilon_x$ を (14.5)、(14.7) で作ると、標準の式 ($a_x = J_xU_0/(2E_0)$、$\bar\Delta_x = \varepsilon_xJ_xU_0/E_0$) に対して

| | $a_x$ の比 | $\bar\Delta_x$ の比 | $\varepsilon_x$ の比 | $J_x$ |
|---|---|---|---|---|
| $\kappa = 0$ | 1.00020 | 0.99988 | 0.99974 | 0.136455 (式 0.136429) |

(残る $2\times10^{-4}$ は $U_0/E_0 = 7.8\times10^{-4}$ のオーダーの補正。) したがって、標準の式からの 2.2% のずれは、すべて RF の kick $\kappa$ の効果である。

### 6.3 1 次の摂動

$\kappa$ を摂動 $V_1 = dT/d\kappa = A\,e_\delta e_z^\top$ とみなす。水平モードの固有値は $z$ の固有値 (縮退した $1$) と離れているので、非縮退の摂動論が使える ($T_0 = A$ の水平モードの右・左固有ベクトル $r_0$、$\ell_0$、固有値 $\lambda_0$、$\ell_0r_0 = 1$)。$(T_0 + \kappa V_1)(r_0 + \kappa r_1) = (\lambda_0 + \kappa\lambda_1)(r_0 + \kappa r_1)$ の $\kappa$ の 1 次から

$$
\lambda_1 = \ell_0V_1r_0,\qquad(T_0 - \lambda_0)\,r_1 = (\lambda_1 - V_1)\,r_0
\tag{14.11}
$$

($\ell_0r_1 = 0$ と選ぶ。$T_0 - \lambda_0$ は特異だが、右辺は $\ell_0$ と直交するので解がある。) 左固有ベクトル $\ell_1$ も同様。

**構造**: $V_1r_0 = A\,e_\delta\,r_{0,z}$ ($r_{0,z}$ は $r_0$ の $z$ 成分)。つまり RF が水平モードの $z$ 成分を $\delta$ の kick に変え、それが分散 ($A$ の $\delta$ の列 $(\eta, R_{56}, 1)$) を通して $x$、$p_x$、$z$ に入る。

- $r_0$ の $z$ 成分: $A$ の $z$ 行は $z' = z + c\cdot x + R_{56}\delta$ ($c = (R_{51}, R_{52})$)、$\delta = 0$ なので、固有値の式 $\lambda_0r_{0,z} = r_{0,z} + c\cdot r_{0,x}$ から $r_{0,z} = c\cdot r_{0,x}/(\lambda_0 - 1)$ (数値で式との比 $1.00000000$)。拡張 Twiss の $z$–$x$ 結合 (`zx` など) に当たる量。
- $r_1$ の $\delta$ 成分: $A$ の $\delta$ 行は $\delta' = \delta$ なので、(14.11) の $\delta$ 行は $(1 - \lambda_0)r_{1,\delta} = -\,r_{0,z}\cdot1$ ($V_1r_0$ の $\delta$ 成分 $= r_{0,z}$ で、$\lambda_1 r_{0,\delta} = 0$)。$r_{1,\delta} = r_{0,z}/(\lambda_0 - 1)$ (数値で比 $1.000000$)。厳密な固有ベクトル ($\kappa = 10^{-4}$) の $\delta$ 成分 $/\kappa$ は $0.06400 - 0.23728i$、1 次の式 $0.06397 - 0.23716i$。
- よって水平モードの $\delta$ 成分は $\kappa\,c\cdot r_{0,x}/(\lambda_0 - 1)^2$: 振幅 $\propto\kappa\propto\nu_s^2$。これが $\varepsilon_x$ の補正が $\nu_s^2$ に比例する理由 (固有ベクトルの $\delta$ 成分が $\kappa$ の 1 次で入り、減衰 (14.5) と励起 (14.7) がその 1 次で変わる)。

**減衰率の変化**: $a_x = -\mathrm{Re}(\ell Dr/\lambda)$ を $\kappa$ で微分する。$D = D'K$ なので $dD/d\kappa = D'e_\delta e_z^\top$ も効く:

$$
\frac{da_x}{d\kappa} = -\mathrm{Re}\left[\frac{\ell_1Dr_0 + \ell_0Dr_1 + \ell_0D'e_\delta e_z^\top r_0}{\lambda_0} - \frac{\ell_0Dr_0}{\lambda_0}\frac{\lambda_1}{\lambda_0}\right]
$$

(第 3 項: RF の kick が $D$ を通して働く分。) 数値: 1 次の摂動 $-1.7162\times10^{-4}$、有限差分 ($\kappa = 10^{-4}$) $-1.7154\times10^{-4}$ (一致)。$a_x$ の相対変化は $\kappa_0$ で $-1.94\%$。

### 6.4 結果と分解

$\kappa$ を $0$ から $4\kappa_0$ まで動かして、(14.5)、(14.7) で $\varepsilon_x$ を求めた (表。基準の標準の式との比):

| $\kappa/\kappa_0$ | $a_x$ の比 | $\bar\Delta_x$ の比 | $\varepsilon_x$ の比 | $(\varepsilon_x\text{比}-1)/\kappa$ [m] |
|---|---|---|---|---|
| 0 | 1.00020 | 0.99988 | 0.99974 | — |
| 0.25 | 0.99538 | 1.00052 | 1.00521 | 3.459 |
| 0.5 | 0.99064 | 1.00121 | 1.01073 | 3.558 |
| 1 | 0.98139 | 1.00279 | 1.02186 | 3.625 |
| 2 | 0.96401 | 1.00677 | 1.04441 | 3.683 |
| 4 | 0.93499 | 1.01871 | 1.08960 | 3.715 |

見方: 最右列が $\kappa$ によらずほぼ一定 (3.5〜3.7 m): $\varepsilon_x$ の補正は $\kappa$ の 1 次で、$\varepsilon_x/\varepsilon_{\rm std} - 1\simeq3.6\,\mathrm{m}\times\kappa$。$\kappa$ を (14.9) で $\nu_s$ に直すと $\kappa = (2\pi\nu_s)^2/|R_{56}^{\rm eff}|$ (小さい $\nu_s$) で、$(\varepsilon_x\text{比}-1)/\nu_s^2 = 3.62\,\mathrm{m}\times(2\pi)^2/16.58\,\mathrm{m} = 8.6$ (測定 8.7)。

**2.2% の内訳** ($\kappa_0$): 減衰率 $a_x$ が $-1.86\%$ (比 $0.98139$)、したがって $J_x$ が $0.13643\to0.13400$ (1.8%)、励起 $\bar\Delta_x$ が $+0.28\%$ (比 $1.00279$)。$\varepsilon_x = \bar\Delta_x/(2a_x)$ なので $1.00279/0.98139 = 1.0218$ ($1.0221$ との差は $\kappa_0$ の 2 次と $U_0/E_0$ のオーダーの項)。つまり 2.2% の大部分は減衰分配 ($J_x$) の変化で、励起の変化は 0.3%。

**まとめ**:
- $J_x$ の 1.8%、$\varepsilon_x$ の 2.2% の差、$\alpha_p$ の 0.2% の差は、$\nu_s = 0.05$ の有限のシンクロトロンチューンによる水平のシンクロ・ベータ結合の補正 ($\propto\nu_s^2$)。
- $J_y = 1.0006$ は別の原因 (§4.2: $\sum a_k$ の $U_0/E_0$ の補正)。
- 係数 (減衰 $-3.1$ m、励起 $+0.46$ m (いずれも $\kappa_0$ のとき)) を、リングの光学関数 ($\eta$、$\beta$、位相、$R_{5j}$) の積分の閉じた式で書く導出は、していない (§9)。上に書いた 1 次の摂動の式 (14.11) と、$r_{0,z}$、$r_{1,\delta}$ の式 (ラティスの行列から計算できる) は、数値の係数と 0.05% 以内で一致する。

## 7. 結合のあるリングの垂直エミッタンス

垂直のエミッタンスは、結合 (スキュー四極など) でモードが混ざると有限になる。(14.7) は結合の有無に依らない ($T$ の 3 組の固有ベクトルについての式) ので、$y$ モードの左固有ベクトル $\ell_y$ に、励起のあるモード ($x$、$z$) の成分が混ざった量で $\bar\Delta_y$ を作る:

$$
\bar\Delta_y = \tfrac12\ell_y'\Delta\ell_y'^\dagger,\qquad\varepsilon_y = \frac{\bar\Delta_y}{1 - e^{-2a_y}}
\tag{14.12}
$$

結合のない $\varepsilon_y$ の励起 ($\ell_yB\ell_y^\dagger = 0$。$B$ は $y$ 方向に成分を持たない) は 0 で、混合 $\ell_y' = \ell_y + \sum_mc_m\ell_m$ が $\bar\Delta_y$ を作る。**混合の 1 次の摂動**: 結合の摂動を $V_c = T - T_{\rm base}$ ($T_{\rm base}$ は結合のない 1 周の行列) とすると、(14.11) と同じ導出で、左固有ベクトルの 1 次の補正は

$$
\ell_y' = \ell_y + \sum_{m\ne y}\ell_m\,\frac{\ell_yV_cr_m}{\lambda_y - \lambda_m}
\tag{14.13}
$$

($m$ は $x$、$z$、$y$ の各モードとその共役の 5 本。) $\bar\Delta_y = \tfrac12\sum_{m,m'}c_mc_{m'}^*\,\ell_mB\ell_{m'}^\dagger$ (対角の項 $|c_m|^2\bar\Delta_m$ と、同じモードとその共役の間などの交差項) で、$c_m\propto k_s$ ($V_c\propto k_s$) なので、すべての項が $k_s^2$ に比例し、$\varepsilon_y\propto k_s^2$。

**数値** ([emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) (F)): 基準のリングの 1 台目の位置に、スキュー四極 (長さ 0.1 m、$K_1 = k_s$、$\pi/4$ の回転) を置く (基準は $k_s = 0$ の同じリング。0.1 m のドリフトが入る)。

| $k_s$ | $\varepsilon_x$ | $\varepsilon_y$ (calc6d) | $\varepsilon_y/\varepsilon_x$ | (14.13) の予想 / calc6d | $J_x$ | $J_y$ |
|---|---|---|---|---|---|---|
| 0.002 | $5.2543\times10^{-5}$ | $1.9687\times10^{-9}$ | $3.747\times10^{-5}$ | 0.9987 | 0.13490 | 1.00062 |
| 0.004 | $5.2547\times10^{-5}$ | $7.8742\times10^{-9}$ | $1.499\times10^{-4}$ | 0.9987 | 0.13490 | 1.00057 |
| 0.008 | $5.2561\times10^{-5}$ | $3.1490\times10^{-8}$ | $5.991\times10^{-4}$ | 0.9989 | 0.13487 | 1.00036 |

見方: $\varepsilon_y/\varepsilon_x$ は $k_s$ を 2 倍にすると 4.000、3.998 倍 ($k_s^2$ に比例)。1 次の摂動 (14.13) の予想は calc6d と 0.13% 以内で一致する。混合の内訳 ($k_s = 0.004$、$\bar\Delta_y = 6.1\times10^{-12}$): 対角の項 $|c_m|^2\bar\Delta_m$ の和は $3.8\times10^{-12}$ で、そのうち $z$ とその共役が $3.4\times10^{-12}$ (89%)、$x$ とその共役が $4.1\times10^{-13}$ (縦の励起 $\bar\Delta_z = 1.8\times10^{-7}$ が $\bar\Delta_x = 5.5\times10^{-9}$ より 30 倍大きいため、**縦モード $z$ との混合が主**)。残りの $2.4\times10^{-12}$ (39%) は交差項である。$J_x$ は結合で $0.1340\to0.1349$ に動く (ベースの 0.1 m のドリフトによる変化を含む)。

(なお、結合の 2 次以降 (混合の振幅 $c_m\sim5\times10^{-3}$) と、共鳴 $\mu_x\pm\mu_y$ の近く ($\lambda_y - \lambda_x\to0$) では、1 次の摂動は使えない。この点は (14.13) の分母から分かるが、数値では確認していない。)

## 8. 確認スクリプト

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `emit_01_equilibrium.py` | (a) `calc6d` のチューン・減衰・エミッタンス、(b) 放射積分の式との比較 ($U_0$、$J_x$、$\varepsilon_x$、$\sigma_\delta$、$\alpha_p$、$\nu_s$)、(c) $\Delta = \Sigma - M\Sigma M^\top$ の $\delta\delta$ 成分と 13 章の $N_\gamma\langle u^2\rangle$、半正定値性、(d) 正規座標の $\Sigma_n$ の対角、(e) (14.7) | [emit_01_equilibrium.out.txt](outputs/emit_01_equilibrium.out.txt) |
| `emit_02_syncbeta.py` | (A) 加速電圧を変えたときの $\varepsilon_x$・$J_x$・$\sum a$、(B) 複素固有ベクトルによるモードの量と `tecalc`、(C) $\kappa$ の掃引・1 次の摂動・RF なしの極限、(D) $\sum a$ と $2U_0/E_0$、$\Delta_{\delta\delta}$ の観測点、バケット高さ、(E) `tsymp` の収束・`tnorm`・`btr`・モード間の相関、(F) スキュー四極つきの $\varepsilon_y$ | [emit_02_syncbeta.out.txt](outputs/emit_02_syncbeta.out.txt) |
| `emit_03_derivations.py` | (sympy) `tsymp` の誤差 $E\to-\tfrac34E^2$、`btr` の対角、和・差・交差の 3×3 行列、縦 2×2 の $\kappa R_{56}$ | [emit_03_derivations.out.txt](outputs/emit_03_derivations.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/emit_01_equilibrium.py` (十数秒)、`emit_02_syncbeta.py` (約 15 秒、`calc6d` を十数回)、`emit_03_derivations.py` (数秒)。

## 9. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV、8 セルの FODO リング (QF = 0.25、QD = −0.30 [1/m]、各 0.3 m、偏向 16 台 × 1.2 m の扇形、周長 41.6 m、$\eta_x = 3.08$ m)、CAVI $V$ = 3 MV、harmonic 64、$\varphi_s = 0.897$ rad、$\nu_s = 0.0505$。偏向の端面なし、6 極なし、結合なし (§7 を除く)、IBS・空間電荷なし。チューンは、RF を含む `calc6d` で $\nu = (1.6157,\ 1.1208)$ (`twiss` の位相)。RF なしの (縦運動のない) 行列では $(1.6165,\ 1.1209)$ で、$x$ の $0.0008$ の差は RF によるチューンのずれ ($\delta\mu_x = 0.814\,\mathrm{m}^{-1}\times\kappa$ の 1 次の摂動 = 0.0049 rad = $0.0008$ 周) である。

**未確認・未導出**:
- `teigenc` の内部 (釣り合い・Hessenberg 化・QR 法の詳細) は、標準の数値線形代数なので導出していない。結果の検査 (固有ベクトルが $T$ を回転ブロックにする) のみ。
- `tnorm` の実固有値の対 (RFSW OFF の $z$ の Jordan 形、不安定面) の枝の導出。Krein 符号と並び替えは §4.1 で説明したが、縮退の処理は未確認。
- `btr` の行・列のスケーリング `ab` の目的は読みで、効果の定量的な確認はしていない。
- $\sum a_k$ と $2U_0/E_0$ の $-1.06\,U_0/E_0$ の差、$\Delta_{\delta\delta}$ と $N_\gamma\langle u^2\rangle$ の $-2.8\,U_0/E_0$ の差の原因 (§3.2、§4.2)。
- 有限の $\nu_s$ の補正の係数 (減衰 $-3.1$ m、励起 $+0.46$ m) を、$\eta$・$\beta$・位相の積分の閉じた式にすること (§6)。固有ベクトルの構造と数値の係数は確認した。$\nu_s$ が大きい ($\ge0.13$) 領域の非線形な依存 (表の $V$ = 12、24 MV) も未導出。
- 結合のあるリングの導出は 1 次の摂動までで、共鳴の近く・強い結合は未確認 (§7)。$\bar\Delta_y$ の交差項を落とした近似 ($\sum|c|^2\bar\Delta$) の精度は評価していない (数値は全体 (14.13) の計算で 0.13%)。
- IBS (`tintraconv`) の反復の収束、空間電荷 (`WSPAC`) との組み合わせ。MINCOUP の下限の数値の効果。
- ソレノイドのある場合の `tmulbs` のうち、$bzh$ の項。
- ヘルプの拡張 Twiss ($V = PBR_6H$) との対応は [12_linear_optics.md](12_linear_optics.md) §8 のとおり未確認 (§6.3 の $r_{0,z}$ は、その $z$–$x$ の量に当たる)。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| equilibrium-beam-envelope | **原文で確認** (TeX 版・平文版とも): $\langle X_iX_j\rangle = M\langle X_iX_j\rangle M^T + \Delta_{ij}$、MINCOUP の式。「$M$ the one-turn transfer matrix including the damping」。(14.2) と同じ。 |
| ヘルプに無いもの | (1) (14.1)〜(14.3) の導出、(2) 正規座標での解法 (14.7)・(14.8)、(3) 減衰率と分配数の定義 (14.6)、(4) 同期運動の式 (14.10)。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

特になし (この章で確認した範囲)。分配数の差 (§6) は放射積分の式との定義の違い (有限の $\nu_s$) で、疑いではない。観察として、(i) $\sum a_k$ が $2U_0/E_0$ に対して $U_0/E_0$ のオーダーで足りない、(ii) $\Delta_{\delta\delta}$ が 13 章の式より $U_0/E_0$ のオーダーで小さい、の 2 点があるが、K64 の誤りと判断する根拠はない (原因は未確認)。

## 付録 C. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `trans(6,12)` | 1:6 列が $T$ (または要素の行列)、7:12 列が放射の微分 $D$ | `temit.f`、`tradke` | |
| `beam(42)` | 21 成分の対称行列 ($\Delta$ の走る和) × 2 | `temit.f` | `st.B` |
| `iaidx(i,j)` | 対称行列の下三角の添字 | `temw` | |
| `r`、`ri`、`rx` | 正規座標の行列、逆行列、`trans(:,1:6)` | `tecalc` | `normal_r`、`normal_ri` |
| `ceig`・`dceig` | 固有値 | `tecalc` | `cd` |
| `cd` | 減衰 + $\mathrm{i}\,\Delta\nu$ (1:3)、$\log\lambda$ (4:6) | `tecalc` | |
| `cc` | $e^{\mathrm{i}\mu}$ | `tecalc` | |
| `rd` | `ri*trans(:,7:12)*r` (減衰の行列) | `tecalc` | |
| `rdb`・`btr` | `rdb(:,1:6)` = `rd`、`rdb(:,7:12)` = 回転 $P(\mu)$。`btr` は (14.8) の $1 - \Lambda'\otimes\Lambda'$ の 21 成分版 | `tecalc` | |
| `emitn`・`beamn` | 正規座標の $\Sigma_n$ | `tecalc` | |
| `sr` | $4/(a_x + a_y + a_z)$ | `tecalc` | |
| `params(ipdampx:...)`・`ipjx` | 減衰率・分配数 | `tecalc` | |
| `taurdx` など | 減衰時間 | `tecalc` | |
| `synchm` | 縦運動が有効 (`rfsw` かつ $\mu_z \ne 0$) | `tecalc` | |
| `alphap`・`omegaz`・`bh` | $\alpha_p$、$\omega_z$、バケット高さ | `tecalc` | `momentum_compaction`、`omegaz` |
| `heff`・`vceff`・`wrfeff`・`dvcacc` | 実効 harmonic、実効電圧、RF の角周波数、電圧の $z$ 微分 | `tmacro` | `heff` |
| `tnorm`・`tsymp`・`teigenc` | 固有ベクトルの規格化・シンプレクティック化・固有値問題 | `tnorm.f`、`toplvl.f`、`teigen.f` | `_tnorm`・`_tsymp`・`teigen` |
