# 16. SAD 物理計算ガイド

最終更新: 2026-11-06 (現状のスナップショット)

式中心。SAD = K64 版 (github KatsOide/SAD)、ソースは `<SAD>/src`、ヘルプは CERN SADHelp。
詳細・ベンチ結果は docs/05 (CALC)、06 (マッチング)、07 (emit)、08 (追跡)、09 (DA)、11 (癖)、work_*.md を参照。行番号は 2026-10 時点の src/ 以下。
凡例: [確認]=ソースを読んだ / [転記]=式は分かるが導出していない / [解釈]=コードからの推論。pysad の欄は 直訳=行単位の移植 / 独自=別実装 / 未=未実装。

## 1. 座標・単位・キーワード

- 1.1 座標 (x, px, y, py, z, δ)。px, py は p0 で規格化した**正準運動量** (tquade.f:39 のコメント)。ドリフト pzi = sqrt((1+δ)^2 - px^2 - py^2)。追跡の内部だけは px/(1+δ) (tconvm.f、§13)。
- 1.2 K1, K2... は**積分強度** (K1*L)。tsolque.f:390 `akk=ak/al`、tbend.f:60、tffs.f:1316。L=0 の QUAD は K1 が薄レンズ強度 (tquad.f:22)。SEXT 等は tthine.f:31 `aki=ak/kord!` (§7)。
- 1.3 BEND の **E1/E2 は ANGLE との比**: `psi1 = E1*ANGLE + AE1` (tffs.f:1261)、tbende.f:448、INTGRL.f:91。逆向き要素では E1<->E2, AE1<->AE2, FB1<->FB2 を入れ替える。
- 1.4 `phib = ANGLE + K0` (磁場の曲げ)、`phi0 = ANGLE` (基準座標)。`rhob = al/phib`, `rho0 = al/phi0`, `drhob = rhob - rho0` (tbende.f:430)。K1≠0 の駆動項 `drhop = 1/rho0 - 1/rhoe` (rhoe = rhob*p)、p=1 なら -K0/L。
- 1.5 キーワード: DISFRIN (内部 FRIN, fringe=(DISFRIN==0))、FRINGE (mfring=nint)、DISRAD (0 で放射)、DISKIN (追跡のドリフトを近軸に)。
- 1.6 速度項 [転記] tturne.f:847-867:
  `dvemit = -(p1+p0)/h1/(p0*h1+p1*h0)*dp + dvfs`、`p1=p0(1+dp)`、dvfs=FSHIFT (RFSW のみ)。要素ごとに `z -= dvemit*L`。(5,6) 成分に `h0/h1^3*L`。
- 1.7 **SAD の `sinc(x) = x cos x - sin x`**、`sinhc(x) = x cosh x - sinh x` (tfloor.f:369-395、小さい x は級数)。sin x/x ではない。qdbend 等の `a12 = rho/phi*sinc(phi)*dphi` はこの定義 (§17)。
- 1.8 FRINGE/DISFRIN の既定: pysad は FRINGE=0, DISFRIN=0 を仮定する。tqfrie は既定で入り、切るのは DISFRIN=1 (A15)。

## 2. 計算の流れ [確認]

```
CALC : tffsa -> tffsmatch -> tffscalc -> qcell1 -> { qcod(qtwiss.f:810) ; qmdiag ; qtwiss(mat=.false.) }   (calc6d なら qcell6d, §18)
EMIT/Emittance[] : temit.f (tcod -> tturne -> tecalc)        RADINT : tffsa.f:617 -> intgrl (INTGRL.f)
TrackParticles : tftrack -> tturn -> tturn1 (tturn.f:286-552)      DA : dapert.n -> tfdapert -> trackd -> tturn
```
要素ルーチンは 3 系統。CALC (qtwiss1、4x5 行列) と 6D (tturne1) は同じ tquade/tbende/tthine を使い、**追跡は別コード** (tquad, tbend, tthin) (A16)。qtwiss1 は IDEAL で要素値 idvalc、通常は各位置の値を読む (qtwiss.f:189)。pqcell/qquad(単体ファイル)/qthin はビルドされない。

## 3. ドリフト

CALC (qtwiss.f:205-224) [確認]: 閉軌道上で線形化。
```
pr=1+d; a=px^2+py^2; dpz=-a/pr/(1+sqrt(1-a/pr^2)); pzi=pr+dpz; ale=L/pzi; alz=ale/pzi^2
a12=ale+px^2 alz  a14=px py alz  a15=-pr px alz  a34=ale+py^2 alz  a35=-pr py alz
x+=px ale; y+=py ale; z += dpz ale - dvemit L
```
追跡 tdrift_free (tdrift.f:10): `al1=L/(1+dpz)`, `z += dpz al1 - dv L`。DISKIN≠0 で近軸 (tturn.f:296)。
**薄い多極 tthine の中のドリフトにも dv 項** (`tdrife0(dvon=.true.)`、`z -= L dv(d)`、A6): これが抜けると TransferMatrices の (5,6) が六極 1 台ごとに 1.3e-9 ずれる。

## 4. 四極

- tquade (tquade.f:1-82): `al=0 -> tthine`, `ak=0 -> tdrife0`, `K1 L<0 -> theta-pi/2, ak=-ak`。順序: tsolrote (x-=DX, y-=DY, 回転) -> tqfrie -> tqlfre -> tsolque -> 出口を逆順。取付誤差は座標変換で入る (feed-down は自動)。
- tsolque (tsolque.f:380-596): `eps=0.1*EPS (EPS=0 なら 0.1)`, `ndiv=min(1000, 1+int(|al hypot(ak,bz)|/eps))`, `akk=ak/al`。tqente(aln/2) -> {本体 -> tqente(aln)}... -> 本体 -> tqente(aln/2)。
  本体 (bz=0): `x'=x dc1 + px s1/w1, px'=-x w1 s1 + px dc1` (y は cosh/sinh)。z 補正は [転記] (tsolque.f:420-425)。
  **z 行はシンプレクティック条件で組む** (tqfrie, tbendebody0, tsolque; A11):
  `T(5,j) = -T(1,j)T(2,6) + T(2,j)T(1,6) - T(3,j)T(4,6) + T(4,j)T(3,6)`。数値微分だけでは一周行列が 6.6e-5 非シンプレクティックになり εy が 2% ずれる。
- 追跡 tquad: tsolqu は `eps=0.2*EPS` (EMIT の 0.1 と違う)、`ndiv=1+int(|L hypot(ak,bz)|/eps)`。**tsolqur (RAD ON, bz≠0) は最後の粒子の bzp を全粒子の最後の半ステップに使う** (tsolqu.f:268,273,519-531): `bzp=bz/(1+d_i)` がスライス内の粒子ごとに更新されるのに最終ステップで再計算されない。結果が粒子の並びに依存 (差 ∝ d_i - d_last、Level 9 で SAD 自身が 1.15e-6 変わる)。放射なしの tsolqu は粒子ループ内で完結。
- tsolqum (放射なし MULT) は z の dv 項が違う: tsolqu `-0.25(...+ap aln0...) - dv aln0`、tsolqum `-0.25(...+(ap-dv) aln0...)` (tsolqu.f:86,250; A2)。
- pysad: 写像は直訳、CALC の行列は数値 Jacobian (eps=1e-7: 1e-6 だと K1=0 BEND の分岐境界で (1,6) に 5.6e-7)。

## 5. 曲げ磁石

- tbende の分岐 (tbende.f:390-423): ANGLE=0 かつ K1=0 -> tsteee、K0=0 -> tquade、他 -> tmulte。phib=0 -> tbdrifte。L=0 -> tbthie (`cod(2) += (phi0-phib) + d phi0`)。通常: tchge -> tblfre(F1) -> tbedge -> 本体 -> tbedge -> tblfre -> tchge。
- K1=0 本体 tbendebody0 (tbende.f:223-317): `xf = rhob(sn pxi - cn dpzi - dcn dp + dpzf) + cn xi - dcn drhob`、`da=asin(..)`、`zf = z - da p rhob - phi0n (drhob + d rhob)`。
- K1≠0 本体 tbendiinit (tbend.f:41):
  `rhosq=rho0 rhoe; drhop=(rhoe-rho0)/rhosq; akxsq=-1/rhosq - akk/p; akysq=akk/p`、
  `dxf = dcxkx(drhop + akxsq xi) + sxkx pxi`、`dpxf = sxkx(drhop + akxsq xi) + dcx pxi` (tbende.f:146)。
  **tbendi (追跡) の `rhoe = rhob (1+d)` は入口で 1 回だけ決め、放射で d が変わっても更新しない** (module bendib、A5)。
- tbendecorr [転記] (tbende.f:61-121): `pxf -= ak1 xf xr (0.5 - xr(2-xr)/12)/p`。スライスは [tbendecorr(平均) + tbendebody] 反復、最後に tbendecorr(akx0/2, alx0/2)。
- 分割 (tbende.f:453-460): `rbc=1-(f1r+f2r)/al0; phic=phi0 rbc; eps=epsbend(1e-3)*EPS; ndiv=1 (K1=0) / 1+int(|phic/eps|)`。
- フリンジ長補正 alcorr [転記] (tbende.f:381): `al = al0 - ((phi0 fb1)^2+(phi0 fb2)^2)/48/al0 * sin(.5(phi0(1-psi1-psi2)-apsi1-apsi2))/sin(.5 phi0)`。同式が pos (LINE["S"]) にも使われる (§20)。
- 追跡 tbendcore (tbend.f:372-527): 入口端面・F1・本体・出口端面を閉じた式で 1 ステップ [転記]。K1≠0 は tbendi (§13)。
- ロール誤差 DROTATE は中点まわりの 3 次元回転 (tbrote/tbgrote、LGEO=L/2, ANGLGEO=ANGLE/2; tfgeo.f:212)。回転の向き `x'=cos x - sin y`。DTHETA は QUAD では ROTATE、BEND では DROTATE。
- pysad: tbende は lattice/emit6d に直訳 (F1・ANGLE=0 分岐は一部)。

## 6. エッジ・フリンジ

- tbedge (tbfrie.f:88-125) は fringe/mfring に無関係に常に呼ばれる: 入口 psi≠0: `tbedgedrift(psi) -> tbedgemaxwell(rhob) -> tbedgebody(rhob,psi)`、出口は逆順 (-rhob)。tbedgemaxwell: `af=(1-(y/rhob)^2/6) y/rhob/pvi`。
- BEND F1 (tblfre, tbende.f:650): `x += dxfr d/pr`, `py += (dyfr - dyfra y^2) y/pr`; `dxfr=fb^2 phibl/24`, `dyfr=fb phibl^2/6`, `dyfra=fringe ? 4 dyfr/fb^2 : 0`、tbendef1 は長さ 0.5 fb rb の tbendebody0 を 2 回。
- tqfrie [転記] (hard-edge): `akk=ak/al/4; a=akk x^2/pr; b=akk y^2/pr; t=(aki(x^2-y^2))^2/6; dx=x(a/3+b+t); dy=-y(a+b/3-t)`。DISFRIN=0 で入る。
- tqlfre (線形): `af1=-ak/al f1|f1|/24; af2=ak/al f2; x -> e^a x + (b/p) px, px -> px e^-a, y -> y e^-a - (b/p) py, py -> py e^a` (a=af1/p, b=af2/p)。mfring=1,3 で入口、2,3 で出口。係数は tsetfringepe (tffs.f:982): K1≠0, L≠0 のときだけ f1in=F1+F1K1F 等。
- 追跡: ttfrin (nord=4/6/一般)、追跡の線形フリンジは `f1=-|akk F1^2|/24, f2=|akk| F2` (K1 の符号によらず f1<0)。hard-edge の条件は mfring>-4 (tquad.f:49,103) で tquade と違う。

## 7. 六極・薄レンズ・MULT

- tthine (CALC/6D, tthine.f:1-93) [確認]: `kord=nord/2-1; aki=ak/kord!`。L≠0: `drift(L/6) -> kick(aki/2) -> drift(2L/3) -> kick(aki/2) -> drift(L/6)`。kick: `cx=(x-iy)^kord; px -= aki Re cx; py -= aki Im cx`。行列は閉軌道上の微分 (feed-down)。enarad=.false.。
- 追跡 tthin (tquad.f:118-322): 4 次積分 `alpha=1/sqrt12, beta=1/6-alpha/2, gamma=1/40-1/(24 sqrt3), alpha1=1/2-alpha`: drift(alpha1 L) -> 半キック -> drift(alpha L) -> 補正キック(beta,gamma) -> drift(alpha L) -> 半キック -> drift。入口出口で ttfrin、放射は最後に tradk(L) 1 回。キックは `aki=K/kord!/(1+d)`。
- MULT: CALC は qmult1 -> tmulte (ANGLE=0)、追跡 tmulti: tsolrot -> ttfrins(n=0..2, dofr) -> [tsolqu/tsolqum + n>=2 キック]xndiv -> tsolqu -> ttfrins -> 戻し。放射は各 tsolqu の後。追跡の MULT にはフリンジ ttfrin が入り、emit の tmulte には入らない (A16)。ANGLE≠0 は tmulta (pysad 未)。

## 8. 空洞

- CALC の qcav は `rfsw = rfsw .and. trpt` で tcave を呼ぶので、**リングの 4D CALC では RFSW に関係なく空洞は長さ L のドリフト** (qcav.f:19; A4)。追跡の CAVI は RFSW に関係なく tcav (RFSW OFF は周回末に z=0 にするだけ、tturn.f:142)。TCAV (横方向) だけが OFF でドリフト。
- tcave (L=0, 6D): `w=omega0 HARM/c`, `phis=trf0 w`, `v=VOLT/amass |charge|`, RFSW ON: `phi=w t + PHI - phis; dh = -va(sin phi + offset1); va = v + v20a x^2 + v02a y^2; v20a = v (w/p0)^2/4`、`offset1 = sin phis` (RADCOD OFF) / 0 (ON)。**PHI の既定は 0** で同期位相は trf0 (`asin(U0/Vc)/w`、R56 の符号で π-asin) として決まる (A3)。長さあり: `ndiv=1+int(min(|wL|, sqrt((v(1/h0+1/h1))^2/3/0.01)))`、tcavfrie が入口出口。
- 追跡 tcav (tcav.f): `tlim=1e4 (RFSW)/0`, `t=clamp(-z h1/p1, +-tlim)`, `dh = max(oneev-h1, V(-2 sin ph cos(ph-dphis) + offset) wsn)`、RAD または TRPT で offset=sin dphis/offset1=0、それ以外 offset=sin dphis - sin phis/offset1=sin phis。

## 9. 閉軌道 qcod (qtwiss.f:810-970) [確認]

```
itmax=10 (最大11 パス); conv = min(1e-10, 1e-18 (1 + (d/1e-3)^2 + nlat))          (qtwiss.f:821,837)
毎パス: cod0 から qtwiss1(mat=.true.) で一周 -> T (4x5); limitcod でクランプ
cx=(T11+T22)/2 (|cx|>1 なら 1/cx) から sx, ax, bx を推定
r = dx^2/bx + bx (dpx + ax/bx dx)^2 + (y 同様)     [= (1+a^2)/b dx^2 + 2a dx dpx + b dpx^2 の Courant-Snyder 不変量]
r <= conv: 収束。Newton: (T-I) cod0_new = T cod0 - cod_end を tsolvg (SVD, 1e-8) で解く
悪化 (r>r0): fact(初期 .5) >= 1e-3 なら cod0=(1-fact) cod00 + fact cod0, fact=|fact|/2; fact<1e-3 なら fact*=16, cod0=(1+fact) cod00 - fact cod0
改善: fact=min(.5, 2 fact)
```
閉じるのは (x,px,y,py) のみ、d 固定、z は閉じない。**返す行列は最後に一周させた出発点 cod0 の行列で収束後に評価し直さない**。ゆえに off-momentum の β は初期値で相対 1e-6 揺らぐ (A13)。初期値は twmov で復元した前回安定解 (先頭 MARK に保存、tffscalc.f:93,118)。Emittance[] の tcod (6D) は z を閉じる d を選ぶ (周長変化で d≠0; 3C で -4.85e-8、A12)。RFSW ON: 先に RF なし 4D 解 -> `z=asin((U0-vcacc)/Vc)/w - trf0` -> 6D Newton、収束 1e-4 (重み codw=(1e-6,1e-7,1e-7,1e-8,1e-5,1e-6))。pysad は割線法/Newton で独自 (物理量は一致)。

## 10. 周期解・Twiss・分散・結合

- qmdiag [転記] (qmdiag.f:7-75): `M=S^-1 Md S, S=(c I, J R^T J; R, c I)`。`trpq=T11+T22-T33-T44; trab2=4[(T31+T24)(T42+T13)+(T14-T32)(T41-T23)] + trpq^2`、`c^2 = (1+|trpq|/sqrt(trab2))/2`、`aa=∓1/(c sqrt trab2)`、`r1=aa(T31+T24), r2=aa(T32-T14), r3=aa(T41-T23), r4=aa(T42+T13)`。trab2<0 は "Sum resonance"。pysad は numpy 固有値。
- |cos mu|>1: 常に 1/cos mu に置き換えて続行。INTRES(>1)/HALFRES(<-1) が立っているときだけ stabx=.false. (qcell.f:209-248; A14、pysad は `lattice._qmdiag`/`_qcell_periodic_sad` で直訳)。`sin mu = sign(sqrt(1-cos^2), a12); alpha=(a11-cos)/sin; beta=a12/sin`。
- **分散は周期解** (qcell.f:267-279): `(a13,a23)=S T(:,5)`, `dcosmu=2 sin^2(mu/2)`, `EX = .5(a13 + a12 a23 - a22 a13)/dcosmu`, `EPX = .5(a23 + a21 a13 - a11 a23)/dcosmu` (x_cod/d の差分ではない; A10)。モード分散で、物理分散 PEX は tgetphysdispu (tgfun.f:239): `PEX=cc EX + r4 EY - r2 EPY, cc=sqrt(1-(r1 r4-r2 r3))` (未)。
- 伝搬 qmat2twiss (qtwiss.f:705): `aa=bx0 T11-ax0 T12; cc=(T12-ax0 aa)/bx0; BX=T11 aa+T12 cc; AX=-(T21 aa+T22 cc); dNX=atan2(T12,aa)`、`EX = T11 ex0 + T12 epx0 + sqrdet T15 - r4 T35 + r2 T45`。pysad の `_transport_twiss` は非結合の標準式 (途中位相は結合なし)。結合素子 (qtwiss.f:541-691、DETR の 1+xyth-detR フラグ) は未。
- Twiss 配列 `twiss(nlat,-1:1,28)` (tffs.f:101): 1-6 AX..NY、7-10 EX..EPY、11-15 R1-R4 DETR、16-21 DX..DDP (閉軌道)、22-28 AZ BZ NZ ZX ZPX ZY ZPY。`twiss(l)` は要素 l の入口。位相は rad、tune=NX/2pi。

## 11. 放射積分 RADINT (INTGRL.f)

RADINT コマンドだけが使う (Emittance[] は INTGRL を使わない; 結合 BEND では 13% 差 [原因未特定。INTGRL は 4D の積分で、bendib の rhob/rho0 を設定せず直前の値 (最後の BEND の L/ANGLE) を全 BEND に使う癖があり (docs/11 A34)、候補の 1 つ])。
**pysad は `pysad/radint.py` が INTGRL.f の直訳で、FFS の RADINT コマンドに接続済み** (ESRF 6 通り、小リング 3 通りが SAD と 7 桁一致。history/work_small_items.md A)。
`COEFF=(elradi/elmass^3/3)1e27`, 対象は BEND/QUAD/SEXT/CAVI。BEND: `IR=min(14,max(1,int(|ANG|/.002))), IRP=2IR`、SYMPS (複合 Simpson)。K0, AE1, F1 は読まない。
被積分: `I2=CURV^2 DL, I3=|CURV|^3 DL, I4 = (CKL0X + 2CKL1X/CURV) CURV^2 (R11 EU + R12 EPU), I5 = |CURV|^3 (GU EU^2 + 2AU EPU EU + BU EPU^2) DL`。エッジ項 TEDGE (CURV tan e)。出力: `DAMPU=1-I4/I2`, `EMITU = CO5 E^2 I5/(I2-I4)`, `SIG = CO4 E sqrt(I3/(2I2+I4u+I4v))`。SEXT の DKL2 は未代入 (寄与実質 0、バグの疑い)。lattice の `_radiation_integrals_bend` は BEND のみ・非結合の別実装。

## 12. 6D エミッタンス (temit.f / photons.f)

- 流れ: tcod -> `irad=12; tinitr12; tsetr0` -> tturne (trans(6,12): 1:6=T, 7:12=D、beam(21)=B) -> tecalc。
- tradke (photons.f:1409): 基準 codr0 から現在の運動量方向の外積 `xpa`, `theta=asinz(xpa)`, `h1=p2h(p0 pr)`, `al1=al-cod(5)+codr0(5)`:
  `anp = anrad h1 theta; uc = cuc h1^3/p0 theta/al1; dg = -cave anp uc (cave=8/(15 sqrt3)); dee = anp uc^2 cuu (cuu=11/27)`。`anrad=5/(2 sqrt3) rcratio, cuc=1.5 r_cl/rcratio` (tfsetparam.f:105)。xpa=0 なら更新なし。RADCOD ON のときだけ cod を更新し `tesetdv(g)`。
  D <- D + tr1 (T+D) (tmuld6)、B <- (I+tr1) B (I+tr1)^T + dbeam (tmulbs; irad<12 なら何もしない)。
- 分割: BEND (enarad) `nrad=1+int(|al0 rbc/epsr1 crad (h0 b)^2|)` (epsrad=1e-6)、`ndiv=max(ndiv, int(nrad emidiv emidib), 1+int(|phib h0 anrad|/epsr1/1e6 ...))`。QUAD `nrad=int(|al crad/epsrad (h0 b)^2|)` (itgetqraddiv)。SEXT 等は放射なし。
- **tecalc の定義** (temit.f:235-510): `rx=T` (RFSW OFF なら `rx(6,1:5)=0, rx(6,6)=1`) -> teigenc -> tnorm (各モードのシンプレクティック積 s=1、z 成分最大を第 3 モード、残りを x 成分で並べる) -> `r=tsymp(r), ri=tinv6(r)`。`rd=ri D r`、対角 2x2 から `cd`: `taurd=-2pi/omega0/Re cd; J=4 Re cd/sum Re cd`。B を ri で正規座標へ、`btr` (21x21) を組み `tsolva(btr, beam, emitn, 1e-8)` (SVD) で `Sigma = M Sigma M^T + B` を解く。`eemx=sqrt|e11 e22-e12^2|, emitp=r emitn r^T, sige=sqrt(emitp(6,6))`。
  synchm=(rfsw and Im cd6≠0): 真 `alphap = -Im(cd6)|Im(cd6)|/(c 2pi/omega0)/(dvcacc/pgev)`、`sigz=sqrt emitp(5,5)`、`bh=sqrt|Vc/(pi|a|h pgev)(2cos phi - (pi-2phi)U0/Vc)|`、`phi=asin(U0/Vc)`。偽: `alphap=-R56/C` (`-rx(5,6)/2pi/c/p0 h0 omega0`)、`sigz=|alphap| sige c p0/h0/omegaz`、z の行列は 0。
- 単位: EnergyLossU0, RfVoltageVc は **eV**。TransferMatrices は先頭からの**累積**行列。
- pysad emit6d.py: 全桁 1.3e-11 (ESRF v03_12、16b/16e)。RADCOD ON も SAD と 2e-11 以内 (tcod の直訳 + 解析ヤコビアン)。INTRA は §19、SynchroBetaEmittance(temits)・POL は未。

## 13. 粒子追跡と放射 (tturn.f, photons.f)

- 内部座標 (tconvm.f): `px <- px/(1+d)`, 粒子ごとの `dv = -d(2+d)/h1/(h1+(1+d)h0) + dvfs`。要素外側・粒子内側 (tturn1)。**tapert**: 最初の要素、以後 15 要素ごと、周末。先に `limitnan` (NaN -> +x、|a|>x は +-x。x,y は 3 LOSSAMPL、px,py は 0.9999、z は 1e10; tturn.f:7-18)、判定は `-alost<x<alost, -alost<y<alost, |z|<=LOSSDZ` を否定形で (NaN も落とす)。落ちた粒子は末尾の生きている粒子と入れ替え (並びが変わる)。RFSW OFF は周末に z=0。
- ndivrad (tffs.f:641): `arad=0.02; aka=(|phi|+hypot(ak,bz) arad)/eps1; ndivrad=int(1+h0/100 aka) (TRPT なら x10)`。K1=0 BEND: `ndiv=min(1024, ndivrad(phir,0,0,eps))`。K1≠0 (tbendi): `max(1+int(|phi0/eps|), ndivrad)`、上限 1000。
- tradkn (FLUC なし) / tradkfn (FLUC あり):
```
theta = asinz(|pp|)  (入口出口の単位運動量の外積);  al1 = al - z + zr0
anp = anrad h1 theta
tradkn : uc = cuc h1^3/p0 theta/al1;  de = cave anp uc
tradkfn: rho = al1/theta; uc = cuc h1^3/p0/rho;  de = dph uc  (dph = sum x_i, tdusrn)
dg = -de/(1+2de);  g = max(gmin, g+dg), gmin=-0.999999
tradkn : ddpx = -.5 dpx dg; x += ddpx al1/3      tradkfn: ddpx = -r1 dpx dg; x += r2 ddpx al1
dv = -g(1+pr)/h2/(h2+p2) + dvfs  (tgetdv の h2+pr h0 とは違うが SAD のまま);  z <- z p2/h2 h1/p
```
  基準点 pxr0, pyr0, zr0 は要素入口と各 tradk の後で更新。tradke (emit) は別式 (dg=-cave anp uc)。
- tdusrn (drndsr.f:71): `p=tran()` を**必ず 1 回**引き (光子 0 個でも)、`p<=pn` になるまで光子を加える (ポアソン累積、打ち切り pn>0.9999 または 10000 個)。光子ごとに drndsr (横谷の逆関数近似) と位置用 `tran()` を 1 回ずつ。NaN の粒子では両条件が偽になり 10000 回ループし約 20001 個の乱数を消費する (gfortran の MAX/MIN は NaN を無視; D2)。損失周の乱数のずれの原因候補 (B6、未解決)。
- 乱数 (random_plugin_sad.c): `seed <- 48828125 seed` (int32 折り返し)、`Random[]=(seed+2^31)/2^32`、GaussRandom は Box-Muller (`rho=sqrt(-2 log(seed/2^32+.5))`, rho sin を次回へ)、`SeedRandom[n]: seed=|n/2| 2+1`、初期 17。FLUC は要素外側・粒子内側、放射ありの BEND/QUAD/MULT はスライスごとに全粒子の本体 -> 全粒子の tradk (B4)。NPARA の既定は 4 (init.n)、9 粒子以上で fork し乱数列が分かれるので比較は NPARA=1。
- RADTAPER [未]、PHOTONS 光子リスト [未]。pysad: track6.py/track6n.py (numba)、Level 6-12・16 で ≤5e-11、乱数状態まで一致。

## 14. マッチング (tffsmatch.f, tffscalc.f, qdcell.f) — 詳細は docs/06

- 残差: BX,BY は `log目標 - log値`、AX,AY は atan、NX 等は差 (qdcell.f tdfun1)。`r = wsum (sum|df wiq|^p/wsum)^(2/p)`、p=ExponentOfResidual。off-momentum は `wi=(OffMomentumWeight/2/sqrt|mfitp|)^2` (ヘルプは sqrt(点数) と違う)。
- 重み wfit (tffscalc.f:442-501) [転記]: EX `dpm/sqrt(emxx bx)`、R1..R4 `coup sqrt(bx/by)` 等。emx が NaN のとき gfortran の MAX/MIN は NaN を無視し `emxx=0.01/1.01 emminv(1e-15)=9.9e-18` となり EX/EPX の重みが 1e6 前後になる (esrf/matching で別解になる原因、work_analytic_deriv)。
- 反復: Newton `qu dval = df` を tsolva(eps=1e-6/fact) で解く、dg>0 なら最急降下、非改善は fmincube の 3 次補間。モード切替は fuzz(`max(0,min(1,(x-a)/(b-a)))`): crate<0 で log10 が NaN -> fuzz=1 -> モードを変えず直線探索 (C7)。収束 `res2r(r) <= CONVERGENCE max(nfcol,1)`、MAXITERATION は上限そのものではない。
- 微分は既定で**解析的** (§17)。FFS$NumericalDerivative=True や SOL があると数値微分。
- pysad: match6.py (Matcher; ループ・fuzz・fmincube・数値微分は直訳)、解析微分は match6d.py。matching.py (Nelder-Mead) は非推奨の旧 API で SAD 非準拠。

## 15. DynamicApertureSurvey (trackd.f, tapert.f) — 詳細は docs/09

- 軸: `Axes` の 1 文字目 ivar1=値リスト (最大 200 点)、2 文字目 ivar2=51 点で走査、残り ivar3 (同じ j で同時に動く)。`a2step=(max-min)/50`。起動順 i=1..n1p, j=51->1 (RAD+FLUC の乱数順に効く)。
- 初期座標 (tinip1): `xa=(x,px,y,py,z,d) x (emx,emx,emx,emx,emz,emz)` -> `r` (temw, 正規->物理) -> `+codin` -> tconv。`emx=sqrt(|EMITX|+|EMITY|)`, `emz=sqrt|EMITZ| (RFSW) / |SIGE|`。PhaseZ=-pi/2 で z=0, d=a。
- 損失窓: 幅 ncons=max(1,nint(DAPWIDTH)) (既定 7) の窓を下へ動かし `nxm(i)` を下げる。Score=sum nsc (nsc は最初に損失した j-1)。判定は tapert (§13)、alost=LOSSAMPL(1)、zlost=LOSSDZ(100)。損失周は `ntloss=完了周数`。trackd は TRPT を .false. にして戻さない。
- Damping: 毎周 `exp(-T0/taurd)` を正規座標で掛ける (tpdamp)、DampEnough は早期に生存扱い [未]。
- pysad: dapert.py (trackd の順序を再現) + track6n.k_apert。NORAD の地図 153 文字一致。

## 16. 要素値と誤差の展開

- 値は 3 層: 要素定義 idval、代表位置 klp (Element[key,"QD"])、各位置 LINE[key,"QD.n"]。`LINE["DX","QD.1"]=v` は CALC 時に tffsadjust (tfsetv.f:95-148) が klp の値を全 QD に複写する。`LINE["DX","QD.5"]` は上書きされる (C1)。1 台だけの誤差は `DELX sigma QD.1`。
- 誤差: DK=K1 相対、DDK=K0 絶対 (errk を CALC で掛ける)、DELX/DELY/DTHETA (QUAD: ROTATE, BEND: DROTATE)。ガウス値は GCUT で棄却、同じコマンドは値を置き換える。
- **部分要素 qfraccomp** (qtwiss.f, Twiss[key,{"Q1",0.3}] や LINE[.,200.5]): `r=rx2-rx1`, `L r, ANGLE r, K0 r, K1 r, VOLT r`、**E1 は rx1=0 のときだけ E1/r、E2 は rx2=1 のときだけ**残り他は 0 (E1 は ANGLE 比なので部分角では 1/r 倍)、FRMD=-f1-2 f2 で端の側だけフリンジ。BEND: FRINGE=0 なら F1=0。pysad は ffscmd.row_frac が線形近似 (未対応: K1 のスケール等)、geo.py の部分要素は qfraccomp どおり。
- IDEAL [未]、expnln の統計誤差チェーン (expln.f:57) [転記・未]。

## 17. 解析微分 (qdcell.f tffsqu, qdtwis.f, qdquad.f, qdbend.f, qdmdia.f) [確認、pysad match6d.py = 直訳]

qu の定義: `qu[cond,var] = df/dv wiq / wvar` (tffsmatch.f:395)。f は比較の形: BX は log b (dtwiss=db/b)、AX は atan (da/(1+a^2))、NX は psi。呼び出し: `tffsqu -> qdcell(k,lp) -> qdtwis -> qddrif/qdbend/qdquad/qdthin/qdtrans(数値)`、周期補正と dgammag。
- 要素: `dtrans = dM_k S_k` (S_k=qtentu の R 変換、5 列は eta_k)。QUAD K1: `k=K1/(L pr), phi=sqrt|k| L, dphi/dK1 = phi/(2 K1)`、`da11=-sin phi dphi, da12=(cos phi - sin phi/phi) dphi/sqrt k, da21=-(cos phi + sin phi/phi) dphi sqrt k` (発散は sinh/cosh)。BEND K1: `a12 = rho/phi sinc(phi) dphi` (sinc は §1.7 の定義; sin x/x と取ると BPI 列が 10-300% ずれる)。端面は `tan(psi1 phi0)` の薄レンズ。
- 輸送: `trans1 = tftmatu(k+1->l) dtrans`、`X = qgettru(k->l)` (k>=l は一周をまたぐ)、`x_ij = X/sqrt(1-detR)`:
  `db_l/b_l = 2(-(dx11 x12 + x11 dx12) a_k + dx11 x11 b_k + dx12 x12 g_k)/b_l`、
  `da_l = (dx11 x22 + x11 dx22 + dx12 x21 + x12 dx21) a_k - (dx11 x21 + x11 dx21) b_k - (dx12 x22 + x12 dx22) g_k`、
  `dpsi = (dx12(b_k x11 - a_k x12) - x12(b_k dx11 - a_k dx12))/(x12^2 + (b_k x11 - a_k x12)^2)`、`deta_l = trans1(1,5)+(R の項)`。
- 周期補正 (CELL): 一周行列の微分 `dM_l = d(M_{k->l}) M_{l->k}`、qdmdia で dr1..dr4、`dmd = dS (M S^-1) + S (dM S^-1) + S (M dS^-1)` から
  `d mu = -(da11+da22)/(2 sin mu)`, `d alpha = (.5(da11-da22) - alpha cos mu d mu)/sin mu`, `db/b = (da12 - b cos mu d mu)/(b sin mu)`。位相: k<l は `dpsi_l = d mu + sin psi'[(cos psi' + a sin psi') db/b - sin psi' da]` (psi'=mu-psi_l)、k>=l は `-sin psi_l[(cos psi_l - a sin psi_l) db/b + sin psi_l da]`。閉軌道 `dx=(I-X)^-1(dM x + dx_k)` (det=4 sin^2(mu/2))。
- 照合: 数値微分と 2e-8 (RING)/4e-7 (ARCA)。esrf/matching の SAD 反復列 (8.7338E-09) と K1 を 12-13 桁再現。未移植: qdthin、qdini、qddtwiss、BMAG、COUPLE、qdtrans (数値で代用)。

## 18. CALC6D の拡張 Twiss (qcell.f:404 qcell6d, tffs.f:2413 tfetwiss, tturne.f:775) [確認、pysad calc6d.py]

`CELL; CALC6D; CALC;` で qcell1 -> qcell6d -> `temit(plot=.true.)` (Emittance[] と同じ)。tcod -> tturne -> tecalc (stab=|Re log lambda|<1e-6) -> 一周を plot 付きで再実行し、要素入口ごとに `tsetetwiss`: `ti = ri inv(trans)` (l=1 のみ temw の ri を掛ける)、`twi = tfetwiss(ti, cod)`、`NX = NX(1) + twi(NX)` (twi<-toln=-0.1 なら +2pi)。NZ は 2pi 補正なし。tfadjustn が NX,NY を単調化。RAD は結果に影響せず (temit は放射減衰込み)、RADCOD は閉軌道に放射損失、RFSW は閉軌道の z と tecalc の rx 6 行目。不安定時 "qcell6d---> Unstable optics" で cell=.false. にし transport 伝搬。
tfetwiss(ri, cod): z ブロック `az=sqrt det ri[5:6,5:6]`, `uz = ri[5:6,5:6]/az`、x,y と z の結合 hx, hy -> H(6x6)、`hi = ri H`、`detm=(det_12 + det_34)/2`、**detm<=xyth(0.375) なら x-y モードが入れ替わりとみなし行を入れ替え**。`axy=sqrt detm`、`R1..R4=R/axy`、`bx=(hi12^2+hi22^2)/axy^2`、`NX=atan2(-hi12, hi22)`、`EX..ZPY = axy h +- R h`、`DETR=R1 R4-R2 R3`、`DX..DDP`=閉軌道、`AZ=bz21 bz22, BZ=crz^2, NZ=atan2(sz,cz)`。逆変換 etwiss2ri (tffs.f:2569)。RFSW OFF では z 固有値が 1 (Jordan) で `r(:,5)=e5, r(:,6)=(eta,0,1)` となり BZ,AZ,NZ は R56 に対する擬似 Twiss (18b: BZ=1.019)。分数位置は qtwissfrac の calc6d 分岐で `trans=inv(etwiss2ri(twiss(l)))` から部分要素 (irad=6) を伝搬。FFS 起動時 calc6d は .false. に固定 (track.f:104)。
照合 (bench_level18a/b/c): ESRF 全 28 成分 1.7e-10。pysad は不安定時の ri を ri0 とする点が SAD と違う可能性あり (未検証)。

## 19. Touschek と IBS (temit.f, tintrb.f, Packages/emit.n) [確認、pysad touschek.py]

INTRA ON の emit は一周追跡の各要素で `tintrb` を呼び、(a) IBS 拡散を beam(1:21) に加え (b) 最終パスで Touschek 損失率表 (touckl/touckm/toucke) を積算する。`temit.f` の converge ループ (`ITMAX=100, RESIB=3e-6, DCMIN=1e-6`、MINCOUP 下限)、`eintrb` で加速、収束後 `caltouck=T` でもう一周。
- 部品: `touckf(x)` (emit.n TouschekF, C(x)/x の近似; x=0 で 1e200)、`fintrbu` (IBS 被積分)、`bint.f` (台形則の細分 + Simpson 外挿)、`teigens33` (thess3+tqr3)、`sols33`、`tmulbs(calint=T)` で beam(22:42) も同じ J で変換、`tradke` は `B2_new = B_new + J (B2 - B) J^T`。表: dp=0.004-0.242 (0.002 刻み, 120)、振幅表 tampl(34,3)、`TouschekTable = [id, [dp, 1/tau], [[nx,ny,nz], mx(34x25), my(34x25)], toucke(nlat x 120)]`。
- `TouschekLifetime[nx,ny,nz]` (emit.n Touschek): spz, spxz, spyz の自然スプライン (spline.f, 範囲外は端区間で外挿)。nx,ny が両方 INF なら dp 受容、片方 INF なら振幅-nz の 2 次元補間、両方有限は NaN。PBUNCH=0 なら IBS も Touschek も 0。
- pysad: 3x3 固有値だけ LAPACK でなく 2 側 Jacobi (eigh は絶対誤差で sigma_y^2 が 1% 狂う)。17a (RADCOD OFF) 寿命 15 ケース 7e-12、17b (RADCOD ON) 1.7e-8。FFS 接続 (run_emit/fn_touschek) は未。

## 20. 要素幾何 GEO (tfgeo.f, tgrot.f, qtwiss.f qfraccomp) [確認、pysad geo.py]

`geo(3,4,k)`: 要素 k 入口のビーム x,y,z 軸の実験室系方向余弦と位置。初期値 geoini: x軸=(0,-1,0), y軸=(0,0,-1), z軸=(1,0,0) (正の ANGLE で GY 増)。角度 (tfchi): `CHI1=atan2(zy,zx), CHI2=asinz(zz), CHI3=atan2(xz,-yz)`。FFS からは `LINE["GX"|"GCHI1"|"GEO"|"OGX"..., loc]` (Twiss["GX"] は K64 ではエラー)。
BEND/MULT (ANGLE=phi≠0): ROTATE=theta (≠0 なら x'=cos x - sin y, y'=sin x + cos y)、`rho=Lgeo/phi, (sp,.,cp,dc)=xsincos(phi), r1=rho sp, r2=-rho dc`、`pos' = pos + r1 z - r2 x'`、`x'' = cp x' + sp z`、`z' = -sp x' + cp z`。ANGLE=0 は直進 (ROTATE=-90 の VD でも)。COORD は回転・並進、MARK の OFFSET は後処理で geo をコピー (整数) か部分要素 (小数)。
`Lgeo = L - ((phi fb1)^2+(phi fb2)^2)/L/48 sin(.5(phi(1-E1-E2)-AE1-AE2))/sin(.5 phi)` (FRINGE≠0 のみ; §5 の alcorr)、これが `pos(i+1)=pos(i)+Lgeo` (LINE["S"]) にも使われる。部分要素は §16 の qfraccomp (E1 は E1/r)。loc が nlat を超えると k-nlat に周回。ORG コマンド: `geo(:,4,:) -= dgo; geo(:,1:4,k) = geo2 geo(:,1:4,k)`、chi は度。getvad は負数を式と読む (`1.5 -2.5` は -1.0)。GEOCAL/GEOFIX フラグ。未: SOL (tsgeo)、MAP (tgmap)、SORG、逆向き要素の部分要素。照合 19A 4.8e-13、19B/19C 1e-13。

## 21. MULT と SOL 領域の追跡 (tmulti.f, tsol.f, tsolqu.f; pysad track6n_ext.py)

- MULT: §7。ANGLE≠0 (tmulta)、F1/FRINGE/VOLT≠0 は未。`akr[n] = K_n cr1^(n+1)` (cr1=theta2 の回転)、ndiv0 は tpara が ampmaxm=0.05 から決め `ndiv=max(ndiv0, ndivrad(|K0|,K1))` (<=1000)。ak=0 で K0 のみの MULT は厳密ドリフトのキック (tdrift_kick)。多極キックは Horner `cx` 複素で n>=2。
- SOL 領域 (tsol.f:510-): BOUND=1 の SOL から次の BOUND=1 の SOL まで。`bzs = charge (BZ+DBZ)/(p0/c)` (tfbzs)。入口 trots -> `tradks`(F1≠0) -> `tsfrin(bzs)` (DISFRIN=0)、出口 `tsfrin(-bzs)`。領域内: DRIFT `tdrift_solenoid` (`px += bzp y/2` で回転)、QUAD `tquad(bz)` -> ttfrin の bz 分岐 + tsolqu/tsolqur の bz≠0 分岐 (tzsetparams; §4 の bzp の癖)、CAVI `tmultiacc`。tapert の la カウンタは領域内で tsol のループ (tsol.f:604) と主ループの両方で数える (損失周の番号付けのみに影響)。
- track6n_ext: track6.py の直訳を numba 化 (演算順序・括弧も同じ)、track6 と 3e-17、FFS 乱数 seed 一致。ソレノイド内四極は EPS を使わない (eps0=0、EPS≠1 で SAD と分割数が違い得る)。

## 22. 未解析・未検証リスト

[転記のみ・導出なし]: dvemit と (5,6) (tturne.f:847)、tsolque の z 補正 (tsolque.f:420)、tbendecorr、tbendcore/tbendi (追跡)、tbedgebody (tbfrie.f:127)、tqfrie の高次項 t、tthin 4 次積分の係数、qmdiag の c,r、qcell の A,B ブロック (qcell.f:144)、tfetwiss の H 構成、tecalc の btr と SVD (temit.f:371)、dbeam の dee 以外の成分、tdusrn の p 打ち切り、wfit、expnln、IBS 被積分 fintrbu、tintrb の Touschek 積分、tsolvg の SVD 細部。

[未検証]: qdthin/qdini/qddtwiss、ttfrin の一般次数の導出、tmulta、tsole1/qsol (SOL の CALC)、tsgeo、tcavfrie 中身、TCAV (ttcav, ttcave)、temits (SynchroBeta)、POL (calpol, spnorm)、INTRA の wspac 分岐、INTGRL の BMAG 状態依存 (bendib)、RADINT と Emittance[] の 13% 差、CHRX の意味、FFS 内 FLUC 初期値、OffMomentumWeight とヘルプの差、DA の Output->0、既定 keyword の個別変更 (`LINE["K1","QD.5"]`)、tinip1 の r 定義、損失周の乱数消費 (B6)、qcell6d 不安定時の ri (§18)、FITP>1 の GO (WSL ビルドで再起動 D1)。

## 23. pysad 対応表

| SAD の計算 | SAD ルーチン | pysad | 状態 |
|---|---|---|---|
| CALC のドリフト | qtwiss1 インライン | lattice `_drift_matrix_at_cod`, `_exact_drift_cod` | 直訳 |
| QUAD 写像 (CALC) | tquade, tsolrote, tqfrie, tqlfre, tsolque, tqente | lattice `_quad_offset_map`, `_tquade_cod`, `_tqfrie_cod`, `_tqlfre_cod`, `_tsolque_cod`, `_tqente_cod` | 直訳 (行列は `_num_jacobian`) |
| BEND K1=0 | tbedge, tbendebody0, tblfre, tbendef1 | lattice `_tbedge_cod`, `_tbendebody0_cod`, `_tblfre_cod`, `_bend_cod_update` | 直訳 |
| BEND K1≠0 | tbende, tbendecorr, tbendebody | emit6d `_tbende`, `_BendIBody`; lattice `_bend_k1_cod_matrix` | 直訳 |
| BEND 誤差・K0 | tchge, tbrote, tbgrote, tdrife1 | lattice `_bend_err_wrap`, `_tbrote_cod`, `_tbgrote_cod`, `_tdrife1_cod` | 直訳 |
| SEXT 等 | tthine | lattice `_tthine_cod`, `_thin_offset_map` | 直訳 |
| MULT (CALC/6D) | tmulte, tmulke, tbfrme, tsoldz | lattice `_mult_sequence`, `_mult_cod_map`, `_tsolconv`, `_tsoldz_cod`; emit6d `_tmulte` | 直訳 (ANGLE=0) |
| 閉軌道 | qcod | lattice `_qcod_periodic`, `_one_turn_cod_matrix` | 直訳 |
| Twiss・分散・結合 | qcell1, qmdiag, qmat2twiss | Lattice.twiss, `_transport_twiss` | 周期解は直訳、伝搬は独自 (非結合)、qmdiag は numpy |
| 6D エミッタンス | temit, tturne, tradke, tecalc, tsolva | emit6d `tturne`, `_State.tradke`, `tecalc`, `_tsolva`, `emittance6d`, `_normal_form`, `_tcave` | 直訳 (閉軌道 tcod は独自) |
| Emittance[] の戻り規則 | emit.n | emit_extra `emittance_rules`, `beam_size_at`, `rules_to_sad_list` | 直訳 |
| CALC6D | qcell6d, tfetwiss, etwiss2ri, tsetetwiss | calc6d `calc6d`, `tfetwiss`, `etwiss2ri`, `twiss_frac` | 直訳 (tcod は独自) |
| Touschek/IBS | tintrb, bint, touckf, spline | touschek `touschek_table`, `touschek_lifetime`, `_tintrb`, `bint`, `touschek_f`, `SadSpline` | 直訳 (固有値のみ Jacobi) |
| 幾何 | tfgeo1, tgrot, tfchi | geo `compute_geo`, `GeoResult`, `tfchi`, `tforbitgeo` | 直訳 |
| マッチング | tffsmatch, tffscalc, tfsolv, tgrad | match6 `Matcher.go`, `_tffscalc`, `_tfsolv`, `_tgrad`, `_fmincube`, `_fuzz` | 直訳 (微分は数値) |
| 解析微分 | tffsqu, qdcell, qdtwis, qdquad, qdbend, qdmdia | match6d `tffsqu`, `qdcell`, `qdtwis`, `qdquad`, `qdbend`, `qdmdia`, `AnalyticMatcher` | 直訳 (qdthin 等は未) |
| 追跡の要素 | tdrift_free, tquad/tsolqu, tthin/ttfrin, tbend*, tcav, tmulti | track6 `_tdrift_free`, `_tquad_ens`, `_tsolqu`, `_tthin_gen`, `_ttfrin_n`, `_tbend_ens`, `_tbendi`, `_tcav`, `_tmulti_ens` / track6n `k_drift`..`k_cav` | 直訳 |
| 追跡の SOL/MULT (numba) | tsol, tmulti | track6n_ext `k_mult`, `k_sol_frin`, `k_sol_drift`, `k_sol_quad`, `ProgExt`, `track_particles_ext` | 直訳 (track6n 本体へ未統合) |
| 追跡の SOL (純 Python) | tsol, tsolqur | track6 `_tsol_region`, `_tsolqu_bz_ens`, `_tsolqu_bz_final`, `_tquad_sol` | 直訳 |
| 放射・乱数 | tradk, tradkn, tradkfn, drndsr, tdusrn | track6 `_tradk`, `_tradk_core`; track6n `tradk_core1`, `tdusrn`; sadrandom `SadRandom` | 直訳 |
| DA | trackd, tapert | dapert `dynamic_aperture`, `rad62`; track6n `k_apert`; ffscmd `fn_dasurvey` | 直訳 (Damping 等は未) |
| RADINT | INTGRL | `radint.compute/report` (FFS `RADINT`) | 直訳 (BEND/QUAD/SEXT、結合、エッジ、SAD の癖込み)。`lattice._radiation_integrals_bend` は BEND のみ・非結合の別実装 |
| 要素値・誤差 | tffsadjust, terror | Lattice.set_element/set_error, ffscmd `set_family`, `_error_cmd` | 独自 |
| 物理分散 PEX, 逆向き要素の mfring, SOL の CALC, qmult(ANGLE≠0), POL, temits | tgetphysdispu 等 | — | 未 |
| matching.py (非推奨の旧 API; Nelder-Mead 等) | — | matching `Match`, `DynamicApertureSurvey` | 独自 (SAD 非準拠) |
