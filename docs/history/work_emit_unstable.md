# 不安定な光学の EMIT / INS + CALC6D / close_z の精度 (2026-10-16)

docs/history/00_RESUME.md の未解決リストの 3 項目 (A, B, C) を、SAD (K64, WSL) の実測とソースの読みで 1 項目ずつ合わせた。
担当: `pysad/emit6d.py`, `pysad/calc6d.py`, `pysad/ffscmd.py`、新規 `pysad/teigen.py`。

| 項目 | 内容 | 結果 | 検証 |
|---|---|---|---|
| A | 不安定な光学での EMIT / Emittance[] (teigen + tnorm の直訳、Deviation 行、Emittance[] の Stable) | **実装** (V 不安定は SAD と行単位で全一致、H 不安定は SAD の T の悪条件の範囲 (1e-5) まで一致) | 28c, 28d, 28f |
| B | INS + CALC6D (qcell6d の cell=.false. の枝) | **実装** (`calc6d.calc6d_ins`) | 30a (≤1.6e-10) |
| C | RFSW ON + RADCOD OFF + FSHIFT≠0 の閉軌道 (tcod の直訳) | **実装** (`_tcod_radcod(radcod=False)`)。旧 close_z の 5e-5 ずれは消えた | 30b (≤1.2e-11) |

回帰: `QUICK=1 run_regression.py` FAIL 0 (28c/28d/28f/30a/30b を追加、48+2 ベンチ)、pytest 164 pass (`--ignore` 3 本、新規 `test_emit_unstable.py` 5 本)、
16b (全桁) は SAD と ≤9e-14、13a/14a/16a/16b/18a/22a/26a〜d/28a/28b は作業前と同じ (差は 1e-17 の雑音の桁のみ)。

---

## A. 不安定な光学での EMIT / Emittance[]

### A.1 SAD の流れ (temit.f)

`temit` は不安定でも止まらない。`tecalc` (temit.f:235) は

1. `r=trans(:,1:6)` (RFSW OFF なら `rx(6,1:5)=0, rx(6,6)=1`)、`teigenc(r,ri,ceig,6,6)` → `tnorm(r,ceig,lfno)` → `r=tsymp(r)`, `ri=tinv6(r)`
2. `cd(4:6)=log(cc(1:5:2))` (Re = 増大率、Im = 位相)、`stab = (|Re cd(4:6)| < 1e-6) .and. fndcod`
3. 以降は stab に関係なく Emittance (EMITX..SIGZ の rsetgl1) まで計算する。stab は `temit` の戻り値 → `Emittance[]` の最初の要素 (`Stable`、不安定なら 0)
   と、INTRA ON のときの `Skip intrabeam because of unstable.` にしか使われない。「Unstable」の文字列は temit では出ない
   (`Vertical Unstable.` は CALC 側の qcell1。`*****qcell6d---> Unstable optics` は CALC6D)
4. temitprint (temit.f:949): 一周行列のシンプレクティック偏差 `so=Σ|rx·tinv6(rx)−1|` が 1e-8 を超えたら
   ` *** Deviation from symplectic matrix = ` + list-directed real*8 (17 桁・指数 3 桁) を `Extended Twiss` の前に出す。
   `Imag.tune:` の行に出るのは Re(cd)/2π (増大率)、`Real tune:` が Im(cd)/2π (チューン): 見出しが逆になっている (SAD の癖)

固有値が実数 (不安定) のとき tnorm は `imag(ceig(i)) == 0` の分岐を通る: 列の大きさを `(sa/sb)^(1/4)` でそろえ、シンプレクティック積 s で規格化、
最後に `r(:,i+1)=r(:,i+1)*r(i,i)−r(i,i+1)*r(:,i)`, `r(:,i)=r(:,i)/r(i,i)`。λ と 1/λ のどちらを先にするかは teigen の返す順で決まる。

### A.2 pysad の以前 (`_normal_form`)

numpy の `eig` で複素固有値のモードだけ取り出し、実固有値は「RF なしの z の Jordan 形」だけ特別扱いして捨てていた。
不安定な面 (実固有値の対) は 2 モードしか見つからず、`BY=3.5e-33`、減衰率 −26 などの別の値になった (28c の旧出力)。

### A.3 実装

- **`pysad/teigen.py` 新規**: teigen.f (tbal の vx, thess の高速 Givens、tqr の二重シフト QR と固有ベクトルの逆代入・Jordan 処理・偶数次元の並べ替え) の直訳。
  1 始まりの配列で書いた。検証: **SAD の `Eigensystem[T]` (SADScript、実行列は teigen を呼ぶ) に pysad の T を渡して並びと固有ベクトルを比べた**
  (28e の一周行列 6×6 で、並び [λx, 1/λx, 複素対, λy, 1/λy]、実固有値のベクトルは全桁一致。複素対のベクトルは複素定数倍だけ違うが tnorm が規格化するので無関係)。
- **`emit6d._tnorm`**: tnorm.f の直訳。`_normal_form(rx)` = teigen + `_tnorm`。旧版は `_normal_form_old` に残し、環境変数 `PYSAD_OLD_NF=1` で切り替えられる (回帰比較用)。
  安定な光学では旧版と同じ r (差は 1e-17 の雑音の桁だけ。13a/14a/16a/16b/18a/22a/26a〜d/28a/28b/3c/4 を新旧で比べた)。
- **`_order_real_pairs`** (SAD の癖への対処、後述 A.4): teigen が返した実固有値の対のうち、2×2 ブロックの対角の差 a−d が雑音 (相対 1e-8 未満) のときだけ増大側 (λ>1) を先にする。
- `emit6d.tecalc` が `rx` を返す。`ffscmd.print_emit` が Deviation 行を出す。Emittance X/Y/Z・Energy spread の書式を `autofg(…,'11.8')` に (旧 11.5。
  固定小数点で出る大きな値 `.02324788` で桁数が違った。指数形では同じなので既存ベンチには影響しない)。
- 新ベンチ: `bench_level28c.sad` (V 不安定: QD ×1.19、既存)、`28d` (H 不安定: QF ×0.9)、`28e` (両方不安定)、`28f` (V と H の `Emittance[]` の戻り値)。
  SAD の結果は `bench_level28{c,d,e,f}.result`。`run_regression.py` の 28c / 28d / 28f。

### A.4 結果と、再現できないもの

| ケース | SAD との差 |
|---|---|
| 28c (V 不安定) | EMIT 報告 37 行が**行単位で一致**。差は 1e-17〜1e-30 の雑音の桁 (閉軌道の Entrance/Exit、R1-R4、EY、Beam tilt の 4e-11) だけ |
| 28f V の `Emittance[]` | `Stable=0`、Emittances・σ・Tunes・DampingRate・α・U0・Vc・OrbitAtExit が一致 (差は雑音のみ) |
| 28d (H 不安定、βx が 2.7e4) | チューン・Twiss・U0・α は一致。減衰率 6e-7、εz 8e-6、σE 4e-6、バンチ長 4e-6 の差。**SAD の T 自身が悪条件** (Deviation = 1.16e-6、pysad は 2.0e-6)。K1 を 1e-12 変えるだけで εx が 3e-6 変わる (実測) ので、これ以上は丸め雑音の再現が必要 = 不可能 |
| 28e (両方不安定) | **再現不能**: SAD の出力自体が NaN と減衰率 1e52 (teigen の並びが雑音で変わり、tnorm が λx と 1/λx を別のモードに並べる)。pysad は teigen の順 (SAD の `Eigensystem` と同じ) で正常な値を返す。回帰には入れない |

SAD の癖: 始点が鏡映対称 (α=0) だと、不安定な面の 2×2 ブロックは a=d で、`a−d` の符号 (丸め雑音) で teigen が λ と 1/λ のどちらを先に返すかが決まる。
28c では SAD は増大側が先 (pysad の自然な雑音は逆だったので、`_order_real_pairs` で増大側を先にした)。28d は自然な順が増大側。28e は SAD が減衰側。つまり SAD の雑音の符号は一定ではなく、
`_order_real_pairs` は「雑音のとき増大側」という仮定 (28c と 28d に合う)。a−d がはっきりあるときは teigen の順のまま。

---

## B. INS + CALC6D

### B.1 SAD の経路

- `qcell` (qcell.f:62): `if(calc6d) call qcell6d(...); return` — CELL か INS かに関係なく 6D 経路に入る。4D には落ちない。
- `qcell6d` (qcell.f:404) の `cell=.false.` の枝 (INS):
  ```
  ri0  = etwiss2ri(twiss(lb,idp,1:28), normali)       ! 先頭 MARK の 28 成分 (tffscalc.f:93 twmov(right)=cmp%value(1:ntwissfun))
  cod0 = twiss(1,idp,mfitdx:mfitddp)                  ! DX DPX DY DPY DZ DDP (DDP は DP0 の別名)
  twiss(lb) = tfetwiss(ri0, cod, .true.)              ! 28 成分を ri0 から作り直す
  beam = tfinibeam(lb)                                ! is==1 なら temw の ri = etwiss2ri(twiss(1)) もここで更新する
  tinitr(trans); tturne0(trans,cod,beam,srot,fbound,iaez,idp,.true.,rt,.true.)
  ```
  tcod も temit も呼ばない (閉軌道は探さない)。
- `tturne0` → `tfsetplot` → `tsetetwiss` (tturne.f:772): lorg=1 のとき `ti = ri·tinv6(T_{1→l})` の `ri` は上の temw の ri (= etwiss2ri(tfetwiss(ri0)) )、
  `norm=normali` (最初の etwiss2ri の normal)。`twi=tfetwiss(ti,cod,norm)`、l=1 は位相 0、以降 `twiss(l0)+twi` (<−toln なら +2π)、最後に `tfadjustn` (NX, NY)。
- 空洞の位相基準 trf0 は直前の temit/tcod が残した tmacro の値。放射は CELL の CALC6D と同じく RAD フラグに依らない。
- INS は MARK に書き戻さない (`if(cell) twmov(..,.false.)` は CELL のみ)。
- 先頭 MARK の 28 成分: `twmov(..,.false.)` は ntwissfun=28 全部を書く。4D の CELL CALC の後の AZ BZ PSIZ ZX ZPX ZY ZPY は直前の 6D の先頭行の値
  (まだなければ AZ=0, BZ=1, 他 0 — SAD の `A_ORG` で確認)。

### B.2 pysad

- `calc6d.calc6d_ins(lattice, p, tw0, rfsw, radcod, trf0, …)`: 上の流れの直訳。`_turn` (CELL の CALC6D と同じ 6D 一周追跡の記録付き) を cod0 から回し、
  `_rows_from_turn(st, ri1, normal0, els, row0=tw1)` (行 1 は `tfetwiss(ri0,cod0,True)`、以降は `ri1·inv(T)`)。以前の「不安定時の枝」とほぼ同じ形だが、
  ri1 を作り直す点と行 1 の扱いを SAD に合わせた。
- `ffscmd.calc`: `not CELL and CALC6D` なら `calc6d_ins` (`_mark_tw28` で MARK の 28 成分を読む。未設定は 0、BX/BY/BZ は 1、DDP は DP0)。trf0 は `self.trf0`。
- `_twmov_to_mark` が縦方向の列 (AZ BZ PSIZ ZX ZPX ZY ZPY) も MARK に書く (4D の CALC のときは `_z_stale`)。

### B.3 検証 `bench_level30a.sad` (ESRF、5 通り、先頭 MARK を含む 8 位置 × 28 列)

A: 4D の CELL CALC → `INS; CALC6D; CALC;`、B: 6D の CELL CALC6D → INS、C: 先頭 MARK に軌道オフセット (DX DPX DY DPY)、
D: BXI/AXI/BYI/AYI と DDP=1e-3、BZ=5、AZ=0.01、E: RFSW RADCOD ON。
**全 32 行 (28 列) が SAD と ≤1.6e-10** (cmp_lines の最大相対差。B の 3e-12 は CELL の周期解の雑音)。先頭 MARK にオフセットを付けた C/D/E も一致。

---

## C. RFSW ON + RADCOD OFF + FSHIFT≠0 の閉軌道

### C.1 SAD の tcod (tcod.f)

RFSW ON: まず `rfsw=.false.` で自分を呼ぶ (RADCOD OFF なら radcod はそのまま)。戻ったら `cod(5)=asinz((u0*pgev−vcacc)/vceff)/wrfeff − trf0`、`im=6`。
本体は 6D の Newton:

- 各反復で `tturne(…, update=.true.)`: RADCOD OFF のとき **毎回 trf0 を更新** する (tturne.f:37-58: `vceff=hypot(vcacc,dvcacc/wrfeff)`、
  `trf0 = asin(u0 pgev/vceff)/wrfeff` (trans(5,6)<0) か `(π−asin)/wrfeff`、`vceff<u0 pgev` なら π/2/wrfeff、`mod` で折り返し)。
- `dtrf0=trf0−trf00`、`r=Σ(dcod1/codw)² + (dtrf0/trw)²` (trw=codw(5)、RFSW OFF の再帰呼びでは codw×0.01)。収束 `r<conv0=1e-4`。
- RADCOD OFF では z を trf0 に移す処理 (`codi(5)−=dz …`) は**ない**。Newton は `tsolvg(trs(1:im,:), dcod1(1:im), dcod)` で `cod=codi; codi=codi+fact*dcod`。
- `rfsw` が偽 (再帰呼び) では `dcod1(5)=0`, `trs(:,5:6)=0` → z も δ も動かさない (4D の閉軌道)。

収束判定が緩い (z は 1e-7 まで) ので、最終値は反復の履歴で決まる。以前の pysad は「z=0 が閉じる δ を割線法 (厳密) で求める」近似で、FSHIFT=1e-6 (δ=−4.3e-3) のとき
SAD の「途中で止まった値」と x で 5e-5 (相対)、px で 6e-9 ずれた。

### C.2 pysad

`_tcod_radcod(els, ph, cod, rfsw, trace, final_full, radcod=True)` に `radcod` 引数を追加 (既存の RADCOD ON の直訳に、RADCOD OFF の枝を足した):

- 各反復の tturne は常に T 付き (`want_db=False`)。直後に `_tturne_update_trf0(st, trf0)` (tturne の update の直訳: vceff は hypot で vc0 へのフォールバックなし、vceff=0 なら trf0 据え置き)
- `r += (dtrf0/trw)²`、`trf0s` の保存 (NaN のとき復帰)、RADCOD ON の z→trf0 の移し替えは `if radcod:` の中だけ
- `closed_orbit6d`: `rfsw and not radcod` で空洞があれば `_tcod_radcod(…, radcod=False)`。収束しなければ従来の割線法に落ちる。
  RFSW OFF の RADCOD OFF (z を閉じる δ の割線法) は変えていない (21 ベンチ以上が一致していて SAD の挙動が別に説明できているため)。
- `calc6d._closed_orbit` と `emittance6d` は `closed_orbit6d` を共用するので CALC6D と Emittance[] の両方が直る。

### C.3 検証 `bench_level30b.sad` (ESRF、`FSHIFT=1e-6`、`RFSW; NORADCOD;`)

FSHIFT = 1e-6 / 3e-7 / −2e-6 / 0 の 4 通りで `Emittance[Orbit->True]` の εx εy εz、U0・Equilibrium dz・Orbit dilation・α・バケット高さ・Vc、チューン、
OrbitAtEntrance/Exit、σE・σz、DampingRate を比較。**最大相対差 1.2e-11** (旧: x が 5e-5、チューン 3e-7、U0 6e-9)。
`bench_level16b` (RFSW、全桁) は SAD と ≤9e-14、`bench_level26a〜d` と `bench_level4` は変わらず一致。

---

## 道具・注意

- 新旧の固有値分解の比較: `PYSAD_OLD_NF=1` を付けて実行すると旧 `_normal_form_old`。
- WSL の `gs` で SADScript の `Eigensystem[T]` を実行すると teigen の結果を直接見られる (T は Python の repr で SAD スクリプトに埋め込む)。
- 28d/28e の `.result` は SAD の悪条件の挙動の記録。28e は回帰に入れない。
- 作業前の基準出力: `base_old/` (旧 `_normal_form`)、`base_new/`・`base_fin/` (新)。

## 未解決として残るもの

1. 28d の 1e-5 までの差、28e の NaN: SAD の T の丸め雑音そのもの (再現不能)。
2. `_order_real_pairs` は「雑音の符号は増大側が先」という仮定 (28c/28d の SAD に合わせた)。別のラティスで SAD が逆の順を返したら合わない。
3. RFSW OFF + RADCOD OFF の閉軌道 (δ の決まり方) は tcod の直訳ではなく近似のまま (SAD では `trs(:,5:6)=0` で δ が動かないはずなのに、実測は z を閉じる δ。
   temit の呼び出し側が cod の δ を決めている可能性 — 未調査)。
4. INS + CALC6D の後の `Twiss[]` 以外 (GO/FIT で INS + CALC6D の条件を使う場合) は未確認。
