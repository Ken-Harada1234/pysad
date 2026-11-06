# CALC6D (6 次元 CALC) の移植 — pysad/calc6d.py 作業ログ

日付: 2026-10-15
対象: SAD K64 (src/qcell.f, temit.f, tturne.f, tffs.f, qtwiss.f, tcod.f) → `pysad/calc6d.py` (新規、既存ファイルは変更なし)
照合: `bench_level18a.sad` (RFSW RAD), `18b` (NORFSW RAD), `18c` (RFSW NORAD) ↔ `bench_level18a_pysad.py`

## 1. SAD の CALC6D の経路 (読み解き)

```
FFS: CELL; CALC6D; CALC;
  tffscalc (tffscalc.f)
    ntfun = ntwissfun (=28)  … calc6d のとき 28 成分すべてを扱う (tffscalc.f:87, 553)
    → qcell → qcell1 (qcell.f:62)   if(calc6d) call qcell6d
  qcell6d (qcell.f:404-458)
    codplt=.true. ; ri0 = etwiss2ri(twiss(lb,:)) ; cod0 = twiss(1, dx:ddp)   (直前の閉軌道が初期値)
    if(cell): call temit(trans,cod,beam,btr, calint1, iaez, plot=.true., params, stab, lfno=0)
              if(.not. stab): '*****qcell6d---> Unstable optics' を出し cell=.false. で再分岐 (下記)
    else   : twiss(lb,:) = tfetwiss(ri0, cod0, .true.)  →  tturne0 (plot=.true.) で transport 線として伝搬
  temit (temit.f:83-)   ※ Emittance[] と同じルーチン
    1. tcod (tcod.f)            6D 閉軌道。RFSW ON なら一旦 RFSW OFF で求め、cod(5) = asin((U0−Vacc)/Veff)/wrfeff − trf0
    2. tturne (plot=.false.)    trans(6,12) = [T | D] (T: シンプレクティック一周行列, D: 放射減衰の摂動)
    3. tecalc (temit.f:339-)    rx = T; RFSW OFF なら rx(6,1:5)=0, rx(6,6)=1
                                teigenc → tnorm → r = tsymp(r), ri = inv(r)
                                stab = |Re log λ_x,y,z| < 1e-6 かつ fndcod
                                params(iptwiss…) = tfetwiss(ri, cod, .true.)
    4. plot かつ postcal なら trans=I, cod=codin からもう一度 tturne(plot=.true.)
  tturne1 (tturne.f:334-)    要素ループの先頭 (要素 l の入口) で tsetr0 → tsetetwiss(trans,cod,beam,lorg=1,l,idp)
  tturne0 (tturne.f:184)     最後の要素 nlat ($$$) は tfsetplot → tsetetwiss(一周の trans, 一周後の cod)
  tsetetwiss (tturne.f:775-826)
      ti = inv(trans) ; lorg<=1 なら ti = ri·ti (temw の ri)
      twi = tfetwiss(ti, cod, normali)
      l=1 : NX=NY=NZ=0
      l>1 : NX = twiss(1,NX) + twi(NX) (+2π if twi(NX) < −toln) ; NY 同様 ; NZ = twiss(1,NZ)+twi(NZ) (2π 補正なし)
  tfadjustn (tturne.f:716)   plot 後に NX, NY を l=2..nlat で単調化 (前より toln 以上減ったら +2π を累積)。NZ はしない
```

- `toln = 0.1` (tffs.f:2316), `xyth = 0.375` (tffs.f:210)。
- 分数位置 `Twiss[key, 200.5]`: qtwissfrac (qtwiss.f:1008 の `if(calc6d)` 分岐)
  `trans = inv(etwiss2ri(twiss(l)))`, cod = twiss(l) の閉軌道、部分要素 (qfracseg(0,fr)) を tturne1 で伝搬 (irad=6, codplt=.false.)、
  `ftwiss = tfetwiss(inv(trans), cod, normal)`、NX/NY/NZ に twiss(l) の位相を足す。
- `Twiss[key, 100]` は twiss(100) (要素 100 の入口) そのまま。

### RAD / RADCOD / RFSW の扱い
- **RAD**: temit の経路 (tturne1 → tbende/tquade/tmulte/tradke) は ffs_flag の `rad` を参照しない
  (tbende の enarad は要素キーワード RAD、tradke は radcod だけ見る)。`rad` を見るのは追跡 (tturn.f, tbendi.f 等) だけ。
  SAD の 18a (RAD) と 18c (NORAD) の出力は全桁一致した。→ CALC6D は RAD に依らず放射減衰込み。
- **RADCOD**: tcod/tradke で閉軌道に放射損失を入れる (RADTAPER と併用)。
- **RFSW**: tcod (閉軌道の z 成分)、tecalc (rx の 6 行目)、tcave (空洞が動くか)。RFSW OFF では z モードの固有値 1 (Jordan)
  で tnorm の実固有値分岐: r(:,5)=e5, r(:,6)=(分散, 0, 1) → BZ, AZ, NZ は「R56 に対する擬似 Twiss」になる
  (18b: BZ=1.019, AZ=0.139, NZ=−0.138 など)。

### 拡張 Twiss 28 成分の定義 (tfetwiss, tffs.f:2413)
ri (物理座標 → 正規モード座標) を受け取り、
1. z ブロックから az = sqrt(det ri[5:6,5:6]), uz = ri[5:6,5:6]/az、hx, hy (x, y と z の結合 2×2) を作る
2. H (6×6) を組み、hi = ri·H で z を分離した 4×4 部分を得る
3. detm = (det 1-2 ブロック + det 3-4 ブロック)/2。detm ≤ xyth なら x-y モードが入れ替わっているとみなし行を入れ替え
4. axy = sqrt(detm), R1..R4 = (R 行列)/axy, βx = (hi12²+hi22²)/axy², αx = …, NX = atan2(−hi12, hi22) (y, z も同様)
5. EX, EPX, EY, EPY, ZX, ZPX, ZY, ZPY = axy·h ± R·h の組合せ、DETR = R1·R4 − R2·R3
6. DX..DDP = 閉軌道、AZ = bz21·bz22, BZ = crz², NZ = atan2(sz, cz)

逆変換 etwiss2ri (tffs.f:2569) は 28 成分から ri を再構成する (normal = DETR < 1)。

### 不安定時
qcell6d は 'Unstable optics' を lfno (=0) に出して cell=.false. にし、直前の twiss(1) から etwiss2ri → tturne0 (plot) で
transport 線として伝搬する。ただし tsetetwiss (lorg=1) は temw モジュールの `ri` (tecalc が最後に作ったもの) を使うので、
twiss(1) を ri0 で置いた直後に l=1 の tsetetwiss が ri で上書きする (ソース上そう読める。SAD で未検証)。
pysad では意図どおり ri0 = etwiss2ri(twiss0) で伝搬する (twiss0 省略時は β=1)。

## 2. pysad/calc6d.py (対応表)

| SAD | pysad.calc6d |
|---|---|
| qcell6d | `calc6d(lattice, momentum_eV, rfsw, radcod, rad, twiss0, out, with_emit)` → `Calc6DResult` |
| tcod | `_closed_orbit()` — emit6d.emittance6d の手順をそのまま使う (tcod の直訳ではない。RFSW/RADCOD の 3 通り) |
| tturne / tturne1 | `_turn(record=True)` — emit6d.tturne と同じ要素分岐。各要素入口で (T_{1→l}, cod_l) を記録、最後に一周後 |
| tsetr0 (要素ごと) | `_turn` 内の `st.setr0()` |
| tecalc (固有値, r, ri, stab) | `emit6d._normal_form` → `_tsymp` → `_tinv6`; stab = |Re log λ| < 1e-6 |
| tecalc (平衡ビーム) | `emit6d.tecalc` を呼んで `result.emit` に入れる (Emittance[] 相当。失敗しても Twiss は返す) |
| tsetetwiss | `_rows_from_turn()` |
| tfadjustn | `_adjust_phase()` |
| tfetwiss | `tfetwiss(ri, cod, normi)` → 長さ 28 の ndarray |
| etwiss2ri | `etwiss2ri(tw)` → (ri, normal) |
| qtwissfrac (calc6d) | `twiss_frac(result, l, fr)` (部分要素は ffscmd.FFS.frac_element を遅延 import) |
| twiss(l,0,1:28) | `result.rows[l-1]` (lattice.TwissRow、列は lattice.py の mfit* 順) |

- SAD は一周を 2 回 (plot=.false. と plot=.true.) 走らせるが同じ初期値・同じ写像なので pysad は 1 回で記録する。
- `_State6d` は emit6d._State の tradke に ON/OFF を付けたもの。既定は ON (SAD の挙動)。`radiation=False` は実験用。
- `rad` 引数は受け取るだけ (SAD と同じく結果に影響しない)。
- TwissRow の `name` は lattice.twiss_table と同じ規約 (rows[0]='*'、rows[i] = 要素 i の出口 = 要素 i+1 の入口なので
  要素 i の名前)。SAD の「要素 l の入口」は rows[l-1]。
- 計算時間: ESRF (3563 要素) で 1 回の calc6d が約 7–10 s (閉軌道探索の数周 + 記録の一周。数値 Jacobian が主)。

## 3. 照合結果 (bench_level18a_pysad.py)

参照: `bench_level18{a,b,c}.sad` = pysad/tests/kek-ls_v03_12.lat + `FFS; $FORM="18.14"; PageWidth=1999; NPARA=1; USE RING; CELL; <flags>; CALC6D; CALC;`
位置: `$$$`, `1`, `QFA1H.3`, `SFA.5`, `MKC.2`, `RFC`, `100`, `200.5` の 28 成分。
(`LINE["LENGTH","$$$"]` は narg エラーになり POS 行は出なかった → 位置は pysad 側で要素名の n 番目から求め、
`LINE["NAME",{100,200,201}]` = {QFA8H.3, DR08.22, SDA.9} と NLAT=3564 で要素列の一致を確認)

| ケース | フラグ | 結果 |
|---|---|---|
| 18a | RFSW RAD | 全位置・全 28 成分で max\|diff\| = 1.7e-10 (NX≈299 rad の 6e-13 相対)。β, α は 1e-11、EX は 1e-13、AZ/BZ/NZ は 1e-12 |
| 18b | NORFSW RAD | 同じく max\|diff\| = 1.7e-10。RFSW OFF の AZ=0.139088, BZ=1.019345, NZ=−0.138201 も 1e-12 で一致 |
| 18c | RFSW NORAD | SAD 出力は 18a と全桁同一。pysad も 18a と同じ計算 → 一致 |

数値ゼロの成分 (EY, R1..R4, DX.. など、SAD で 1e-16 以下) は丸め誤差レベルの差のみ。
出力は `bench_level18{a,b,c}_pysad.out` に保存。
`result.emit` (tecalc) の εx=1.549e-10, σδ=6.95e-4, σz=3.77 mm (RFSW) は emit6d と同じ値。

## 4. ffscmd.py への接続案 (他担当向け)

`FFS.calc()` (ffscmd.py:410) で CELL かつ `self.flags['CALC6D']` のとき:
```python
from .calc6d import calc6d, twiss_frac
if self.flags['CELL'] and self.flags['CALC6D']:
    prev = self.rows[0] if self.rows else None          # SAD: 直前の twiss(1) が不安定時の初期値
    res = calc6d(self.lat, self.momentum(), rfsw=self.flags['RFSW'], radcod=self.flags['RADCOD'],
                 twiss0=prev, out=self)                  # out.write に '*****qcell6d---> Unstable optics'
    self.rows = res.rows
    self.calc6d_result = res                             # twiss_frac と Emittance 相当 (res.emit) に使う
    self.trf0 = res.trf0
    return self.rows
```
- `row_frac(k, fr)` (ffscmd.py:525) は CALC6D のとき `twiss_frac(self.calc6d_result, k, fr)` を返す。
- `twiss_value` の `_TWKEY` は既に AZ/BZ/NZ/ZX/ZPX/ZY/ZPY を含むので Twiss[] はそのまま動く。
- `_lat.RFSW` の切替は calc6d 内で保存・復元している (ffscmd 側で触らなくてよい)。
- 不安定時 SAD は optstat%stabx/staby=.false. にする (マッチングの UnstableLevel 等)。`res.stable` を使う。
- INS (CELL OFF) + CALC6D は未実装 (qcell6d の else 分岐 = transport 伝搬。`calc6d` の不安定分岐と同じ手順で
  twiss0 を先頭 MARK の値にすれば作れる)。
- 対応要素は emit6d.tturne と同じ (DRIFT, BEND, QUAD, SEXT/OCT…, CAVI, MULT (ANGLE=0), 長さ 0 の MARK/MONI)。SOL 等は NotImplementedError。
