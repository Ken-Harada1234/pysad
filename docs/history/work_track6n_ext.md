# 作業ログ: track6n_ext.py (MULT とソレノイド領域の numba カーネル) — 2026-10-15

## 1. 目的と成果物

track6n.py (track6.py の numba 版) が未対応だった MULT (tmulti) とソレノイド領域 (tsol) を、
track6n.py と同じ流儀の numba カーネルにして **新規ファイル `pysad/track6n_ext.py`** に置いた
(track6n.py は別担当が編集中のため変更していない。他の既存ファイルも変更なし)。

- `pysad/track6n_ext.py` (約 600 行): カーネル `k_mult`, `k_sol_frin`, `k_sol_drift`, `k_sol_quad` と補助関数、
  `ProgExt` (track6n.Prog の派生)、`turn_ext` (track6n.turn + 新しい kind)、`track_particles_ext`
  (track6n.track_particles_fast と同じ仕様、フォールバックなし)。
- 照合: track6.track_particles (純 Python) と NORAD / RAD / FLUC で 20 粒子 × 5 周、座標差 ≤ 4.3e-17、FLUC の RNG.seed 一致 (§4)。
- 本書: 対応表 (§3)、照合結果 (§4)、track6n.py への組み込み案 (§5)、気付いた点 (§6)。

使い方 (暫定、組み込み前):
```python
from pysad.track6n_ext import track_particles_ext
out = track_particles_ext(elements, coords, nturn, momentum_eV, trf0=..., rad=True, fluc=True, seed=12345)
```

## 2. 方針

- **track6.py を直訳する** (SAD → track6 は Level 8/9/11 で SAD と ≤1.5e-10 の一致を確認済み)。演算の順序・括弧も
  track6 と同じにし、`x**2` は `(x*x)` に置き換える (gfortran と同様に pow(x,2)=x·x)。
- 配列は track6n と同じ S(n,10) = [x, px, y, py, z, g, dv, pxr0, pyr0, zr0]、`@njit(**NJ)` (error_model='numpy', cache)、
  `fmax/fmin` (NAN_MODE)、乱数は `rs` (int64[1]) を `tradk_all` に渡す。
- track6n の既存ヘルパーを import して再利用: `sqrt1, pxy2dpz, asinz, xsincos, xsincosh, drift1, rot1, setr0_1,
  tradk_core1, tradk_all, tzsetparam0, tsolqu1, tsolqur_all, ttfrin_n_1, ndivrad, k_apert`。
- 粒子外側/スライス内側の順序 (乱数の消費順、tsolqur の bzp の持ち越し) は track6 と同じ構造で書く。

## 3. 対応表

### 3.1 MULT (track6 `_tmulti_ens` ← SAD tmulti.f tmulti / tffs.f tpara)

| SAD | track6.py | track6n_ext.py | 備考 |
|---|---|---|---|
| tpara (akr, theta2, ndiv, dofr) | `_lat._mult_params`, `_tmulti_ens` 冒頭 | `compile_mult(el)` (compile 時) | akr は complex128 配列、dofr は bool 配列。ANGLE/F1/FRINGE/VOLT ≠ 0 は NotImplementedError (track6 と同じ)。nmmax<0 は K_DRIFT |
| tsolrot (ent/exit, bz=0) | `c[0]-=dx; _rot(c,ct,st)` | `k_mult` 内 `rot1` | |
| ndivrad | `_ndivrad` | `track6n.ndivrad` | krad のとき `max(ndiv, ndivrad(|akr0|, akr1, 0, EPS, h0))`, `min(1000, ndiv)` |
| ttfrins (複素 ak のフリンジ) | `_ttfrins` | `ttfrins1(i,S,nord,akre,akim,al)` | 主軸へ回して `track6n.ttfrin_n_1` (nord=2,4,6,一般)、戻す |
| tsolqu / tsolqum (K0 付き) | `_tsolqu(…, ak0x, ak0y, qum)` | `tsolqu_m1(i,S,al,ak,eps0,ak0x,ak0y,qum)` | track6n.tsolqu1 に dx0=ak0x/ak, dy0=ak0y/ak と ak=0 分岐を足したもの。qum (放射なし) は z の dv 項が −0.25(ap−dv)aln0 |
| tdrift (bz=0, ak0≠0; tsolconv) | `_tdrift_kick`, `_lat._tsolconv` | `tdrift_kick1`, `tsolconv` | ak=0 で K0 だけの MULT (厳密ドリフト)。放射は tradk で別に掛けるので enarad=False のみ |
| 多極キック (Horner) | `_tmulti_ens` 内の complex ループ | `mult_kick1(i,S,akrm,nmmax,nmin)` | nmin=2: 長さあり (n≥2, 最後に ·cx1/(1+δ))、nmin=1: L=0 の薄い多極 (+akrm[0]) |
| tradk (各分割の後と最後) | `_tradk` 全粒子 | `tradk_all` | 粒子ループの外で呼ぶので乱数の順は SAD/track6 と同じ |
| 全体 | `_tmulti_ens(parts, el, krad)` | `k_mult(S,n,al,nmmax,akr,dofr,ndiv0,dx,dy,ct,st,fringe,eps_rad,krad,fluc,rs,p0,h0,anrad,cuc)` | |

### 3.2 ソレノイド領域 (track6 `_tsol_region` ← SAD tsol.f tsol)

| SAD | track6.py | track6n_ext.py | 備考 |
|---|---|---|---|
| tfbzs | `_tfbzs(el, ph)` | compile 時に `_t6._tfbzs` | bzs = charge·(BZ+DBZ)/brhoz [1/m] |
| 領域の走査 (BOUND≠0 → 次の BOUND≠0) | `_tsol_region` | `ProgExt._compile_sol_region` | 始点 SOL → K_SOL_FRIN(bzs)、終点 SOL → K_SOL_FRIN(0−bzs)、途中の SOL → K_SOL_FRIN(bz1−bzs)。DISFRIN≠0 なら dbz=0。DX/DY/DZ/CHI/F1 ≠ 0、BOUND=0 の SOL、領域内の BEND/MULT/CAVI は NotImplementedError (track6 と同じ) |
| tsfrin | `_tsfrin(p, dbz)` | `tsfrin1`, `k_sol_frin(S,n,dbz)` | |
| tdrift_solenoid (+tradk1) | `_tdrift_solenoid(p, al, bz, enarad)` | `tdrift_solenoid1`, `k_sol_drift(S,n,al,bz,krad,…)` | 放射の基準は (pxi, pyi, zi) を S[i,7..9] に置いて `tradk_core1` (setr0 はしない: track6 と同じ)。bz=0 は drift1 |
| tquad (bz≠0): tsolrot 回転 → ttfrin (bz) → tsolqu/tsolqur → 戻し | `_tquad_sol` | `k_sol_quad(S,n,L,ak,ct,st,fringe,bz,krad,…)` | DX/DY/F1/FRINGE ≠ 0, L=0, K1=0 は NotImplementedError。bz=0 なら track6n の `tsolqur_all` / `tsolqu1` (eps0=0: track6 と同じく EPS は使わない) |
| ttfrin (nord=4, bz≠0) | `_ttfrin4_bz` | `ttfrin4_bz1` | |
| tzsetparams + tsolqu 本体 (bz≠0) | `_TZ.__init__`, `_tsolqu_bz_ens` のループ | `tz_step1(i,S,aln0,akk,bz,first)` → bzp を返す | _TZ の 30 個ほどの係数は関数内のローカル変数 |
| tztafs | `_TZ.tztafs` | `tztafs(…, 13 個の係数)` | |
| 最後の半ステップ | `_tsolqu_bz_final` | `tsolqu_bz_final1(i,S,aln0,bzp)` | **tsolqur (krad) では `bzp_last` = 最後のスライスの最後の粒子の bzp を全粒子に使う** (SAD tsolqu.f:519-531 の癖)。放射なしは粒子ごとの bz/(1+δ) |
| 全体 (bz≠0) | `_tsolqu_bz_ens(parts, al, ak, bz, krad)` | `tsolqu_bz_all(S,n,al,ak,bz,eps0,krad,…)` | スライス外側・粒子内側、各スライスの後に `tradk_all` |

### 3.3 kind 番号と params

| kind | 値 | params | sol_mode |
|---|---|---|---|
| K_MULT | 10 | (al, nmmax, akr[complex128], dofr[bool], ndiv0, dx, dy, ct, st, fringe, EPS, radok) | 0 |
| K_SOL_FRIN | 11 | (dbz,) | 1 (始点) / 2 (途中の SOL) / 3 (終点) |
| K_SOL_DRIFT | 12 | (al, bzs, radok) | 2 / 3 |
| K_SOL_QUAD | 13 | (L, ak, ct, st, fringe, bzs, radok) | 2 / 3 |

`ProgExt.sol_mode` (SM_NONE=0, SM_START=1, SM_IN=2, SM_END=3) は tapert の la の扱い (§5.3) のために持つ。

## 4. 照合結果

テストラティス (Python で組む、`scratchpad/test_ext.py`): Level 11 の cell (2.5 GeV) をもとに 1 cell だけ。
- 共通: MARK, CAVI (L=0, HARM=12), MULT MQF (L=0.3, K1L=0.54), DRIFT, MULT MSF (L=0.2, K0=2e-4, K2L=0.6), DRIFT, BEND (0.8 m, 9°, E1=E2), MULT MK0 (L=0.1, K0=1e-4: ak=0 → tdrift_kick1), **MULT MTH (L=0, K0+K2: 薄い多極)**, **MULT MS2 (L=0.2, K2L=6: tpara の ndiv=2)**, DRIFT
- 'mult' 版 (18 要素): + DRIFT, MULT MQD (L=0.3, K1L=−0.51 (akang 90°), **K3L=50 (dofr n=3 のフリンジ)**, DX=1e-4), DRIFT, MULT MSD (K2L=−1.2, DY=1e-4, ROTATION=1e-3), DRIFT, BEND, DRIFT
- 'sol' 版 (20 要素): + **SOL SB (BZ=2 T, BOUND=1), DRIFT 0.1, QUAD QD (K1=−1.7), SOL SE (BOUND=1)**, DRIFT, MULT MSD, DRIFT, BEND, DRIFT

20 粒子 (x,y ~ 3e-4, px,py ~ 1e-4, z ~ 1e-3, δ ~ 5e-4 の正規乱数) を 5 周、trf0=1e-4、seed=12345:

| ラティス | モード | max\|Δ\| (track6 − ext) | RNG.seed (track6 / ext) |
|---|---|---|---|
| mult | NORAD | 3.1e-17 | 12345 / 12345 |
| mult | RAD | 4.0e-17 | 12345 / 12345 (乱数を引かない) |
| mult | FLUC | 4.3e-17 | −2025222455 / −2025222455 一致 |
| sol | NORAD | 2.0e-17 | 12345 / 12345 |
| sol | RAD | 3.1e-17 | 12345 / 12345 |
| sol | FLUC | 3.0e-17 | −1861977247 / −1861977247 一致 |

- いずれも 1e-13 以下 (ほぼ丸め誤差の範囲: 複素数の積・商や numba の最適化の違い)。ext の 2 回実行は bit 一致。
- tsolqur の癖の確認 (`scratchpad/test_quirk.py`): 'sol' 版で粒子の並びを逆にすると RAD では結果が 5.3e-7 変わり、NORAD では 0 (track6 と同じ挙動、Level 9 で SAD と確認済みのもの)。
- apert 付き `turn_ext` (§5.3 の la 処理込み) の動作確認: 損失なしで n 不変、x=5 m の粒子を 1 つ入れると n が 1 減り lost_turn に周番号が入る。
- 速度 (2000 粒子 × 50 周、i7 相当): mult 版 NORAD 86 ns/粒子/要素、RAD 287、FLUC 407; sol 版 99 / 306 / 356。
  初回の JIT コンパイルは約 3 s (cache=True で 2 回目以降は 0.05-0.4 s)。

## 5. track6n.py への組み込み案

### 5.1 kind 番号

```python
K_DRIFT, K_QUAD, K_THIN, K_BEND, K_CAV, K_STEER0, K_NONE = 0, 1, 2, 3, 4, 5, 6
K_MULT, K_SOL_FRIN, K_SOL_DRIFT, K_SOL_QUAD = 7, 8, 9, 10      # track6n_ext では 10..13 (衝突回避のため)
```
組み込み時は track6n_ext のカーネル群 (tsolconv … k_sol_quad) をそのまま track6n.py に移し、`from .track6n import …` を外す。

### 5.2 Prog の差分案

```python
 class Prog:
     def __init__(self, elements, ph):
-        for el in elements:
-            k, p = self._compile(el, ph)
-            ...
+        self.sol_mode = []
+        i = 0
+        while i < len(elements):
+            el = elements[i]
+            if type(el).__name__ == 'Solenoid':           # ソレノイド領域 (BOUND≠0 … 次の BOUND≠0)
+                if _t6._g(el, 'BOUND') == 0.0:
+                    raise NotImplementedError('track6n: BOUND=0 の SOL (領域外)')
+                i = self._compile_sol_region(elements, i, ph) + 1   # track6n_ext.ProgExt._compile_sol_region
+                continue
+            k, p = self._compile(el, ph)
+            ... (従来どおり; la_count は MARK/MONI と DRIFT 以外 True)
+            self.sol_mode.append(0)
+            i += 1

     @staticmethod
     def _compile(el, ph):
         ...
+        if tn == 'Multipole':
+            return compile_mult(el)                        # track6n_ext.compile_mult (nmmax<0 は K_DRIFT)
         if tn in ('Marker', 'Monitor') or el.L == 0.0:
             return K_NONE, ()
```
注意: 現在の `_compile` は `el.L == 0.0` の要素を黙って K_NONE にするので、**L=0 の MULT (薄い多極) と SOL は
NotImplementedError にならず無視される** (track6 へのフォールバックが効かない: §6.2)。Multipole/Solenoid の判定は
この行より前に置くこと。

### 5.3 turn() の差分案

```python
 def turn(prog, S, n, ph, trf0, rfsw, rad, fluc, rs, apert=None):
     ...
-    for kind, pr, cnt in zip(prog.kinds, prog.params, prog.la_count):
-        if apert is not None:
-            if la <= 0:
-                n = k_apert(...); if n <= 0: return 0
-                la = 15
-            if cnt:
-                la -= 1
+    nreg = 0
+    for kind, pr, cnt, sm in zip(prog.kinds, prog.params, prog.la_count, prog.sol_mode):
+        if apert is not None:
+            if sm == SM_NONE or sm == SM_START:            # tturn1 の主ループ
+                if la <= 0:
+                    n = k_apert(...); if n <= 0: return 0
+                    la = 15
+                if cnt and sm == SM_NONE:
+                    la -= 1
+                if sm == SM_START:
+                    nreg = 0
+            else:                                          # tsol 内のループ (tsol.f:604-615)
+                nreg += 1
+                if la <= 0:
+                    n = k_apert(...); if n <= 0: return 0
+                    la = 15
+                else:
+                    la -= 1
         if kind == K_DRIFT:
             ...
+        elif kind == K_MULT:
+            (al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad, radok) = pr
+            k_mult(S, n, al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad,
+                   rad and radok, fluc, rs, p0, h0, anrad, cuc)
+        elif kind == K_SOL_FRIN:
+            k_sol_frin(S, n, pr[0])
+        elif kind == K_SOL_DRIFT:
+            (al, bz, radok) = pr
+            k_sol_drift(S, n, al, bz, rad and radok, fluc, rs, p0, anrad, cuc)
+        elif kind == K_SOL_QUAD:
+            (L, ak, ct, st, fringe, bz, radok) = pr
+            k_sol_quad(S, n, L, ak, ct, st, fringe, bz, rad and radok, fluc, rs, p0, h0, anrad, cuc)
+        if apert is not None and sm == SM_END:             # 主ループの後処理 (tturn.f:240-256, 553)
+            la -= 1                                        # 始点 SOL の 1020
+            for _ in range(nreg):                          # l=k+1..ke: sol=l<ke; go to 1020
+                if la <= 0:
+                    n = k_apert(...); if n <= 0: return 0
+                    la = 15
+                la -= 1
```
la の扱いの根拠 (SAD tturn.f / tsol.f):
- 主ループ (tturn1) は始点 SOL k の前で `la<=0` なら tapert、`case(icSOL)` で tsol を呼んだ後 1020 で la−1。
  その後 l=k+1..ke も主ループを通り (`if(sol) … sol=l<ke; go to 1020`)、各 l の先頭で `la<=0` なら tapert (座標は領域通過後のもの)、1020 で la−1。
- tsol の内側ループ (l=k+1..ke) は各要素の前で「`la<=0` なら tapert, la=15、そうでなければ la−1」(tsol.f:604-615)。
  要素の種類によらず数える (主ループの DRIFT/MARK の例外はない)。
- `track6n_ext.turn_ext` はこのとおりに実装してある (座標には影響しない。損失周の番号付けのみ)。

### 5.4 track_particles_fast

`Prog(elements, ph)` が MULT/SOL で NotImplementedError を出さなくなるので、フォールバックは ANGLE≠0 の MULT、SOL の DX/DY/CHI/F1、
BOUND=0 の SOL、領域内の BEND/MULT/CAVI などに限られる (これらは track6 も未対応なので、実質フォールバックは
L=0 四極や ANGLE=0 K1≠0 BEND 等の既存の分だけ)。

## 6. 気付いた点

### 6.1 track6.py: L=0 の MULT で n>2 の K_n があると ZeroDivisionError
`_tmulti_ens` の `dofr = [... abs(akr[n]) / al * 0.05**(n+1) ...]` は al=0 で n>2 のとき 0 除算になる
(dofr は nz のときしか使わないが、内包表記は先に評価される)。track6n_ext の `compile_mult` では al=0 なら dofr を全部 False にした。
track6.py 側は `n <= 2 or al == 0.0 or ...` を足せば直る (別担当の範囲なので未修正)。テストでは薄い MULT を K0+K2 にして回避した。

### 6.2 track6n.py: L=0 の要素を黙って K_NONE にする (`_compile` 末尾)
`if tn in ('Marker','Monitor') or el.L == 0.0: return K_NONE` のため、L=0 の MULT (薄い多極キック) と SOL (L=0) は
エラーにならず無視され、`track_particles_fast` は track6 にフォールバックせずに **MULT のキックなし・ソレノイドなし** で追跡する
(確認: MULT を含まない SOL 付きの cell (QUAD/SEXT/BEND/DRIFT + SOL BZ=2 T) を track_particles_fast に渡すと
track6/ext と 6.2e-4 ずれる = ソレノイドなしの結果。MULT (L≠0) が一緒にあるときだけ NotImplementedError で track6 に落ちて正しくなる)。§5.2 のとおり Multipole/Solenoid の判定をこの行より前に置くか、
この行を `tn in ('Marker','Monitor')` だけにして L=0 の未知要素は NotImplementedError にするのがよい。

### 6.3 ソレノイド内の四極の EPS
track6 `_tquad_sol` は `_tsolqur_ens(parts, L, ak)` / `_tsolqu_bz_ens(…)` を eps0=0 で呼ぶ (SAD の tquad は cmp の EPS を渡す)。
ext も track6 のまま eps0=0 にした。EPS≠1 のソレノイド内四極は SAD と分割数が変わり得る (track6 と同じ制限、Level 9 は EPS 既定なので影響なし)。

### 6.4 複素数演算の一致
`mult_kick1` の `cx / float(n)` は Python の complex/int (Smith の除算で虚部 0 → 成分ごとの除算と同じ) と一致、
複素積は両者とも (ac−bd, ad+bc)。差は 1e-17 程度に収まっている (§4)。

## 7. ファイル
- `pysad/track6n_ext.py` — 本体 (新規)
- `docs/history/work_track6n_ext.md` — 本書
- 検証スクリプト (scratchpad、成果物ではない): `test_ext.py` (照合)、`test_quirk.py` (並び依存)、`test_perf.py` (速度・apert)
