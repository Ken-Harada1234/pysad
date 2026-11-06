# 作業ログ: 要素幾何 (GX/GY/GZ/CHI1-3) の移植 — pysad/geo.py (2026-10-15)

対象: SAD K64 `src/tfgeo.f` (tfgeo, tfgeo1), `src/tgrot.f` (module geolib: geoini, tfchi, tfchitogeo, tfrotgeo,
tfderotgeo, tgrot, tforbitgeo, tfgeofrac), `src/tfloor.f` (xsincos, asinz), `src/qtwiss.f` (qfraccomp),
`src/tftwiss.f` (tfline の 'GX'..'GCHI3', 'GEO', 'OGX'.., ktfgeol), `src/tffsa.f` (ORG コマンド)。
新規ファイルは `pysad/geo.py` だけ。既存ファイル (ffscmd.py, parser.py, lattice.py) は触っていない。

## 1. SAD 側の構造

- `geo(3,4,nlat)` (ffs_pointer): 要素 k の**入口**の幾何。`geo(:,1,k)` = ビーム座標 x 軸、`geo(:,2,k)` = y 軸、
  `geo(:,3,k)` = z 軸 (進行方向) の実験室系での単位ベクトル (方向余弦)、`geo(:,4,k)` = 入口の位置 (GX, GY, GZ)。
  `pos(nlat)` = 入口の経路長 s。nlat = 要素数 + 1 で末尾が `$$$`。
- 初期値 `geoini` (tgrot.f):  x 軸 = (0,−1,0), y 軸 = (0,0,−1), z 軸 = (1,0,0), 位置 = 0。
  つまり実験室系 X が進行方向、ビームの x は −Y、y は −Z。ESRF では正の ANGLE で GY が増える。
- 角度 (tgrot.f `tfchi(geo,i)`): z 軸 = (zx, zy, zz), x 軸の第 3 成分 xz, y 軸の第 3 成分 yz として
  | キー | 定義 | 備考 |
  |---|---|---|
  | CHI1 | atan2(zy, zx) | zy = 0 なら 0。水平方位角。ESRF では 1 セルで π/10 |
  | CHI2 | asinz(zz) | 仰角。asinz は [−1,1] に切った asin (小さい値は級数) |
  | CHI3 | atan2(xz, −yz) | ロール角。xz = 0 なら 0 |
- **FFS からの参照は `LINE[key, loc]`**: key = `"GX" "GY" "GZ" "GCHI1" "GCHI2" "GCHI3"` (スカラー)、`"GEO"`
  ({{gx,gy,gz},{chi1,chi2,chi3}})、閉軌道を乗せた `"OGX".."OCHI3"`, `"OGEO"` (tforbitgeo)。
  `Twiss["GX",loc]` は K64 では `???General::wrongval: #1 (GX) is must be to be name of optical function`
  になる (tftwiss.f の nlist 検索は mfito までで GX は mfitgx > mfito; Twiss 配列も ntwissfun 個しかない)。
  "CHI1" 等の名前は FIT 条件 (`CHI1 0.1` など、tdgeo.f で微分を作る) と `ORG` の引数、DISP 表示でだけ使われる。
- loc が nlat を超える番号 (例 4000, nlat=3564) は SAD では **k − nlat** (= 436) の要素になる (tflinestk の周回)。
  小数付き (200.5) は要素 200 の入口から 0.5 の部分要素 (tfgeofrac → qfracseg/qfraccomp → tfgeo1s)。

## 2. 計算 (tfgeo1) と geo.py の対応表

| SAD (Fortran) | geo.py | 内容 |
|---|---|---|
| `geoini` | `GEOINI` | 初期幾何 3×4 (列優先 reshape を展開) |
| `xsincos(x,s,xs,c,dc)` | `xsincos(x) -> (s, xs, c, dc)` | 小角でも桁落ちしない sin/cos (しきい値 xthcs1/2 も同じ式) |
| `asinz` | `asinz` | |
| `tfchi(geo,i)` | `tfchi(geo,i)` | CHI1..3 |
| `tfchitogeo / tfchitorot` | 同名 | chi → 軸 / 回転行列 |
| `tfrotgeo / tfderotgeo` | 同名 | 軸を chi で回す / 戻す (ORG, SOL の GEO で使用) |
| `tgrot(geo0,geo1)` | 同名 | geo0 から見た geo1 の (chi1,chi2,chi3) (SOL の CHI1-3 の書き戻しで使用) |
| `tforbitgeo(geo,cod)` | 同名 | OGX.. 用 |
| `tfgeo1` の要素ループ | `compute_geo` 内ループ + `_advance` + `_geo_length` | 下表 |
| `tfgeo` の MARK OFFSET 処理 | `compute_geo` 後半 | `xp = clamp(offset + i, 1, nlat)`、整数ならコピー、小数なら部分要素 |
| `tfgeo` の ORG 変換 | `compute_geo` 末尾 | `dgo`, `geo2` (§4) |
| `tfgeofrac(l,fr)` | `GeoResult.geo_frac(k,fr)` | 部分要素の幾何 (§3) |
| tfline 'GX'..'GEO' | `GeoResult.line_value(key,k,fr)` | FFS の `LINE[]` の値 |
| tfline 'OGX'..'OGEO' | `GeoResult.orbit_value(key,k,cod,fr)` | 閉軌道付き |
| `tsgeo` (SOL の BOUND/GEO) | **未移植** (直進扱い、警告) | ソレノイド内の軌道追跡 (tsgeo1: tdrift/tquad/tmulte) が要る |
| `tgmap` (MAP 要素) | 未移植 (直進扱い) | |

要素ごとの進め方 (`_advance`; geo_i = 入口, g1 = 出口):

| 要素 | 式 (tfgeo1 の直訳) |
|---|---|
| ANGLE キーを持つ型 (BEND, MULT) で ANGLE=φ≠0 | θ=ROTATE。θ≠0 なら x' = cosθ x − sinθ y, y' = sinθ x + cosθ y、θ=0 なら x'=x。ρ=L_geo/φ, (sp,·,cp,dc)=xsincos(φ), r1=ρ sp, r2=−ρ dc。位置' = 位置 + r1 z − r2 x'、x'' = cp x' + sp z、z' = −sp x' + cp z。θ≠0: y_new = −sinθ x'' + cosθ y', x_new = cosθ x'' + sinθ y'。θ=0: y_new = y |
| 同型で ANGLE=0 (ESRF の HD/VD: ROTATE=−90° でも) | 直進: 位置' = 位置 + z·L_geo、軸不変 |
| COORD | CHI1..3, DX/DY/DZ, DIR の回転・並進 (`_coord`; 行列 r11..r33 をそのまま写した) |
| MARK (OFFSET≠0) | 本ループでは直進 (L=0)、後処理で `geo(:,:,i)=geo(:,:,xp)` |
| その他 (DRIFT, QUAD, SEXT, CAVI, MONI, MARK, ...) | 直進 |

`L_geo` (`_geo_length`) = L。ただし BEND/MULT で `FRINGE≠0` かつ (fb1,fb2)≠0 のとき
`L −= ((φ fb1)² + (φ fb2)²)/L/48 · sin(½(φ(1−E1−E2) − AE1 − AE2)) / sin(½φ)`
(BEND: fb1 = F1+FB1, fb2 = F1+FB2; MULT: FB1, FB2)。この L_geo が `pos(i+1) = pos(i) + L_geo` にも使われる
(= `LINE["S"]` / `LINE["LENG"]`)。E1/E2 は SAD の意味 (ANGLE 比)。

要素値の取り方 (`_val(el, key)`): SAD のキーワード名で `el._raw` (parser が残した元の数値) → 属性の順。
parser が ANGLE=0 の BEND を Drift にしたものは `_sad_type='BEND'` で元の型を見る。
pysad の `Bend.E1/E2` は rad 値 (E1·ANGLE+AE1) なので、`_raw` に E1/E2 がない場合だけ ANGLE で割って比に戻す。

## 3. 部分要素 (tfgeofrac / qfraccomp, rx1=0, rx2=fr → f1=1, f2=0, r=fr)

- L·r, ANGLE·r, **E1 → E1/r** (比なので絶対端角は不変), E2 → 0, AE1/AE2 はそのまま。
- BEND: FRINGE=0 なら F1=0、FRINGE≠0 なら FRINGE=−f1−2f2=−1。FB1/FB2 不変。
- MULT (ANGLE≠0): FRINGE=0 なら FB1=FB2=0、≠0 なら FRINGE=−1。
- ROTATE は不変。L=0 かつ ANGLE=0 の要素は変化なし (chg=.false. → 入口の値)。
- 逆向き要素 (cmp%ori=.false.) では E1·f2/r=0, E2·f1/r=E2/r になるが、pysad の要素は向きを持たないので未対応
  (parser は逆向き BEND の E1/E2 属性を入れ替えるが `_raw` は元のまま)。
- 最初は E1 をそのままにしていて、19B の B4 (FRINGE=1 F1=0.3 E1=0.2 E2=0.3 AE1=0.01) の 11.37 で 8e-6 ずれた。
  qfraccomp を読み直して E1/r に直したら 1e-15 で一致。

## 4. ORG (tffsa.f `ORI_GIN`/`ORG`) と tfgeo 末尾の変換

`ORG 要素 gx gy gz chi1 chi2 chi3` (chi は度; `scale(mfitchi1..3)=π/180`)。要素名なしの `ORG` は
iorgr=1, geo0=geoini にリセット。`geo0(:,1:3) = tfrotgeo(geoini(:,1:3), −chi0·π/180)`、`geo0(:,4) = (gx,gy,gz)`。
tfgeo の末尾 (calgeo のとき):

```
dgo = geo(:,4,iorg) + gx0·z_iorg − gy0·x_iorg − gz0·y_iorg
geo(:,4,:) −= dgo
geo2 = geo0(:,1:3) · geo(:,1:3,iorg)^T          (iorg の軸系 → geo0 の軸系)
geo(:,1:4,k) = geo2 · geo(:,1:4,k)   for all k
```
→ 既定 (iorg=1, chi0=0, g0=0) では恒等。chi0=0 なら要素 iorg の位置は (−gx, −gy, −gz) になる
(原点の方を要素から見た位置で与える形)。`SORG` フラグなら pos も pos(iorg) を引く (未実装; pos は常に先頭基準)。
`compute_geo(..., gx0, gy0, gz0, chi0=(rad), iorg=k)` に対応。

**SAD の getvad の注意**: `ORG B3 1.5 -2.5 0.5 10 20 30` は `1.5 -2.5` を式 1.5−2.5=−1.0 と読むので
gx=−1.0, gy=0.5, gz=10, chi1=20, chi2=30, chi3 欠落 → `?ORG-Missing parameters` が出る
(欠落分は 0 になった; 直前の ORG の値は残らなかった)。負の数は `0-2.5` のように続けて書くか別行にすること。

## 5. フラグ GEOCAL / GEOFIX (tffsa.f, tffscalc.f)

- 既定 GEOCAL (flags(24)=geocal、pysad ffscmd の `_FLAGS` にも 'GEOCAL': True がある)。
- `tfgeo(calgeo0)`: GEOFIX (geocal=.false.) のときは p0 が変わった場合に gammab だけ更新して戻る (幾何は触らない)。
  USE 直後 (chguse) は GEOFIX でも一度だけ GEOCAL 扱いで計算する (tffsa.f:114-118)。
- CALC (tffscalc.f:79): `tfgeo(geomet .or. .not. fitflg)` — CALC (fitflg=.false.) または GX..CHI3 の FIT 条件がある
  (geomet) ときだけ軸を進め直す (calgeo)。それ以外 (GO) は pos と gammab のみ。ただし SOL (idtype 20) が 1 つでもあれば
  常に calgeo。STOP/EXIT/QUIT でも `tfgeo(.false.)`。EMIT 系は `tfgeo(.true.)`。
- pysad では幾何は要素値の関数なので「CALC のたびに compute_geo し直す (GEOCAL) / USE 時の結果を保持 (GEOFIX)」で足りる。

## 6. 照合結果

### 6.1 Level 19A — ESRF (bench_level19a.sad → bench_level19a.result, bench_level19a_pysad.py → bench_level19a_pysad.result)

pysad/tests/kek-ls_v03_12.lat + `USE RING; CELL; CALC;` で `LINE[{GX,GY,GZ,GCHI1,GCHI2,GCHI3}, loc]` を 21 箇所
(`$$$`, QFA1H.3, BPI6.3, VD.1, HD.5, SFA.5, MKC.2, RFC, QDA2.7, BID5.20, BPI7.13, 2, 100, 500, 1000, 2000, 3000, 4000,
200.5, 1500.25, `^^^`)、`LINE["GEO",...]` 2 箇所、`LINE["S",locs]`。nlat=3564。

| 項目 | 最大 |差| |
|---|---|
| GX, GY, GZ (21 箇所) | 4.8e-13 (値 ~190 m) |
| GCHI1..3 | 2.6e-13 |
| LINE["GEO"] | 4.2e-13 |
| LINE["S"] (pos) | 1.1e-13 |
| **全体** | **4.8e-13 (目標 1e-10 達成)** |

周回の確認: `$$$` で GX=−1.513e-7, GY=−5e-14, CHI1=−1.58e-9 (SAD) と同じ値 (角度の和が 2π から 1.6e-9 ずれている
ラティスそのものの性質を再現)。loc=4000 は SAD の周回規則 (k−nlat=436, QFA8H.9) に合わせた。ESRF は平面リングなので
GZ, CHI2, CHI3 は全て 0 → 19B で分岐を別途確認。

### 6.2 Level 19B — 分岐テスト (bench_level19b.sad / .result, bench_level19b_pysad.py / .result)

小さな線 `M0 D1 B1 D2 QF B2 M1 D1 B3 D2 B4 M2 D1 B5 QF D2 ME`:
B1 (ANGLE 0.3), B2 (ROTATE 30°), B3 (ANGLE −0.15, ROTATE 90° = 垂直偏向), B4 (FRINGE=1 F1=0.3 E1=0.2 E2=0.3 AE1=0.01),
B5 (ROTATE −60°, FRINGE=1 FB1=0.2 FB2=0.1 F1=0.05), M1 (OFFSET=2.3), M2 (OFFSET=−1.5)。
整数位置 1..18 と 1.37..17.37 の 35 箇所 × 6 値 + pos。

| 項目 | 最大 |差| |
|---|---|
| 整数位置 (ROTATE, 垂直偏向, OFFSET 含む) | 5.3e-15 |
| 小数位置 (部分要素; B4/B5 のフリンジ補正を含む) | 5.3e-15 |
| pos (フリンジ長補正を含む) | 5.3e-15 |
| ORG B3 (gx=−1, gy=0.5, gz=10, chi=(20°,30°,0)) | 8.0e-14 |

### 6.3 Level 19C — ORG の各引数 (bench_level19c.sad / .result, bench_level19c_pysad.py / .result)

gx/gy/gz のみ、chi1 のみ (10°)、chi2 のみ (20°)、chi3 のみ (30°)、全部、別要素 (D2.1) の 6 通り × 18 点: 最大 |差| 9.1e-14。

## 7. ffscmd.fn_twiss / fn_line への接続案 (既存ファイルは未変更)

1. `FFS.use()` の最後 (または `calc()` の中) で幾何を作り保持する:
   ```python
   from . import geo as _geo
   self.geo = None                       # use() で None に
   def ensure_geo(self):
       if self.geo is None or self.flag('GEOCAL'):      # GEOFIX なら USE 時の結果を保持
           self.geo = _geo.compute_geo(self.els, gx0=self.org_g[0], gy0=self.org_g[1], gz0=self.org_g[2],
                                       chi0=self.org_chi, iorg=self.org_k)
       return self.geo
   ```
   `invalidate()` (要素値の変更) で `self.geo = None`。ORG コマンドは `_exec_word` に
   `ORI_GIN`/`ORG`: 要素名 → `org_k` (なければリセット: org_k=1, org_g=0, org_chi=0)、続く 6 値 (chi は度 → rad) を追加。
2. `line_value()` に分岐を追加:
   ```python
   if key in ('GX', 'GY', 'GZ', 'GCHI1', 'GCHI2', 'GCHI3', 'GEO'):
       return self.ensure_geo().line_value(key, k, fr)
   if key in ('OGX', 'OGY', 'OGZ', 'OCHI1', 'OCHI2', 'OCHI3', 'OGEO'):
       r = self.row_frac(k, fr) if fr else self.row(k)
       return self.ensure_geo().orbit_value(key, k, (r.dx, r.dpx, r.dy, r.dpy), fr)
   ```
   `LINE["S"/"LENG"]` も `ensure_geo().pos_frac(k, fr)` にすると SAD のフリンジ長補正付き pos と一致する
   (今の実装は TwissRow.s の線形内挿)。
3. `twiss_value()` の `'GX'`/`'GY'` は現在 γ=(1+α²)/β を返しているが、SAD では γ は `GMX`/`GMY` で、
   `Twiss["GX"]` 自体はエラー。SAD に合わせるなら 'GMX','GMY','GMZ' を γ にし、'GX'..'CHI3' は
   `SADError('... is must be to be name of optical function')` にするか、互換のため geo の値を返す (どちらにするかは要判断)。
4. FIT 条件 `GX/GY/GZ/CHI1/CHI2/CHI3 値` (tdgeo.f の微分、重み tffscalc.f:488-493: GX..GZ は
   1e-3·sqrt(max(NX,NY)/em/L)、CHI は 1e3) は match6 側の仕事で、今回は対象外。
5. `Element["GEO","SOL名"]`, SOL の CHI1-3 書き戻し (tsgeo) は未移植。SOL を含むラティスでは警告を出す。

## 8. 未移植・注意点のまとめ

- `tsgeo` (SOL の BOUND/GEO: ソレノイド中の軌道追跡で幾何と SOL の DX/DY/DZ/CHI1-3 を決める) と `tgmap` (MAP) は未移植。
- `SORG` フラグ (pos を ORG 要素基準にする) は未実装。
- 逆向き要素の部分要素 (E1/E2 の入れ替え) は未対応 (pysad の要素に向き情報がない)。
- `ORG` の負数は getvad の式解釈に注意 (§4)。
- `LINE[..., 4000]` のような範囲外は呼び出し側 (ffscmd._resolve_loc) で k−nlat に周回させる必要がある。

## 9. ファイル

- `pysad/geo.py` — 新規 (本体)
- `bench_level19a.sad` / `.result` / `bench_level19a_pysad.py` / `bench_level19a_pysad.result` — ESRF 照合
- `bench_level19b.sad` / `.result` / `bench_level19b_pysad.py` / `bench_level19b_pysad.result` — 分岐テスト
- `bench_level19c.sad` / `.result` / `bench_level19c_pysad.py` / `bench_level19c_pysad.result` — ORG
- SAD 実行: `wsl -d Ubuntu -- bash -lc 'cd "<repo>" && <SAD>/bin/gs bench_level19a.sad > bench_level19a.result'`
- pysad 実行: `PYTHONIOENCODING=utf-8 python -u bench_level19a_pysad.py`

## 10. 作業ログ

1. tfgeo.f / tsgeo.f / tgrot.f / tfloor.f (xsincos, asinz) / tftwiss.f (tfline の GX..) / tffsa.f (ORG) / tffscalc.f / qtwiss.f (qfraccomp) を読み、
   `Twiss["GX"]` が K64 では使えないことを SAD 実行で確認 → ベンチは `LINE[]` に変更。
2. geo.py を書き、19A (ESRF) で 4.8e-13 一致。loc=4000 は SAD の周回 (k−nlat) と判明。
3. 19B で ROTATE/垂直偏向/フリンジ長補正/OFFSET/部分要素を確認。部分要素の E1 の扱い (E1/r) を qfraccomp に合わせて修正。
4. ORG の引数解釈 (式として読まれる負数、欠落 chi3=0) を 19C で切り分け、6 通り全て 1e-13 で一致。
