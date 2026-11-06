# 作業ログ: BEND の K0 (閉軌道・分散・追跡・Emittance) と LINE["E1"] (2026-10-15)

## 1. 結論

pf2 スクリプト (bump2/disp_response/disp_kick_effect) で応答が 0 だった原因は 1 つではなく 4 つだった。

1. **ANGLE=0 の BEND (K1..K5, KK1 など) の COD 写像が無かった** (最大の原因)。`lattice._one_turn_cod_matrix` の K1=0 の BEND の分岐は `_bend_cod_update` に行き、ANGLE=0 だと COD を変えずに返し、行列は数値微分で恒等になる。K0 を無視するだけでなく**軌道があるとドリフトの長さも落ちていた**。SAD は `tbende.f` の phi0==0 の分岐 (tsteee / tthine / tquade / tmulte)。→ `_tsteee_cod_matrix` を追加 (解析行列 `_tsteee_body`)、emit6d に `_tsteee`、追跡に `track6._tsteer_ens`。
2. **`LINE["K0","BPI6.1"]=v` の意味**。SAD は読み出しでは 1 台だけ変わるが、次の CALC/Emittance[]/TrackParticles で**ファミリー全体が代表 (.1) の値に揃う** (`.3` だけ設定すると CALC 後に 0 に戻る。bench で確認)。pysad は 1 台だけだった。→ `FFS._spread_pending` / `spread_pending()` (calc, fn_emittance, fn_track, fn_dasurvey の先頭)。ドット無しの `LINE["K0","BPI6"]`、FFS の `BPI6 K0 v` は即座にファミリー全体。`BPI6.1 K0 v;` は SAD では式として評価され何も起きない。
3. **`LINE["E1"]`** が rad だった。→ `FFS._e_ratio` (ANGLE 比。ANGLE=0 は `_raw` の値、逆向き要素は `_rev` で E1/E2 を戻す)。`LINE["ANGLE",..]=` で ANGLE を変えても比を保つ。`LINE["E1",..]=` は属性 (rad) と `_raw` を更新。
4. **追跡**: `TrackParticles[{1,beam},1,1,N]` の N 周を無視していた (`ffscmd.fn_track` は args[1] を周回数と読んでいた) → nbegin..nend から mt=nend−nbegin+1。ANGLE=0 の薄いステアリング (L=0) の符号が逆 (SAD: `px −= K0/(1+δ)`、pysad は +)。L>0 の tsteer は未実装 (NotImplementedError)。

ANGLE≠0 の BEND の K0 (phib=ANGLE+K0) は CALC・emit・追跡とも元から正しかった (bench で確認)。

## 2. 見つけて直した付随バグ

- **`lattice._tsolque_body_cod` の z 写像の符号**: SAD の `xsincos` は xs = x − sin x (`xsinch` は x − sinh x) なのに `xs1 = s1 − phi1`, `xsh2 = sh2 − phi1` と逆だった。x²/y² の項が半分になり、軌道がずれた四極の z (経路長) が理想値 (RK 積分、track6 と一致) の約半分だった。K0 で軌道が mm 級になる/δ が大きい場合の Emittance[] の δ・U0・αc に効く (23c: U0 の差 3e-7→7e-11)。
- `$$$` 行の COD: SAD は最後に評価した一周の出口 COD (入口 cod0 と conv 分ずれる: 軌道が 1e-5 で相対 1e-4)。pysad は入口だった → `_LAST_COD_EXIT`、`twiss_table` の最終行の x,px,y,py を出口にした (23d/23f/23h の `$$$` が 1e-10 で一致)。
- `emit6d/calc6d` の z を閉じる δ の割線法の反復を 4→30 回 (K0 の総和が 0 でないとき δ が 1e-3 動く。実際の原因は上の spread と z の符号で、反復回数の影響は小さい)。
- `ffscmd.get_value` の E1/E2、`parser._rev` の追加、docs/02 の「K0 が属性に渡らない」は誤り (訂正)。

## 3. 参照ベンチ (SAD K64 で実行し .result を保存、pysad は `python run_sad.py bench_level23X.sad`)

比較は `tmp23/cmp3.py` (Print の数値行を順に突き合わせ)。差は絶対差。

| ベンチ | 内容 | 結果 |
|---|---|---|
| 23a | esrf `LINE["K0","BPI6.1"]=1e-4` (K1≠0 BEND)、`BPI2.1`=1.5e-4 (K1=0)、E1/ANGLE/K0 の読み出し (ファミリー展開) | DX,DPX,EX,EPX,NX,BX が 3e-9 以内 |
| 23b | FFS の要素値コマンド `BPI6 K0 2e-4`、`BPI2 K0 -1e-4` | 4e-8 以内 (BX≈10 で相対 3e-9、SAD の COD の収束の緩さ) |
| 23c | K0 を 2 ファミリー + Emittance[] (RFSW NORAD) | U0 7e-11 相対、δ_COD 1e-11、EMT/TUN 1e-9 |
| 23c2 | 同上で K0 の総和が 0 でない (δ≈−1e-3) | U0 5e-11 相対、COD の δ 1e-11 |
| 23d | pf2 v30_10 (USE RING): BMS/BMA の K0、LINE["E1"] (=.5 を返す)、ANGLE 変更、emit | CALC 全点 2e-11 ($$$ DX は 4.7e-9→出口 COD 化後 1e-10)、TUN 1e-9、U0 8e-9 |
| 23e | esrf: K0 を 2 ファミリーに入れ TrackParticles NORAD/RAD NOFLUC 1,20 周 + Emittance | 追跡 3e-9 (20 周)、1 周 1e-14、emit U0 7e-11 |
| 23f | esrf に ANGLE=0 の BEND ST (L=.15) と SL (L=0) を各セルに挿入 (`ST.1`=1e-5, `SL.1`=-2e-5): CALC、Emittance、追跡 1,3 周 | CALC 1e-10、追跡 N1 1e-16 / N3 1e-11 / R3 2e-13、U0 7e-10、TUN 3e-8 (δ の収束差) |
| 23g | pf2 v9_25 (USE RING): K1..K5 (ANGLE=0, L=.15, E1=.5) の K0 を 1 つずつ + 2 つ同時、emit | 非 $$$ 行 2e-11、EMT 3e-7 相対 (εy)、TUN 2e-8 |
| 23h | ファイルで定義した K0 (ST, SL, BPI2, BPI6) の CALC/emit/追跡 | `$$$` 含め CALC 2e-10、追跡 1e-16〜3e-13 |

<work-pf2> (実行は <work-pf2> で `python ../run_sad.py <lattice> <script>`):
- disp_response_r (v9_25): 534 行の PM 応答表が SAD と 1e-8 µm 以内 (以前は DX 応答が全て 0)。
- disp_kick_effect_r (v9_16_inj): INJ の DX 表 (11.4926477282 など) が SAD と 10 桁一致。
- bump2_r (v9_18): K0 の fit (`FREE BK1 K0`) が働き、MKSTRSS の DX が 0→2e-8 まで追い込まれる。SAD も同じく 2e-9 の残差で止まる (DX = −2.0e-8、pysad は +1.7e-8。どちらも収束の下限)。値は一致しないが、応答が 0 だった以前とは別物。

## 4. 回帰

- 13a/13b/14a/16a/18a/19a/22a: 実行前の保存が無かったので、10/14〜10/15 午前の既存出力と diff。13a/13b/14a は AX の 1e-14 の丸め (表示)、16a は実行時間・警告の行番号、18a は時間、19a は完全一致、22a は INS の数値が 1e-10 変わる (`$$$` を出口 COD にしたため。SAD の値 5.373644720983 に対し 168→070 で同程度)。
- `pytest -q pysad/tests --ignore=…` → 57 passed。
- `QUICK=1 python run_regression.py`: 21 ベンチ中 FAIL 1、他は PASS / 既知 (KNOWN)。
  - **FAIL: 3c の `Emittance[] EMITY`** 差 5.3e-16 (許容 3e-16、相対 −1.95e-6)。四極 z 写像の符号修正 (§2) が原因: 戻すと +3.4e-7 (許容内) に戻る。
    修正後の z 写像は RK 積分・SAD の解析行列 (x のみ、y のみの δ 微分) と一致し (`tmp23/dbg6.py`)、23c/23e の U0 は 3e-7→7e-11 に改善するので正しいと判断して残した。
    3c (四極ロール・結合) の εy だけが 2e-6 動くのは、px,py≠0 の組み合わせで数値微分の (5,6) が SAD の解析式と簿記が違うため (未調査)。許容を厳しく見るなら `_tsolque_body_cod` の xs1/xsh2 の符号を 1 行戻せば元通り。
  - KNOWN: 3c の EX/EPX (元からの差。z の符号修正を戻しても、`$$$` の出口 COD を戻しても同値)、14a の EMIT dz (SAD 表示桁)。
  - SKIP (QUICK=1): 13a, 16a, 7, 12t (13a/16a は上の diff で確認済み)。

## 5. 実装メモ

- tsteee (CALC/emit): `tbende(…,−phib)` で呼ばれ、tsteee の phib 引数 = −K0。本体 rhob = L/phib_arg = −L/K0、エッジは `tbedge(L, −phib_arg=K0, apsi)` で rhob_e = +L/K0 (本体と逆符号。tsteer の追跡でも同じ向きの Maxwell 項)。放射は ndiv 分割 + tradke、F1 フリンジは未対応 (NotImplementedError)。
- tthine(nord=2): L=0 は `px −= K0` ('1/(1+δ)' なし、正準運動量)。追跡は内部の px/(1+δ) なので `px −= K0/(1+δ)`。
- K1≠0 かつ ANGLE=0 かつ K0≠0 の BEND は `_bend_as_mult` で Multipole 相当にして tmulte 経路 (CALC, emit)。**SAD との照合は未実施**。追跡は未対応 (NotImplementedError)。
- track6n (numba): L=0 の thin は符号を修正して継続、L>0 のステアリングと、ずれ/ロール付きの L=0 は NotImplementedError → `track_particles_fast` が track6 へフォールバック (純 Python)。numba カーネル化は未実施。

## 6. SAD の挙動メモ

- tsteer (追跡) の出口の tanp2 は fb2≠0 のときだけ入る (fb2=0 なら 0)。CALC の tsteee は両端で tbedge(psi=apsi) を使う。ANGLE=0, AE1=AE2=0 では差が出ない。
- `$$$` の Twiss は最後の評価の出口 COD。`^^^` と先頭要素の COD は収束点 (入口)。
- Emittance[] の δ (RFSW ON でも path-length で決まる) は SAD の tcod の収束 (codw: δ で ~1e-8) の分だけ揺れる。

## 7. 未解決

- ROTATE=−90 の垂直ステアリング (VD/PV/VS) の K0: パーサが Drift にするため K0 が効かない。emit6d は BEND の ROTATE 自体を扱わない。
- ANGLE=0, K1≠0, K0≠0 の BEND: 追跡 (NotImplementedError)、SAD 照合なし。
- tsteer/tsteee の F1/FB1/FB2 フリンジ。
- track6n の tsteer カーネル (現在は track6 へフォールバック、遅い)。
- bump2 の fit の最終残差・K0 値が SAD と一致しない (どちらも COD 収束の下限で止まる)。
- 3c の EX 2.8e-7 (元からの差)。

## 8. 追記 (2026-10-16): §7 の ROTATE=−90 の垂直ステアリング

解決した。docs/history/work_rotate.md 参照 (parser が Drift にしていたのをやめ Bend のまま保持、CALC・emit・追跡・GEO が SAD と一致、track6n の K_STEER カーネルも追加)。
ANGLE=0, K1≠0, K0≠0 の BEND は SAD K64 が FFS 初期化で SIGSEGV するため照合不可と確定。
