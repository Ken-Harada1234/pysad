# -*- coding: utf-8 -*-
"""SAD のラティス (リング) を elegant の lte に変換する (pysad で展開した要素の並びから書く)。

使い方 (リポジトリ直下から):
  PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/sad2ele.py astrid2|kekls [出力のフォルダ]
出力: <名前>.lte (elegant のラティス)、<名前>.map.json (要素の番号 → 名前・s・種類、pysad の位置 loc との対応)。

要素の対応 (SAD → elegant):
  DRIFT  → DRIF
  QUAD   → KQUAD (K1 [1/m^2]、SYNCH_RAD=1、ISR=0、N_KICKS=nk_q、端のフリンジなし)
  SEXT   → KSEXT (K2 [1/m^3]、SYNCH_RAD=1、ISR=0)
  OCT    → KOCT  (K3 [1/m^4]、SYNCH_RAD=1、ISR=0)
  BEND   → CSBEND (L は弧長、ANGLE、K1 [1/m^2] (複合機能)、K2 = −K1·ANGLE/L (K1≠0 のとき。SAD の勾配の扱いに合わせる)、E1・E2 = SAD の E1・E2 × ANGLE [rad]、SYNCH_RAD=1、ISR=0、N_KICKS=nk_b、TILT=ROTATE)
           L = 0 かつ ANGLE = 0 の BEND (座標回転だけの要素) は MARK (恒等)
  CAVI   → RFCA (VOLT、FREQ、PHASE = 180 − asin(U0/V) [度]、L = 0。U0 は --u0 で渡す)
  MARK・MONI → MARK
SAD の K1 は 1/m^2 (pysad の要素の K1 属性が既に L で割った値)。FSE・誤差・マルチポール・SOL は無い (どちらのリングも)。
LMA の損失判定: 各要素の出口に MAXAMP (±aper、矩形) を置く (pysad の check_every=1 と同じ判定の位置)。
(elegant の momentum_aperture は、RF バケットから出ただけでは損失にしない。CLEAN で |t| を切る方法は、t が周回で増える絶対時間なので使えなかった)
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)

CASES = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), ring='RING'),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), ring='RING'),
}


def convert(name, outdir, u0_eV, aper=0.010, nk_b=40, nk_q=8, nk_s=4, nk_o=4, radiate=True, maxamp=True, sad_k1=True, ap_sel=None, tag='', am_unique=False):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    cs = CASES[name]
    src = open(cs['path'], encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, cs['ring'], rad=True, rfsw=True, alost=aper)
    sr = 1 if radiate else 0
    defs = {}          # 名前 → 定義文字列 (同じ名前で中身が違えば連番で別名にする)
    order = []         # 並び (定義の名前)
    smap = []
    s = 0.0
    nap = 0
    for i, e in enumerate(st.els):
        t = type(e).__name__
        L = float(e.L)
        if t == 'Drift':
            nm, d = 'D_%s' % e.name, 'DRIF, L=%.12g' % L
        elif t == 'Quadrupole':
            assert abs(getattr(e, 'K1S', 0.0)) == 0 and e.ROTATION == 0 and e.DX == 0 and e.DY == 0 and e.ROTATE == 0
            nm, d = 'Q_%s' % e.name, 'KQUAD, L=%.12g, K1=%.15g, N_KICKS=%d, SYNCH_RAD=%d, ISR=0, INTEGRATION_ORDER=4' % (L, e.K1, nk_q, sr)
        elif t == 'Sextupole':
            assert e.ROTATE == 0 and e.ROTATION == 0
            nm, d = 'S_%s' % e.name, 'KSEXT, L=%.12g, K2=%.15g, N_KICKS=%d, SYNCH_RAD=%d, ISR=0, INTEGRATION_ORDER=4' % (L, e.K2, nk_s, sr)
        elif t == 'Octupole':
            assert e.ROTATE == 0 and e.ROTATION == 0
            nm, d = 'O_%s' % e.name, 'KOCT, L=%.12g, K3=%.15g, N_KICKS=%d, SYNCH_RAD=%d, ISR=0, INTEGRATION_ORDER=4' % (L, e.K3, nk_o, sr)
        elif t == 'Bend':
            if L == 0.0 and e.ANGLE == 0.0:
                nm, d = 'M_%s' % e.name, 'MARK'
            else:
                assert getattr(e, 'DX', 0) == 0 and e.ROTATION == 0 and e.K0 == 0
                e1, e2 = e.E1 * e.ANGLE, e.E2 * e.ANGLE
                nm = 'B_%s' % e.name
                # K1 != 0 の BEND: SAD (K64、gs と pysad は 1e-12 で同じ) の本体は、勾配の項に曲率の因子 (1 + h x) が掛からない。
                # elegant の CSBEND (厳密な曲線座標) に K2 = -h*K1 (h = ANGLE/L) を足すと SAD と同じ δ 依存の光学 (チューンの δ 依存) になる
                # (一様な磁場指数のリングで 4 桁一致: gs_bend_k1.py・bend_k1_elegant.py)。sad_k1=False なら K2 を足さない (厳密な勾配)。
                # EDGE_ORDER=2: 端の効果に傾き x' の 1 次の項 (オフモメンタムの軌道の端での垂直方向の集束) を入れる (既定の 1 では SAD と ξy が違う)
                k2b = -e.K1 * e.ANGLE / L if (sad_k1 and L != 0.0) else 0.0
                d = ('CSBEND, L=%.12g, ANGLE=%.15g, K1=%.15g, K2=%.15g, E1=%.15g, E2=%.15g, TILT=%.15g, HGAP=%.6g, FINT=%.6g, '
                     'N_KICKS=%d, NONLINEAR=1, EDGE_ORDER=2, SYNCH_RAD=%d, ISR=0, INTEGRATION_ORDER=4' % (
                         L, e.ANGLE, e.K1, k2b, e1, e2, e.ROTATE, e.HGAP, e.FINT, nk_b, sr))
        elif t == 'Cavity':
            phase = 180.0 - math.degrees(math.asin(u0_eV / e.VOLT))
            nm, d = 'RF_%s' % e.name, 'RFCA, L=0, VOLT=%.12g, FREQ=%.12g, PHASE=%.12g, N_KICKS=0' % (e.VOLT, e.FREQ, phase)
        elif t in ('Marker', 'Monitor'):
            nm, d = 'M_%s' % e.name, 'MARK'
        else:
            raise ValueError('未対応の要素 %s (%s)' % (t, e.name))
        key = nm
        k = 1
        while key in defs and defs[key] != d:
            k += 1
            key = '%s_%d' % (nm, k)
        defs[key] = d
        order.append(key)
        s += L
        smap.append(dict(i=i + 1, name=str(e.name), type=t, s_in=s - L, s_out=s, ele=key))
        if maxamp:
            # 要素 i+1 の出口の MAXAMP。ap_sel (要素の番号の集合) の位置は名前 AM (momentum_aperture の対象)、ほかは AP
            if ap_sel is not None and (i + 1) in ap_sel:
                if am_unique:       # 位置ごとに別の名前 AM001、AM002、… (elegant の include_name_pattern で 1 位置だけ選べる。skip_elements は要素の数え方が違った)
                    nm_ap = 'AM%03d' % (sorted(ap_sel).index(i + 1) + 1)
                    defs[nm_ap] = 'MAXAMP, X_MAX=%.6g, Y_MAX=%.6g, ELLIPTICAL=0' % (aper, aper)
                else:
                    nm_ap = 'AM'
            else:
                nm_ap = 'AP'
            order.append(nm_ap)
            nap += 1
    defs['AP'] = 'MAXAMP, X_MAX=%.6g, Y_MAX=%.6g, ELLIPTICAL=0' % (aper, aper)
    defs['AM'] = defs['AP']
    lines = ['! SAD -> elegant (sad2ele.py): %s' % os.path.basename(cs['path'])]
    for k, d in defs.items():
        lines.append('%s: %s' % (k, d))
    # LINE (行の長さに上限があるので折り返す)
    body = ', '.join(order)
    chunks, cur = [], ''
    for tok in body.split(', '):
        if len(cur) + len(tok) > 400:
            chunks.append(cur)
            cur = ''
        cur += (tok + ', ')
    chunks.append(cur.rstrip(', '))
    lines.append('RING: LINE=(' + ' &\n'.join(chunks) + ')')
    path = os.path.join(outdir, '%s%s.lte' % (name, tag))
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(lines) + '\n')
    with open(os.path.join(outdir, '%s%s.map.json' % (name, tag)), 'w') as f:
        json.dump(dict(circumference=s, nlat=st.nlat, elements=smap), f)
    return path, st, ffs


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('outdir', nargs='?', default=HERE)
    ap.add_argument('--u0', type=float, default=None, help='1 周のエネルギー損失 [eV] (RF の位相用。省略なら pysad の値)')
    ap.add_argument('--aper', type=float, default=0.010)
    ap.add_argument('--no-rad', action='store_true')
    ap.add_argument('--no-maxamp', action='store_true')
    ap.add_argument('--exact-k1', action='store_true', help='K1≠0 の BEND に K2 = -h K1 を足さない (厳密な勾配の曲線座標)')
    a = ap.parse_args()
    u0 = a.u0
    if u0 is None:
        from pysad import compat
        from pysad.trksetup import TrackSetup
        compat.set_mode('v2')
        src = open(CASES[a.name]['path'], encoding='utf-8').read()
        u0 = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True)[0].conditions()['U0_eV']
    p, st, ffs = convert(a.name, a.outdir, u0, aper=a.aper, radiate=not a.no_rad, maxamp=not a.no_maxamp, sad_k1=not a.exact_k1)
    print('wrote', p, ' U0 for RF phase = %.3f eV' % u0)
