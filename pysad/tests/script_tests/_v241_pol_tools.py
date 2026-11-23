# -*- coding: utf-8 -*-
"""v2.4.1: 偏光 (POL) の FLUC つきスピン追跡の gs・pysad 比較の道具 (台本の生成と出力の解析)。

台本 (`scripts/v241_pol_*.sad`) は gs と pysad (k64) に同じものを渡す。生成は `py -3.12 pysad/tests/script_tests/_v241_pol_tools.py gen ...`
(閉軌道・n0 (k64 の 1 周のスピン行列)・N0 (v2 の spin_orbit) を pysad で計算してリテラルで台本に埋める。数分かかる)。
試験は台本と gs の出力 (`.gs.txt`) だけを使うので、生成の部分は試験では走らない。
リング: fodo6 / fodo24 = 12 セル FODO (4.2 GeV、RF 6 / 24 MV、HARM 64、四極 2 台に縦ずれ ±3 mm)、
pf = PF (2.5 GeV、example/pf/lattice/v9_16_asym.sad の複製を読むだけ。四極 78 台に縦ずれ rms 0.4 mm、乱数の種 5)、
astrid2 = ASTRID2 (0.58 GeV、example/astrid/lattice/ASTRID2_v01_01.sad。四極 24 台に縦ずれ rms 2 mm、種 5。k64 の閉軌道が 2.5 mm で収束しない)。
"""
import io
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SP = os.path.join(HERE, 'scripts')
CACHE = os.environ.get('V241_POL_CACHE', os.path.join(HERE, 'scripts'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def fodo_lat(E, volt, dy=3e-3):
    return f'''MOMENTUM = {E}GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3=(L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QFE = (L=0.3 K1=0.7 DY={dy}) QDE = (L=0.3 K1=-0.6 DY=-{dy});
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT={volt}e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QFE D2 B D2 QD D1);
LINE CELL3 = (QF D2 B D2 QDE D1);
LINE RING = (ST 2*CELL1 CELL2 4*CELL1 CELL3 4*CELL1 RFC);
'''


def read(p):
    return open(os.path.join(ROOT, p), encoding='utf-8').read()


def lattice(name):
    """戻り: (格子の文字列, ライン名, DY 乱数の rms [m] か None, 台数)"""
    if name == 'fodo6':
        return fodo_lat(4.2, 6), 'RING', None
    if name == 'fodo24':
        return fodo_lat(4.2, 24), 'RING', None
    if name == 'pf':
        return read('example/pf/lattice/v9_16_asym.sad'), 'RING', 0.4e-3
    if name == 'astrid2':
        return read('example/astrid/lattice/ASTRID2_v01_01.sad'), 'RING', 2.0e-3
    raise KeyError(name)


def dy_list(nq, rms, seed=5):
    rng = np.random.default_rng(seed)
    return [float(x) for x in np.round(rng.normal(0, rms, nq), 9)]


def sadlist(v):
    return '{' + ','.join('%.12g' % x for x in v) + '}'


HEAD = '''FFS;
$FORM = "14.8";
PageWidth = 1999;
NPARA = 1;
USE {ring};
CELL;
{dyset}CALC;
'''


def setup_text(name, nq=None):
    lat, ring, rms = lattice(name)
    dyset = ''
    if rms is not None:
        dyset = 'LINE["DY","Q*"] = %s;\n' % sadlist(dy_list(nq, rms))
    return lat + HEAD.format(ring=ring, dyset=dyset), ring


NQ = {'fodo6': None, 'fodo24': None, 'pf': 78, 'astrid2': 24}


def prep(name, want_n0=True, want_N0=True):
    """cod (k64 の Emittance[] の入口の閉軌道)・n0 (1 周のスピン行列の固有ベクトル)・N0 (v2 の spin_orbit) を pysad で作る (キャッシュ)"""
    fn = os.path.join(CACHE, 'v241_pol_prep_%s.json' % name)
    if os.path.exists(fn):
        return json.load(open(fn))
    from pysad import compat, emit6d, track6, ffscmd
    txt, ring = setup_text(name, NQ[name])
    res = {}
    compat.set_mode('k64')
    ffs, _ = ffscmd.run_sad_string(txt + 'RFSW; RAD; NOFLUC; POL;\nExit;\n', out=io.StringIO())
    p = ffs.momentum()
    e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True, pol=True)
    sp_ = e['spin']
    els = ffs.lat.elements
    circ = emit6d._emit_circ(els)
    trf0 = e.get('trf0', 0.0)
    cod = np.array(e['orbit_entrance'])
    M = np.zeros((3, 3))
    for j in range(3):
        sj = np.zeros(3); sj[j] = 1
        _o, s_ = track6.track_particles(els, np.array([cod]), 1, p, trf0, True, True, False, spins=[sj], radpol=False, circ=circ)
        M[:, j] = s_[0]
    w, v = np.linalg.eig(M)
    k = np.argmin(abs(w - 1)); n0 = np.real(v[:, k]); n0 /= np.linalg.norm(n0)
    if n0[1] < 0: n0 = -n0
    res['cod'] = cod.tolist(); res['n0'] = n0.tolist(); res['M'] = M.tolist()
    res['k64'] = {k_: (v_ if np.ndim(v_) == 0 else np.asarray(v_).tolist()) for k_, v_ in sp_.items() if k_ in ('nup', 'smu', 'sdamp', 'taup', 'sdamp_rel', 'equpol')}
    res['momentum'] = p
    res['circ'] = circ
    compat.set_mode('v2')
    ffs, _ = ffscmd.run_sad_string(txt + 'RFSW; RAD; NOFLUC; POL;\nExit;\n', out=io.StringIO())
    e2 = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    s2 = e2['spin']
    res['N0'] = np.asarray(s2['spin_orbit']).tolist() if s2.get('spin_orbit') is not None else None
    res['v2'] = {k_: (v_ if np.ndim(v_) == 0 else np.asarray(v_).tolist()) for k_, v_ in s2.items() if k_ in ('nup', 'taup', 'sdamp', 'sdamp_rel', 'equpol', 'depol_rate', 'depol_rate_linear', 'taudep', 'spin_spread', 'spin_spread_sb', 'mod_spin', 'mod_chrom', 'equpol_k64')}
    res['v2cod'] = np.asarray(e2['orbit_entrance']).tolist()
    for key in ('emittance_x', 'emittance_y', 'momentum_spread', 'bunch_length', 'damping_rate'):
        res['v2_' + key] = np.asarray(e2[key]).tolist() if np.ndim(e2[key]) else float(e2[key])
        res['k64_' + key] = np.asarray(e[key]).tolist() if np.ndim(e[key]) else float(e[key])
    json.dump(res, open(fn, 'w'), indent=1)
    return res


BODY = r'''
RFSW; RAD; NOFLUC; NORADPOL;
POL;
em = Emittance[];
Print["EM ", {PolarizationTime, SpinTune, NominalSpinTune, Polarization, PolarizationVector}/.em];
Print["EMB ", {MomentumSpread, EquilibriumPosition}/.em, " ", Emittances/.em];
n0 = %(n0)s;
N0 = %(N0)s;
cod = %(cod)s;
np = %(np)d;
one = Table[1, {i, np}];
%(init)s
SeedRandom[%(seed)d];
FLUC; NOPOL;
If[%(neq)d > 0, beam = TrackParticles[{1, beam}, 1, 1, %(neq)d][[2]]];
rms[v_] := Sqrt[Plus@@(v^2)/Length[v]];
Print["EQ turns=", %(neq)d, " rms x px y py z dp ", rms[beam[[1]]-cod[[1]]], " ", rms[beam[[2]]-cod[[2]]], " ", rms[beam[[3]]-cod[[3]]], " ", rms[beam[[4]]-cod[[4]]], " ", rms[beam[[5]]-cod[[5]]], " ", rms[beam[[6]]-cod[[6]]]];
u = Table[{beam[[1,i]]-cod[[1]], beam[[2,i]]-cod[[2]], beam[[3,i]]-cod[[3]], beam[[4,i]]-cod[[4]], beam[[5,i]]-cod[[5]], beam[[6,i]]-cod[[6]]}, {i, np}];
s = Table[n0 + N0.u[[i]], {i, np}];
s = Table[s[[i]]/Sqrt[s[[i]].s[[i]]], {i, np}];
ang = Table[ArcTan[s[[i,1]], s[[i,3]]], {i, np}];
sy = Table[s[[i,2]], {i, np}];
beam = {beam[[1]], beam[[2]], beam[[3]], beam[[4]], beam[[5]], beam[[6]], ang, sy, one, one, one};
sn[b_] := Table[Sqrt[1-b[[8,i]]^2]*Cos[b[[7,i]]]*n0[[1]] + b[[8,i]]*n0[[2]] + Sqrt[1-b[[8,i]]^2]*Sin[b[[7,i]]]*n0[[3]], {i, Length[b[[1]]]}];
snl[b_] := Table[Module[{uu, nn, ss},
  uu = {b[[1,i]]-cod[[1]], b[[2,i]]-cod[[2]], b[[3,i]]-cod[[3]], b[[4,i]]-cod[[4]], b[[5,i]]-cod[[5]], b[[6,i]]-cod[[6]]};
  nn = n0 + N0.uu; nn = nn/Sqrt[nn.nn];
  ss = {Sqrt[1-b[[8,i]]^2]*Cos[b[[7,i]]], b[[8,i]], Sqrt[1-b[[8,i]]^2]*Sin[b[[7,i]]]};
  nn.ss], {i, Length[b[[1]]]}];
Print["S0 1-<S.n0> = ", 1 - Plus@@sn[beam]/np, " L= ", 1 - Plus@@snl[beam]/np];
POL;
Do[
  FFS["FLUC"];
  beam = TrackParticles[{1, beam}, 1, 1, %(n1)d][[2]];
  bd = beam;
  a = sn[beam];
  Print["STAGE ", k, " turns ", k*%(n1)d, " A= ", 1 - Plus@@a/np, " sdA= ", Sqrt[(Plus@@(a^2)/np - (Plus@@a/np)^2)/np], " L= ", 1 - Plus@@snl[beam]/np,
        " rms x y dp z ", rms[beam[[1]]-cod[[1]]], " ", rms[beam[[3]]-cod[[3]]], " ", rms[beam[[6]]-cod[[6]]], " ", rms[beam[[5]]-cod[[5]]]];
  FFS["NOFLUC"];
  bd = TrackParticles[{1, bd}, 1, 1, %(nd)d][[2]];
  d = sn[bd];
  Print["DAMP ", k, " turns ", k*%(n1)d, " D= ", 1 - Plus@@d/np, " sdD= ", Sqrt[(Plus@@(d^2)/np - (Plus@@d/np)^2)/np],
        " orb ", Max[Abs[bd[[1]]-cod[[1]]]], " ", Max[Abs[bd[[3]]-cod[[3]]]]];
  Print["DV ", k, " ", 1-d];
  Print["AV ", k, " ", 1-a];
  Print["LV ", k, " ", 1-snl[beam]],
  {k, 1, %(nstage)d}];
Print["SPART ", Table[1-sn[beam][[i]], {i, Min[np, 8]}]];
Exit;
'''


def sigma(name):
    fn = os.path.join(CACHE, 'v241_pol_sigma_%s.json' % name)
    if os.path.exists(fn):
        return json.load(open(fn))
    from pysad import compat, emit6d, ffscmd
    compat.set_mode('k64')
    txt, ring = setup_text(name, NQ[name])
    ffs, _ = ffscmd.run_sad_string(txt + 'RFSW; RAD; NOFLUC;\nExit;\n', out=io.StringIO())
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=False)
    S = np.asarray(e['beam_matrix']).tolist()
    json.dump(S, open(fn, 'w'))
    return S


def gauss_init(name, np_, seed):
    r = prep(name)
    S = np.array(sigma(name))
    rng = np.random.default_rng(1000 + seed)
    x = rng.multivariate_normal(np.zeros(6), S, size=np_) + np.array(r['cod'])
    rows = ['{' + ','.join('%.12g' % v for v in x[:, j]) + '}' for j in range(6)]
    return 'beam = {' + ','.join(rows) + ', one};'


def script(name, np_, neq, n1, nd, nstage, seed, init='cod'):
    r = prep(name)
    txt, ring = setup_text(name, NQ[name])
    n0 = sadlist(r['n0'])
    N0s = '{' + ','.join(sadlist(row) for row in r['N0']) + '}'
    if init == 'cod':
        ini = 'beam = {cod[[1]]*one, cod[[2]]*one, cod[[3]]*one, cod[[4]]*one, cod[[5]]*one, cod[[6]]*one, one};'
    else:
        ini = gauss_init(name, np_, seed)
    return txt + BODY % dict(n0=n0, N0=N0s, cod=sadlist(r['cod']), np=np_, seed=seed, neq=neq, n1=n1, nd=nd, nstage=nstage, init=ini)


TUNE_BODY = r'''
RFSW; RAD; NOFLUC; NORADPOL;
POL;
em = Emittance[];
Print["EM ", {PolarizationTime, SpinTune, NominalSpinTune, Polarization, PolarizationVector}/.em];
cod = %(cod)s;
b = {{cod[[1]]},{cod[[2]]},{cod[[3]]},{cod[[4]]},{cod[[5]]},{cod[[6]]},{%(ang)s},{%(sy)s},{1},{1},{1}};
Do[
  b = TrackParticles[{1, b}, 1, 1, 1][[2]];
  Print["T ", k, " ", Sqrt[1-b[[8,1]]^2]*Cos[b[[7,1]]], " ", b[[8,1]], " ", Sqrt[1-b[[8,1]]^2]*Sin[b[[7,1]]], " ", Table[b[[j,1]]-cod[[j]], {j, 6}]],
  {k, 1, %(nt)d}];
Exit;
'''


def tune_script(name, nt=100):
    r = prep(name)
    n0 = np.array(r['n0'])
    ex = np.array([1.0, 0, 0])
    m = ex - ex.dot(n0) * n0
    m /= np.linalg.norm(m)
    ang = float(np.arctan2(m[2], m[0]))
    sy = float(m[1])
    txt, ring = setup_text(name, NQ[name])
    return txt + TUNE_BODY % dict(cod=sadlist(r['cod']), ang='%.15g' % ang, sy='%.15g' % sy, nt=nt)


num = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?'


def read_stat(fn):
    txt = open(fn, encoding='utf-8', errors='replace').read()
    out = {}
    m = re.search(r'^EM\s*\{(.*?)\}\s*\n', txt, re.M | re.S)
    out['EM'] = [float(x) for x in re.findall(num, m.group(1).replace('{', ' '))] if m else None
    m = re.search(r'^EMB\s*\{(.*?)\}\s*\{(.*?)\}', txt, re.M | re.S)
    if m:
        out['EMB'] = [float(x) for x in re.findall(num, m.group(1))] + [float(x) for x in re.findall(num, m.group(2))]
    m = re.search(r'^EQ turns=\s*(\S+) rms x px y py z dp\s+(.*)$', txt, re.M)
    if m:
        out['EQ'] = [float(x) for x in re.findall(num, m.group(2))]
    for tag in ('DV', 'AV', 'LV'):
        out[tag] = {}
        for m in re.finditer(r'^%s\s+(\S+)\s*\{(.*?)\}' % tag, txt, re.M | re.S):
            out[tag][int(round(float(m.group(1))))] = np.array([float(x) for x in re.findall(num, m.group(2))])
    stages = []
    for m in re.finditer(r'^STAGE\s+(\S+) turns\s+(\S+) A=\s*(\S+) sdA=\s*(\S+) L=\s*(\S+) rms x y dp z\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', txt, re.M):
        stages.append(dict(k=int(float(m.group(1))), turns=float(m.group(2)), A=float(m.group(3)), sdA=float(m.group(4)), L=float(m.group(5)),
                           rms_x=float(m.group(6)), rms_y=float(m.group(7)), rms_dp=float(m.group(8)), rms_z=float(m.group(9))))
    dmp = []
    for m in re.finditer(r'^DAMP\s+(\S+) turns\s+(\S+) D=\s*(\S+) sdD=\s*(\S+) orb\s+(\S+)\s+(\S+)', txt, re.M):
        dmp.append(dict(k=int(float(m.group(1))), turns=float(m.group(2)), D=float(m.group(3)), sdD=float(m.group(4)), orbx=float(m.group(5)), orby=float(m.group(6))))
    out['stages'] = stages
    out['damp'] = dmp
    m = re.search(r'^S0 1-<S.n0> =\s*(\S+) L=\s*(\S+)', txt, re.M)
    if m:
        out['S0'] = (float(m.group(1)), float(m.group(2)))
    return out


def rate(out, tag='DV', n1=None):
    """粒子ごとの時間への直線当てはめ (切片つき) の傾きの平均と標準誤差 [/周]。段数が 1 のときは (値 − S0) / turns"""
    stg = {s['k']: s['turns'] for s in out['stages']}
    ks = sorted(out[tag])
    if len(ks) < 2:
        return None
    t = np.array([stg[k] for k in ks], float)
    Y = np.array([out[tag][k] for k in ks])      # (nk, np)
    tm = t.mean()
    sl = ((t - tm)[:, None] * (Y - Y.mean(0))).sum(0) / ((t - tm) ** 2).sum()
    N = Y.shape[1]
    return float(sl.mean()), float(sl.std(ddof=1) / np.sqrt(N)), float(Y.mean(1)[0] - sl.mean() * t[0]), N


def read_tune(fn):
    em = None
    S = []
    for l in open(fn, encoding='utf-8', errors='replace'):
        if l.startswith('EM '):
            em = [float(x) for x in re.findall(num, l[3:])]
        elif l.startswith('T '):
            v = [float(x) for x in re.findall(num, l[2:])]
            S.append(v[1:4])
    return em, np.array(S)


def tune_from_spin(S, n0):
    """S: 毎周のスピン (x,y,z)、n0: 軸。n0 まわりの 1 周の回転角の平均 (符号つき、1 周あたり /2π)"""
    n0 = np.asarray(n0, float)
    m = np.array([1.0, 0, 0]) - n0[0] * n0
    m /= np.linalg.norm(m)
    l = np.cross(n0, m)
    c = S @ m + 1j * (S @ l)
    ph = np.unwrap(np.angle(c))
    k = np.arange(len(ph))
    p = np.polyfit(k, ph, 1)
    res = ph - np.polyval(p, k)
    return p[0] / (2 * np.pi), float(np.abs(res).max()), float(np.abs(c).mean()), float(np.abs(S @ n0).mean())


if __name__ == '__main__':
    # stat RING NP NEQ N1 ND NSTAGE SEED OUT [cod|gauss]    統計・同一性の台本
    # tune RING NT OUT                                      スピンチューンの台本
    cmd = sys.argv[1]
    if cmd == 'stat':
        name = sys.argv[2]
        np_, neq, n1, nd, nstage, seed = [int(x) for x in sys.argv[3:9]]
        out = sys.argv[9]
        init = sys.argv[10] if len(sys.argv) > 10 else 'cod'
        with open(out, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(script(name, np_, neq, n1, nd, nstage, seed, init))
    elif cmd == 'tune':
        with open(sys.argv[4], 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(tune_script(sys.argv[2], int(sys.argv[3])))
