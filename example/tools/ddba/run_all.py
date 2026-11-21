# -*- coding: utf-8 -*-
"""全 DDBA_NN で cell calc emit を pysad で走らせ、結果を JSON に保存する。
使い方: DDBA_SRC=<元の台本のフォルダ> py -3.12 example/tools/ddba/run_all.py <出力 json> [k64|v2] [lattice|src]
lattice = example/ddba/lattice の整えた台本、src = 元の台本 (整える前の確認用)"""
import sys, os, json, re, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ddba_lib as L
sys.path.insert(0, os.path.join(L.EX, 'tools', 'pf_history'))
def emit_table(txt):
    import pf_lib
    return pf_lib.emit_table(txt)
def parse(out):
    g = lambda p: re.search(p, out)
    d = {}
    m = g(r'circumference \[m\] =\s*([-\d.E+]+)'); d['C'] = float(m.group(1)) if m else None
    m = g(r'tune nx ny =\s*([-\d.E+]+)\s+([-\d.E+]+)'); d['nx'], d['ny'] = (float(m.group(1)), float(m.group(2))) if m else (None, None)
    m = g(r'MOMENTUM \[eV\] =\s*([-\d.E+]+)'); d['P'] = float(m.group(1)) if m else None
    m = g(r'bends =\s*([\d.]+)\s+sum of ANGLE \[rad\] =\s*([-\d.E+]+)\s+small-angle bends =\s*([\d.]+)\s+quads =\s*([\d.]+)\s+sexts =\s*([\d.]+)\s+octs =\s*([\d.]+)\s+cavities =\s*([\d.]+)')
    if m: d['nb'], d['sa'], d['nr'], d['nq'], d['ns'], d['no'], d['nc'] = [float(x) for x in m.groups()]
    m = g(r'CAVI total VOLT \[V\] FREQ \[Hz\] HARM =\s*([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)')
    if m: d['V'], d['f'], d['h'] = [float(x) for x in m.groups()]
    try: d.update(emit_table(out))
    except Exception as e: d['emit_err'] = 'no emit table'
    return d
if __name__ == '__main__':
    outp = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else 'k64'; which = sys.argv[3] if len(sys.argv) > 3 else 'lattice'
    res = {}; os.makedirs(os.path.join(os.path.dirname(outp), 'raw'), exist_ok=True)
    for nn, srcname, _ in L.MAP:
        t = time.time()
        p = L.lat_path(nn) if which == 'lattice' else L.find_src(srcname)
        try:
            src = open(p, encoding='utf-8').read().replace('\r', '')
            out = L.run(src, mode)
            open(os.path.join(os.path.dirname(outp), 'raw', 'DDBA_%02d_%s_%s.txt' % (nn, mode, which)), 'w', encoding='utf-8').write(out)
            d = parse(out); d['sec'] = round(time.time() - t, 1)
            res[nn] = d
        except Exception as e:
            res[nn] = {'error': repr(e)[:300]}
        print(nn, res[nn].get('C'), res[nn].get('nx'), res[nn].get('ny'), res[nn].get('ex'), res[nn].get('error', ''), flush=True)
    json.dump(res, open(outp, 'w'), indent=1)
