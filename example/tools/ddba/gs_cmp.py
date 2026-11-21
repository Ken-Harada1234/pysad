# -*- coding: utf-8 -*-
"""代表の DDBA_NN を gs (WSL) と pysad (k64) で走らせて比べる。gs は 1 本ずつ。
使い方: py -3.12 example/tools/ddba/gs_cmp.py 1 15 25 28   (出力: example/ddba/check/DDBA_NN.gs.txt / .pysad.txt と一致表)"""
import sys, os, subprocess, tempfile, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ddba_lib as L, run_all as R
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import sadenv
GS = sadenv.gs_wsl()
OUT = os.path.join(L.EX, 'ddba', 'check'); os.makedirs(OUT, exist_ok=True)
def gs_run(nn, timeout=900):
    src = open(L.lat_path(nn), encoding='utf-8').read()
    script = src + '\nON LOG;\nON RAD FLUC COD RFSW RADCOD;\n' + L.CALC
    tmp = tempfile.gettempdir(); name = '_ddba_gs_%02d' % nn
    open(os.path.join(tmp, name + '.sad'), 'w', encoding='utf-8', newline='\n').write(script)
    t = tmp.replace(chr(92), '/').replace('C:', '/mnt/c')
    env = dict(os.environ, MSYS_NO_PATHCONV='1')
    r = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c', 'cd %s && timeout %d %s %s.sad < /dev/null' % (t, timeout, GS, name)], capture_output=True, env=env)
    return (r.stdout + r.stderr).decode('utf-8', 'replace')
if __name__ == '__main__':
    res = {}
    for nn in map(int, sys.argv[1:]):
        g = gs_run(nn); open(os.path.join(OUT, 'DDBA_%02d.gs.txt' % nn), 'w', encoding='utf-8').write(g)
        p = L.run(open(L.lat_path(nn), encoding='utf-8').read(), 'k64'); open(os.path.join(OUT, 'DDBA_%02d.pysad.txt' % nn), 'w', encoding='utf-8').write(p)
        a, b = R.parse(g), R.parse(p); res[nn] = (a, b)
        print(nn, 'gs', a.get('C'), a.get('nx'), a.get('ny'), a.get('ex'), '| pysad', b.get('C'), b.get('nx'), b.get('ny'), b.get('ex'), flush=True)
    json.dump(res, open(os.path.join(OUT, 'gs_vs_pysad.json'), 'w'), indent=1)
