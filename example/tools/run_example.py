# -*- coding: utf-8 -*-
"""example/NAME.sad を gs と pysad で走らせ、NAME.gs.txt と NAME.pysad.txt に画面表示を保存する。
使い方: py -3.12 example/tools/run_example.py example/FFS_CALCULATE.sad [timeout]"""
import os, subprocess, sys, shutil, tempfile
GS = '/mnt/d/work/SAD-master/SAD-master/bin/gs'
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
def run(path, to=120):
    path = os.path.abspath(path); base = path[:-4]; name = os.path.basename(base)
    tmp = tempfile.gettempdir()
    shutil.copy(path, os.path.join(tmp, '_ex_%s.sad' % name))
    t = tmp.replace(chr(92), '/').replace('C:', '/mnt/c')
    env = dict(os.environ, MSYS_NO_PATHCONV='1')
    r = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c',
                        'cd %s && timeout %d %s _ex_%s.sad < /dev/null' % (t, to, GS, name)],
                       capture_output=True, env=env)
    open(base + '.gs.txt', 'wb').write(r.stdout + r.stderr)
    env = dict(os.environ, PYTHONPATH=ROOT, PYTHONIOENCODING='utf-8',
               PYSAD_SAD_PACKAGES='D:/work/SAD-master/SAD-master/Packages')
    r = subprocess.run(['py', '-3.12', '-c',
                        "from pysad.ffscmd import run_sad_string;import sys;run_sad_string(open(sys.argv[1],encoding='utf-8').read())", path],
                       capture_output=True, env=env, timeout=to * 5)
    open(base + '.pysad.txt', 'wb').write(r.stdout + r.stderr)
if __name__ == '__main__':
    for p in sys.argv[1:]:
        run(p)
