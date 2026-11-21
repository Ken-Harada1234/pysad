# -*- coding: utf-8 -*-
"""example/ddba/v07_14 の .sad を gs と pysad で走らせ、NAME.gs.txt と NAME.pysad.txt に画面表示を保存する。

使い方 (リポジトリ直下から):
    py -3.12 example/tools/run_ddba_v07_14.py example/ddba/v07_14/design/V0714_match_tune.sad [--gs-timeout 1800] [--py-timeout 3600] [--only gs|py]

- `example/ddba/v07_14` の木を作業用の場所 (環境変数 V0714_WORK、既定は一時フォルダ) に複製し、その中の台本のフォルダを
  カレントにして走らせる (`Get["../modules/..."]` の相対パスが効く。台本が書くファイル (datas.txt、SAVE ALL の出力) は複製側に出る)。
- gs は WSL。走らせる前に `sad1.exe` が動いていないことを確認する (他の gs と同時には走らせない)。終了後に残った sad1.exe は PID で止める。
- pysad は k64 モード (gs と同じ結果が要るので)。`PYTHONUTF8=1`。
- 台本が書き出した PNG (`NAME*.png`) と `NAME*.dat` は、リポジトリ側の台本のフォルダへ戻す。
- 所要時間は標準出力と `NAME.time.txt` (作業側) に出す。
"""
import os, re, shutil, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SRC = os.environ.get('V0714_SRC') or os.path.join(ROOT, 'example', 'ddba', 'v07_14')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import sadenv
GS = sadenv.gs_wsl()
PKG = sadenv.packages()


def _work():
    return os.environ.get('V0714_WORK') or os.path.join(tempfile.gettempdir(), 'v0714_run')


def _wsl_path(p):
    p = os.path.abspath(p).replace('\\', '/')
    return '/mnt/' + p[0].lower() + p[2:]


def _wsl(cmd, timeout=None):
    env = dict(os.environ, MSYS_NO_PATHCONV='1')
    return subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c', cmd], capture_output=True, env=env, timeout=timeout)


def sad1_pids():
    r = _wsl("ps -eo pid,comm | grep -i sad1 | grep -v grep | awk '{print $1}'")
    return set(r.stdout.decode().split())


def sync_tree(tag):
    dst = os.path.join(_work(), tag, 'tree')
    os.makedirs(dst, exist_ok=True)
    for dp, dn, fn in os.walk(SRC):
        rel = os.path.relpath(dp, SRC)
        for f in fn:
            if f.endswith(('.sad', '.n')) and not f.endswith(('.gs.txt', '.pysad.txt')):
                os.makedirs(os.path.join(dst, rel), exist_ok=True)
                shutil.copy2(os.path.join(dp, f), os.path.join(dst, rel, f))
    return dst


def run_gs(sad, timeout):
    tree = sync_tree('gs')
    rel = os.path.relpath(os.path.dirname(os.path.abspath(sad)), SRC)
    name = os.path.basename(sad)
    wd = os.path.normpath(os.path.join(tree, rel))
    while sad1_pids():
        print('  (sad1.exe が動いている。60 s 待つ)', flush=True)
        time.sleep(60)
    before = sad1_pids()
    t = time.time()
    try:
        r = _wsl('cd %s && timeout %d %s %s < /dev/null' % (_wsl_path(wd), timeout, GS, name), timeout + 120)
        out = r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        out = b'(python timeout)'
    dt = time.time() - t
    left = sad1_pids() - before
    for p in left:
        _wsl('kill %s' % p)
    if left:
        print('  残った sad1.exe を kill:', sorted(left))
    base = os.path.abspath(sad)[:-4]
    open(base + '.gs.txt', 'wb').write(out)
    return dt


def _clean(out):
    # numpy の RuntimeWarning (作業フォルダの絶対パスが入る) の行と、続くソースの 1 行を除く
    ls = out.decode('utf-8', errors='replace').split(chr(10))
    r = []
    i = 0
    while i < len(ls):
        if 'Warning:' in ls[i] and '.py:' in ls[i]:
            i += 2 if i + 1 < len(ls) and ls[i + 1].startswith('  ') else 1
            continue
        r.append(ls[i])
        i += 1
    return chr(10).join(r).encode('utf-8')


def run_py(sad, timeout):
    tree = sync_tree('py')
    rel = os.path.relpath(os.path.dirname(os.path.abspath(sad)), SRC)
    name = os.path.basename(sad)
    wd = os.path.normpath(os.path.join(tree, rel))
    env = dict(os.environ, PYTHONPATH=ROOT, PYTHONIOENCODING='utf-8', PYTHONUTF8='1',
               PYSAD_COMPAT='k64', PYSAD_SAD_PACKAGES=PKG, MPLBACKEND='Agg')
    t = time.time()
    try:
        r = subprocess.run(['py', '-3.12', '-c',
                            "import sys;from pysad.ffscmd import run_sad_string;"
                            "run_sad_string(open(sys.argv[1],encoding='utf-8').read(),cwd=sys.argv[2])", name, wd],
                           capture_output=True, env=env, timeout=timeout, cwd=wd)
        out = r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        out = b'(python timeout)'
    dt = time.time() - t
    base = os.path.abspath(sad)[:-4]
    out = _clean(out)
    open(base + '.pysad.txt', 'wb').write(out)
    bn = os.path.basename(base)
    for f in os.listdir(wd):
        if f.startswith(bn) and f.endswith(('.png', '.dat')):
            shutil.copy2(os.path.join(wd, f), os.path.join(os.path.dirname(os.path.abspath(sad)), f))
    return dt


if __name__ == '__main__':
    a = sys.argv[1:]
    g = int(a[a.index('--gs-timeout') + 1]) if '--gs-timeout' in a else 900
    p = int(a[a.index('--py-timeout') + 1]) if '--py-timeout' in a else 1800
    only = a[a.index('--only') + 1] if '--only' in a else None
    for s in [x for x in a if x.endswith('.sad')]:
        res = []
        if only in (None, 'gs'):
            res.append('gs %.0f s' % run_gs(s, g))
        if only in (None, 'py'):
            res.append('pysad %.0f s' % run_py(s, p))
        line = '%s  %s' % (os.path.basename(s), ', '.join(res))
        print(line, flush=True)
        open(os.path.abspath(s)[:-4] + '.time.txt', 'a').write(time.strftime('%Y-%m-%d %H:%M ') + line + '\n')
