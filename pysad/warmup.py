# -*- coding: utf-8 -*-
"""Precompilation of numba kernels (`pysad-warmup`), run once after installation so later calls load cached machine code.
numba のカーネルの事前コンパイル (`pysad-warmup`)。

pysad の計算の核 (光学・追跡・放射・エミッタンスなど約 200 本) は numba の `@njit(cache=True)` で書いてある。
最初に呼ばれたときに機械語へ変換する (CPU 依存、全体で 1 分ほど)。変換した結果はディスクにキャッシュされ、2 回目以降は読み込むだけ (1 秒以下)。
`pysad-warmup` は、インストールの直後に 1 回、小さなラティスで全種類の計算を走らせて、このキャッシュを先に作る。

キャッシュの場所は numba が決める。
  1. 環境変数 `NUMBA_CACHE_DIR` があればその下。
  2. なければ `pysad/__pycache__` (パッケージのフォルダが書き込めるとき)。
  3. 書き込めなければ (site-packages が読み取り専用など) ユーザーのキャッシュのフォルダ (Windows は %LOCALAPPDATA%、Linux は ~/.cache)。
実行時も同じ規則で探すので、`pysad-warmup` と台本の実行で `NUMBA_CACHE_DIR` を揃えること。キャッシュは CPU の種類 (命令セット) と numba の版ごと。別の CPU・版では使われず、再変換される。

使い方 (コマンドの `pysad-warmup` は `py -m pysad.warmup` でもよい):
  pysad-warmup           全カーネルを変換してキャッシュを作る (作ってあるものは読み込むだけで速い)
  pysad-warmup --check   キャッシュが有効か調べる (有効なら終了コード 0)。変換はしない
  pysad-warmup --clear   pysad のキャッシュを消す
  pysad-warmup --quiet   進み具合を表示しない (最後の要約だけ)
"""
from __future__ import annotations

import importlib
import io
import json
import os
import re
import sys
import tempfile
import time

_PKG = os.path.dirname(os.path.abspath(__file__))
RECORD_NAME = 'pysad_warmup.json'

# 全種類のカーネルを呼ぶ最小の入力。ラティスは台本の中に埋め込んである (example/ や外部のファイルに依存しない)。
# 3 GeV、12 セルの FODO (周長 42 m) を土台に、要素ごとに小さな台本に分ける (どれかが失敗しても残りは走る)。
_HEAD = """MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.5235987755982988 E1=0.2617993877991494 E2=0.2617993877991494);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=1.5E6 HARM=48 PHI=0);
"""
_BEAM = "b = {{1e-4,2e-4},{0,0},{1e-4,2e-4},{0,0},{0,0},{0,1e-3},{1,1}};\n"
# RF なし・放射なし、RF あり・放射あり (ゆらぎつき) の順に追跡し、Emittance[] を呼ぶ
_TRACK1 = "NORFSW; NORAD; NOFLUC;\nr = TrackParticles[{1, b}, 1, 1, 3];\nRFSW; RAD; FLUC;\nr = TrackParticles[{1, b}, 1, 1, 3];\n"
_TRACK = _TRACK1 + "RAD; NOFLUC;\ne = Emittance[];\n"
_RING = "LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);\nLINE CELLX = (%s);\nLINE RING = (ORIGIN RF CELLX 11*CELL);\n"


def _ring(defs, cellx, body=_TRACK):
    """12 セルの FODO のうち 1 セルを cellx に替えたリングで、CALC → body を走らせる台本"""
    return _HEAD + defs + _RING % cellx + "FFS;\nUSE RING; CELL; CALC;\n" + _BEAM + body + "Exit;\n"


_WAKE = ('WakeFunction[Longitudinal, "%s"] = {{-1, 1e9}, {0, 1e9}, {1e-9, 0}, {1, 0}, {2, 0}};\n'
         'WakeFunction[Transverse, "%s"] = {{-1, 1e12}, {0, 1e12}, {1e-9, 0}, {1, 0}, {2, 0}};\n')
_SOLDEF = "SOL SA = (L=0 BZ=1 BOUND=1 GEO=1%s) SB = (L=0 BZ=0 BOUND=1%s);\n"
_B30 = "E1=0.2617993877991494 E2=0.2617993877991494"

SCRIPTS = [
    ('光学・エミッタンス・追跡', _HEAD + """LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);
LINE RING = (ORIGIN RF 12*CELL);
FFS;
USE RING; CELL; CALCULATE;
RFSW; RAD;
e = Emittance[];
EMIT;
""" + _BEAM + _TRACK1 + """LOSSMAP; LOSSAMPL = 0.05;
r = TrackParticles[{1, b}, 1, 1, 3];
Exit;
"""),
    ('偏向磁石 (K1・フリンジ)',
     _ring("BEND BK = (L=1.2 ANGLE=0.5235987755982988 K1=0.05 " + _B30 + " F1=0.04 FRINGE=3);\n", "SF DS QF D2 BK D2 QD DS SD D1")),
    ('ステアリング (K0)',
     _ring("BEND SA = (L=0.2 ANGLE=0 K0=0.001 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.02 AE2=-0.01);\n"
           "BEND ST = (L=0 ANGLE=0 K0=1e-5);\n", "SF DS QF D2 B D2 QD DS SD SA ST D1")),
    ('MULT・OCT・APERT',
     _ring("MULT MM = (L=0.1 K1=0 K2=0.5 K3=10);\nOCT OC = (L=0.1 K3=1);\nAPERT AP = (AX=0.05 AY=0.05);\n"
           "MULT MB = (L=1.2 ANGLE=0.5235987755982988 " + _B30 + " K1=0.01);\n", "SF DS QF D2 MB D2 QD DS SD MM OC AP D1")),
    ('COORD・TCAVI',
     _ring("COORD CI = (DX=1e-4 DY=-2e-4 CHI3=0.001 DIR=0) CO = (DX=1e-4 DY=-2e-4 CHI3=0.001 DIR=1);\n"
           "TCAVI TC = (L=0.2 K0=3e-4 FREQ=4.5e8 PHI=0);\n", "SF DS QF D2 CI B CO D2 QD DS SD D1 TC")),
    ('UND・PHSROT',
     _ring("UND U1 = (L=1 FBX=1.0 POLE=20 SLICE=200) U2 = (L=1 POLE=20 SLICE=200);\n"
           "PHSROT PR = (AX=0 BX=4 PX=0.001 AY=0 BY=2 PY=0.001 AZ=0 BZ=10 PZ=0.001 D11=1e-3 D22=1e-3 D33=1e-3 D44=1e-3 "
           "D55=1e-3 D66=1e-3 B11=1e-11 B22=1e-11 B33=1e-13 B44=1e-13 B55=1e-8 B66=1e-8);\n",
           "SF DS QF D2 B D2 QD DS SD D1 U1 U2 PR")),
    ('CAVI の横方向電圧・DX',
     _ring("CAVI RFL = (L=0.3 VOLT=1.5E5 HARM=48 PHI=0 DX=1e-4 V1=1e3 V20=1e3 V11=1e3 V02=1e3);\n"
           "BEND SU = (L=0.2 ANGLE=0 K0=1e-5);\n", "SF DS QF D2 B D2 QD DS SD D1 RFL SU")),
    ('BEAMBEAM・WAKE',
     _ring("BEAMBEAM BB = (BX=1 BY=1 EMIX=1E-8 EMIY=1E-8 SIGZ=1E-3 DP=1E-3 SLICE=2 NP=1E9);\nMARK WK = ();\n"
           "CAVI RFW = (L=0.3 VOLT=1.5E6 HARM=48 PHI=0);\n", "BB WK RFW SF DS QF D2 B D2 QD DS SD D1",
           "NORAD; NOFLUC; PBUNCH = 1e10; NBUNCH = 1;\n" + _WAKE % ('WK', 'WK') + _WAKE % ('RFW', 'RFW') + """np = 20;
x0 = Table[1e-4 + 1e-5*Sin[0.37*i], {i, np}];
z0 = Table[0.01*Sin[0.13*i], {i, np}];
zero = Table[0, {i, np}];
beam = {x0, zero, zero, zero, z0, zero, Table[1, {i, np}]};
LWAKE; TWAKE; RFSW;
r = TrackParticles[{1, beam}, 1, 1, 2][[2]];
TRPT; INS;
r = TrackParticles[{1, beam}, 1, 1, 2][[2]];
""")),
    ('ソレノイド領域',
     _ring(_SOLDEF % (' F1=0.05 DX=1e-4 CHI3=0.001', ' F1=0.05') + "MULT MQ = (L=0.3 K1=0.01 K2=0.5);\nMULT MK = (L=0.2 K0=1e-5);\n"
           "BEND SU = (L=0.2 ANGLE=0 K0=1e-5);\n"
           "BEND SV = (L=0.2 ANGLE=0 K0=1e-5 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.02 AE2=-0.01);\n"
           "BEND B2 = (L=1.2 ANGLE=0.5235987755982988 " + _B30 + " F1=0.04 FRINGE=3);\n", "QF D2 SA B2 D2 SU MK SV QD DS D1 SB DS")),
    ('加速 (TRPT)・RF つき MULT・SOL 内の空洞', """MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3) M1=();
DRIFT D1 = (L=1.0) DS=(L=0.3);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25) QS=(L=0.2 K1=0.1);
CAVI CA = (L=0.5 VOLT=2E8 FREQ=1.3E9 PHI=-1.2);
MULT MA = (L=0.4 VOLT=1E8 FREQ=1.3E9 PHI=-1.0 K1=0.05);
""" + _SOLDEF % ('', '') + """LINE TL = (M0 QF D1 DS CA DS QS DS MA DS D1 QD D1 SA DS CA DS MA DS SB D1 M1);
FFS;
USE TL;
TRPT; RFSW; INS; CALC;
""" + _BEAM + """r = TrackParticles[{1, b}, 1, 1, 1];
RAD; FLUC;
r = TrackParticles[{1, b}, 1, 1, 1];
NORAD; NOFLUC; PBUNCH = 1e10; NBUNCH = 1;
""" + _WAKE % ('CA', 'CA') + """LWAKE; TWAKE;
r = TrackParticles[{1, b}, 1, 1, 1];
Exit;
"""),
]


def kernel_modules():
    """Return the names of the pysad modules that contain numba kernels (those using njit).
    numba のカーネルを持つモジュール名 (pysad/*.py のうち njit を使うもの)"""
    from . import nbcache
    out = []
    for f in sorted(os.listdir(_PKG)):
        if f.endswith('.py') and f[:-3] not in ('nbcache', 'cli', 'warmup'):
            if nbcache._kernel_module(f[:-3]):
                out.append(f[:-3])
    return out


def dispatchers():
    """Return [(module name, function name, numba Dispatcher)] for kernels visible at module level.
    [(モジュール名, 関数名, numba の Dispatcher)]。Python から見える (モジュール直下の) カーネルだけ"""
    from numba.core.registry import CPUDispatcher
    out = []
    for m in kernel_modules():
        try:
            mod = importlib.import_module('pysad.' + m)
        except Exception:
            continue
        for k, v in sorted(vars(mod).items()):
            if isinstance(v, CPUDispatcher) and getattr(v.py_func, '__module__', '') == 'pysad.' + m:
                out.append((m, k, v))
    return out


def _magic():
    from numba.core.registry import cpu_target
    return cpu_target.target_context.codegen().magic_tuple()


def _cache_file(disp):
    c = getattr(disp, '_cache', None)
    return getattr(c, '_cache_file', None)


def kernel_cached(disp, magic=None):
    """Tell whether a usable cache exists on disk for this CPU, numba version and source content.
    この CPU・この版の numba・この内容のソースで使えるキャッシュが、ディスクにあるか"""
    cf = _cache_file(disp)
    if cf is None:
        return False
    magic = magic or _magic()
    try:
        idx = cf._load_index()
    except Exception:
        return False
    for key, fname in idx.items():
        if key[1] == magic and os.path.exists(os.path.join(cf._cache_path, fname)):
            return True
    return False


def cache_dirs(ds=None):
    """Return the numba cache folders of the pysad kernels (a list of paths).
    pysad のカーネルの numba キャッシュのフォルダ (パスの並び) を返す。"""
    ds = dispatchers() if ds is None else ds
    out = []
    for _, _, d in ds:
        cf = _cache_file(d)
        if cf is not None and cf._cache_path not in out:
            out.append(cf._cache_path)
    return out


def _cache_files(dirs):
    """dirs の中の pysad のカーネルのキャッシュ (.nbi・.nbc) の絶対パス"""
    mods = kernel_modules()
    pat = re.compile(r'^(%s)\.[^/\\]*\.(nbi|nbc)$' % '|'.join(re.escape(m) for m in mods))
    out = []
    for d in dirs:
        if os.path.isdir(d):
            out += [os.path.join(d, f) for f in sorted(os.listdir(d)) if pat.match(f)]
    return out


def cache_size(dirs):
    """Return the total size in bytes of the cache files in `dirs`.
    `dirs` の中のキャッシュのファイルの合計の大きさ (バイト) を返す。"""
    return sum(os.path.getsize(f) for f in _cache_files(dirs))


def _fmt_size(n):
    return '%.1f MB' % (n / 1e6)


def _record_path(dirs):
    return os.path.join(dirs[0], RECORD_NAME) if dirs else None


def _read_record(dirs):
    for d in dirs:
        try:
            with open(os.path.join(d, RECORD_NAME), encoding='utf-8') as fh:
                return json.load(fh)
        except (OSError, ValueError):
            continue
    return None


def run_scripts(scripts=None, verbose=True, out=None):
    """Run each embedded script once; returns [(name, seconds, error line or None)].
    埋め込みの台本を 1 回ずつ走らせる。戻り値: [(名前, 秒, エラーの行 or None)]"""
    from .ffscmd import run_sad_string
    out = out or sys.stdout
    res = []
    old = os.getcwd()
    with tempfile.TemporaryDirectory(prefix='pysad_warmup_') as tmp:
        os.chdir(tmp)
        try:
            for name, txt in (SCRIPTS if scripts is None else scripts):
                t = time.perf_counter()
                buf = io.StringIO()
                err = None
                try:
                    run_sad_string(txt, out=buf, cwd=tmp)
                except BaseException as ex:      # どの台本が失敗しても、残りは走らせる
                    err = '%s: %s' % (type(ex).__name__, str(ex)[:100])
                dt = time.perf_counter() - t
                res.append((name, dt, err))
                if verbose:
                    print('  %-28s %6.1f 秒%s' % (name, dt, ('  失敗: ' + err) if err else ''), file=out, flush=True)
        finally:
            os.chdir(old)
    return res


def status(ds=None):
    """Return {module: (cached kernels, kernels)} and the overall (cached, total) counts.
    {モジュール: (キャッシュのあるカーネル数, カーネル数)} と、全体の (あり, 全体)"""
    ds = dispatchers() if ds is None else ds
    magic = _magic()
    per = {}
    for m, _, d in ds:
        a, n = per.get(m, (0, 0))
        per[m] = (a + (1 if kernel_cached(d, magic) else 0), n + 1)
    return per


def do_clear(out=None):
    """Delete the pysad numba cache files and their records (the warmup record and the dependency record), and print what was removed. Returns 0.
    pysad の numba キャッシュのファイルと記録 (変換の記録・依存の記録) を消し、消した数と大きさを表示する。0 を返す。"""
    out = out or sys.stdout
    ds = dispatchers()
    dirs = cache_dirs(ds)
    files = _cache_files(dirs)
    size = sum(os.path.getsize(f) for f in files)
    for d in dirs:                                   # 記録 (変換の記録・依存の記録) も消す
        for extra in (RECORD_NAME, 'nbdeps.json'):
            p = os.path.join(d, extra)
            if os.path.exists(p):
                files.append(p)
    n = 0
    for f in files:
        try:
            os.remove(f)
            n += 1
        except OSError as ex:
            print('消せません: %s (%s)' % (f, ex), file=out)
    print('pysad の numba キャッシュを消した: %d ファイル (%s)' % (n, _fmt_size(size)), file=out)
    for d in dirs:
        print('  場所: ' + d, file=out)
    return 0


def do_check(out=None):
    """Return 0 if the cache is valid, 1 otherwise.
    キャッシュが有効なら 0、そうでなければ 1"""
    out = out or sys.stdout
    ds = dispatchers()
    dirs = cache_dirs(ds)
    per = status(ds)
    tot = sum(n for _, n in per.values())
    have = sum(a for a, _ in per.values())
    rec = _read_record(dirs)
    print('numba のキャッシュ (pysad のカーネル %d 本、%d モジュール)' % (tot, len(per)), file=out)
    for d in dirs:
        print('  場所: %s' % d, file=out)
    print('  大きさ: %s、キャッシュのあるカーネル: %d / %d' % (_fmt_size(cache_size(dirs)), have, tot), file=out)
    from . import nbcache
    try:
        stale = nbcache.stale_modules()
    except Exception:
        stale = {}
    for m, dep in sorted(stale.items()):
        print('  依存先が変わって古い: %s (%s)' % (m, ', '.join(dep)), file=out)
    if rec is None:
        print('  事前コンパイルの記録がありません: `pysad-warmup` を実行してください', file=out)
        return 1
    names = set(rec.get('kernels', []))
    magic = _magic()
    missing = sorted('%s.%s' % (m, k) for m, k, d in ds if '%s.%s' % (m, k) in names and not kernel_cached(d, magic))
    new = sorted('%s.%s' % (m, k) for m, k, d in ds if '%s.%s' % (m, k) not in names)
    import numba
    if rec.get('numba') != numba.__version__:
        print('  numba の版が違います (記録 %s、現在 %s)' % (rec.get('numba'), numba.__version__), file=out)
    ok = not missing and not stale and rec.get('numba') == numba.__version__
    if missing:
        print('  無効 (CPU・numba・ソースが変わった、または消された) のカーネル: %d 本 (例 %s)' % (len(missing), ', '.join(missing[:3])), file=out)
    print('  記録した事前コンパイル: %d 本 (%s)。記録後に増えたカーネル: %d 本' % (len(names), rec.get('time', '?'), len(new)), file=out)
    print('  判定: ' + ('有効' if ok else '無効 (`pysad-warmup` を実行してください)'), file=out)
    return 0 if ok else 1


def do_warmup(verbose=True, out=None):
    """Run every pysad numba kernel once on small embedded lattices to create the cache (an existing cache is only loaded), write the record used by `--check`, and print a summary. Returns 0, or 1 if some script failed.
    全カーネルを小さな埋め込みのラティスで 1 回ずつ走らせて numba のキャッシュを作り (あるキャッシュは読み込むだけ)、`--check` が使う記録を書いて要約を表示する。失敗した台本があれば 1、なければ 0 を返す。"""
    out = out or sys.stdout
    import numba
    t0 = time.perf_counter()
    ds0 = dispatchers()
    dirs = cache_dirs(ds0)
    before = status(ds0)
    nb0 = sum(a for a, _ in before.values())
    if verbose:
        print('pysad の numba カーネルを事前コンパイルする (カーネル %d 本。作ってあるキャッシュは読み込むだけ)' % len(ds0), file=out, flush=True)
        if os.environ.get('NUMBA_CACHE_DIR'):
            print('  NUMBA_CACHE_DIR = ' + os.environ['NUMBA_CACHE_DIR'], file=out)
    try:
        from . import nbcache
        nbcache.purge()                              # 依存先が変わって古いキャッシュを先に捨てる (書き込めないときは何もしない)
    except OSError:
        pass
    res = run_scripts(verbose=verbose, out=out)
    magic = _magic()
    covered = sorted('%s.%s' % (m, k) for m, k, d in ds0 if kernel_cached(d, magic))
    uncovered = sorted('%s.%s' % (m, k) for m, k, d in ds0 if not kernel_cached(d, magic))
    dirs = cache_dirs(ds0)
    try:
        from . import nbcache
        nbcache.purge()                              # 依存先のハッシュを記録する
    except OSError:
        pass
    wrote = False
    if dirs:
        try:
            with open(_record_path(dirs), 'w', encoding='utf-8') as fh:
                json.dump({'numba': numba.__version__, 'kernels': covered, 'time': time.strftime('%Y-%m-%d %H:%M:%S'),
                           'seconds': round(time.perf_counter() - t0, 1)}, fh, ensure_ascii=False, indent=0)
            wrote = True
        except OSError:
            wrote = False
    dt = time.perf_counter() - t0
    fails = [r for r in res if r[2]]
    print('完了: %.1f 秒。キャッシュのあるカーネル %d / %d (実行前は %d)' % (dt, len(covered), len(ds0), nb0), file=out)
    print('  キャッシュの大きさ: %s' % _fmt_size(cache_size(dirs)), file=out)
    for d in dirs:
        print('  場所: ' + d, file=out)
    if uncovered:
        print('  この入力では呼ばれない (最初に使うときに変換される) カーネル: %d 本 (%s)' % (len(uncovered), ', '.join(uncovered[:6]) + (' ...' if len(uncovered) > 6 else '')), file=out)
    if fails:
        print('  失敗した台本: ' + ', '.join(n for n, _, _ in fails), file=out)
    if not wrote:
        print('  (記録を書けませんでした: `--check` は使えません)', file=out)
    return 1 if fails else 0


def main(argv=None):
    """Command `pysad-warmup [--check | --clear] [--quiet]`: precompile the numba kernels (create the cache), check it, or clear it.
    コマンド `pysad-warmup [--check | --clear] [--quiet]`: numba のカーネルを事前に変換してキャッシュを作る。`--check` で有効か調べ、`--clear` で消す。"""
    args = list(sys.argv[1:] if argv is None else argv)
    if '-h' in args or '--help' in args:
        print('usage: pysad-warmup [--check | --clear] [--quiet]\n'
              '  (none)   run every kernel once on a small lattice to create the numba cache\n'
              '  --check  check whether the cache is valid (exit code 0 if valid); no compilation\n'
              '  --clear  delete the pysad numba cache\n'
              '  --quiet  do not show progress\n'
              '使い方: pysad-warmup [--check | --clear] [--quiet]\n'
              '  (なし)   全カーネルを小さなラティスで 1 回ずつ走らせ、numba のキャッシュを作る\n'
              '  --check  キャッシュが有効か調べる (有効なら終了コード 0)。変換はしない\n'
              '  --clear  pysad の numba キャッシュを消す\n'
              '  --quiet  進み具合を表示しない')
        return 0
    known = {'--check', '--clear', '--quiet'}
    bad = [a for a in args if a not in known]
    if bad or ('--check' in args and '--clear' in args):
        print('引数が正しくありません: %s (--help を見てください)' % ' '.join(bad or ['--check と --clear は同時に使えません']), file=sys.stderr)
        return 2
    if '--clear' in args:
        return do_clear()
    if '--check' in args:
        return do_check()
    return do_warmup(verbose='--quiet' not in args)


if __name__ == '__main__':
    sys.exit(main())
