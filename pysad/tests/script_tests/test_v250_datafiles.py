# -*- coding: utf-8 -*-
"""pysad に同梱した小さなデータファイル (pysad/data/、pysad/datafiles.py) の試験。

KEK-LS (3 GeV、20 セル) のラティス `kek-ls_v03_12.lat` が、パッケージの中にあり (pip で入る)、`pysad/tests/` のものと同じ内容で、
pip では、パッケージの中 (pysad/data) にだけ入る (例 (example/) は pip に入れないので、例の相対パスの位置には置かない)。
"""
import filecmp
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
NAME = 'kek-ls_v03_12.lat'


def test_data_path_and_list():
    from pysad import datafiles
    p = datafiles.data_path(NAME)
    assert os.path.isfile(p) and NAME in datafiles.list_data()
    assert os.path.abspath(p).replace(chr(92), '/').endswith('pysad/data/' + NAME)
    try:
        datafiles.data_path('no_such_file.lat')
    except FileNotFoundError as e:
        assert NAME in str(e)
    else:
        raise AssertionError('FileNotFoundError にならない')


def test_bundled_lattice_is_the_same_as_the_tests_copy():
    from pysad import datafiles
    assert filecmp.cmp(datafiles.data_path(NAME), os.path.join(ROOT, 'pysad', 'tests', NAME), shallow=False)


def test_bundled_lattice_loads_and_has_cells():
    from pysad import datafiles, load_sad_file
    r = load_sad_file(datafiles.data_path(NAME))
    assert len(r.element_names()) > 30
    assert 'KEK-LS' in open(datafiles.data_path(NAME), encoding='utf-8', errors='replace').read()


def test_packaging_includes_the_data_file():
    toml = open(os.path.join(ROOT, 'pyproject.toml'), encoding='utf-8').read()
    assert 'package-data' in toml and '"data/*"' in toml, 'package-data に data/* が無い'
    man = open(os.path.join(ROOT, 'MANIFEST.in'), encoding='utf-8').read()
    assert 'graft pysad/data' in man


def test_examples_refer_to_the_lattice_by_the_documented_relative_path():
    """例が参照するのは pysad/tests/kek-ls_v03_12.lat (リポジトリの根、または <prefix>/share/pysad をカレントにして走らせる)"""
    ex = os.path.join(ROOT, 'example', '02_optics_cod', 'fn_MEASURE_ring_a.sad')
    t = open(ex, encoding='utf-8', errors='replace').read()
    assert 'pysad/tests/' + NAME in t
    assert os.path.isfile(os.path.join(ROOT, 'pysad', 'tests', NAME))
