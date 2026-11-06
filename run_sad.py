#!/usr/bin/env python3
"""run_sad.py — SAD スクリプト (MAIN 部 + FFS 部) を pysad で実行する。
使い方: python run_sad.py lattice_file [script_file ...]   (SAD の `gs` に複数ファイルを連結して渡すのと同じ)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pysad.ffscmd import run_sad

if __name__ == '__main__':
    run_sad(sys.argv[1:])
