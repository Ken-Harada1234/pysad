# -*- coding: utf-8 -*-
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'example', 'tools'))
import run_example
if __name__ == '__main__':
    to = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 120
    t = time.time(); run_example.run(sys.argv[1], to); print(sys.argv[1], 'sec', round(time.time() - t))
