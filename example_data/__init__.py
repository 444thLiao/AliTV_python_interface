import sys
from os.path import dirname, join

example_dir = dirname(__file__)
lib_dir = dirname(dirname(__file__))
sys.path.insert(0, lib_dir)
