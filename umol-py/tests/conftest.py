import sys


collect_ignore = []
if sys.version_info < (3, 10):
    collect_ignore.append("test_structural_pattern_matching.py")

