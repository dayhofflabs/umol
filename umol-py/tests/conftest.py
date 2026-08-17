import sys


collect_ignore = []
if sys.version_info < (3, 10):
    collect_ignore.append("test_class_constructor_signatures.py")
    collect_ignore.append("test_structural_pattern_matching.py")
