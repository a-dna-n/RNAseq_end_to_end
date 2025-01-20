import os
import sys
def extend_path():
    file_dir = os.path.dirname(os.path.realpath(__file__))
    parent_dir = os.path.dirname(file_dir)
    sys.path.append(parent_dir)
