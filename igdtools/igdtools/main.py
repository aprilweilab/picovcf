import subprocess
import sys
import os

THISDIR = os.path.dirname(os.path.realpath(__file__))

def main():
    try:
        subprocess.check_call([os.path.join(THISDIR, "igdtools")] + sys.argv[1:])
    except subprocess.CalledProcessError as e:
        return e.returncode
    return 0
