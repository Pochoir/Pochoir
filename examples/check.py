import subprocess
import sys
import re
import glob
import os
import re

def green(s):
    return '\033[1;32m%s\033[m' % s

def red(s):
    return '\033[1;31m%s\033[m' % s

def log(*m):
    print >> sys.stderr, " ".join(m)

def check(name):
    if check_existence(name) and check_correctness(name):
        log(green("PASS"), "%s correctness check" % name)

def check_existence(name):
    try:
        with open(name):
            return True
    except IOError:
        log(red("FAIL"), "%s not found" % name)
        return False

def check_correctness(name):
    return True

if __name__ == "__main__":
    stencils = [os.path.splitext(filename)[0] for filename in glob.glob('*.cpp') if "tb_" in filename]
    stencils = set([re.sub("_pochoir", "", stencil) for stencil in stencils])
    stencils = set([re.sub("_kernel_info", "", stencil) for stencil in stencils])
    stencils = set([re.sub("tb_", "", stencil) for stencil in stencils])

    for stencil in sorted(stencils):
        check(stencil)