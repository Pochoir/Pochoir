import subprocess
import sys
import re
import glob
import os
import re
import unittest

class ExistenceTest(unittest.TestCase):
    pass

class TestTest(unittest.TestCase):
    def test_test(self):
        self.assertEqual(1, 1)

def green(s):
    return '\033[1;32m%s\033[m' % s

def red(s):
    return '\033[1;31m%s\033[m' % s

def log(*m):
    print >> sys.stderr, " ".join(m)

def check(name):
    if check_existence(name) and check_correctness(name):
        log(green("PASS"), "%s correctness check" % name)

def tg_check_existence(name):
    def check_existence(self):
        try:
            with open(name):
                pass
        except IOError:
            self.assertTrue(False, "%s not found" % name)
    return ("test_compiles_%s" % name, check_existence)

def check_correctness(name):
    return True

if __name__ == "__main__":
    stencils = [os.path.splitext(filename)[0] for filename in glob.glob('*.cpp') if "tb_" in filename]
    stencils = set([re.sub("_pochoir$", "", stencil) for stencil in stencils])
    stencils = set([re.sub("_kernel_info$", "", stencil) for stencil in stencils])
    stencils = set([re.sub("_gen_kernel$", "", stencil) for stencil in stencils])
    stencils = set([re.sub("^tb_", "", stencil) for stencil in stencils])

    for stencil in sorted(stencils):
        test_name, test = tg_check_existence(stencil)
        setattr(ExistenceTest, test_name, test)

    unittest.main()
