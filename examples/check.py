import subprocess
import sys
import re
import glob
import os
import re
import unittest
import inspect

def get_benchmarks():
    stencils = [os.path.splitext(filename)[0] for filename in glob.glob('*.cpp') if "tb_" in filename]
    stencils = set([re.sub("_pochoir$", "", stencil) for stencil in stencils])
    stencils = set([re.sub("_kernel_info$", "", stencil) for stencil in stencils])
    stencils = set([re.sub("_gen_kernel$", "", stencil) for stencil in stencils])
    stencils = set([re.sub("^tb_", "", stencil) for stencil in stencils])
    return stencils

def get_testcases():
    return [classname for (classname, c) in inspect.getmembers(sys.modules[__name__], inspect.isclass)]

def file_exists(filename):
    try:
        with open(filename):
            return True
    except IOError:
        return False

# CONSOLE PRINTING UTILITY METHODS
def green(s):
    return '\033[1;32m%s\033[m' % s

def red(s):
    return '\033[1;31m%s\033[m' % s

def log(*m):
    print >> sys.stderr, " ".join(m)

# META TESTS
class Test_Existence_Test(unittest.TestCase):
    pass

def tg_test_test_existence(name):
    """
    Generates a test that ensures that a test suite class exists
    for a given benchmark.
    """
    def test_test_existence(self):
        for case in get_testcases():
            if name.lower() in case.lower():
                return
        self.assertTrue(False, "no unit tests for %s" % name)
    return ("test_test_existence_%s" % name, test_test_existence)

# GENERAL TEST CLASSES
class Existence_Test:
    def test_existence(self):
        self.assertTrue(file_exists(self.name), "%s not found" % self.name)

# BENCHMARK TEST CLASSES
class Fib_Test(unittest.TestCase, Existence_Test):
    name = "fib"

class Life_Test(unittest.TestCase, Existence_Test):
    name = "life"

class Heat_1D_NP_Test(unittest.TestCase, Existence_Test):
    name = "heat_1D_NP"

class Heat_2D_NP_Test(unittest.TestCase, Existence_Test):
    name = "heat_2D_NP"

class Heat_2D_NP_Zero_Test(unittest.TestCase, Existence_Test):
    name = "heat_2D_NP_zero"

class Heat_2D_NP_Ref_Test(unittest.TestCase, Existence_Test):
    name = "heat_2D_NP_ref"

class Heat_2D_P_Test(unittest.TestCase, Existence_Test):
    name = "heat_2D_P"

class Heat_2D_P_Diff_Slope_Test(unittest.TestCase, Existence_Test):
    name = "heat_2D_P_diff_slope"

class Heat_2D_P_Macro_H_Test(unittest.TestCase, Existence_Test):
    name = "heat_2D_P_macro_h"

class Heat_3D_NP_Test(unittest.TestCase, Existence_Test):
    name = "heat_3D_NP"

class Overlap_Tile_2D_Test(unittest.TestCase, Existence_Test):
    name = "overlap_tile_2D"

class Overlap_Tile_3D_Test(unittest.TestCase, Existence_Test):
    name = "overlap_tile_3D"

if __name__ == "__main__":
    # dynamically generate test cases for all benchmarks
    # to ensure that each benchmark in the examples directory 
    # has a corresponding test case
    for benchmark in get_benchmarks():
        test_name, test = tg_test_test_existence(benchmark)
        setattr(Test_Existence_Test, test_name, test)

    unittest.main()
