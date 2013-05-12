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

# CONSOLE PRINTING UTILITY METHODS
def green(s):
    return '\033[1;32m%s\033[m' % s

def red(s):
    return '\033[1;31m%s\033[m' % s

def log(*m):
    print >> sys.stderr, " ".join(m)

# META TESTS
class Unit_Test_Existence_Test(unittest.TestCase):
    pass

def tg_test_unit_test_existence(name):
    """
    Generates a test that checks whether a unit test class exists
    for a given benchmark.
    """
    def test_unit_test_existence(self):
        for case in get_testcases():
            if name.lower() in case.lower():
                return
        self.assertTrue(False, "no unit tests for %s" % name)
    return ("test_%s" % name, test_unit_test_existence)

# GENERAL TEST CLASSES
class Existence_Test:
    name = None

    def file_exists(self):
        try:
            with open(self.name):
                return True
        except IOError:
            return False

    def test_existence(self):
        self.assertTrue(self.name != None, "must specify name field!")
        self.assertTrue(self.file_exists(), "%s binary not found" % self.name)

class Standard_Test:
    name = None
    inputs = None

    def run_stencil(self, inp):
        proc = subprocess.Popen(["./%s" % self.name] + [str(i) for i in inp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = proc.communicate()[0]

        if proc.returncode != 0:
            return False

        if "failed" in out.lower():
            return False

        return True

    def test_benchmark(self):
        self.assertTrue(self.name != None, "must specify name field!")
        self.assertTrue(self.inputs != None, "must specify inputs field!")

        for inp in self.inputs:
            self.assertTrue(
                self.run_stencil(inp), 
                "execution failure for %s %s" % (self.name, " ".join([str(i) for i in inp]))
            )

# BENCHMARK TEST CLASSES
class Fib_Test(unittest.TestCase, Existence_Test, Standard_Test):
    name = "fib"
    inputs = [(0,), (10,), (25,)]

class Life_Test(unittest.TestCase, Existence_Test, Standard_Test):
    name = "life"
    inputs = [(0, 0), (10, 0), (10, 10), (100, 10), (100, 100)]

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
        test_name, test = tg_test_unit_test_existence(benchmark)
        setattr(Unit_Test_Existence_Test, test_name, test)

    unittest.main()
