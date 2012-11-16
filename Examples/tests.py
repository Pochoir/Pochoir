# Tests
import subprocess
import sys
import re

def green(s):
    return '\033[1;32m%s\033[m' % s

def red(s):
    return '\033[1;31m%s\033[m' % s

def log(*m):
    print >> sys.stderr, " ".join(m)

def check_tn_stencil(name, cmd):
    test_pass = False
    for N_SIZE, T_SIZE in ((100, 0), (100, 1), (100, 10), (100, 100)):
        f = open("%s_%s_%s.out" % (name, N_SIZE, T_SIZE), "w+")
        p = subprocess.Popen([cmd, str(N_SIZE), str(T_SIZE)], stdout=f)
        retcode = p.wait()
        f.flush()       
        f.seek(0)

        out = f.read()
        if re.search("PASSED", out):
            test_pass = True

        if re.search("FAILED", out) or retcode != 0:
            test_pass = False
            log(red("FAIL"), "%s correctness check for N_SIZE=%s T_SIZE=%s" % (name, N_SIZE, T_SIZE))

        f.close()

    if test_pass:
        log(green("PASS"), "%s correctness check" % name)
    else:
        log(red("FAIL"), "%s correctness check" % name)


if __name__ == "__main__":
    tn_stencils = ["heat_1D_NP",
                   "heat_2D_NP",
                   "heat_2D_NP_zero",
                   "heat_2D_P",
                   "heat_3D_NP",
                   "life"]

    for stencil in tn_stencils:
        check_tn_stencil(stencil, "./tb_%s_pochoir" % stencil)