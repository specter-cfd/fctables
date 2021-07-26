"""
Python3 convenience script for creating several FC-Gram
tables using FC-Tables.

2021 - Mauro Fontana - University of Buenos Aires
"""
import os
import signal
import itertools

_DIGITS = 1024            # Decimal digits of precision to use
_BLEND  = True            # True = also generate blend to zero operators (A)


ds = list(range(3,10))    # Values of d to consider
Cs = list(range(15, 35))  # Values of C to consider

Neumann  = True           # True = also build Q1n
Neumann2 = True           # True = also build Q2n

params = itertools.product(reversed(ds), reversed(Cs))

_DIGITS = str(_DIGITS)
if _BLEND:
    _BLEND = 'yes'
else:
    _BLEND = 'no'

def build_fctables(digits=256, blend=False):
    """
    Builds FC-Tables with the appropriate digits and blend to zero vars.
    All the rest of the Makefile variables are retained.
    Current `Makefile.in` is saved as `.Makefile.in` and restored afterwards.
    """
    from shutil import copy2 as copy
    import re
    import os
    import subprocess

    print("Building FC-Tables...", end=" ", flush=True)
    copy('Makefile.in', '.Makefile.in')
    with open("Makefile.in", 'r') as fh:
        makefile = fh.read()
        new_makefile = re.sub("^DIGITS.*", "DIGITS="+_DIGITS, makefile,
                              flags=re.MULTILINE)
        new_makefile = re.sub("^BLEND.*", "BLEND="+_BLEND, new_makefile,
                              flags=re.MULTILINE)

    try:
        with open('Makefile.in','w') as fh:
            fh.write(new_makefile)


        make_result = subprocess.run('make', stdout=subprocess.DEVNULL,
                                     stderr=subprocess.PIPE, text=True)

        if make_result.returncode == 0:
            print("Sucessfully built FC-Tables")
        else:
            print("Building FC-Tables failed with error:")
            print(make_result.stderr)

    finally:
        copy('.Makefile.in', 'Makefile.in')
        os.remove('.Makefile.in')

def backup_parameters():
    """
    Backups the current `parameter.inp` to `.parameter.inp`.
    """
    from shutil import copy2 as copy
    copy('parameter.inp', '.parameter.inp')


def restore_parameters(*args):
    """
    Restores `.parameter.inp` to `parameter.inp`.
    *args are optional arguments added to make the function API
    compliant with `signal.signal`, but are unused.
    """
    import os
    from shutil import copy2 as copy
    copy('.parameter.inp', 'parameter.inp')
    os.remove('.parameter.inp')


def set_parameters(**kwargs):
    """
    Sets the desired value for kwarg in `parameter.inp`.
    """
    from shutil import copy2 as copy
    import re

    with open("parameter.inp", 'r') as fh:
        parameter = fh.read()

        for key, value in kwargs.items():
            search  = "^{}\s*=.*".format(key)
            replace = "{}={}".format(key, value)
            parameter = re.sub(search, replace, parameter,
                               flags=re.MULTILINE | re.IGNORECASE)

        with open('parameter.inp','w') as fh:
            fh.write(parameter)


def compute_table():
    """
    Runs FC-Tables to compute the desired table. It uses a pipe to capture
    both the retval and the SVD error in case the blend-to-zero operator was
    computed (which is the most numerically challenging task).
    """
    import subprocess
    import re

    err_pattern = 'Error in SVD decomposition:(.+)'
    pipe = subprocess.Popen('./fc_tables', stdout=subprocess.PIPE, text=True)

    error = None
    for line in pipe.stdout:
        capture = re.search(err_pattern, line)
        if capture is not None:
            error = float(capture.group(1))

    pipe.wait()
    ret = True if pipe.returncode == 0 else False

    return ret, error


def create_report():
    """
    Helper function to write report header
    """
    with open("batch_compute.out", 'w') as fh:
        fh.write("  Table  | Status |   Error   \n")
        fh.write("------------------------------\n")


def write_report_row(der, stat, d, C=None, error=None):
    """
    Helper function to write report row
    """
    stat = "  OK  " if stat else "FAILED"
    with open("batch_compute.out", 'a') as fh:
        if C is not None:
            fh.write("  A{}-{}  | {} |  {:.2e}\n".format(C, d, stat, error))
        else:
            ind = "{}  ".format(d) if der == 0 else "{}n{}".format(der,d)
            fh.write("  Q{}   | {} |\n".format(ind, stat))

# Build FC-Tables with desired options and backup current parameters.
build_fctables(digits=_DIGITS, blend=_BLEND)
backup_parameters()

# Bind restore_parameters to SIGINT and SIGTERM
# in case the job gets (politely) killed
signal.signal(signal.SIGINT, restore_parameters)
signal.signal(signal.SIGTERM, restore_parameters)

# If output info file doesn't exists, create it
if not os.path.isfile("batch_compute.out"):
    create_report()
try:
    for d, C in params:
        # Compute A and Q
        set_parameters(d=d, c=C, tkind=0)
        ret = compute_table()
        if _BLEND == "yes":
            write_report_row(0, ret[0], d, C, ret[1])
        write_report_row(0, ret[0], d)

        # Compute Q1n
        if Neumann:
            set_parameters(tkind=1)
            ret = compute_table()
            write_report_row(1, ret[0], d)

        # Compute Q2n
        if Neumann2:
            set_parameters(tkind=2)
            ret = compute_table()
            write_report_row(2, ret[0], d)

finally:
    restore_parameters()
