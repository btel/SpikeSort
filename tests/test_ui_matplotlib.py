from nose.tools import ok_, raises
from nose import with_setup
import os
import shutil

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.testing.decorators import compare_images
from functools import wraps
import matplotlib.pyplot as plt
import glob

test_dir = os.path.join(os.path.dirname(__file__),'mpl_testimgs')
baseline_dir = os.path.join(os.path.dirname(__file__), 'mpl_baseline')
log_dir = os.path.join(os.path.dirname(__file__), 'mpl_errors')

def setup():
    os.mkdir(test_dir)
def teardown():
    shutil.rmtree(test_dir)


def mpl_compare(baseline=None, tol=1e-6):
    def wrapper(func):
        if not baseline:
            _baseline = func.func_name
        else:
            _baseline = baseline

        @wraps(func)
        def _compare():
            generated = os.path.join(test_dir, _baseline +'.png')
            original = os.path.join(baseline_dir, _baseline +'.png')
            func()
            plt.savefig(generated)
            try:
                err = compare_images(original, generated, tol=tol)
            except IOError:
                err = ("Baseline image %s does not exist."
                       "If you run this test for the first time copy the image from %s"  % (_baseline, log_dir))
            if err is None:
                return
            else:
                file_list = glob.glob(os.path.join(test_dir, "*"+_baseline+"*"))
                if not os.path.exists(log_dir):
                    os.mkdir(log_dir)
                for f in file_list:
                    shutil.copy(f, log_dir)
                raise AssertionError, err
        return _compare
    return wrapper

import spike_sort.ui.plotting as splt

@with_setup(setup, teardown)
@mpl_compare()
def test_mpl_compare():
    fig = plt.figure()
    plt.plot([1,2])

@with_setup(setup, teardown)
@mpl_compare()
def test_legend():
    splt.legend([1,2,3])
