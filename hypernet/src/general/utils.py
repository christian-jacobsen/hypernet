import os
import sys
import time
import signal
import inspect
import datetime
import numpy as np

from varname import nameof
from functools import wraps
from hypernet.config import headers


# Printing functions
###############################################################################
def print_main(text, start='\n', verbose=True, *args, **kwargs):
    if verbose:
        print(start + headers['main'] + text, *args, **kwargs)

def print_submain(text, start='', verbose=True, *args, **kwargs):
    if verbose:
        print(start + headers['main'] + '  ' + text, *args, **kwargs)

def input_main(text, start='', *args, **kwargs):
    return input(start + headers['main'] + text, *args, **kwargs)

def input_submain(text, start='', *args, **kwargs):
    return input(start + headers['main'] + '  ' + text, *args, **kwargs)

def warning(text, start='', *args, **kwargs):
    print(start + headers['warning'] + text, *args, **kwargs)

def raise_value_err(text):
    message = ''.join([headers['val_err'], text])
    raise ValueError(message)

# Decorators
###############################################################################
def main_decorator(f):
    """Decorator for the `main` function."""
    @wraps(f)
    def wrapper(*args, **kwargs):
        print_main('START >>>')
        result = f(*args, **kwargs)
        print_main('<<< END')
        return result
    return wrapper

def timing(f):
    """Decorator for measuring the execution time of methods."""
    @wraps(f)
    def wrapper(*args, **kwargs):
        ts = time.time()
        result = f(*args, **kwargs)
        te = time.time()
        delta = datetime.timedelta(seconds=(te - ts))
        print_submain(">> Function %r took %s" % (f.__name__, delta))
        sys.stdout.flush()
        return result
    return wrapper

# Decorator for applications --------------------------------------------------
LEN = 98

def compose(text='', start='', between=' '):
    end = ''.join(reversed(start))
    l = len(text)
    return start + between*int((LEN-len(start)*2-l)/2) \
        + text + between*int((LEN-len(start)*2-l)/2) + end

def app_decorator(name=None):
    def dec(fun):
        @wraps(fun)
        def wrapper(name=name, *args, **kwargs):
            ts = time.time()
            start = '# '
            print(
                "\n" + compose(start='# ', between='*') +
                "\n" + compose(text='PrODE ', start='# ') +
                "\n" + compose(start='# ', between='-') +
                "\n" + compose(
                    text='Machine-Learning-Based Approximators', start='# '
                ) +
                "\n" + compose(
                    text='for Ordinary Differential Equations (ODEs)',
                    start='# '
                ) +
                "\n" + compose(start='# ') +
                "\n" + compose(start='# ', between='*')
            )
            print_main('App `{}` starts\n'.format(name) + "="*LEN + "\n")
            result = fun(*args, **kwargs)
            te = time.time()
            delta = datetime.timedelta(seconds=(te - ts))
            if not name:
                name = fun.__name__
            print(app_epilog(name=name, delta=delta))
            # print("\n" + "=" * LEN)
            # print_main('App `{}` terminates'.format(name), start='')
            # print_submain(">> It took {}\n".format(delta))
            return result
        return wrapper
    return dec

def app_epilog(name=None, delta=None):
    s = '\n'+'='*LEN+'\n'+headers['main']+'App `{}` terminates\n'.format(name)
    if delta:
        return s + headers['main'] + '  ' + '>> It took {}\n'.format(delta)
    else:
        return s + ' '

# Handling data types
###############################################################################
def get_class(module, name):
    """Return a class object given its name and the module which belongs to."""
    for name_i, obj_i in inspect.getmembers(module, inspect.isclass):
        if name_i == name:
            return obj_i
    raise_value_err(
        "Class `{}` not found in module `{}`.".format(name, module.__name__)
    )

def get_class_name(cls):
    """Return the class name."""
    return cls.__class__.__name__

def make_dict(keys, values):
    """Convert two lists or two variables into a dictionary."""
    if isinstance(keys, (list, tuple)):
        if len(keys) != len(values):
            raise_value_err("Keys and values have different length.")
        return dict(zip(keys, values))
    return {keys: values}

def list_to_str(nums, precision=2):
    """Convert a list/tuple/array of numbers into a string."""
    if nums is None:
        return ""
    if not isinstance(nums, (list, tuple, np.ndarray)):
        return "{:.{}e}".format(nums, precision)
    return "[{:s}]".format(
        ", ".join(["{:.{}e}".format(x, precision) for x in nums])
    )

def get_num_args(func):
    """Get the number of arguments of a Python function."""
    if sys.version_info[0] == 2:
        return len(inspect.getargspec(func).args)
    sig = inspect.signature(func)
    return len(sig.parameters)

def convert_to_array(arg):
    """Try to convert the argument into an array."""
    if not isinstance(arg, np.ndarray):
        try:
            arg = np.array(arg)
            if len(arg.shape) == 0:
                arg = np.expand_dims(arg, 0)
        except:
            raise_value_err(
                "`{}` can't be converted into an array.".format(nameof(x))
            )
    return arg

def search_string(file, string):
    with open(file, 'r') as f:
        for i, line in enumerate(f):
            if string+' ' in line:
                return i, line

def check_format(data_type, string):
    if data_type == float and 'd' in string:
        return string.replace("d", "e")
    else:
        return string

# Others
###############################################################################
def kill_process():
    """Send terminating signal to job `pid`."""
    os.kill(os.getpid(), signal.SIGTERM)
