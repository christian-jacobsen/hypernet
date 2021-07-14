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
def print_main(text, *args, **kwargs):
    print(headers['main'] + text, *args, **kwargs)

def print_submain(text, *args, **kwargs):
    print(headers['submain'] + text, *args, **kwargs)

def input_main(text, *args, **kwargs):
    return input(headers['main'] + text, *args, **kwargs)

def input_submain(text, *args, **kwargs):
    return input(headers['submain'] + text, *args, **kwargs)

def warning(text, *args, **kwargs):
    print(headers['warning'] + text, *args, **kwargs)

def raise_value_err(text):
    message = ''.join([headers['val_err'], text])
    raise ValueError(message)


# Decorators
###############################################################################
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


# Handling data types
###############################################################################
def get_class(module, name):
    """Return a class object given its name and the module which belongs to."""
    for name_i, obj_i in inspect.getmembers(module, inspect.isclass):
        if name_i == name:
            return obj_i
    raise_value_err("Class `{}` not found in module `{}`.".format(name, module.__name__))

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
    return "[{:s}]".format(", ".join(["{:.{}e}".format(x, precision) for x in nums]))

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


# Others
###############################################################################
def kill_process():
    """Send terminating signal to job `pid`."""
    os.kill(os.getpid(), signal.SIGTERM)
