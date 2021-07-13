import numpy as np

def convert_to_array(x):        
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    return x
