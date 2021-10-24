from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# TensorFlow 
###############################################################################
import os
import hypernet
from hypernet.config import tf_setup
os.environ['TF_CPP_MIN_LOG_LEVEL'] = tf_setup['TF_CPP_MIN_LOG_LEVEL']
import tensorflow as tf
tf.get_logger().setLevel('ERROR')
tf.keras.backend.set_floatx(tf_setup['DTYPE'])
tf.keras.backend.set_epsilon(tf_setup['EPSILON'])

# HyperNet
###############################################################################
__author__  = "Ivan Zanardi"
__email__   = "zanardi3@illinois.edu"
__url__     = "https://github.com/ivanZanardi/hypernet"
__license__ = "Apache-2.0"
__version__ = "0.0.1"

# Import all modules
from hypernet import apps
from hypernet import config
from hypernet import database
from hypernet import src

__all__ = [
    "apps",
    "config",
    "database",
    "src"
]
