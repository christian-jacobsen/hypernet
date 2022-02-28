# HyperNet
**Machine Learning-Based library for modeling multi-component non-equilibrium thermochemical processes.**

[![Build Status](https://travis-ci.org/ivanZanardi/hypernet.svg?branch=main)](https://travis-ci.org/github/ivanZanardi/hypernet)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://github.com/ivanZanardi/prode/hypernet/main/LICENSE)

HyperNet has been designed to 
- solve forward and inverse ordinary differential equations (ODEs) via physics-informed neural network (PINN);
- approximate nonlinear operators via deep operator network (DeepONet);
- approximate functions from a dataset.

## Features

HyperNet supports

- two types of neural networks: unstacked fully connected neural network, and residual neural network;
- many different losses, metrics, optimizers, learning rate schedules, initializations, regularizations, etc.;
- useful techniques, such as dropout and batch normalization;
- callbacks to monitor the internal states and statistics of the model during training;
- enables the user code to be compact, resembling closely the mathematical formulation.

All the components of PrODE are loosely coupled, and thus PrODE is well-structured and highly configurable. It is easy to customize PrODE to meet new demands.

## Installation

To install PrODE simply type the following commands in your Unix shell:

```
$ git clone https://github.com/ivanZanardi/prode.git
$ cd prode
$ ./install -a install
```

- Dependencies

  - [Matplotlib](https://matplotlib.org/)
  - [NumPy](http://www.numpy.org/)
  - [pandas](https://pandas.pydata.org/)
  - [SciPy](https://www.scipy.org/)
  - [TensorFlow](https://www.tensorflow.org/)>=2.4.1


## Explore more

- [Examples](https://github.com/ivanZanardi/HyperNet/tree/main/examples)

## License

[Apache license 2.0](https://github.com/ivanZanardi/hypernet/blob/main/LICENSE)
