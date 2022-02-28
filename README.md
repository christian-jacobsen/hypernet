# PrODE
**Machine-Learning-Based Approximators for Ordinary Differential Equations (ODEs)**

[![Build Status](https://travis-ci.org/ivanZanardi/prode.svg?branch=main)](https://travis-ci.org/github/ivanZanardi/prode)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://github.com/ivanZanardi/prode/blob/main/LICENSE)

PrODE is a deep learning library written on top of the deep learning API [Keras](https://keras.io/), which in turn runs on top of the machine learning platform [TensorFlow](https://www.tensorflow.org/).
PrODE has been designed to 
- solve forward and inverse ordinary differential equations (ODEs) via physics-informed neural network (PINN);
- approximate nonlinear operators via deep operator network (DeepONet);
- approximate functions from a dataset.

## Features

PrODE supports

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
  - [Netron](https://netron.app/)
  - [NumPy](http://www.numpy.org/)
  - [pandas](https://pandas.pydata.org/)
  - [pyDOE](https://pythonhosted.org/pyDOE/)
  - [scikit-learn](https://scikit-learn.org)
  - [SciPy](https://www.scipy.org/)
  - [TensorFlow](https://www.tensorflow.org/)>=2.4.1
  - [TensorFlow Addons](https://www.tensorflow.org/addons/)>=0.12.1
  - [tqdm](https://tqdm.github.io/)

Install also external Ubuntu packages:
- [Graphviz](https://graphviz.org/)
- [FFmpeg](https://www.ffmpeg.org/)

```
$ sudo apt-get install graphviz ffmpeg
```

## Use with Docker

If you have downloaded a [TensorFlow Docker image](https://www.tensorflow.org/install/docker), you can build an image of PrODE by following the instructions given in `docker/build_image.sh` and then typing the following command:

```
$ ./docker/build_image.sh
```

Afterwards, you can run your `main.py` script by running the new built image. To do that, follow the instructions given in `docker/run_script.sh`.

## Explore more

- [Examples](https://github.com/ivanZanardi/PrODE/tree/main/examples)

## License

[Apache license 2.0](https://github.com/ivanZanardi/prode/blob/main/LICENSE)