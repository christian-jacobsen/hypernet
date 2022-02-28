# HyperNet
**Machine Learning-Based library for modeling multi-component non-equilibrium thermochemical processes**

[![Build Status](https://travis-ci.org/ivanZanardi/hypernet.svg?branch=main)](https://travis-ci.org/github/ivanZanardi/hypernet)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://github.com/ivanZanardi/prode/hypernet/main/LICENSE)

Written in Python, HyperNet has been designed to model multi-component non-equilibrium thermochemical processes with common numerical techniques and/or ML-based surrogate models.

## Workflow

1) save a copy of the level-to-bin mapping file in `hypernet/database/air/grouping` [from USER + before installation]
2) create bin-avaraged rates from state-to-state one at fixed temperatures (use app `avarageRates` [MISSING])
3) fit bin-avaraged rates as a function of temperature (use app `fitRates`)
4) run the 0D master equations for the reduced system and the state-to-state one (use app `box`)

## Features

HyperNet (currently) supports

- [x] symmetric 3-atomic systems [O3 (tested) and N3]
- [x] isothermal 0D chemical reactor
- [ ] adiabatic 0D chemical reactor
- [x] one-temperature model (only translational temperature)
- [ ] multi-temperature model (translational + internal bins temperatures)
- [ ] state-to-state system

Coming soon:
1) app `avarageRates`
2) adiabatic 0D chemical reactor
3) multi-temperature model

## Installation

To install HyperNet simply type the following commands in your Linux/Unix shell:

```
$ git clone https://github.com/ivanZanardi/hypernet.git
$ cd hypernet
$ ./installer -a install
```
For more details, type:

```
$ ./installer -h
```

- Dependencies

  - [Matplotlib](https://matplotlib.org/)
  - [NumPy](http://www.numpy.org/)
  - [pandas](https://pandas.pydata.org/)
  - [SciPy](https://www.scipy.org/)
  - [TensorFlow](https://www.tensorflow.org/)>=2.4.1


## Explore more

- [Examples](https://github.com/ivanZanardi/HyperNet/tree/main/examples)

To run an example, just type the name of the app you want to use in your Linux/Unix shell inside the working folder. For more details, type:
```
$ <app-name> -h
```

## License

[Apache license 2.0](https://github.com/ivanZanardi/hypernet/blob/main/LICENSE)
