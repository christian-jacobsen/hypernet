import io
from setuptools import setup
from setuptools import find_packages

with io.open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
    name="hypernet",
    version="0.0.1",
    description="Machine-Learning-Based library for multi-component non-equilibrium thermochemical processes modeling.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Ivan Zanardi",
    author_email="zanardi3@illinois.edu",
    url="https://github.com/ivanZanardi/hypernet",
    license="Apache-2.0",
    install_requires=install_requires,
    entry_points ={
        'console_scripts': [
            'boxNet = hypernet.apps.box:main'
        ]
    },
    classifiers=[
        "Development Status :: 1 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache-2.0 License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    keywords=[
        "Deep Learning",
        "Machine Learning",
        "Neural Networks",
        "Scientific computing",
        "Differential equations",
        "ODE solver",
        "Non-equilibrium thermodynamics",
        "Non-equilibrium chemistry",
        "Hypersonic flows"
    ],
    packages=find_packages(),
    include_package_data=True,
    python_requires=">=3.6",
)
