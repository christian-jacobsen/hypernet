#!/bin/bash

find . -type d -name __pycache__ -exec rm -r {} \;
find . -type f -name .DS_Store   -exec rm -r {} \;
find . -type f -name *.pyc       -exec rm -r {} \;
find . -type f -name *.bak       -exec rm -r {} \;
