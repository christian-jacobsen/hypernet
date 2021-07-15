#!/bin/bash

# Cleaning the run folder
if [[ -d "./__pycache__" ]]; then

    rm -rf output __pycache__

fi
