#!/bin/bash

# Cleaning the run folder
if [[ -d "./__pycache__" ]]; then
    rm -rf __pycache__
fi
