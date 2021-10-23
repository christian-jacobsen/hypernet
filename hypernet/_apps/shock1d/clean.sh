#!/bin/bash

# Cleaning the run folder
if [[ -d "./postprocessing" ]]; then

    rm -rf postprocessing

    read -p "Do you want to delete the 'chemSurrogateModel' folder? [y/n] " asw
    if [[ $asw = y* ]]; then
		rm -rf chemSurrogateModel
	fi

fi
