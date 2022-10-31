#!/bin/bash

source ~/.bashrc

cd GaiaDR3/

conda activate py3

python Gaia_DR3_download.py > out_gaia.txt
