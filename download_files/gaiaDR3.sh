#!/bin/bash


# Session Name
session="gaiaDR3"

# Start New Session with our name
tmux new-session -d -s $session

tmux send-keys "exec ./gaia_DR3_base.sh" 'C-m'
