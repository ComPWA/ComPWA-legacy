#! /bin/bash
# make the script verbose and exit as soon as one command fails
set -ev

# just put all test application commands here
./bin/gendalitzApp
./bin/dalitzfitApp

exit 0
