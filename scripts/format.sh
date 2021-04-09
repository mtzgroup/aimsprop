#!bin/sh -e

set -x

black scripts aimsprop tests
isort scripts aimsprop tests
flake8 tests aimsprop/fms90.py scripts