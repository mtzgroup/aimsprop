#!bin/sh -e

set -x

black scripts aimsprop tests examples
isort scripts aimsprop tests examples
flake8 tests aimsprop/fms90.py scripts