#!/bin/sh -e

set -x

autoflake --remove-all-unused-imports --recursive --remove-unused-variables --in-place scripts aimsprop tests examples --exclude=__init__.py
black scripts aimsprop tests examples
isort scripts aimsprop tests examples
