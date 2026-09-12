#!/bin/sh

dir=$(pwd)

# Remove autoconf cashe.
rm -rf autom4te.cache
rm -rf aclocal.m4

if [ -f Makefile ]; then
    make distclean || echo "make distclean failed; continuing with a clean bootstrap"
fi

mkdir -p config

# Every remaining step must succeed.
set -e

./bootstrap
./configure --prefix="$dir/projects" "$@"

make
