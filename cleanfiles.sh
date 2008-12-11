#!/bin/bash
echo Removing redundant files
rm -vr debian/python-*
rm -vr debian/pycompat
rm -vr debian/compat
rm -vr debian/files
rm -vr debian/stamp-makefile-build
rm -vr python-build-stamp-* 

rm -vr test/tmp
rm -vr test/output

rm -v MANIFEST
rm -vr dist
rm -vr build

rm -vr ext/build

