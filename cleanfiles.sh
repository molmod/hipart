#!/bin/bash
echo Removing redundant files
for i in `find hipart scripts  | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak|\.so$"` ; do rm -v ${i}; done
rm -vr python-build-stamp-* 
rm -vr HiPart.egg-info

rm -v MANIFEST
rm -vr dist
rm -vr build

