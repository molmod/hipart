#!/bin/bash
echo Cleaning python code in $PWD and subdirectories
for file in `find * | egrep "(\.py$)|(\.f90$)|(^scripts/tr-)"`; do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
done
for i in `find * | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$"` ; do rm -v ${i}; done


