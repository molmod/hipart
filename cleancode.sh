#!/bin/bash
echo Cleaning python code in $PWD and subdirectories
for file in `find lib scripts | egrep "(\.py$)|(\.f90$)|(^scripts/tr-)"`; do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
done

