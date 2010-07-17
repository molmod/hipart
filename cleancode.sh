#!/bin/bash
echo Cleaning python code in $PWD and subdirectories
for file in $(find hipart scripts | egrep "(\.py$)|(\.f90$)|(\.c$)|(\.h$)|(\.pyf$)|(^scripts/hi-)"); do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]*$//' ${file}
  sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' ${file}
done

