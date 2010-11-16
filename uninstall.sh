#!/bin/bash
# This is a very simplistic uninstall scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  rm -vr /usr/local/lib/python*/site-packages/hipart
  rm -vr /usr/local/lib/python*/site-packages/Hipart*
  rm -vr /usr/local/bin/hi-*.py
else
  rm -vr $HOME/lib/python/hipart-*
  rm -vr $HOME/lib/python/HiPart*
  rm -vr $HOME/bin/hi-*.py
fi
