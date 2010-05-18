#!/bin/bash
# This is a very simplistic uninstall scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  rm -v /usr/local/bin/hi-*
  rm -vr /usr/local/lib/python*/site-packages/hipart
else
  rm -v $HOME/bin/hi-*
  rm -vr $HOME/lib/python/hipart
fi
