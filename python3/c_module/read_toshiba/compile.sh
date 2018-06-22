#!/bin/bash

gcc -c -fPIC read_toshiba.c
gcc -shared -fPIC -o read_toshiba.so read_toshiba.o
