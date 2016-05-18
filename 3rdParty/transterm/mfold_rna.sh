#!/bin/sh

if [ ${#} -lt 2 ] ; then
    echo "usage: $0 program_name  [program options]" > /dev/stderr
    exit 3
fi

prog=$1
shift

$prog --gc=-2.1 --au=-1.3 --gu=-0.8 --mm=3.5 --gap=3.9 \
    --loop-penalty=4.1,4.9,4.4,4.7,5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.8,5.9,6 $*
    
