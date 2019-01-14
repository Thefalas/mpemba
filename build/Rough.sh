#!/bin/sh
module load intel

# usage: "./Rough.sh 44" where 44 is the simulation number
prefix="000"
if (($1 > 9 ))
then
    prefix="00"
fi
if [ $1 -gt 99 ]
then
    prefix="0"
fi
if [ $1 -gt 999 ]
then
    prefix=""
fi
nombre="sRWN${prefix}${1}"
cd /home/malopez/mpemba/build
./RWN > $nombre
