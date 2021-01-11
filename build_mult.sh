#!/bin/bash

######################################################################################################################
# this script uses build.sh to just calculate multiple CONTCAR files from XDATCARS
#
# >>> ./build_mult.sh $/file/path/to/XDATCAR $start $step $finish $global_step
#
# or $finish can be left out and the script will just iterate until the end of the XDATCAR file is reached
#
#
# note t=0 is the initial POSCAR not the first set of coords in the XDATCAR.
#
# mode = 'm' is automaticall selected when calling on build.sh.
######################################################################################################################

# checking input arguments
args=$(echo $#)

if [ $args -lt 4 ] ; then
    echo "Did not supply enough input arguments."
    echo "/data/home/sshepar2/bin/build_mult.sh /path/to/XDATCAR start step [finish] global_step"
    exit 1
elif [ $args -lt 5 ] ; then
    echo "It has been assumed that the [finish] input argument has been left out. The XDATCAR will be parsed until the end"
    global=$4
    finish=$(grep -c "Direct" $1)
else
    finish=$4
    global=$5
fi

i=$2
tglobal=$(echo $2 + $global| bc)
while [ $i -le $finish ] ; do
    echo "Building file for t=$i [tg=$tglobal]"
    /data/home/sshepar2/bin/build.sh $1 s $i $global
    i=$(echo $i + $3 | bc)
    tglobal=$(echo $tglobal + $3 | bc)
done

echo "done"
