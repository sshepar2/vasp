#!/bin/bash

########################################################################################################################
# For converting the coordinates in a given AIMD step in the XDATCAR file into a CONTCAR and .xyz file.
# For script to work must not concatenate XDATCAR files, or specifically the lines 'Direct configuration: ###'
# in the XDATCAR file must have unique number labels.
# With the assumption that the there will be multiple XDATCAR files for the trajectory an input option is
# given so that the output file has the proper global time step in its name. i.e. XDATCAR_1 XDATCAR_2 each
# have 500 time steps, for step 500 of XDATAR_2 we want to specifiy the 500th step, but the output CONTCAR to
# be labeled CONTAR_1000, since it is the 1000th step globally.
#
# The .xyz file is used in the program Nanodcal, which can take in columns specifying the type of LCAO bases.
# If they differ from atom type to atom type.
# An option is available in this script to add that column. A declaration is used to assign either SZP or
# DZP to a given atom type. (see code below). The option can be toggled on or off below as well.
#
# In 's' (single) mode this wil build a single CONTCAR style file by specifying the directory containing the XDATCAR
# you are interested in followed by the time step in that XDATCAR file.
# In 'm' (mult) mode this will do the same as 's' but specify $start, $step, $finish. One can leave out $finish and
# will be done until out of time steps in the XDATCAR file. First time step is assumed to be 1 (not zero)
#
# In 's' mode
# >>> build.sh $file/path/to/XDATCAR $mode $which-time-step $globalstarttimestep (general example)
# >>> build.sh ../1/XDATCAR s 143 1500 (real example) (output file will be denoted "CONTCAR_1643")
#
# In 'm' mode
# >>> build.sh $file/path/to/XDATCAR $mode $start $step $finish $globalstarttimestep
# >>> build.sh $file/path/to/XDATCAR m 1 50 450 1500 (output files will be denoted "CONTCAR_1501", CONTCAR_1551,...)
#
########################################################################################################################


# getting preliminary info

# build header file

head -7 $1 > header

# getting atom type from file 'header'
tail -2 header | head -1 > temp
cat temp | tr ' ' '\n' > temp2 && mv temp2 temp
sed -i '/^$/d' temp
mapfile -t atomtypes < temp

# getting atom numbers from file 'header'
tail -1 header > temp                            # getting line of interest
cat temp | tr ' ' '\n' > temp2 && mv temp2 temp  # turning into a single column
sed -i '/^$/d' temp                              # removing empty lines
mapfile -t atomnums < temp                       # saving lines as an array

numel=$(echo "${#atomnums[@]}")
atomtot=0
maxindex=$(echo $numel - 1 | bc)

for i in $(seq 0 $maxindex) ; do
    atomtot=$(echo $atomtot+${atomnums[$i]} | bc)
done

# if statements for mode type (IGNORE MODE TYPE NOT USING)
if [ $2 == "s" ] ; then
    length=$(echo ${#3})
    global_t=$(echo $3 + $4 | bc)
    glength=$(echo ${#global_t})
    # if statements to get proper grep text for 's' mode
    if [ $length -eq 1 ] ; then
        grep -A $atomtot "Direct configuration=     $3" $1 > coords
        # if statements to get correct numbering for output file
        if [ $glength -eq 1 ] ; then
            prefix="0000$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 2 ] ; then
            prefix="000$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 3 ] ; then
            prefix="00$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 4 ] ; then
            prefix="0$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 5 ] ; then
            output="CONTCAR_$global_t"
            cat header coords > $output
        else
            echo "Global number labeling only setup for numbers less than 100000"
        fi
    elif [ $length -eq 2 ] ; then
        grep -A $atomtot "Direct configuration=    $3" $1 > coords
        # if statements to get correct numbering for output file
        if [ $glength -eq 1 ] ; then
            prefix="0000$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 2 ] ; then
            prefix="000$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 3 ] ; then
            prefix="00$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 4 ] ; then
            prefix="0$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 5 ] ; then
            output="CONTCAR_$global_t"
            cat header coords > $output
        else
            echo "Global number labeling only setup for numbers less than 100000"
        fi
    elif [ $length -eq 3 ] ; then
        grep -A $atomtot "Direct configuration=   $3" $1 > coords
        # if statements to get correct numbering for output file
        if [ $glength -eq 1 ] ; then
            prefix="0000$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 2 ] ; then
            prefix="000$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 3 ] ; then
            prefix="00$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 4 ] ; then
            prefix="0$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 5 ] ; then
            output="CONTCAR_$global_t"
            cat header coords > $output
        else
            echo "Global number labeling only setup for numbers less than 100000"
        fi
    elif [ length $3 -eq 4 ] ; then
        grep -A $atomtot "Direct configuration=  $3" $1 > coords
        # if statements to get correct numbering for output file
        if [ $glength -eq 1 ] ; then
            prefix="0000$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 2 ] ; then
            prefix="000$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 3 ] ; then
            prefix="00$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 4 ] ; then
            prefix="0$global_t"
            output="CONTCAR_$prefix"
            cat header coords > $output
        elif [ $glength -eq 5 ] ; then
            output="CONTCAR_$global_t"
            cat header coords > $output
        else
            echo "Global number labeling only setup for numbers less than 100000"
        fi
    else
        echo "Step choice must be smaller than 10000"
    fi
#elif [ '$2' = 'm']
else
    echo "Choose either 's' or 'm' mode option"
fi

rm coords temp

# Here will convert the output to cartesian .xyz in nanodcal format. you can add options to include other columns
xyz="$output.xyz"
echo "$atomtot" > $xyz

# making long array of atomtypes
touch temp
for i in $(seq 0 $maxindex) ; do
    for j in $(seq 1 ${atomnums[$i]}) ; do
        printf "%s\n" "${atomtypes[$i]}" >> temp
    done
done

mapfile -t atomtypelist < temp
rm temp

# getting unit cell info to convert to cartesian
ax=$(head -3 header | tail -1 | awk '{print $1}')
ay=$(head -3 header | tail -1 | awk '{print $2}')
az=$(head -3 header | tail -1 | awk '{print $3}')

bx=$(head -4 header | tail -1 | awk '{print $1}')
by=$(head -4 header | tail -1 | awk '{print $2}')
bz=$(head -4 header | tail -1 | awk '{print $3}')

cx=$(head -5 header | tail -1 | awk '{print $1}')
cy=$(head -5 header | tail -1 | awk '{print $2}')
cz=$(head -5 header | tail -1 | awk '{print $3}')

orbs=1 # 1 will include orbital info. Add details yourseld as needed

if [ $orbs -eq 0 ] ; then
    printf "%8s%3s%15s%15s\n" "AtomType" "X" "Y" "Z" >> $xyz # or below to specify orbital type
    for i in $(seq 1 $atomtot) ; do
        line=$(echo 8 + $i | bc)
        da=$(head -$line $output | tail -1 | awk '{print $1}')
        db=$(head -$line $output | tail -1 | awk '{print $2}')
        dc=$(head -$line $output | tail -1 | awk '{print $3}')
        x=$(echo "scale= 10; $da*$ax + $db*$bx + $dc*$cx" | bc)
        y=$(echo "scale= 10; $da*$ay + $db*$by + $dc*$cy" | bc)
        z=$(echo "scale= 10; $da*$az + $db*$bz + $dc*$cz" | bc)
        atid=$(echo $i - 1 | bc)
        atom=${atomtypelist[$atid]}
        printf "%2s%20.10f%15.10f%15.10f\n" "$atom" "$x" "$y" "$z" >> $xyz
    done
elif [ $orbs -eq 1 ] ; then
    printf "%8s%3s%15s%15s%24s\n" "AtomType" "X" "Y" "Z" "OrbitalType" >> $xyz
    #specify which atom types will be DZP, the rest will be SZP
    declare -A dzp=([N]=1 [C]=1 [H]=1 [O]=1)
    # building long string of OrbitalTypes
    for i in $(seq 1 $atomtot) ; do
        atid=$(echo $i - 1 | bc)
        atom=${atomtypelist[$atid]}
        if [ -n "${dzp[$atom]}" ] ; then
            echo "DZP" >> temp
        else
            echo "SZP" >> temp
        fi
    done
    mapfile -t orbitaltypelist < temp
    rm temp
    for i in $(seq 1 $atomtot) ; do
        line=$(echo 8 + $i | bc)
        da=$(head -$line $output | tail -1 | awk '{print $1}')
        db=$(head -$line $output | tail -1 | awk '{print $2}')
        dc=$(head -$line $output | tail -1 | awk '{print $3}')
        x=$(echo "scale= 10; $da*$ax + $db*$bx + $dc*$cx" | bc)
        y=$(echo "scale= 10; $da*$ay + $db*$by + $dc*$cy" | bc)
        z=$(echo "scale= 10; $da*$az + $db*$bz + $dc*$cz" | bc)
        atid=$(echo $i - 1 | bc)
        atom=${atomtypelist[$atid]}
        orbital=${orbitaltypelist[$atid]}
        printf "%2s%20.10f%15.10f%15.10f%5s\n" "$atom" "$x" "$y" "$z" "$orbital">> $xyz
    done
fi

rm header
