#!/bin/tcsh -f

# Usage:  Demo.csh  N [10]  fb [0.1]  freeze [no]  dt [2^-5]

set N = 10
if ($#argv > 0) set N = $1
set f = 0
if ($#argv > 1) set f = $2
alias FREEZE cat
if ($#argv > 2) alias FREEZE freeze
set D = -5
if ($#argv > 3) set D = $4

echo ""
echo Demonstration of xstarplot for a $N-body system, binary fraction $f
echo ""
set echo

mkplummer -n $N -i | mkmass -l 0.5 -u 2 | mksecondary -f $f | scale -s \
    | mkbinary -l 100 -u 100 | FREEZE \
    | (kira -t 100 -d 100 -D $D \
    | xstarplot -b -l 1.5 -s 1.5 -m -u -t -p 1 -P 0.9) >& /dev/null

