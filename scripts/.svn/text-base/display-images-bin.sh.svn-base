#!/bin/sh
## display-images.sh for  in /Users/pouchet
##
## Made by Louis-Noel Pouchet
## Contact: <pouchet@cs.ucla.edu>
##
## Started on  Sat Apr 26 19:12:02 2014 Louis-Noel Pouchet
## Last update Sun Jan 25 21:53:07 2015 Louis-Noel Pouchet
##


echo "** Displaying $1 as pdf file of $2 slices of size $3 x $4 **";
myx="x";
str="$3$myx$4";
c=`convert -depth 32 -size $str -define quantum:format=integer gray:$1 $1.pdf`;
open $1.pdf;
