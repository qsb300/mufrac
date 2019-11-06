#!/bin/bash

#input virial coeffcients based on sigma instead of v as in doi:10.1135/cccc2009113, 
#THE EFFECT OF TRUNCATION AND SHIFT ON VIRIAL COEFFICIENTSOF LENNARD–JONES POTENTIALS

b2=`echo $1 |awk '{print $1/(3.1415926/6.0)}'`
b3=`echo $2 |awk '{print $1/((3.1415926/6.0)**2)}'`
b4=`echo $3 |awk '{print $1/((3.1415926/6.0)**3)}'`
b5=`echo $4 |awk '{print $1/((3.1415926/6.0)**4)}'`
set -x
python ../src/AllowedFraction.py  0 50 0.01 $b2 $b3 $b4 $b5
