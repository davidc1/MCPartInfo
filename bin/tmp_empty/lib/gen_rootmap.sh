#!/bin/bash

temp=make_rootmap.$$.temp

rootlibmap() {
 ROOTMAP=$1
 SOFILE=$2
 LINKDEF=$3
 shift 3
 DEPS=$*
 if [[ -e $SOFILE && -e $LINKDEF ]]; then
     rlibmap -f -o $ROOTMAP -l $SOFILE -d $DEPS -c $LINKDEF 2>> $temp
     rm -f $temp
 fi
}

######################################################
# Package_Name
rootlibmap libMCPartInfo_Package_Name.rootmap libMCPartInfo_Package_Name.so $LARLITE_USERDEVDIR/MCPartInfo/Package_Name/LinkDef.h 
