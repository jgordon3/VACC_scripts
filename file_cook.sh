#!/bin/bash

$1=file1
$2=file2
$3=name

awk 'FNR==NR {data[FNR]=$0; next} {if(FNR%4==1) $0=data[FNR]}1' $file1 $file2 > $name

