#!/bin/bash
file=$1 #you need to have ped and map file in the folder
plink --file $file  --r2 inter-chr --ld-window-r2 0 --make-bed --make-founders
