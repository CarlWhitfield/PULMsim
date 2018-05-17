#!/bin/sh

#  repeatrun.sh
#  
#
#  Created by Carl Whitfield on 07/01/2018.
#  
exec_file=PULMsim_linux

for filename in *.lung
do
 dir=$(echo $filename | cut -d'.' -f1)
 mkdir $dir
 cp $exec_file $dir
 mv $filename $dir
 cd $dir
	./$exec_file $filename
 cd ../
done
