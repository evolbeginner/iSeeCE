#! /bin/bash


# check whether Perl libs are installed.


###################################################
declare -a libs


###################################################
while [ $# -gt 0 ]; do
	case $1 in
		--lib)
			for i in `echo $2 | awk -F',' '{for(i=1;i<=NF;i++){print $i}}'`; do
				libs=(${libs[@]} $i)
			done
			shift
			;;
	esac
	shift
done


###################################################
for lib in ${libs[@]}; do
	perl -M$lib -e 1 >/dev/null 2>&1
	if [ $? == 0 ]; then
		echo -e "${lib} is detected."
	else
		echo -e "\e[0;1;31m${lib}\e[0m has NOT been installed."
	fi
done


