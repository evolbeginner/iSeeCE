#! /bin/bash


#################################################################
CWD=`dirname $0`
cd $CWD >/dev/null
CWD=`pwd`
cd - >/dev/null

run_orthomcl=$CWD/scripts/run_orthomcl.sh
cpu=4


#################################################################
while [ $# -gt 0 ]; do
	case $1 in
		--indir)
			indir=$2
			shift
			;;
		--cpu|--blast_CPU|--CPU)
			cpu=$2
			shift
			;;
	esac
	shift
done


#################################################################
for pep_file in $indir/*pep.fas; do
	b=`basename $pep_file`
	c=${b%%.*}
	a=(${a[@]} "-i $pep_file,$c")
done


#################################################################
echo "bash $run_orthomcl -s -b --blast_CPU $cpu ${a[@]}"


