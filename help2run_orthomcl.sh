#! /bin/bash


#################################################################
CWD=`dirname $0`
cd $CWD >/dev/null
CWD=`pwd`
cd - >/dev/null

run_orthomcl=$CWD/scripts/run_orthomcl.sh
cpu=4
is_diamond=false
blast_prog="blastp"


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
		-p)
			passwd=$2
			shift
			;;
		--orthomcl_config)
			orthomcl_config_file=$2
			shift
			;;
		--blast_prog)
			blast_prog=$2
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
echo "bash $run_orthomcl --orthomcl_config $orthomcl_config_file --blast_prog $blast_prog -s -b --blast_CPU $cpu ${a[@]}"


