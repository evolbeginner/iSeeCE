#! /bin/bash


# check the requirements for iSeeGC


########################################################
RED='\e[0;31m'
BOLD='\e[0;1m'
NC='\e[0m' # No Color


########################################################
dir=`dirname $0`
lib_dir=$dir/lib

run_orthomcl=$dir/scripts/run_orthomcl.sh


########################################################
check_which=$dir/check_which.sh
check_ruby_lib=$dir/check_ruby_lib.rb
check_perl_lib=$dir/check_perl_lib.sh

source $check_which


########################################################
echo -e "${BOLD}Checking bc ......${NC}"
if ! which bc 2>/dev/null 1>&2; then
	echo "bc needs to be installed. Please find the information at https://www.thelinuxfaq.com/159-bc-command-not-found-in-centos-rhel-fedora-ubuntu."
fi
echo ''


echo -e "${BOLD}Checking Ruby and Perl ......${NC}"
software=(ruby perl)
check_which ${software[@]}
v=`perl -v | grep -o 'v[0-9]\(\.[0-9]\+\)\+' | grep -o '[0-9].\+'`
a=`echo $v | awk -F '.' '{if($1>=5 && $2>=10){;}else{print "The version of Perl has to be no lower than 5.10 (currently "$0")";}}'`
if [ ! -z "$a" ]; then
	echo $a
	exit 1
fi


echo ''


echo -e "${BOLD}Checking Ruby libs ......${NC}"
ruby $check_ruby_lib --lib 'getoptlong,find,bio,parallel,Dir' --lib_dir $lib_dir
echo ''


echo -e "${BOLD}Checking Perl libs ......${NC}"
bash $check_perl_lib --lib 'Getopt::Long,Carp,List::Util,Bio::SeqIO,Bio::Perl,Bio::Tools::CodonTable,autovivification' --lib_dir $lib_dir
echo ''


echo -e "${BOLD}Checking alignment and phylogenetic reconstruction software ......${NC}"
phylo_software=(mafft FastTree raxmlHPC-PTHREADS)
check_which ${phylo_software[@]}
echo ''


echo -e "${BOLD}Checking requirements for OrthoMCL ......${NC}"
orthomcl_progs=(blastp makeblastdb diamond mcl orthomclAdjustFasta \
	orthomclBlastParser orthomclLoadBlast orthomclPairs orthomclDumpPairsFiles orthomclMclToGroups)
check_which ${orthomcl_progs[@]}

for line in `grep '^\(orthomcl_config_file\|install_schema\)' $run_orthomcl`; do
#orthomcl_config_file=/home/sswang/software/sequence_analysis/orthomclSoftware-v2.0.4/test1/orthomcl.config
#install_schema_log_file=/home/sswang/software/sequence_analysis/orthomclSoftware-v2.0.4/test1/install_schema.log
	file_name=`echo $line | cut -f 1 -d =`
	path=`echo $line | cut -f 2 -d =`
	if [ ! -f $path ]; then
		echo "The path to $file_name ($path) cannot be found!"
	fi
done
echo ''


