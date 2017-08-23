#! /bin/bash


# check the requirements for iSeeGC


########################################################
RED='\e[0;31m'
BOLD='\e[0;1m'
NC='\e[0m' # No Color


########################################################
dir=`dirname $0`
lib_dir=$dir/lib


########################################################
check_ruby_lib=$dir/check_ruby_lib.rb
check_perl_lib=$dir/check_perl_lib.sh


########################################################
echo -e "${BOLD}Checking Ruby libs ......${NC}"
ruby $check_ruby_lib --lib 'getoptlong,find,bio,Dir' --lib_dir $lib_dir
echo ''

echo -e "${BOLD}Checking Perl libs ......${NC}"
bash $check_perl_lib --lib 'Getopt::Long,Carp,List::Util,Bio::SeqIO,Bio::Perl,Bio::Tools::CodonTable' --lib_dir $lib_dir
echo ''


