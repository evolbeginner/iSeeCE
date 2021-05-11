#! /bin/bash


# run_orthomcl:	run orthomcl automatically
# last updated on 2014-04-08
# To see its usage, please type 'run_orthomcl.sh --h' to see its usage.
# Author: Sishuo Wang (sishuowang@hotmail.com; wangsishuo@yeah.net) from the department of Botany, the University of British Columbia


##########################################################################
# fill in the path of the following two files here
#orthomcl_config_file=/home/sswang/software/sequence_analysis/orthomclSoftware-v2.0.4/test1/orthomcl.config
#install_schema_log_file=/home/sswang/software/sequence_analysis/orthomclSoftware-v2.0.4/test1/install_schema.log
install_schema_log_file=install_schema.log
blast_prog="blastp"
blastp=blastp
#user_argu="-uroot"
user_argu=''


##########################################################################
read_para(){
[ -z "$1" ] && show_help
while [ $# -gt 0 ]
do
	case $1 in
	-b|--b|--blast)
		echo "blast will be performed" 
		blast='true';
		;;
	--blast_evalue)
		echo "blast e_value is $2"
		blast_evalue=$2
		shift
		;;
	--blast_CPU|--cpu|--CPU)
		echo "CPU of blast is $2"
		blast_CPU=$2
		shift
		;;
	--blast_prog|--prog)
		echo "blast_prog is $2"
		blast_prog=$2
		shift
		;;
	-i|--i|--in|--input)
		local tmp=$2
		local input_tmp=${tmp%,*}
		local abbr_tmp=${tmp##*,}
		echo $input_tmp, $abbr_tmp
		input=(${input[@]} $input_tmp)
		abbr=(${abbr[@]} $abbr_tmp)
		shift
		;;
	--orthomcl_config)
		orthomcl_config_file=$2
		shift
		;;
	-s|--s|--silent)
		silent='true'
		;;
	-u)
		user_argu='-u'$2
		shift
		;;
	-p)
		passwd_argu="-p$2"
		shift
		;;
	-h|--h|--help)
		show_help
		;;
	*)
		echo "unkonw argument $1"
		show_help
        	;;
        esac
	shift
done
echo

[ -z $blast_CPU ] && blast_CPU=2
[ -z $blast_evalue ] && blast_evalue=1e-10
if [ -z $blast ]; then
	if ! ls all_VS_all.out.tab 1>/dev/null 2>/dev/null; then
		echo -e "all_VS_all.out.tab cannot not found!\nExiting ......" >&2
		exit -1
	fi
fi
}


function read_orthomcl_config(){
	local config_file=$1
	user_argu="-u"`grep '^dbLogin' $config_file | cut -f 2 -d =`
	passwd_argu="-p"`grep '^dbPassword' $config_file | cut -f 2 -d =`
	mysql_db=`grep ^dbConnectString $config_file | awk -F':' '{print $NF}'`
}


AdjustFasta(){
	echo "AdjustFasta ......"
	for i in ${!input[@]}; do
		echo $i, ${input[$i]}
		orthomclAdjustFasta ${abbr[$i]} ${input[$i]} 1 > /dev/null;
	done
	echo "AdjustFasta finished"
	echo
}


do_blast(){
[ -e compliantFasta ] && rm -rf compliantFasta
mkdir compliantFasta
for i in ${abbr[@]}; do
	echo $i
	ln -s `pwd`/$i.fasta compliantFasta/$i;
done
cat compliantFasta/* > compliantFasta/all.fasta
echo $*
local blast;
blast=$1;
blast_programme=$2;
echo -ne "BLASTing\t";
if [ "$blast" == 'true' ]; then
	type="prot"
	echo "makeblastdb"
	cd compliantFasta
	local input_basename=`basename all.fasta`
	[ $blast_programme == "blastn" ] && type="nucl"
	if [ $blast_programme == "diamond" ]; then
		diamond  makedb --in $input_basename -d $input_basename -p 2
		blast_command="diamond blastp -q all.fasta -o all_VS_all.out.tab -d $input_basename -p $blast_CPU -e 1e-3"
	else
		makeblastdb -in $input_basename -dbtype $type -out $input_basename
		blast_command="$blastp -query all.fasta -db $input_basename -out all_VS_all.out.tab -outfmt 6 -evalue 1e-3 -num_threads $blast_CPU"
	fi
	echo $blast_command;
	$blast_command;
	mv all_VS_all.out.tab ../
	cd -
fi
}


prepare_compliant(){
content=`ls compliantFasta`;
if [ "$silent" == 'true' ]
then
	echo 
	cd compliantFasta; ls | grep -v '.*.fasta$' | xargs rm; cd ..;
else
	echo "Y/N";
	read input; 
	if  [[ $input =~ Y|y ]]
	then
		cd compliantFasta; ls | grep -v '.*.fasta$' | xargs rm; cd ..;
	fi
fi
}


orthomcl_2(){
	orthomclBlastParser all_VS_all.out.tab compliantFasta > similarSequences.txt
	awk 'BEGIN{FS="\t"}{pairs[$1"\t"$2]+=1; all[$1"\t"$2]=$0}END{for(i in all){print all[i]}}' similarSequences.txt | sponge similarSequences.txt
	orthomclLoadBlast $orthomcl_config_file similarSequences.txt
	if [ -e pairs ]; then rm pairs -rf; fi;
	orthomclPairs $orthomcl_config_file orthomcl_pairs.log cleanup=no
	orthomclDumpPairsFiles $orthomcl_config_file
	mcl mclInput --abc -I 1.5 -o mclOutput
	orthomclMclToGroups GF_ 1 < mclOutput > groups.txt
}


clear_files(){
	for i in ${abbr[@]}; do
		rm $i.fasta
	done
}


show_help(){
	echo "`basename $0` <--orthomcl_config orthomcl_config> <--i|--input infile,abbr>:
		abbr refers to the abbreviation of the name representing the infile.
		abbr should be seperated from the name of input file by a comma.
		If abbr is not specified, its value will be the same as the input file." >&2
	echo "Options:
		[--b|--blast]:           to do blast or not
			                 default off
		[--blast_CPU CPU]:       the number of threads for BLAST search
		                         default 2
		[--blast_evalue evalue]: default 1e-10
		[--s|--silent]
		[--h|--help]" >&2
        exit 1
}


###################################################################################
###################################################################################
read_para $*;

read_orthomcl_config $orthomcl_config_file

#mysql $user_argu $passwd_argu -e 'drop database if exists orthomcl; create database orthomcl;';
mysql $user_argu $passwd_argu -e "drop database if exists $mysql_db; create database $mysql_db;";


###	-------------------------------------------------------------------	###
orthomclInstallSchema $orthomcl_config_file $install_schema_log_file;

AdjustFasta ${input[@]}

do_blast $blast $blast_prog;

prepare_compliant

orthomcl_2

clear_files


