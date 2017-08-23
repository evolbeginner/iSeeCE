#! /bin/bash


# run_orthomcl:	run orthomcl automatically
# last updated on 2014-04-08
# To see its usage, please type 'run_orthomcl.sh --h' to see its usage.
# Author: Sishuo Wang (sishuowang@hotmail.ca, wangsishuo@yeah.net, tomassonwss@gmail.com) from the department of Botany, the University of British Columbia


##########################################################################
# fill in the path of the following two files here
orthomcl_config_file=/home/sswang/software/sequence_analysis/orthomclSoftware-v2.0.4/test1/orthomcl.config
install_schema_log_file=/home/sswang/software/sequence_analysis/orthomclSoftware-v2.0.4/test1/install_schema.log
blast_prog="blastp"
blastp=~/bin/blastp


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
	-s|--s|--silent)
		silent='true'
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
		echo -e "all_VS_all.out.tab cannot not found!\nExiting ......"
		exit -1
	fi
fi
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


blast(){
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
if [ "$blast" == 'true' ]; then
	type="prot"
	echo "makeblastdb"
	cd compliantFasta
	local input_basename=`basename all.fasta`
	[ $blast_programme == "blastn" ] && type="nucl"
	makeblastdb -in $input_basename -dbtype $type -out $input_basename
	#formatdb -i $input_basename -p $type;
	echo -ne "BLASTing\t";
	blast_command="$blastp -query all.fasta -db $input_basename -out all_VS_all.out.tab -outfmt 6 -evalue 1e-3 -num_threads $blast_CPU"
	#blast_command="blastall -p $blast_programme -i all.fasta -d $input_basename -e $blast_evalue -o all_VS_all.out.tab -a $blast_CPU -m8";
	#blast_command="blastall -p $blast_programme -i all.fasta -d $input_basename -e $blast_evalue -o all_VS_all.out.tab -a $blast_CPU -m8";
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
	echo "`basename $0` <--i|--input infile,abbr>:
		abbr refers to the abbreviation of the name representing the infile.
		abbr should be seperated from the name of input file by a comma.
		If abbr is not specified, its value will be the same as the input file."
	echo "Options:
		[--b|--blast]:           to do blast or not
			                 default off
		[--blast_CPU CPU]:       the number of threads for BLAST search
		                         default 2
		[--blast_evalue evalue]: default 1e-10
		[--s|--silent]
		[--h|--help]"
        exit 1
}


###################################################################################
###################################################################################
read_para $*;

mysql -uroot -e 'drop database orthomcl; create database orthomcl;';

###	-------------------------------------------------------------------	###
orthomclInstallSchema $orthomcl_config_file $install_schema_log_file;

AdjustFasta ${input[@]}

blast $blast $blast_prog;

prepare_compliant

orthomcl_2

clear_files

mysql -uroot -e 'drop database orthomcl; create database orthomcl;';

