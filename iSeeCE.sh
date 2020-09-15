#! /bin/bash

##############################################################
# A wrapper of scripts to identify recurrent gene conversion using phylogenomic methods
# Please type "bash iSeeCE.sh" for help message.
# Last updated: 2017/11/21


##############################################################
CWD=`dirname $0`
cd $CWD >/dev/null
CWD=`pwd`
cd - >/dev/null
export RUBYLIB=$RUBYLIB:$CWD

source "$CWD/scripts/processbar.sh"
source "$CWD/word_limit_per_line.sh"


##############################################################
# Path to Mauve
Mauve=''
Mauve_jar=''
ruby=ruby


##############################################################
RED='\e[0;31m'
BOLD='\e[0;1m'
NC='\e[0m' # No Color


##############################################################
# Variables
gc_count_min=5
ortholog_count_min=2
is_force=false
blast_cpu=2
mafft_cpu=2
bootstrap_min=0
is_mauve=false


##############################################################
function parse_genbank_files(){
	local indir=$1
	local genbank2cds=$2
	local count=0
	local total=0

	echo -e "Parsing Genbank files ......"
	total=`ls -1 $indir/*gb|wc -l`

	for i in $indir/*gb; do 
		bash $genbank2cds --gb $i --out_prefix ${i%.gb};
		count=`expr $count + 1`
		processbar $count $total
	done
	echo ''
}


function run_orthomcl(){
	local indir=$1
	local orthomcl_dir=$2
	local blast_cpu=$3

	echo "Running OrthoMCL (which usually takes a long time) ......"
	run_orthomcl_cmd=`bash $help2run_orthomcl --indir $indir --cpu $blast_cpu`
	cd $orthomcl_dir >/dev/null
	$run_orthomcl_cmd 1>/dev/null 2>/dev/null

	processbar 1 4
	count=1
	for i in 6 4 2; do
		mcl mclInput --abc -I $i -o mclOutput$i 2>/dev/null
		count=`expr $count + 1`
		processbar $count 4
	done
	echo ''

	cp mclOutput mclOutput1.5
	cd - >/dev/null
}


function run_mauve(){
	local mauve_dir=$1
	local Mauve=$2
	local Mauve_jar=$3

	echo "Doing genome alignment with Mauve ......"
	cd $mauve_dir >/dev/null
	$Mauve --output=mauve ../sequences/*gb
	java -cp $Mauve_jar org.gel.mauve.analysis.OneToOneOrthologExporter -f mauve -o mauve.ortho
	[ ! -s mauve.ortho ] && echo "The file size of mauve.ortho equals 0! Exiting ......" && exit 1
	cd - >/dev/null
}


function run_mafft(){
	local outdir=$1
	local aln_dir=$2
	local selectOrthoFromMauveOrthomcl=$3
	local check_aln_group=$4
	local mafft_cpu=$5
	local count=0

	cd $outdir >/dev/null
	echo "Running MAFFT (which usually takes a long time) ......"
	for i in 1.5 2 4 6; do
		$ruby $selectOrthoFromMauveOrthomcl --mauve mauve/mauve.ortho --orthomcl OrthoMCL/mclOutput$i --seq_indir sequences --seq_file_suffix cds.fas --outdir aln/aln$i --force --group_size 4,200 --mauve_size 0,200 --cpu $mafft_cpu >/dev/null;
		count=`expr $count + 1`
		processbar $count 4
		continue
	done
	echo '' 
	cd - >/dev/null

	$ruby $check_aln_group --indir $aln_dir/aln1.5/ --indir $aln_dir/aln2/ --indir $aln_dir/aln4/ --indir $aln_dir/aln6 --suffix aln --outdir $aln_dir/combined --force
}


function run_FastTree(){
	local aln_dir=$1
	local count=0
	local total=0

	echo "Running FastTree ......"
	cd $aln_dir/combined/
	local total=`ls -1 | wc -l`

	for i in `ls`; do
		FastTree -quiet -nt < $i > ${i/aln/FastTree.tre} 2>/dev/null
		count=`expr $count + 1`
		processbar $count $total
	done
	echo ''

	mv *tre ../../FastTree/trees/;
	cd - >/dev/null
}


function parse_FastTree(){
	local fasttree_dir=$1
	local mauve_dir=$2
	local ortholog_count_min=$3
	local bootstrap_min=$4
	local is_mauve=$5

	local count=0

	if [ $is_mauve == true ]; then
		echo "Parsing FastTree result with Mauve ......"
		[ -f $fasttree_dir/trees.mauve.FastTree_result ] && rm $fasttree_dir/trees.mauve.FastTree_result
		local total=`ls -1 $fasttree_dir/trees/* | wc -l`
		for i in $fasttree_dir/trees/*; do
			count=`expr $count + 1`
			processbar $count $total
			$ruby $putativeGC_from_seqSimilarity --tree $i --mauve $mauve_dir/mauve.ortho -b 0 >> $fasttree_dir/trees.mauve.FastTree_result 2>/dev/null
		done
		$ruby $convertMauveResultToRawResult -i $fasttree_dir/trees.mauve.FastTree_result --ortho_count_min $ortholog_count_min > $fasttree_dir/GC.mauve.raw_result
		echo
	fi

	echo "Parsing FastTree result ......"
	bash $newFastTreeWrapper --n1 1 --n2 6 --min 3 -b 0.9 -o $fasttree_dir/trees.FastTree_result --indir $fasttree_dir/trees --o2 $fasttree_dir/trees.summary
}


function findBestReciprocalBlast(){
	local full_length_dir=$1
	local is_mauve=$2
	
	#[ $is_mauve == true ] && return

	echo "Finding potential converted genes by gene synteny based on "
	echo "best reciprocal BLAST hits ......"
	cd $full_length_dir >/dev/null
	$ruby $findBestReciprocalBlast -i ../OrthoMCL/all_VS_all.out.tab --summary ../FastTree/trees.summary -o ../FastTree/GC.raw_result

	if [ $? == 0 ]; then
		for i in `seq 10`; do processbar $i 10; sleep 0.1; done
	fi
	echo ''

	cd - >/dev/null
}


function getCandidateGC(){
	local integrateTD2GC=$1
	local full_length_dir=$2
	local gc_count_min=$3
	local ortholog_count_min=$4
	local bootstrap_min=$5
	local is_mauve=$6

	cd $full_length_dir >/dev/null
	echo "Generating final output ......"

	if [ $is_mauve == true ]; then
		$ruby $getCandidateGC -i ../FastTree/GC.mauve.raw_result --gc_count_min $gc_count_min --ortho_count_min $ortholog_count_min -b $bootstrap_min > ../FastTree/GC.mauve.result
	fi

	$ruby $getCandidateGC -i ../FastTree/GC.raw_result --gc_count_min $gc_count_min --ortho_count_min $ortholog_count_min -b $bootstrap_min > ../FastTree/GC.result

	ruby $integrateTD2GC --td ../TD/TD.list -i ../FastTree/GC.result 

	if [ $? == 0 ]; then
		for i in `seq 10`; do processbar $i 10; sleep 0.1; done;
		echo -e "\nDone!";
	fi

	cd - >/dev/null
}


function generate_TD(){
	local orthomcl_dir=$1;
	local outdir=$2;
	local get_TD=$3;
	local TD_dir=$4;
	echo "Identitying tandem duplicates ......"
	for i in $outdir/sequences/*gff; do
		gff_arg="$gff_arg -g $i"
	done
	python $get_TD $gff_arg -i $orthomcl_dir/all_VS_all.out.tab --type gff3 --gene_regexp ".+\|(.+)" --attr ID -n 5 -d 20000 > $TD_dir/TD.list
	for i in `seq 10`; do processbar $i 10; sleep 0.1; done;
	echo
}


function usage(){
	if which figlet 2> /dev/null 1>&2; then
		figlet "iSeeCE" -k
	elif which banner 2>/dev/null 1>&2; then
		banner "iSeeCE"
	else
		echo -e "${BOLD}iSeeCE$NC";
	fi
	echo "================================================="

	word_limit_per_line "iSeeCE: a tool to identify recurrent concerted evolution (CE) based on gene synteny and phylogeny" 78
	echo

	echo -e "${BOLD}USAGE$NC"
	echo -e "Required arguments:"
	echo -ne "--indir\t\t\t\t\t"
	word_limit_per_line "The directory containing input files in the format of Genbank" 38 40
	echo

	echo -ne "--outdir\t\t\t\t"
	word_limit_per_line "The output directory (cannot be the current directory)" 38 40
	echo

	echo "Optional arguments:"

	echo -ne "--gc_count_min|--GC_count_min\t\t"
	word_limit_per_line "The minimum of the number of converted genes" 38
	echo -e "\t\t\t\t\tdefault: 5"
	echo

	echo -ne "--ortho_count_min|--ortholog_count_min\t"
	word_limit_per_line "The minimum of the number of syntenic orthologs in the flanking region of a converted gene" 38 40
	echo -e "\t\t\t\t\tdefault: 2"
	echo

	echo -ne "--blast_cpu\t\t\t\t"
	echo "The number of threads for BLASTp"
	echo -e "\t\t\t\t\tdefault: 2"
	echo

	echo -ne "--mafft_cpu\t\t\t\t"
	echo "The number of threads for MAFFT"
	echo -e "\t\t\t\t\tdefault: 2"
	echo

	echo -ne "-b|--bootstrap\t\t\t\t"
	 word_limit_per_line "The minimum support value in FastTree for a node to be considered as converted genes. It should range from 0 to 1. Note that it is different from the traditional bootstrap value. For more details, please visit http://www.microbesonline.org/fasttree/" 38 40
	echo -e "\t\t\t\t\tdefault: 0"
	echo

	echo -ne "--is_mauve|--is_Mauve\t\t\t"
	word_limit_per_line "Identification of syntenic orthologs will be performed with Mauve." 38 40
	echo -e "\t\t\t\t\tdefault: disabled"
	echo

	echo -ne "--mauve|--Mauve\t\t\t\t"
	echo "The path to progressiveMauve"
	echo -e "\t\t\t\t\tdefault: disabled"
	echo

	echo -ne "--mauve_jar|--Mauve_jar\t\t\t"
	echo "The path to Mauve.jar"
	echo -e "\t\t\t\t\tdefault: disabled"
	echo

	echo -e "-h|--h|--help\t\t\t\thelp message"
	echo

	echo -e "Feel free to use, edit and deliver the codes for any use."
	echo -e "Please write to sishuowang@hotmail.ca for bug report or any question. Your help is highly appreciated!"
	echo
	exit 1
}


function determine_YN(){
	local outdir=$1

	echo -e "Do you want to ignore this issue (the outdir $outdir will not be overwritten but new files will still be output to $outdir)? ${BOLD}[Y/N]${NC}"
	read -r -n1 char
	echo ''
	case $char in
		y|Y)
			echo "Continue ......"
			;;
		n|N)
			echo "Please use --force if you want to remove the outdir $outdir."
			echo "Exiting ......"
			usage
			;;
		*)
			echo -e "Input has to be ${BOLD}Y/N${NC} (case insensitive)."
			determine_YN $outdir
			;;
	esac
}


##############################################################
[ $# == 0 ] && usage

while [ $# -gt 0 ]; do
	case $1 in
		--indir)
			cd $2 >/dev/null 2>&1
			indir=$PWD
			cd - >/dev/null
			shift
			;;
		--outdir)
			outdir=$2
			shift
			;;
		--gc_count_min|--GC_count_min)
			gc_count_min=$2
			shift
			;;
		--ortho_count_min|--ortholog_count_min)
			ortholog_count_min=$2
			shift
			;;
		--blast_cpu)
			blast_cpu=$2
			shift
			;;
		--mafft_cpu)
			mafft_cpu=$2
			shift
			;;
		-b|--bootstrap)
			bootstrap_min=$2
			shift
			;;
		--is_mauve)
			is_mauve=true
			;;
		--mauve)
			Mauve=$2
			shift
			;;
		--mauve_jar)
			Mauve_jar=$2
			shift
			;;
		--force)
			is_force=true
			;;
		-h|--h|--help)
			usage
			;;
		*)
			echo -e "Unknown argument ${BOLD}$1${NC}. Exiting ......"
			usage
			;;
	esac
	shift
done


##############################################################
# check arguments
if [ -z $indir ]; then
	echo -e "${BOLD}indir$NC has to be given by '--indir'!"
	usage
elif [ -z $outdir ]; then
	echo -e "${BOLD}outdir$NC has to be given by '--outdir'!"
	usage
elif [ $is_mauve == true ]; then
	if [ -z "$Mauve" -o ! -f "$Mauve" ]; then
		echo "Mauve has to be given by '--mauve' if '--is_mauve' is specified."
		usage
	elif [ -z "$Mauve_jar" -o ! -f "$Mauve_jar" ]; then
		echo "Mauve_jar has to be given by '--mauve_jar' if '--is_mauve' is specified."
		usage
	fi
fi


##############################################################
if [ -d $outdir ]; then
	if [ $outdir == '.' -o $outdir == './' ]; then
		echo "outdir cannot be the current directory! Exiting ......"
		usage
	fi

	if [ -d $outdir ]; then
		if [ $is_force == true ]; then
			rm -rf $outdir
		else
			echo "The outdir $outdir has already existed."
			determine_YN $outdir
			#usage
		fi
	fi
	mkdir -p $outdir
fi


cd $indir >/dev/null
indir=`pwd`
cd - >/dev/null


##############################################################
# Path to other scripts
genbank2cds=$CWD/genbank2cds.sh
selectOrthoFromMauveOrthomcl=$CWD/selectOrthoFromMauveOrthomcl.rb
check_aln_group=$CWD/check_aln_group.rb
help2run_orthomcl=$CWD/help2run_orthomcl.sh
newFastTreeWrapper=$CWD/newFastTreeWrapper.sh
convertMauveResultToRawResult=$CWD/convertMauveResultToRawResult.rb
putativeGC_from_seqSimilarity=$CWD/putativeGC_from_seqSimilarity.rb
filterFastTree=$CWD/filterFastTree.rb
findBestReciprocalBlast=$CWD/bestHit.rb
getCandidateGC=$CWD/getCandidateGC.rb
get_TD=$CWD/scripts/get_TD.py
integrateTD2GC=$CWD/scripts/integrateTD2GC.rb


##############################################################
TD_dir=$outdir/TD
aln_dir=$outdir/aln
mauve_dir=$outdir/mauve
orthomcl_dir=$outdir/OrthoMCL
fasttree_dir=$outdir/FastTree
full_length_dir=$outdir/full_length

mkdir -p $TD_dir
mkdir -p $aln_dir
mkdir -p $mauve_dir
mkdir -p $orthomcl_dir
mkdir -p $fasttree_dir/trees
mkdir -p $full_length_dir/flanking2

ln -s $indir $outdir/sequences 2>/dev/null


##############################################################
parse_genbank_files $indir $genbank2cds

[ $is_mauve == true ] && run_mauve $mauve_dir $Mauve $Mauve_jar

run_orthomcl $indir $orthomcl_dir $blast_cpu

run_mafft $outdir $aln_dir $selectOrthoFromMauveOrthomcl $check_aln_group $mafft_cpu

run_FastTree $aln_dir

generate_TD $orthomcl_dir $outdir $get_TD $TD_dir

parse_FastTree $fasttree_dir $mauve_dir $ortholog_count_min $bootstrap_min $is_mauve

findBestReciprocalBlast $full_length_dir $is_mauve

getCandidateGC $integrateTD2GC $full_length_dir $gc_count_min $ortholog_count_min $bootstrap_min $is_mauve


##############################################################
#$ruby $CWD/makeAlnConcise.rb --aln_dir aln/combined/ --mauve mauve/mauve.ortho --indir pre_raxml --tree_result $fasttree_dir/trees.FastTree_result | sort > $fasttree_dir/trees.FastTree_result.info


