#! /bin/bash


##############################################################
Mauve=/mnt/bay3/sswang/software/genome/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve
Mauve_jar=/mnt/bay3/sswang/software/genome/mauve_snapshot_2015-02-13/Mauve.jar


##############################################################
while [ $# -gt 0 ]; do
	case $1 in
		--indir)
			cd $2 >/dev/null
			indir=$PWD
			cd - >/dev/null
			shift
			;;
		--outdir)
			outdir=$2
			shift
			;;
		*)
			echo "Unknown argument $1. Exiting ......"
			exit 1
			;;
	esac
	shift
done


##############################################################
if [ -d $outdir ]; then
	if [ $outdir != '.' -a $outdir != './' ]; then
		echo "haha!"
		#echo "The outdir $outdir has already existed. Exiting ......"
		#exit 1
	fi
else
	mkdir -p $outdir
fi

CWD=`dirname $0`
cd $CWD >/dev/null
CWD=`pwd`
cd - >/dev/null

cd $indir >/dev/null
indir=`pwd`
cd - >/dev/null


##############################################################
mauve_dir=$outdir/mauve
orthomcl_dir=$outdir/OrthoMCL
fasttree_dir=$outdir/FastTree/8_4/
mkdir $mauve_dir
mkdir $orthomcl_dir
mkdir -p $fasttree_dir

ln -s $indir $outdir/sequences


##############################################################
echo "Parsing Genbank files in $indir ......"
for i in $indir/*gb; do 
	bash $CWD/genbank2cds.sh --gb $i --out_prefix ${i%.gb};
done

echo "Doing genome alignment with Mauve ......"
cd $mauve_dir >/dev/null
#$Mauve --output=mauve ../sequences/*gb
java -cp $Mauve_jar org.gel.mauve.analysis.OneToOneOrthologExporter -f mauve -o mauve.ortho
[ ! -s mauve.ortho ] && echo "The file size of mauve.ortho equals 0! Exiting ......" && exit 1
cd - >/dev/null

echo "Running OrthoMCL ......"
run_orthomcl_cmd=`bash $CWD/help2run_orthomcl.sh --indir $indir`
cd $orthomcl_dir >/dev/null
#$run_orthomcl_cmd >/dev/null
for i in 6 4 2; do mcl mclInput --abc -I $i -o mclOutput$i; done
cp mclOutput mclOutput1.5
cd - >/dev/null


##############################################################
cd $outdir >/dev/null
echo "Running MAFFT ......"
for i in 1.5 2 4 6; do
	#ruby $CWD/selectOrthoFromMauveOrthomcl.rb --mauve mauve/mauve.ortho --orthomcl OrthoMCL/mclOutput$i --seq_indir sequences --seq_file_suffix cds.fas --outdir aln/8_4/aln$i --force --group_size 4,200 --mauve_size 4,200;
	echo ""
done

ruby2.1 $CWD/check_aln_group.rb --indir aln/8_4/aln1.5/ --indir aln/8_4/aln2/ --indir aln/8_4/aln4/ --indir aln/8_4/aln6 --suffix aln --outdir aln/8_4/combined --force


echo "Running FastTree ......"
cd aln/8_4/combined/
#for i in `ls`; do FastTree -quiet -nt < $i > ${i/aln/FastTree.tre} 2>/dev/null; done;
mv *tre ../../../FastTree/8_4/;
cd - >/dev/null


echo "Parsing FastTree result ......"
for i in FastTree/8_4/*; do
	ruby $CWD/putativeGC_from_seqSimilarity.rb --tree $i --mauve mauve/mauve.ortho -b 0
done > FastTree/8_4.FastTree_result 2>/dev/null
ruby $CWD/makeAlnConcise.rb --aln_dir aln/8_4/combined/ --mauve mauve/mauve.ortho --indir pre_raxml --tree_result FastTree/8_4.FastTree_result | sort > FastTree/8_4.FastTree_result.info


