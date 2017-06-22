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

CWD=`pwd`


##############################################################
mauve_dir=$outdir/mauve
mkdir $mauve_dir


##############################################################
echo "Parsing Genbank files in $indir ......"
for i in $indir/*gb; do 
	bash genbank2cds.sh --gb $i --out_prefix ${i%.gb};
done

echo "Doing genome alignment with Mauve"
cd $mauve_dir >/dev/null
#$Mauve --output=mauve $indir/*gb
java -cp $Mauve_jar org.gel.mauve.analysis.OneToOneOrthologExporter -f mauve -o mauve.ortho
cd - >/dev/null


