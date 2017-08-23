#! /bin/bash


######################################################################
CWD=`dirname $0`
cd $CWD >/dev/null
CWD=`pwd`
cd - >/dev/null
for i in ~/tools/self_bao_cun/packages/bash/processbar.sh $CWD/scripts/processbar.sh; do
	[ -f $i ] && source $i
done


######################################################################
dir=`dirname $0`

getNumOfTaxa=$dir/getNumOfTaxa.rb
putativeGC_from_seqSimilarity=$dir/putativeGC_from_seqSimilarity.rb
filterFastTree=$dir/filterFastTree.rb

#new_FastTree_result=FastTree/new.FastTree_result2


######################################################################
bootstrap=0.9
count_min=5
n1=1
n2=6
suffix=Fasttree.tre
summary_outfile=''


######################################################################
while [ $# -gt 0 ]; do
	case $1 in
		--indir)
			indir=$2
			shift
			;;
		-o|--o1)
			new_FastTree_result=$2
			shift
			;;
		--o2)
			summary_outfile=$2
			shift
			;;
		--n1)
			n1=$2
			shift
			;;
		--n2)
			n2=$2
			shift
			;;
		--min)
			count_min=$2
			shift
			;;
		-b|--bootstrap)
			bootstrap=$2
			shift
			;;
		--suffix)
			suffix=$2
			shift
			;;
	esac
	shift
done


if [ -z $summary_outfile ]; then
	summary_outfile='-'
fi


######################################################################
if [ ! -d $indir ]; then
	echo "echo FastTree/8_4 does not exist. Exiting ......"
	exit 1
fi

num_of_taxa=`ruby $getNumOfTaxa --indir $indir`

otu_min=`expr $num_of_taxa \* $n1`
otu_max=`expr $num_of_taxa \* $n2`


######################################################################
[ -f $new_FastTree_result ] && rm $new_FastTree_result
total=`ls -1 $indir/*FastTree.tre | wc -l`
count=0
for i in $indir/*FastTree.tre; do
	if [ `declare -f processbar > /dev/null; echo $?` == 0 ]; then
		count=`expr $count + 1`
		processbar $count $total
	fi
	ruby2.1 $putativeGC_from_seqSimilarity --tree $i -b 0 --otu_minmax $otu_min,$otu_max >> $new_FastTree_result
done
echo ''

#ruby $filterFastTree -i $new_FastTree_result --count_min 3 --tree_dir full_length/phylo/ --aln_dir full_length/aln/ --bootstrap 0.9
#ruby $filterFastTree -i $new_FastTree_result --count_min 5 --tree_dir full_length/phylo/ --aln_dir full_length/aln/ --bootstrap 0.9
ruby $filterFastTree -i $new_FastTree_result --count_min $count_min --bootstrap $bootstrap -o $summary_outfile


