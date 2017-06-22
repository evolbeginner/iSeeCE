#! /bin/bash


#####################################################################
genbank2gff=scripts/genbank2gff.rb
extractGenomeFromGenbank=scripts/extractGenomeFromGenbank.rb
gff2fasta=scripts/gff2fasta.pl
translate_seq=scripts/translate_seq.rb
#genbank2fasta=~/tools/self_bao_cun/basic_process_mini/genbank2fasta.py


#####################################################################
while [ $# -gt 0 ]; do
	case $1 in
		--gb)
			gb_file=$2
			shift
			;;
		--out_prefix)
			out_prefix=$2
			shift
			;;
		*)
			echo "Unkonw argument $1! Exiting ......"
			exit
			;;
	esac
	shift
done


if [ -z $out_prefix ]; then
	echo "out_prefix has to be specified! Exiting ......"
	exit
fi

mkdir -p `dirname $out_prefix`
out_genome=$out_prefix.genome.fa
out_gff=$out_prefix.gff
out_cds=$out_prefix.cds.fas
out_pep=$out_prefix.pep.fas


#####################################################################
ruby $extractGenomeFromGenbank -i $gb_file --no_plasmid > $out_genome

ruby $genbank2gff --seq $gb_file --feature gene,CDS > $out_gff

perl $gff2fasta --gff $out_gff -d --fasta $out_genome --feature CDS --length 1 --out $out_cds --attribute ID --no_process_seqid >/dev/null
sed -i '/^>/s/|.\+//' $out_cds

ruby ~/tools/self_bao_cun/basic_process_mini/translate_seq.rb -i $out_cds -o $out_pep
#python $genbank2fasta -i $gb_file -o $out_pep


