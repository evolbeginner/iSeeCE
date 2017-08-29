#! /bin/bash


###############################################################
gene=$1
n=$2
is_locus=$3

if [ -z $n ]; then n=100; fi



###############################################################
if [ -z $is_locus ]; then
	grep -P $1 ../sequences/*gb -A $n -B $n | awk 'BEGIN{ORS="";FS="product="}{if($0~/product/){FS="product="; if($0!~/"$/){a=1;print $2" "; next}else{print $2"\n"}} if(a==1){print $0"\n"; a=0}}'  | sed 's!\.\./sequences[^ ]\+ \+!!' | sed 's/^ //; s/^"//; s/"$//'
	#grep -P $1 ../sequences/*gb -A $n -B $n | grep product #awk 'BEGIN{ORS=" "}{if($0~/product/){FS="product="; if($0!~/"$/){a=1;print $2" "; next}else{print $2"\n"}} if(a==1){print $0"\n"; a=0}}'  | sed 's!\.\./sequences[^ ]\+ \+!!' | sed 's/  \+/ /g'
else
	grep -P $1 ../sequences/*gb -A $n -B $n | grep 'product\|locus_tag' | grep -v old | uniq
	#grep -P $1 ../sequences/*gb -A $n -B $n | awk 'BEGIN{ORS="";FS="product="}{if($0~/product\locus_tag/){FS="="; if($0!~/"$/){a=1;print $2" "; next}else{print $2"\n"}} if(a==1){print $0"\n"; a=0}}'  | sed 's!\.\./sequences[^ ]\+ \+!!' | sed 's/^ //; s/^"//; s/"$//'
fi


