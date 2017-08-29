#! /bin/bash


function word_limit_per_line(){
	local input=$1
	local word_limit=$2
	local num_of_space=$3

	local total_word_count=0
	local num_of_line=0
	local sentence
	local is_pass=false
	local arr=()

	for i in `echo $input | grep -o "\([^ ]\+\)"`; do
		arr=(${arr[@]} $i);
	done

	for index in ${!arr[@]}; do
		i=${arr[$index]}
		word_count=`expr length $i`
		total_word_count=`expr $word_count + $total_word_count`
		#echo $word_count
		#echo $total_word_count

		if [ $total_word_count -gt $word_limit ]; then
			((num_of_line++))
			if [ $num_of_line -gt 1 ]; then
				echo_space $num_of_space
			fi
			echo $sentence
			sentence=''
			total_word_count=$word_count
		fi

		sentence=$sentence" "$i
	done

	if [ ! -z "$sentence" ]; then
		echo_space $num_of_space
		echo $sentence
	fi
}


function echo_space(){
	local num_of_space=$1
	if [ ! -z $num_of_space ]; then
		for j in `seq $num_of_space`; do
			echo -ne " "
		done
	fi
}


