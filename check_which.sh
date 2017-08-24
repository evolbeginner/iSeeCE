function check_which(){
	for i in $@; do
		which $i 1>/dev/null 2>&1;
		if [ $? == 0 ]; then
			echo -e "${i} is detected."
		else
			echo -e "\e[0;1;31m${i}\e[0m has NOT been installed."
		fi
	done
}
