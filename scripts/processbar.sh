function processbar() {
  local current=$1; local total=$2;
  local maxlen=80; local barlen=66; local perclen=14;
  local format="%-${barlen}s%$((maxlen-barlen))s"
  #local perc="[$current/$total]"
  local perc=`echo "scale=2;$current*100/$total" | bc`"%"
  local progress=$((current*barlen/total))
  local prog=$(for i in `seq 0 $progress`; do printf '#'; done)
  printf "\r$format" $prog $perc
}


