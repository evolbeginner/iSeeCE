#! /bin/env ruby


require 'getoptlong'


###################################################################
infile = nil
td_file = nil

final_lines = Array.new


###################################################################
def read_td_file(td_file)
  tds = Hash.new()
  in_fh = File.open(td_file, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    line_arr.each_with_index do |i, index1|
      line_arr.each_with_index do |j, index2|
        next if index1 == index2
        tds[[i,j].sort] = ''
      end
    end
  end
  in_fh.close
  return(tds)
end


###################################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--TD', '--td', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^--(TD|td)$/
      td_file = value
  end
end


###################################################################
tds = read_td_file(td_file)


###################################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  final_lines << line and next if line =~ /^#/
  line.chomp!
  #270	MA_RS20650,MA_RS20670:1.0	MBAR_RS03690,MBAR_RS03705:1.0	MSTHT_RS08150,MSTHT_RS08165:0.964
  line_arr = line.split("\t")
  id = line_arr[0]
  strs = line_arr[1,line_arr.size-1]
  strs.each_with_index do |str, index|
    str_arr = str.split(':')
    genes = str_arr[0].split(',')
    if tds.include?(genes.sort)
      strs[index] = [str, '*'].join('')
    end
  end
  final_lines << [id, strs].flatten.join("\t")
end
in_fh.close


###################################################################
out_fh = File.open(infile, 'w')
final_lines.each do |line|
  out_fh.puts line
end
out_fh.close


