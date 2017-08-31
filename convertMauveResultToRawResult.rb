#! /bin/env ruby


require 'getoptlong'


######################################################
infile = nil
ortholog_count_min = 2


######################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--ortholog_count_min', '--ortho_count_min', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^--ortholog_count_min$|^--ortho_count_min$/
      ortholog_count_min = value.to_i
  end
end


######################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  #FROM
  #trees/193.FastTree.tre	Mace|MA_RS00755	Mace|MA_RS15580	MBAR_RS19280	MBAR_RS04540,MSTHT_RS10190	0.999
  #TO
  #200	MA_RS20650	3	MA_RS20670	4	0.996
  line_arr = line.split("\t")
  genes = line_arr[1,2]
  genes.map!{|i|i.split('|')[1]}
  orthologs = line_arr[3,2]
  bootstrap = line_arr[-1].to_f
  file_name = line_arr[0]
  b = File.basename(file_name)
  b =~ /(\d+)\./
  corename = b
  next if orthologs.select{|i|i.split(',').size >= ortholog_count_min}.size < 2
  puts [b, genes[0], orthologs[0].split(',').size, genes[1], orthologs[1].split(',').size, bootstrap].join("\t")
end


