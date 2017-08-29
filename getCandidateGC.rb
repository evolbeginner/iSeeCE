# /bash/env ruby


#ruby getCandidateGC.rb -i haha/FastTree/GC.raw_result  --gc_count_min 3 --ortholog_count_min 4
require 'getoptlong'


##############################################
infile = nil
gc_count_min = 5
ortholog_count_min = 2

gc_info = Hash.new{|h,k|h[k]={}}


##############################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--gc_count_min', GetoptLong::REQUIRED_ARGUMENT],
  ['--ortho_count_min', '--ortholog_count_min', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--gc_count_min'
      gc_count_min = value.to_i
    when '--ortho_count_min', '--ortholog_count_min'
      ortholog_count_min = value.to_i
  end
end


##############################################
in_fh = (infile.nil? or infile == '-') ? STDIN : File.open(infile, 'r')
in_fh.each_line do |line|
  #200	MA_RS20650	3	MA_RS20670	4
  line.chomp!
  line_arr = line.split("\t")
  locus = line_arr[0]
  genes = line_arr.values_at(1,3)
  if line_arr.values_at(2,4).select{|i|i.to_i >= ortholog_count_min}.size == 2
    gc_info[locus][genes] = ''
  end
end
in_fh.close unless infile.nil?


gc_info.each_pair do |locus, v|
  if v.size >= gc_count_min
    puts [locus, v.keys.map{|i|i.join(',')}].join("\t")
  end
end


