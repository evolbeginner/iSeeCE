# /bash/env ruby


#ruby getCandidateGC.rb -i haha/FastTree/GC.raw_result  --gc_count_min 3 --ortholog_count_min 4
require 'getoptlong'


##############################################
infile = nil
bootstrap_min = 0
gc_count_min = 5
ortholog_count_min = 2

gc_info = Hash.new{|h,k|h[k]={}}


##############################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['-b', '--bootstrap', GetoptLong::REQUIRED_ARGUMENT],
  ['--gc_count_min', GetoptLong::REQUIRED_ARGUMENT],
  ['--ortho_count_min', '--ortholog_count_min', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '-b', '--bootstrap'
      bootstrap_min = value.to_f
    when '--gc_count_min'
      gc_count_min = value.to_i
    when '--ortho_count_min', '--ortholog_count_min'
      ortholog_count_min = value.to_i
  end
end


##############################################
in_fh = (infile.nil? or infile == '-') ? STDIN : File.open(infile, 'r')
in_fh.each_line do |line|
  #200	MA_RS20650	3	MA_RS20670	4 0.956
  line.chomp!
  line_arr = line.split("\t")
  locus = line_arr[0]
  genes = line_arr.values_at(1,3)
  bootstrap = line_arr[-1].to_f
  next if bootstrap < bootstrap_min
  if line_arr.values_at(2,4).select{|i|i.to_i >= ortholog_count_min}.size == 2
    gc_info[locus][genes] = bootstrap
  end
end
in_fh.close unless infile.nil?


##############################################
puts %w[#id gene1,gene2:bootstrap].join("\t")

gc_info.each_pair do |locus, v|
  if v.size >= gc_count_min
    puts [locus, v.keys.map{|i|i.join(',') + ':' + v[i].to_s}].join("\t")
  end
end


