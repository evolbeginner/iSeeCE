#! /bin/env ruby


require "getoptlong"
require "bio"
require 'parallel'

require "Dir"


#######################################################################
orthomcl_file = nil
mauve_file = nil
pattern = nil
seq_files = Array.new
seq_indir = nil
seq_file_suffix = nil
prefixes = Array.new
group_size_min = nil
group_size_max = nil
mauve_size_min = nil
mauve_size_max = nil
outdir = nil
cpu = 1
is_force = false
is_output = true

seq_objs = Hash.new
prefix_rela = Hash.new
mauve_rela = Hash.new


#######################################################################
def read_mauve(mauve_file)
  mauve_rela = Hash.new
  File.open(mauve_file, 'r').each_line do |line|
    #0:SACI_RS03420:575816-576373	2:STK_RS01680:318927-319484
    line.chomp!
    line_arr = line.split("\t")
    genes = line_arr.map{|item|item.split(':')[1]}
    line_arr.each do |item|
      item_arr = item.split(':')
      gene = item_arr[1]
      mauve_rela[gene] = genes.select{|i|i if i != gene}
    end
  end
  return(mauve_rela)
end


def get_orthomcl(orthomcl_file, mauve_rela, group_size_min, group_size_max, mauve_size_min, mauve_size_max)
  orthomcl_groups = Array.new
  in_fh = orthomcl_file == '-' ? STDIN : File.open(orthomcl_file, 'r')
  in_fh.each_line do |line|
    temp_rela = Hash.new
    line.chomp!
    line_arr = line.split("\t")
    # criteria
    next if not (group_size_min .. group_size_max).include?(line_arr.size)

    orgns = line_arr.map{|i|i.split('|')[0]}
    genes = line_arr.map{|i|i.split('|')[1]}
    next if orgns.inject(Hash.new(0)) {|hash,word| hash[word] += 1; hash}.values.count{|i|i>1} < 2

    genes.each do |gene|
      if ! mauve_rela.empty? and ! mauve_rela[gene].empty? and mauve_rela.include?(gene)
        temp_rela[gene] = genes.select{|i|i!=gene}.map{|i|i if mauve_rela[gene].include?(i)}.compact
      end
    end
    # criteria
    if (mauve_size_min .. mauve_size_max).include?(temp_rela.size)
      orthomcl_groups << genes
    end
  end
  in_fh.close if orthomcl_file != '-'
  return(orthomcl_groups)
end


def read_seq(seq_file, prefix)
  seq_objs = Hash.new
  prefix_rela = Hash.new
  Bio::FlatFile.open(seq_file).each_entry do |f|
    seq_name = f.definition
    prefix_rela[seq_name] = prefix
    seq_objs[seq_name] = f
  end
  return([seq_objs, prefix_rela])
end


def output_seq(orthomcl_groups, seq_objs, prefix_rela, outdir, cpu)
  #orthomcl_groups.each_with_index do |orthomcl_group, index|
  Parallel.map(orthomcl_groups.each_with_index, in_processes: cpu) do |orthomcl_group, index|
    outfile = File.join([outdir, (index+1).to_s+".fas"])
    out_aln = File.join([outdir, (index+1).to_s+".aln"])
    out_fh = File.open(outfile, 'w')
    orthomcl_group.each do |gene|
      puts gene if not prefix_rela.include?(gene)
      out_fh.puts ">" + prefix_rela[gene] + '|' + gene
      out_fh.puts seq_objs[gene].seq
    end
    out_fh.close
    `mafft --thread #{cpu} --quiet #{outfile} > #{out_aln}`
  end
end


#######################################################################
opts = GetoptLong.new(
  ['--orthomcl', GetoptLong::REQUIRED_ARGUMENT],
  ['--mauve', GetoptLong::REQUIRED_ARGUMENT],
  ['--pattern', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_file_suffix', GetoptLong::REQUIRED_ARGUMENT],
  ['--prefix', GetoptLong::REQUIRED_ARGUMENT],
  ['--group_size', GetoptLong::REQUIRED_ARGUMENT],
  ['--mauve_size', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--no_output', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when "--orthomcl"
      orthomcl_file = value
    when "--mauve"
      mauve_file = value
    when "--pattern"
      pattern = value
    when "--seq"
      seq_files << value.split(',')
    when "--seq_indir"
      seq_indir = File::expand_path(value)
    when '--seq_file_suffix'
      seq_file_suffix = value
    when "--prefix"
      prefixes << value.split(',')
    when "--group_size"
      group_size_min, group_size_max = value.split(',').map{|i|i.to_i}
    when "--mauve_size"
      mauve_size_min, mauve_size_max = value.split(',').map{|i|i.to_i}
    when "--outdir"
      outdir = value
    when "--cpu"
      cpu = value.to_i
    when "--force"
      is_force = true
    when "--no_output"
      is_output = false
  end
end


if not seq_indir.nil? and seq_files.empty? and not seq_file_suffix.nil?
  Dir.foreach(seq_indir).each do |file_basename|
    next if file_basename =~ /^\./
    next if file_basename !~ /#{seq_file_suffix}/
    seq_files << File.join([seq_indir, file_basename])
  end
end


seq_files.flatten!
if prefixes.empty?
  seq_files.each do |seq_file|
    basename = File.basename(seq_file)
    basename =~ /^([^.]+)/
    prefix = $1
    prefixes << prefix
  end
else
    prefixes.flatten!
end


#######################################################################
if seq_files.empty?
  STDERR.puts "seq_files empty! Exiting ......"
  exit
end

if prefixes.empty?
  STDERR.puts "prefix not given! Exiting ......"
  exit
end

if not [group_size_min, group_size_max, mauve_size_min, mauve_size_max].select{|i|i.nil?}.empty?
  STDERR.puts "group_size and mauve_size have to be given! Exiting ......"
  exit
end 

mkdir_with_force(outdir, is_force) if is_output


#######################################################################
seq_files.each_with_index do |seq_file, index|
  prefix = prefixes[index]
  arr = read_seq(seq_file, prefix)
  seq_objs.merge!(arr[0])
  prefix_rela.merge!(arr[1])
end

mauve_rela = read_mauve(mauve_file) if ! mauve_file.nil? and File.exists?(mauve_file)

orthomcl_groups = get_orthomcl(orthomcl_file, mauve_rela, group_size_min, group_size_max, mauve_size_min, mauve_size_max)
puts orthomcl_groups.size

output_seq(orthomcl_groups, seq_objs, prefix_rela, outdir, cpu) if is_output


