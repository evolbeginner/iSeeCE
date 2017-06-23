#! /bin/env ruby


require 'getoptlong'
require 'find'

require 'bio'


############################################################################
indir = nil
mauve_file = nil
aln_dir = nil
tree_result_file = nil

pairs = Hash.new{|h,k|h[k]=[]}
gap_props = Hash.new{|h,k|h[k]={}}


############################################################################
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


def read_aln(aln_dir)
  seq_num_in_aln = Hash.new
  Dir.foreach(aln_dir) do |file_basename|
    next if file_basename =~ /^\./
    file_basename =~ /^(\d+)/
    file_fullname = File.join([aln_dir, file_basename])
    corename = $1
    seq_num_in_aln[corename] = `grep '^>' #{file_fullname} | wc -l`.chomp
  end
  return(seq_num_in_aln)
end


############################################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--mauve', GetoptLong::REQUIRED_ARGUMENT],
  ['--aln_dir', GetoptLong::REQUIRED_ARGUMENT],
  ['--tree_result', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
    when '--mauve'
      mauve_file = value
    when '--aln_dir'
      aln_dir = value
    when '--tree_result'
      tree_result_file = value
  end
end


############################################################################
mauve_rela = read_mauve(mauve_file)

seq_num_in_aln = read_aln(aln_dir)


############################################################################
if not tree_result_file.nil?
  File.open(tree_result_file, 'r').each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    pair = line_arr.values_at(1,2).join("\t")
    File.basename(line_arr[0]) =~ /^([^.]+)/
    corename = $1
    pairs[pair] << corename
  end
elsif not indir.nil?
  Find.find(indir) do |path|
    b = File.basename(path)
    d = File.dirname(path)
    b2 = File.basename(d)
    if b =~ /^[^.]+\.gc_range$/
      in_fh = File.open(path, 'r')
      in_fh.each_line do |line|
        line.chomp!
        line_arr = line.split("\t")
        if line_arr.size == 2
          pair = line
          pairs[pair] << b2
        end
      end
    end
  end
end


pairs.each do |pair, v|
  basic_info_item = v.sort_by{|i|seq_num_in_aln[i].to_i}.map{|i|i+'('+seq_num_in_aln[i].to_s+')'}.join("\t")
  genes = pair.split("\t")
  orgn1, corename1 = genes[0].split('|')
  orgn2, corename2 = genes[1].split('|')
  ortho_info_item = [mauve_rela[corename1].join(','), mauve_rela[corename2].join(',')].flatten.join("\t")
  puts [pair, ortho_info_item, basic_info_item].flatten.join("\t")
end


