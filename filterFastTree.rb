#! /bin/env ruby


require 'getoptlong'
require 'bio'


##########################################################
infile = nil
outfile = nil
count_min = 0
tree_dir = nil
aln_dir = nil
bootstrap_min = 0.0

file_orgn_info = Hash.new{|h1,k1|h1[k1]=Hash.new{|h2,k2|h2[k2]={}}}
tree_basenames = Array.new
genes_in_aln = Hash.new


##########################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['-o', GetoptLong::REQUIRED_ARGUMENT],
  ['--count_min', GetoptLong::REQUIRED_ARGUMENT],
  ['--tree_dir', GetoptLong::REQUIRED_ARGUMENT],
  ['--aln_dir', GetoptLong::REQUIRED_ARGUMENT],
  ['--bootstrap', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '-o'
      outfile = value
    when '--count_min'
      count_min = value.to_i
    when '--tree_dir'
      tree_dir = value
    when '--aln_dir'
      aln_dir = value
    when '--bootstrap'
      bootstrap_min = value.to_f
  end
end


##########################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  file, gene1, gene2 = line_arr
  bootstrap = line_arr[-1]
  bootstrap = bootstrap.to_f
  orgn = gene1.split('|')[0]
  file_orgn_info[file][orgn][[gene1, gene2]] = bootstrap
end
in_fh.close


if ! tree_dir.nil?
  Dir.foreach(tree_dir) do |b|
    next if b =~ /^\./
    tree_basenames << b
  end
end


if not aln_dir.nil?
  tree_basenames.each do |b|
    aln_file = File.join(aln_dir, b+'.aln')
    next if not File.exist?(aln_file)
    Bio::FlatFile.open(aln_file, 'r').each_entry do |f|
      genes_in_aln[f.definition] = ''
    end
  end
end


##########################################################
# output results
out_fh = (! outfile.nil? or outfile == '-') ? File.open(outfile, 'w') : STDOUT
file_orgn_info.each_pair do |file, v|
  next if v.size <  count_min
  b = File.basename(file)
  c = b.split('.')[0]
  next if tree_basenames.include?(c)
  
  passed_taxa = Hash.new
  count = 0
  pair_count = 0
  passed_pair_count = 0
  v.each_pair do |taxon, genes_arr|
    genes_arr.each_pair do |genes, bootstrap|
      pair_count += 1
      genes.each do |gene|
        if ! genes_in_aln.empty? and genes_in_aln.include?(gene)
          count += 1
        end
      end
      if bootstrap >= bootstrap_min
        passed_taxa[taxon] = ''
        passed_pair_count += 1
      end
    end
  end
  
  out_fh.puts [file, v.size, passed_taxa.size, pair_count, passed_pair_count, v.values.map{|genes_arr|genes_arr.keys.map{|genes|[genes.join(','),genes_arr[genes]].join(':')}}].join("\t") if count <= 1
end


