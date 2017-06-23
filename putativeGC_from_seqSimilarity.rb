#! /bin/env ruby


require "getoptlong"
require "bio"
require "Dir"


####################################################################
infiles = Array.new
geneconv_files = Array.new
tree_files = Array.new
gc_length_min = 50
pvalue_bc_max = 0.05
mauve_file = nil
mauve_size_min = 1
is_de_prefix = false
length_min = 30
bootstrap_min = nil
is_no_bootstrap = false

geneconv_info = Hash.new


####################################################################
class GENECONV
  attr_reader :pvalue, :pvalue_bc, :start, :stop, :length
  def initialize(pvalue, pvalue_bc, start, stop, length)
    @pvalue = pvalue
    @pvalue_bc = pvalue_bc
    @start = start
    @stop = stop
    @length = length
  end
end


####################################################################
def parse_geneconv_file(geneconv_file, geneconv_info, mauve_rela, mauve_size_min)
  file_corename = get_corename_of_file(geneconv_file)
  geneconv_info.each_pair do |gene1, v1|
    orgn1, corename1 = gene1.split('|')
    v1.each_pair do |gene2, v2|
      orgn2, corename2 = gene2.split('|')
      if orgn2 == orgn1
        v2.each do |i|
          next if not (mauve_rela.include?(corename1) and mauve_rela.include?(corename2) and mauve_rela[corename1].size >= mauve_size_min and mauve_rela[corename2].size >= mauve_size_min)
          #puts [gene1, i.pvalue, i.pvalue_bc, i.start, i.stop, i.length].join("\t")
          puts [geneconv_file, [gene1, gene2].sort, mauve_rela[corename1].join(','), mauve_rela[corename2].join(','), i.pvalue, i.pvalue_bc, i.start, i.stop, i.length].flatten.join("\t")
        end
      end
    end
  end
end


def parse_infile(infile, seq_objs, length_min, mauve_rela, mauve_size_min)
  similarities = Hash.new{|h,k|h[k]={}}
  lengths = Hash.new{|h,k|h[k]={}}
  posi_pairs = Hash.new

  file_corename = get_corename_of_file(infile)
  seq_objs.each_pair do |gene1, f|
    seq_objs.keys.each do |gene2|
      next if gene2 == gene1
      genes = [gene1, gene2].sort
      sim, length1, sim2, length2 = calculate_pairwise_similarity(seq_objs[gene1].seq, seq_objs[gene2].seq)
      next if similarities.include?(genes)
      similarities[gene1][gene2] = sim2
      lengths[gene1][gene2] = length2
    end
  end

  similarities.each_pair do |gene1, homolog_info|
    sim1 = nil
    sim2 = nil
    gene2 = nil
    corename1 = nil
    corename2 = nil

    prefix = gene1.split('|')[0]
    gene2, sim2 = parse_sim(homolog_info, prefix)
    if not sim2.nil? and not sim2 == 0 and sim2 == homolog_info.values.compact.max 
      pair = [gene1, gene2].sort.join("\t")
      next if posi_pairs.include?(pair)
      paralog_paralog_sim = similarities[gene2][gene1]
      next if paralog_paralog_sim != similarities[gene2].values.compact.max
      corename1 = gene1.split('|')[1]
      corename2 = gene2.split('|')[1]
      next if not (mauve_rela.include?(corename1) and mauve_rela.include?(corename2) and mauve_rela[corename1].size >= mauve_size_min and mauve_rela[corename2].size >= mauve_size_min)
      next if lengths[gene1][gene2] < length_min
      puts [infile, [gene1, gene2].sort, mauve_rela[corename1].join(','), mauve_rela[corename2].join(','), ((sim2+paralog_paralog_sim)/2).to_s, lengths[gene1][gene2]].flatten.join("\t")
      posi_pairs[pair] = ""
    end
  end
end


def parse_tree_file(tree_file, bootstrap_min, mauve_rela, mauve_size_min, is_no_bootstrap)
  posi_pairs = Hash.new
  treeio = Bio::FlatFile.open(Bio::Newick, tree_file)
  tree = treeio.next_entry.tree
  bootstraps = tree.nodes.map{|node|node.bootstrap}.compact
  bootstrap_style = nil

  file_corename = get_corename_of_file(tree_file)
  tree.nodes.each do |node|
    next if node.name !~ /\w/
    node.name.gsub!(" ", "_")
    parent_node =  tree.parent(node)
    neighbor_nodes = tree.children(parent_node).select{|i|i.name != node.name}
    next if neighbor_nodes.size > 3
    parent_node.bootstrap=0 and bootstrap_min=0 if is_no_bootstrap
    if parent_node.bootstrap.nil?
      next if tree.adjacent_nodes(parent_node).select{|i|i unless i.bootstrap.nil?}.empty?
      parent_node_bootstrap = tree.adjacent_nodes(parent_node).select{|i|i unless i.bootstrap.nil?}[0].bootstrap
      parent_node.bootstrap = parent_node_bootstrap
    end

    if not bootstrap_min.nil?
      ;
    else
      if not is_no_bootstrap
        if bootstraps.max <= 1
          bootstrap_style = "FastTree"
        else
          bootstrap_style = "Newick"
        end
      end
      if bootstrap_style == "FastTree"
        bootstrap_min = 0.9
      elsif bootstrap_style == "Newick"
        bootstrap_min = 70
      else
        puts "bootstrap error! Exiting ......"
        exit
      end
    end
    next if parent_node.bootstrap < bootstrap_min

    orgn1, corename1 = node.name.split('|')
    neighbor_nodes.each do |neighbor_node|
      next if neighbor_node.name !~ /\w/
      orgn2, corename2 = neighbor_node.name.split('|')
      next if orgn1 != orgn2
      next if not (mauve_rela.include?(corename1) and mauve_rela.include?(corename2) and mauve_rela[corename1].size >= mauve_size_min and mauve_rela[corename2].size >= mauve_size_min)
      gene1, gene2 = node.name, neighbor_node.name
      pair = [gene1, gene2].sort.join("\t")
      next if posi_pairs.include?(pair)
      puts [tree_file, [gene1, gene2].sort, mauve_rela[corename1].join(','), mauve_rela[corename2].join(','), parent_node.bootstrap].flatten.join("\t")
      posi_pairs[pair] = ""
    end
  end
end


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


def parse_sim(homolog_info, prefix)
  paralog = nil
  paralog_sim = nil
  homolog_info.each_pair do |homolog, sim|
    homolog_prefix = homolog.split('|')[0]
    if prefix == homolog_prefix
      paralog = homolog
      paralog_sim = sim
    end
  end
  return([paralog, paralog_sim])
end


def calculate_pairwise_similarity(seq1, seq2)
  identity, identity2 = nil, nil
  identical_counter = 0
  both_gap_counter = 0
  total_length = seq1.size
  arr1 = seq1.split("")
  arr2 = seq2.split("")
  arr1.each_with_index do |nucl, index|
    if nucl == arr2[index]
      if nucl == '-'
        both_gap_counter += 1
      else
        identical_counter += 1
      end
    end
  end

  shortest_length_no_gap = [seq1.gsub('-',"").size, seq2.gsub('-',"").size].min
  if seq1.count('-') >= total_length
    identity = nil
  else
    identity = identical_counter/(total_length-both_gap_counter).to_f
  end

  if shortest_length_no_gap == 0
    identity2 = nil
  else
    identity2 = identical_counter/(shortest_length_no_gap).to_f
  end

  return([identity, total_length-both_gap_counter, identity2, shortest_length_no_gap])
end


def read_seq(seq_file, is_de_prefix=false)
  seq_objs = Hash.new
  Bio::FlatFile.open(seq_file).each_entry do |f|
    seq_title = nil
    seq_title = f.definition
    if is_de_prefix
      seq_title.sub!(/.+\|/, "")
    end
    seq_objs[seq_title] = f
  end
  return(seq_objs)
end


def read_geneconv_file(geneconv_file, gc_length_min, pvalue_bc_max)
  geneconv_info = Hash.new{|h1,k1|h1[k1]=Hash.new{|h2,k2|h2[k2]=[]}}
  File.open(geneconv_file, 'r').each_line do |line|
  # #   Names                                   Pvalue  Pvalue   Begin  End   Len  Poly Dif  Difs Pen.
  # GI  sto|STK_RS10565;sisM164|M164_RS03400    0.0319  0.28175  1168   1184   17    12   0  125  None
    line.chomp!
    next if line !~ /^GI/
    line_arr = line.split(/\s+/)
    gene_pair = line_arr[1]
    genes = gene_pair.split(';')
    pvalue, pvalue_bc, start, stop, length = line_arr.values_at(2,3,4,5,6).map{|i|i.to_f}
    next if length < gc_length_min
    next if pvalue_bc > pvalue_bc_max
      geneconv_obj = GENECONV.new(pvalue, pvalue_bc, start, stop, length)
      geneconv_info[genes[0]][genes[1]] << geneconv_obj
    end
  return(geneconv_info)
end


def get_corename_of_file(file)
  basename = File.basename(file)
  corename = basename.split(".")[0]
  return(corename)
end


####################################################################
opts = GetoptLong.new(
  ["--mauve", GetoptLong::REQUIRED_ARGUMENT],
  ["-i", "--aln", GetoptLong::REQUIRED_ARGUMENT],
  ["--gc", "--geneconv", GetoptLong::REQUIRED_ARGUMENT],
  ["--tree", GetoptLong::REQUIRED_ARGUMENT],
  ["--gc_length_min", GetoptLong::REQUIRED_ARGUMENT],
  ["--pvalue_bc_max", GetoptLong::REQUIRED_ARGUMENT],
  ["--mauve_size_min", GetoptLong::REQUIRED_ARGUMENT],
  ["--de_prefix", GetoptLong::NO_ARGUMENT],
  ["--length_min", GetoptLong::REQUIRED_ARGUMENT],
  ['-b', "--bootstrap", "--bootstrap_min", GetoptLong::REQUIRED_ARGUMENT],
  ["--no_bootstrap", GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--mauve'
      mauve_file = value
    when "-i", "--aln"
      infiles << value.split(',')
    when "--gc", "--geneconv"
      geneconv_files << value.split(',')
    when "--tree"
      tree_files << value.split(',')
    when "--gc_length_min"
      gc_length_min = value.to_f
    when "--pvalue_bc_max"
      pvalue_bc_max = value.to_f
    when '--mauve_size_min'
      mauve_size_min = value.to_i
    when '--de_prefix'
      is_de_prefix = true
    when '--length_min'
      length_min = value.to_i
    when '-b', '--bootstrap', '--bootstrap_min'
      bootstrap_min = value.to_f
    when '--no_bootstrap'
      is_no_bootstrap = true
  end
end


infiles.flatten!
geneconv_files.flatten!
tree_files.flatten!


####################################################################
mauve_rela = read_mauve(mauve_file)

if not geneconv_files.empty?
  geneconv_files.each do |geneconv_file|
    geneconv_info = read_geneconv_file(geneconv_file, gc_length_min, pvalue_bc_max)
    parse_geneconv_file(geneconv_file, geneconv_info, mauve_rela, mauve_size_min)
  end
elsif not infiles.empty?
  infiles.each do |infile|
    seq_objs = read_seq(infile, is_de_prefix)
    parse_infile(infile, seq_objs, length_min, mauve_rela, mauve_size_min)
  end
elsif not tree_files.empty?
  tree_files.each do |tree_file|
    parse_tree_file(tree_file, bootstrap_min, mauve_rela, mauve_size_min, is_no_bootstrap)
  end
else
  puts "geneconv_file, seqfile or tree_file has to be given! Exiting ......"
  exit
end


