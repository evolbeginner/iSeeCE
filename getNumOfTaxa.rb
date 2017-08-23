#! /bin/env ruby


require 'getoptlong'
require 'bio'


############################################################
indir = nil

infiles = Array.new
taxa = Hash.new


############################################################
def get_taxa(tree_file)
  taxa = Hash.new
  treeio = Bio::FlatFile.open(Bio::Newick, tree_file)
  tree = treeio.next_entry.tree
  tree.nodes.each do |node|
    next if node.name !~ /\w/
    taxon = node.name.split('|')[0]
    taxa[taxon] = ''
  end
  return(taxa)
end


############################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
  end
end


############################################################
Dir.foreach(indir) do |b|
  next if b =~ /^\./
  infile = File.join(indir, b)
  infiles << infile
end


infiles.each do |infile|
  taxa = taxa.merge!(get_taxa(infile))
end

p taxa.size


