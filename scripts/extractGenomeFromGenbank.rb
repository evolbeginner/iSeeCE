#! /bin/env ruby


require "getoptlong"
require "bio"


####################################################################
infile = nil
is_no_plasmid = false


####################################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--no_plasmid", GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i"
      infile = value
    when "--no_plasmid"
      is_no_plasmid = true
  end
end


####################################################################
#in_fh = File.open(infile, "r")
#first_line = in_fh.readline
#locus_name = first_line.split(/\s+/)[1]
#in_fh.close

Bio::FlatFile.open(infile).each_entry do |f|
  next if f.definition =~ /plasmid/ if is_no_plasmid
  next if f.entry_id.nil?
  puts '>'+f.entry_id
  puts f.seq
end


