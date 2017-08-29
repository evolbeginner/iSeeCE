#! /bin/env ruby


require "getoptlong"
require "bio"
require "Dir"


##########################################################
indirs = Array.new
suffix = nil
outdir = nil
is_force = false

seq_titles = Hash.new{|h,k|h[k]=0}
source_rela = Hash.new


##########################################################
opts = GetoptLong.new(
  ["--indir", GetoptLong::REQUIRED_ARGUMENT],
  ["--suffix", GetoptLong::REQUIRED_ARGUMENT],
  ["--outdir", GetoptLong::REQUIRED_ARGUMENT],
  ["--force", GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "--indir"
      indirs << value.split(',')
    when "--suffix"
      suffix = value
    when "--outdir"
      outdir = value
    when "--force"
      is_force = true
  end
end

indirs.flatten!

mkdir_with_force(outdir, is_force)


##########################################################
indirs.each do |indir|
  Dir.foreach(indir) do |file|
    next if file =~ /^\./
    next if file !~ /#{suffix}$/ if not suffix.nil?
    temp_seq_titles = Array.new
    full_name = File.join([indir, file])
    File.open(full_name,'r').each_line do |line|
      line.chomp!
      next if line !~ /^>(.+)/
      temp_seq_titles << $1
    end
    seq_titles[temp_seq_titles.sort!] += 1
    source_rela[temp_seq_titles] = full_name
  end
end


##########################################################
counter = 0
seq_titles.each_pair do |seq_title_arr, num_of_share|
  counter += 1
  source = source_rela[seq_title_arr]
  `cp #{source} #{outdir}/#{counter}.aln`
end


