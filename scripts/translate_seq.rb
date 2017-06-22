#! /bin/env ruby

# Objective:
#   1. trim sequences based on the number of nucleotides of the CDS sequences
#   2. remove stop codons in CDS sequences and translate CDS sequences
# Written by Sishuo Wang from the department of Botany, the University of British Columbia (sishuowang@hotmail.ca)


###########################################################################
require 'bio'
require 'getoptlong'
require 'Dir'

cds_file=nil
codon_table=1
frames=Array.new
is_remove_star=false
outfile = nil
outdir='./'
force=false
cds_file_basename=nil


###########################################################################
def translate_cds_seq(seq=nil, frames=[1], codon_table=1)
  aaSeqs = Array.new
  naSeq_obj = Bio::Sequence::NA.new(seq)
  frames.each do |frame|
    aaSeq=naSeq_obj.translate(frame,codon_table)
    aaSeqs.push aaSeq
    #puts aaSeq
  end
  return aaSeqs
end


def remove_star(aaSeq, naSeq)
  star_starts = Array.new
  aaSeq.scan('*') do |star|
    star_starts.push($`.size)
  end
  star_starts.each_with_index do |ele,index|
    star_starts[index]=ele-1 if index >= 1
  end
  star_starts.each do |start|
    naSeq.slice!(start*3,3)
  end
  aaSeq.gsub!(/\*/,"")
  return ([aaSeq,naSeq])
end


def usage
  script_basename=File.basename $0
  puts "Usage of #{script_basename}: ruby #{script_basename} arguments"
  print <<EOF
Arguments
Mandatory arguments:
-i|--cds_file     input file (CDS sequences)
Optional arguments:
--frame|--frames  the frame(s) used (separated by ',')
                  default: 1
--codon_table     the number of the corresponding codon table
                  only functional when '--translate' is specified
                  default: 1
-h|--help         see usage
EOF
  puts
  print <<EOF
Codon table
1. "Standard (Eukaryote)"
2. "Vertebrate Mitochondrial"
3. "Yeast Mitochondorial"
4. "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma"
5. "Invertebrate Mitochondrial"
6. "Ciliate Macronuclear and Dasycladacean"
9. "Echinoderm Mitochondrial"
10. "Euplotid Nuclear"
11. "Bacteria"
12. "Alternative Yeast Nuclear"
13. "Ascidian Mitochondrial"
14. "Flatworm Mitochondrial"
15. "Blepharisma Macronuclear"
16. "Chlorophycean Mitochondrial"
21. "Trematode Mitochondrial"
22. "Scenedesmus obliquus mitochondrial"
23. "Thraustochytrium Mitochondrial"
EOF
  puts "\nPlease write e-mail to sishuowang@hotmail.ca if you have any question and/or suggestion."
  puts
  exit
end


###########################################################################
opts=GetoptLong.new(
  ['-i', '--cds_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--codon_table', GetoptLong::REQUIRED_ARGUMENT],
  ['--frame','--frames', GetoptLong::REQUIRED_ARGUMENT],
  ['--remove_star','--is_remove_star', GetoptLong::NO_ARGUMENT],
  ['-o', '--out', '--outfile', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['-h','--help', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '-i', '--cds_file'
      cds_file=value
      cds_file_basename = File.basename(cds_file)
    when '--frame', '--frames'
      value.split(',').each do |frame| frames.push(frame.to_i) end
    when '--codon_table'
      codon_table=value.to_i
    when '--remove_star', '--is_remove_star'
      is_remove_star=true
    when '-o', '--out', '--outfile'
      outfile = value
    when '--outdir'
      outdir=value
    when '--force'
      force=true
    when '-h', '--help'
      usage
  end
end


if frames.empty?
  frames=[1]
end


if outdir != './'
  mkdir_with_force(outdir,force)
end




###########################################################################
outfile_name_core = nil
if cds_file_basename =~ /(.+) \.[^.]+$/x
  outfile_name_core = $1
end

pep_fhs = Hash.new
cds_fhs = Hash.new

if outfile.nil?
  frames.each_with_index do |frame,index|
    file_name = File.join(outdir,[outfile_name_core, '_pep', '_'+frames[index].to_s + '.fas'].join(''))
    pep_fh = File.open(file_name, 'w')
    pep_fhs[index] = pep_fh
    file_name = File.join(outdir,[outfile_name_core, '_cds', '_'+frames[index].to_s + '.fas'].join(''))
    cds_fh = File.open(file_name, 'w')
    cds_fhs[index] = cds_fh
  end
else
  pep_fh = File.open(outfile, 'w')
  pep_fhs[0] = pep_fh
end


fh = Bio::FlatFile.open(cds_file)
fh.each_entry do |f|
  naSeq=f.seq
  aaSeqs=translate_cds_seq(f.seq, frames, codon_table)
  aaSeqs.each_with_index do |aaSeq,index|
    pep_fh = pep_fhs[index]
    if is_remove_star
      cds_fh = cds_fhs[index]
      (aaSeq,naSeq) = remove_star(aaSeq, naSeq)
      if ! cds_fh.nil?
        cds_fh.puts ['>', f.definition].join('')
        cds_fh.puts naSeq
      end
    end
    pep_fh.puts ['>', f.definition].join('')
    pep_fh.puts aaSeq
  end
end

pep_fhs.each_value do |pep_fh| pep_fh.close end
cds_fhs.each_value do |cds_fh| cds_fh.close end if ! cds_fhs.empty?


