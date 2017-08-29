#! /bin/env ruby


# alt name: findBestReciprocalBlast.rb


##########################################################
require 'getoptlong'


##########################################################
dir = File.dirname($0)
$examineFlanking = File.join(dir, 'examineFlanking.sh')


##########################################################
$flanking_dir = 'flanking2'
$sequence_dir = '../sequences'


##########################################################
infile = nil
outfile = nil
evalue_cutoff = 1e-10
n = 10
target = nil
locus = nil
summary_file = nil
num_positional_orthologs_min = 2
num_flanking_genes_with_positional_orthologs_min = 2


evalue_info = Hash.new{|h,k|h[k]={}}
best_info = Hash.new
best_reciprocal_info = Hash.new{|h,k|h[k]={}}
flanking_genes_with_positional_orthologs = Hash.new


##########################################################
def analyze_flanking_orthologs(target, locus, n, best_reciprocal_info, num_positional_orthologs_min, num_flanking_genes_with_positional_orthologs_min)
  flanking_genes_with_positional_orthologs = Hash.new
  n2 = n + 1

  upstream_genes=`grep #{target} -B300 #{$sequence_dir}/*gb | grep '/locus_tag=' | uniq | tail -#{n2} | head -n #{n}`.split("\n").map{|i|i.gsub(/.+\/locus_tag=/, "").gsub('"', '')}
  downstream_genes=`grep #{target} -A300 #{$sequence_dir}/*gb | grep '/locus_tag=' | uniq | head -#{n2} | tail -n #{n}`.split("\n").map{|i|i.gsub(/.+\/locus_tag=/, "").gsub('"', '')}

  [upstream_genes, downstream_genes].each do |arr|
    arr.each do |gene|
      positional_orthologs = Array.new
      
      best_reciprocal_info[gene].each_key do |ortholog|
        a = Array.new
        a = `grep #{ortholog} #{$flanking_dir}/#{locus}.locus | uniq`.split("\n").map{|i|i.gsub(/.+\/locus_tag=/, "").gsub('"', '')}
        positional_orthologs << a if not a.empty?
      end
      positional_orthologs.flatten!

      if positional_orthologs.size >= num_positional_orthologs_min
        flanking_genes_with_positional_orthologs[gene] = ''
      end
    end
  end

  if flanking_genes_with_positional_orthologs.size >= num_flanking_genes_with_positional_orthologs_min
    return flanking_genes_with_positional_orthologs.size
  else
    return 0
  end
end


def get_evalue_info(infile, evalue_cutoff)
  evalue_info = Hash.new{|h,k|h[k]={}}
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    query, subject, evalue = line_arr.values_at(0,1,-2)
    evalue = evalue.to_f
    orgns = query.split('|')[0], subject.split('|')[0]
    genes = query.split('|')[1], subject.split('|')[1]
    next if orgns[0] == orgns[1]
    next if evalue > evalue_cutoff
    #evalue_info[query][subject] = evalue
    evalue_info[genes[0]][genes[1]] = evalue
  end
  in_fh.close
  return(evalue_info)
end


def get_best_info(evalue_info)
  best_info = Hash.new
  evalue_info.each_pair do |gene, v|
    best_orthologs = Array.new
    v.each_pair do |ortholog, evalue|
      if best_orthologs.empty?
        best_orthologs << ortholog
      else
        if evalue < v[best_orthologs[0]]
          best_orthologs = [ortholog]
        elsif evalue == v[best_orthologs[0]]
          best_orthologs << ortholog
        end
      end
    end
    best_info[gene] = best_orthologs
  end
  return(best_info)
end


def get_best_reciprocal_info(best_info)
  best_reciprocal_info = Hash.new{|h,k|h[k]={}}
  best_info.each_pair do |gene, best_orthologs|
    best_orthologs.each do |best_ortholog|
      next if not best_info.include?(best_ortholog)
      if best_info[best_ortholog].include?(gene)
        best_reciprocal_info[gene][best_ortholog] = ''
      end
    end
  end
  return(best_reciprocal_info)
end


def read_summary_file(summary_file)
  pair_info = Hash.new{|h,k|h[k]={}}
  in_fh = File.open(summary_file, 'r')
  in_fh.each_line do |line|
    #FastTree/8_4/660.FastTree.tre	3	3	3	3	Ncom|AAW31_RS06815,Ncom|AAW31_RS10580:0.942	Neut|NEUT_RS11925,Neut|NEUT_RS11975:1.0	Nmul|NMUL_RS06280,Nmul|NMUL_RS06290:0.999
    line.chomp!
    line_arr = line.split("\t")
    n1, n2, n3, n4 = line_arr[1,4].map{|i|i.to_i}
    file_name = line_arr[0]
    b = File.basename(file_name)
    b =~ /^(\d+)\.FastTree\.tre$/
    locus = $1
    pair_infos = line_arr[5,line_arr.size-5]
    pair_infos.each do |i|
      genes = i.split(':')[0].split(',').map{|i|i.split('|')[1]}
      bootstrap = i.split(':')[1].to_i
      pair_info[locus][genes] = ''
    end
  end
  in_fh.close
  
  pair_info.each_pair do |locus, v|
    flanking_gene_file = File.join($flanking_dir, locus + '.locus')
    #next if File.exists?(flanking_gene_file)
    out_fh = File.open(flanking_gene_file, 'w')
    v.each_key do |genes|
      genes.each do |gene|
        `bash #{$examineFlanking} #{gene} 240 1 >> #{flanking_gene_file}`
      end
    end
    out_fh.close
  end

  return(pair_info)
end


def output_results(outfile, results)
  out_fh = (outfile.nil? or outfile == '-') ? STDOUT : File.open(outfile, 'w')
  results.each_with_index do |arr1, index|
    if index % 2 == 0
      if index+1 <= results.size-1
        arr2 = results[index+1]
        #if [arr1, arr2].select{|arr|arr[2] != false and arr[2]>=2}.size == 2
        locus = arr1[0]
        out_fh.puts [locus, arr1[1,2], arr2[1,2]].flatten.join("\t")
        #end
      end
    else
      next
    end
  end
  out_fh.close if ! outfile.nil?
end


##########################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['-o', GetoptLong::REQUIRED_ARGUMENT],
  ['-e', '--evalue', GetoptLong::REQUIRED_ARGUMENT],
  ['-n', GetoptLong::REQUIRED_ARGUMENT],
  ['--gene', GetoptLong::REQUIRED_ARGUMENT],
  ['--locus', GetoptLong::REQUIRED_ARGUMENT],
  ['--summary', GetoptLong::REQUIRED_ARGUMENT],
  ['--flanking_dir', GetoptLong::REQUIRED_ARGUMENT],
  ['--sequence_dir', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^-o$/
      outfile = value
    when /^(--evalue|-e)$/
      evalue_cutoff = value.to_f
    when /^-n$/
      n = value.to_i
    when /^--gene$/
      target = value
    when /^--locus$/
      locus = value
    when /^--summary$/
      summary_file = value
    when /^--flanking_dir$/
      $flanking_dir = value
    when /^--sequence_dir$/
      $sequence_dir = value
  end
end


##########################################################
unless summary_file.nil?
  pair_info = read_summary_file(summary_file)
end

evalue_info = get_evalue_info(infile, evalue_cutoff)

best_info = get_best_info(evalue_info)

best_reciprocal_info = get_best_reciprocal_info(best_info)

##########################################################
results = Array.new
if not target.nil?
  results << [locus, target, analyze_flanking_orthologs(target, locus, n, best_reciprocal_info, num_positional_orthologs_min, num_flanking_genes_with_positional_orthologs_min)]
elsif not pair_info.empty?
  pair_info.each_pair do |locus, gene_sh|
    gene_sh.each_key do |genes|
      genes.each do |gene|
        results << [locus, gene, analyze_flanking_orthologs(gene, locus, n, best_reciprocal_info, num_positional_orthologs_min, num_flanking_genes_with_positional_orthologs_min)]
      end
    end
  end
end


##########################################################
output_results(outfile, results)


