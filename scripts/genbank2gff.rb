#! /bin/env ruby


require 'bio'
require 'getoptlong'


###########################################################
seq_file=nil
features_included = Hash.new
products_included = Hash.new
translation_output=false


###########################################################
def read_seq_file(seq_file=nil,features_included={},qualifiers_included={},products_included={},translation_output=nil)
  if translation_output then
    translation_output_fh = File.open(translation_output,'w')
  end

  fh = Bio::FlatFile.open(seq_file)
  fh.each_entry do |f|
    f.features.each do |feature_obj|
      qualifiers=Hash.new 
      qualifiers_content=nil
      parent_feature=feature_obj.feature

      (next if not features_included[parent_feature]) if ! features_included.empty?
      strand = feature_obj.position=~/complement/ ? '-' : '+'
      parent_start, parent_stop = $1, $2 if feature_obj.position =~ /(<?\d+).+((?<=\.\.)>?\d+)/

      counter=0
      feature_obj.position.scan(/(<?\d+)\.\.(>?\d+)/) do
        counter+=1
        #feature = feature_obj.feature == 'CDS' ? 'exon' : parent_feature
        feature = parent_feature
        start,stop = $1,$2
        if counter <= 1 then
          feature_obj.qualifiers.each do |i|
            (next if not qualifiers_included[feature].include? i.qualifier) if ! qualifiers_included.empty?
            if not products_included.empty?
              next if not products_included.include?(i.value)
            end
            #next if not products_included.include?(i.qualifier.product)
            if i.qualifier == 'translation' then
              if translation_output_fh then
                translation_output_fh.puts(i.value)
              end
            else
              qualifier = qualifiers_included[feature][i.qualifier]
              qualifiers[qualifier]=i.value
            end
          end
        end
        qualifiers_content = qualifiers.keys.map{|k|[k,qualifiers[k]].join('=')}.join(';')

        if features_included.include?(feature_obj.feature)
          puts [f.entry_id,'genbank',feature,start,stop,1,strand,'.',qualifiers_content].join("\t")
        end
      end

      #if parent_feature == 'CDS' then puts parent_start or exit if parent_stop =~ />/ end
      #next if parent_feature == 'CDS'
      #feature, start, stop = parent_feature, parent_start, parent_stop
      #puts [f.entry_id,'genbank',feature,start,stop,1,strand,'.',qualifiers_content].join("\t")
    end
    #[:locus, :entry_id, :circular, :division, :date, :strand, :natype, :each_cds, :each_gene, :basecount, :seq, :naseq, :nalen, :seq_len, :date_modified, :classification, :strandedness, :to_biosequence, :definition, :accessions, :versions, :acc_version, :accession, :version, :gi, :nid, :keywords, :segment, :source, :common_name, :vernacular_name, :organism, :taxonomy, :references, :comment, :features, :origin, :tags, :exists?, :get, :fetch]
  end
  translation_output_fh.close if translation_output_fh
end


def build_qualifiers_included(hash={})
  qualifiers_included=Hash.new{|h,k|h[k]=Hash.new}
  if not hash.empty? then
    hash.each_pair do |k,v|
      v.each do |i|
        value = i[1].nil? ? i[0] : i[1]
        qualifiers_included[k][i[0]]=value
      end
    end
  end
  return qualifiers_included
end


###########################################################
opts = GetoptLong.new(
  ['--genbank','--seq',GetoptLong::REQUIRED_ARGUMENT],
  ['--feature',GetoptLong::REQUIRED_ARGUMENT],
  ['--product',GetoptLong::REQUIRED_ARGUMENT],
  ['--output_translation',GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--genbank', '--seq'
      seq_file=value
    when '--feature'
      value.split(',').each{|i|features_included[i]=1}
    when '--product'
      value.split(',').each{|i|products_included[i]=1}
    when '--translation_output'
      translation_output=value
  end
end


###########################################################
qualifiers_included = build_qualifiers_included(
                          {
                           'CDS'  =>  [%w[locus_tag ID],%w[protein_id],%w[translation]],
                           'gene' =>  [%w[locus_tag ID],%w[protein_id]],
                           'exon' =>  [%w[locus_tag Parent],%w[protein_id]],
                           'intron'=> [%w[locus_tag Parent],%w[protein_id]],
                           'rRNA'=> [%w[product],%w[protein_id]],
                          }
)

read_seq_file(seq_file,features_included,qualifiers_included,products_included,translation_output)


