require "bio"
require "rest-client"
##require "net/http"

#gene_file = ARGV[0]

#gene_names = File.open(gene_file)
#gene_lines = gene_names.readlines()

g_l = "AT1g05207"
g_l2 = "AT5g19120"
g_l3 = "AT4g27030"
g_l4 = "AT5g54270"
g_l5 = "AT4g05180"

triallist = [g_l,g_l2,g_l3,g_l4,g_l5]

def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
  
  rescue RestClient::ExceptionWithResponse => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
end

def exam_sequences(gene_id)
  @genes = Hash.new()
#for gene_id_n in gene_lines
  #gene_id = gene_id_n.split("\n")
  address = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_id.to_s.tr("[]","")}"
  res = fetch(address);
  if res
    embl = Bio::EMBL.new(res.body)
    @genes[gene_id]=embl
  else
    puts "the Web call failed - see STDERR for details..."
  end
  return @genes
end

def match_seq(seq)
    search = Bio::Sequence::NA.new("cttctt")
    re = Regexp.new(search.to_re)
    match = seq.seq.match(re)
    if match
        start = seq =~ /cttctt/
        ending = start + 5
        #puts "cttctt starts at #{start} and ends at #{ending}"
        return start,ending
    else
      puts "----NO MATCH FOR #{seq}"
      return nil
    end
end

def new_feature(start,ending,strand)
  
end

def generate_gff3(seq,source="BioRuby",type="direct_repeat",score=".",phase=".")
  #File.open("gff3_gene_report.gff3","w+") do |line|
   # line.puts"##gff-version 3"
    @embl.features.each {|cod,seq|
      puts cod
      puts seq
      #next unless feature.feature == "match_in_exon"
      puts feature.locations.first
      puts feature.assoc["strand"]
      puts cod
      puts seq
      }
  #end
end


triallist.each {|gl|
  exam_sequences(gl).each {|id,embl|
    puts "-----------"
    #seq = Bio::Sequence::NA.new("cttttaaagagagagcttagatcgtatagc")
    seq = embl.to_biosequence
    embl.features do |feature|
      #puts feature
      puts feature.feature
      #puts feature.class
      puts feature.feature.class
      puts feature.position
      puts feature.assoc
      next if feature.feature != "exon"
      #puts feature.locations
      forw = /\A\d+\.\.\d+/ # To look for exons in the forward strand
      rev = /complement\(\d+\.\.\d+\)/ # To look for exons in the reverse strand
      if feature.locations.to_s.match(forw) # For exons in the forward strand
          ini = feature.locations.to_s.split("..")[0].to_i # Save initial coordinate
          fin = feature.locations.to_s.split("..")[1].to_i # Save ending coordinate
          #puts ini
          #puts fin
          #puts seq
          exonseq = seq.subseq(ini,fin)
          #puts seq.subseq(ini,fin)
          #puts "match with forw"
          puts embl.features.length
          unless match_seq(exonseq).nil?
            start,ending = match_seq(exonseq) # Look for the "cttctt" sequence in the exon of the forward strand
            found = Bio::Feature.new("match_in_exon","#{start}..#{ending}")
            found.append(Bio::Feature::Qualifier.new("motif","cttctt"))
            #found.append(Bio::Feature::Qualifier.new("function","insertion site"))
            found.append(Bio::Feature::Qualifier.new("strand","+"))
            embl.features << found
            puts embl.features.length
          end
      elsif feature.locations.to_s.match(rev) # For exons in the reverse strand
          ini = feature.locations.to_s.split("(")[1].split("..")[0].to_i 
          fin = feature.locations.to_s.split("(")[1].split("..")[1].to_i
          #puts ini
          #puts fin
          #puts seq
          #puts seq.length
          seq = seq.reverse_complement
          #puts seq
          #puts seq.length
          #puts seq.subseq(ini,fin)
          exonseq = seq.subseq(ini,fin)
          puts embl.features.length
          unless match_seq(exonseq).nil?  
            start,ending = match_seq(exonseq,ini,fin,st)
            found = Bio::Feature.new("match_in_exon","Complement:#{start}..#{ending}")
            found.append(Bio::Feature::Qualifier.new("motif","cttctt"))
            #found.append(Bio::Feature::Qualifier.new("function","insertion site"))
            found.append(Bio::Feature::Qualifier.new("strand","-"))
            embl.features << found
            puts embl.features.length
            #puts "match with rev"
          end
      end
    end}
  }
puts "................."

triallist.each {|g_l|  
  exam_sequences(g_l).each {|seq,embl|
    embl.features do |feature|
      puts feature.feature
    end
    #generate_gff3(embl)
    }}
