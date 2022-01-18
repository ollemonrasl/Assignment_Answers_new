require "bio"
require "rest-client"

gene_file = ARGV[0]
gene_names = File.open(gene_file)
gene_lines = gene_names.readlines()

gendict = Hash.new # HASH WHERE WE'LL STORE GENE INFORMATION

# FUNCTION TO CONNECT TO DATABASE
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

# FUNCTION TO RETRIEVE INFORMATION FROM THE DATABASE FOR EACH GENE OF THE FILE
def exam_sequences(gene_id)
  address = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_id.to_s.tr("[]","")}"
  res = fetch(address);
  if res
    embl = Bio::EMBL.new(res.body)
  else
    puts "the Web call failed - see STDERR for details..."
  end
  return gene_id,embl
end

# FUNCTION TO LOOK FOR CTTCTT REPEATS IN THE PROVIDED EXON SEQUENCE. IT RETURNS A STRING WITH PAIRS OF COORDINATES AS STRINGS "START::END"
def match_seq(seq)
    posit = Array.new # ARRAY WHERE ALL PAIRS OF COORDINATES WILL BE STORED AS SINGLE STRINGS
    search = Bio::Sequence::NA.new("cttctt")
    re = Regexp.new(search.to_re)
    match = seq.seq.match(re)
    if match
      positions = seq.enum_for(:scan, /(?=(cttctt))/i).map { Regexp.last_match.begin(0) } # RETURNS A STRING WITH STARTING POSITIONS OF EACH FOUND REPEAT
      positions.each {|st|
        en = st + 5 # THE END IS 5 POSITIONS FURTHER
        stenpairs = st.to_s + "::" + en.to_s # STORE IT AS A SINGLE STRING INSIDE THE ARRAY
        posit << stenpairs}
      return posit
    else
      return nil
    end
end

# FUNCTION TO CREATE NEW FEATURE WITH POSITIONS OF EACH REPEAT ON THE GENE AND ON THE WHOLE CHROMOSOME
def new_feature(st,en,chr,gst,gen,strand,note)
  found = Bio::Feature.new("match_in_exon","#{st}..#{en}")
  found.append(Bio::Feature::Qualifier.new("motif","cttctt"))
  found.append(Bio::Feature::Qualifier.new("start",st))
  found.append(Bio::Feature::Qualifier.new("ending",en))
  found.append(Bio::Feature::Qualifier.new("chr",chr))
  found.append(Bio::Feature::Qualifier.new("genestart",gst))
  found.append(Bio::Feature::Qualifier.new("geneend",gen))
  found.append(Bio::Feature::Qualifier.new("strand",strand))
  found.append(Bio::Feature::Qualifier.new("attrib",note))
  return found
end

# STORING RETRIEVED INFO INTO GENDICT HASH
puts "RETRIEVING INFO FROM CHOSEN DATABASE..."
gene_lines.each {|line|
  gl = line.strip()
  gen_id,embl = exam_sequences(gl)
  gendict[gen_id]=embl}

# TRANSFORMING GENE SEQUENCES INTO EXONS AND SCANNING EXON SEQUENCES LOOKING FOR CTTCTT REPEATS 
puts "SEARCHING FOR CTTCTT REPEATS..."
norepfilelist = Array.new 
  gendict.each {|id,embl|
    seq = embl.to_biosequence
    embl.features do |feature|
      next if feature.feature != "exon"
      forw = /\A\d+\.\.\d+/ # TO LOOK FOR EXONS IN THE FORWARD STRAND
      rev = /complement\(\d+\.\.\d+\)/ # TO LOOK FOR EXONS IN THE REVERSE STRAND
      if feature.locations.to_s.match(forw) # FOR EXONS IN THE FORWARD STRAND
          ini = feature.locations.to_s.split("..")[0].to_i # SAVE INITIAL COORDINATE
          fin = feature.locations.to_s.split("..")[1].to_i # SAVE ENDING COORDINATE
          exonseq = seq.subseq(ini,fin) # CUT SEQUENCE INTO EXON
          unless match_seq(exonseq).nil? # IF CTTCTT REPEATS ARE FOUND
            note = feature.assoc.values[0] 
            posit = match_seq(exonseq) # FUNCTION RETRIEVING LIST WITH START::END COORDINATES IN EXON OF REPETITIONS CTTCTT
            posit.each {|stenpair|
              start = stenpair.split("::")[0] # STORING START OF REP IN EXON
              ending = stenpair.split("::")[1] # STORING END OF REP IN REGION
              genst = embl.sv().split(":")[3].to_i + start.to_i # STORING START OF REP IN CHROM
              genen = embl.sv().split(":")[3].to_i + ending.to_i # STORING END OF REP IN CHROM
              newf = new_feature(start,ending,embl.sv().split(":")[2],genst,genen,"+",note) # ADDING FEATURE WITH DATA
              embl.features << newf # INCLUDING NEW FEATURE
              }
          else
            norepfilelist << id unless norepfilelist.include?(id)
          end
      elsif feature.locations.to_s.match(rev) # FOR EXONS IN THE REVERSE STRAND
          ini = feature.locations.to_s.split("(")[1].split("..")[0].to_i
          fin = feature.locations.to_s.split("(")[1].split("..")[1].to_i 
          seq = seq.reverse_complement
          exonseq = seq.subseq(ini,fin)
          unless match_seq(exonseq).nil?
            note = feature.assoc.values[0]
            posit = match_seq(exonseq)
            posit.each {|stenpair|
              start = stenpair.split("::")[0]
              ending = stenpair.split("::")[1]
              genst = embl.sv().split(":")[3].to_i + start.to_i
              genen = embl.sv().split(":")[3].to_i + ending.to_i
              newf = new_feature(start,ending,embl.sv().split(":")[2],genst,genen,"-",note)
              embl.features << newf}
          else
            norepfilelist << id unless norepfilelist.include?(id)
          end
      end
    end}

# TO WRITE THE GFF3 FILE
puts "WRITING .GFF3 FILE..."
gff3_file = File.open("gene_report.gff3","a+")
gff3_file.puts "##gff-version 3"
source = "BioRuby"
type = "cttctt_repeat"
score = "."
phase = "."
gendict.each {|seq,embl|
  embl.features {|feature|
    next unless feature.feature == "match_in_exon"
    seqid = seq
    start = feature.assoc["start"]
    ending = feature.assoc["ending"]
    strand = feature.assoc["strand"]
    attrib = feature.assoc["attrib"]
    gff3_file.puts "#{seqid}\t#{source}\t#{type}\t#{start}\t#{ending}\t#{score}\t#{strand}\t#{phase}\t#{attrib}"}}

# TO WRITE TXT FILE WITH GENES WITH NO CTTCTT MATCHES
puts "WIRITING .TXT FILE WITH GENES THAT DO NOT CONTAIN CTTCTT REPEATS..."
norepfile = File.open("no_rep_file.txt","a+")
norepfile.puts "Genes of the list that DO NOT contain CTTCTT repeats:"
norepfilelist.each {|gene|
  norepfile.puts gene}

# TO WRITE THE GFF3 WITH CHROMOSOMIC COORDINATES 
puts "WRITTING .GFF3 FILE WITH CHROMOSOMIC COORDINATES..."
gff3_file = File.open("chrom_report.gff3","a+")
gff3_file.puts "##gff-version 3"
source = "BioRuby"
type = "cttctt_repeat"
score = "."
phase = "."
gendict.each {|seq,embl|
  embl.features {|feature|
    next unless feature.feature == "match_in_exon"
    seqid = feature.assoc["chr"]
    genestart = feature.assoc["genestart"].to_i
    geneend = feature.assoc["geneend"].to_i
    strand = feature.assoc["strand"]
    attrib = feature.assoc["attrib"]
    gff3_file.puts "#{seqid}\t#{source}\t#{type}\t#{genestart}\t#{geneend}\t#{score}\t#{strand}\t#{phase}\t#{attrib}"}}