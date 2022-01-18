require "bio"
require "stringio"


def ret_type(ff)
  seq_type = Bio::Sequence.auto(ff.next_entry.to_s).guess
  if seq_type.to_s.split("::")[2] == "NA"
    type = "nucl"
  elsif seq_type.to_s.split("::")[2] == "AA"
    type = "prot"
  end
  return type
end

def create_db(ff,f)
  type = ret_type(ff)
  out = f.to_s.split(".")[0]
  puts "\nCreating local database..."
  system("makeblastdb -in #{f} -dbtype #{type} -out #{out}")
end

def bl_type(a,b)
  if ret_type(a) == "prot" and ret_type(b) == "prot"
    dbtype = "blastp"
  elsif ret_type(a) == "nucl" and ret_type(b) == "nucl"
    dbtype = "blastn"
  elsif ret_type(a) == "nucl" and ret_type(b) == "prot"
    dbtype = "tblastn"
  elsif ret_type(a) == "prot" and ret_type(b) == "nucl"
    dbtype = "blastx"
  end
  return dbtype
end

arabid = ARGV[0]
spombe = ARGV[1]

arabidf = Bio::FlatFile.auto(arabid)
spombef = Bio::FlatFile.auto(spombe)

create_db(arabidf,arabid)
create_db(spombef,spombe)

# from https://doi.org/10.1371/journal.pone.0101850: evalue = 1e-6 ("-e 1e-6")
# from https://doi.org/10.1093/bioinformatics/btm585: filtering = "m S" ("-F 'm S'")

def ortologues(f1,f2,ff1,ff2)
  eval = "-e 1e-6"
  filt = "-F 'm S'"
  besthit = Hash.new
  orth = Hash.new
  ldb1 = f1.to_s.split(".")[0]
  ldb2 = f2.to_s.split(".")[0]
  dbtype1 = bl_type(ff1,ff2)
  dbtype2 = bl_type(ff2,ff1)
  ff2.each_entry do |q1|
    id1 = q1.entry_id
    factory1 = Bio::Blast.local("#{dbtype1}","#{ldb1}","'#{filt} #{eval}'")
    results1 = factory1.query(q1)
    if results1.hits[0]
      val = results1.hits[0].definition.split("|")[0].strip
      besthit[id1] = val
    end
  puts "Besthit dict fin"
  end
  ff1.each_entry do |q2|
    id2 = q2.entry_id
    next unless besthit.has_value?(id2)
      factory2 = Bio::Blast.local("#{dbtype2}","#{ldb2}","'#{filt} #{eval}'")
      results2 = factory2.query(q2)
      if results2.hits[0]
       reciprocal = results2.hits[0].definition.split("|")[0].strip
       if besthit[reciprocal] == id2
         orth[reciprocal] = id2
       end
      end
    end
  puts "Reciprocal hits fin"
  return orth
end
  
pairs = ortologues(arabid,spombe,arabidf,spombef)

File.open("orthologues.txt","w+") do |line|
  line.puts "List of putative orthologue pairs between Arabidopsis thaliana and Schizosaccharomyces pombe:\n"
  for v,k in pairs
    line.puts "#{v} and #{k}"
  end
end