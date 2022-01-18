require "./HybridCross.rb"
require "./SeedStock.rb"

class Gene
    
    attr_accessor :gene_id
    attr_accessor :gene_name
    attr_accessor :mutant_phen
    attr_accessor :linked_to
    @@genes = []
    
    def initialize(params = {})
      #testing gene_id format
      id = params.fetch(:gene_id, "AT0G00000")
      match_id = Regexp.new(/A[Tt]\d[Gg]\d\d\d\d\d/)
      unless match_id.match(id)
        puts "Gene identificator #{gene_id} not valid"
        @gene_id = "AT0G00000"
      else
        @gene_id = id
      end
      @gene_name = params.fetch(:gene_name, "XXX")
      @mutant_phen = params.fetch(:mutant_phen, "xxxxxxxxx")
      @linked_to = params.fetch(:linked_to, "XXXX")
      @@genes << self
    end
    
    def Gene.get_genes(pathG)#reading gene_information.tsv and putting each variable in its corresponding attribute.
        gene_file = File.open(pathG, "r")
        File.readlines(gene_file).drop(1).each do |line|
            gi, gn, mp = line.strip.split("\t")
            Gene.new(gene_id: gi, gene_name: gn, mutant_phen: mp)
        end
        gene_file.close
    end
    
    def Gene.linkages#retrieving data from Hybrid.Cross.rb in order to include the gene of linkage if it is the case. 
        puts "\nFinal Report:\n"
        HybridCross.ret_crosses.each do |cross|#iterating over the HybridCross class attributes to retrieve linked genes.
            cross.linked_genes.each do |k,v|#extracting linked genes from hash. 
                @@genes.each {|g|
                    if g.gene_name == k
                        g.linked_to = v#updating attribute with liked gene
                        puts "#{k} is linked to #{v}"
                    end}
            end
        end
    end
    
    def self.ret_genes#in order to retrieve information in other ruby scripts
        return @@genes
    end
   
end