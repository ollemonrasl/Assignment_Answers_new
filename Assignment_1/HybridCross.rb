require "./SeedStock.rb"

class HybridCross
    attr_accessor :p1
    attr_accessor :p2
    attr_accessor :f2_wild
    attr_accessor :f2_p1
    attr_accessor :f2_p2
    attr_accessor :f2_p1p2
    attr_accessor :linked_genes
    @@crosses = []
    
    def initialize(params = {})
      @p1 = params.fetch(:p1, "X0000")
      @p2 = params.fetch(:p2, "X0000")
      @f2_wild = params.fetch(:f2_wild, 000).to_f
      @f2_p1 = params.fetch(:f2_p1, 000).to_f
      @f2_p2 = params.fetch(:f2_p2, 000).to_f
      @f2_p1p2 = params.fetch(:f2_p1p2, 000).to_f
      @linked_genes = Hash.new
      @@crosses << self
    end
    
    def HybridCross.get_crosses(pathC)#reading cross_data.tsv and putting each variable in its corresponding attribute.
        cross_file = open(pathC, "r")
        File.readlines(cross_file).drop(1).each do |line|
            p1, p2, f2w, f2p1, f2p2, f2p1p2 = line.strip.split("\t")
            HybridCross.new(p1: p1, p2: p2, f2_wild: f2w, f2_p1: f2p1, f2_p2: f2p2, f2_p1p2: f2p1p2)
        end
        cross_file.close
    end
    
    def store_linked_genes=(lkg)#storing linked genes into hash to later retrieve them in the Gene.rb script.
        @linked_genes[lkg[0]] = lkg[1]
        @linked_genes[lkg[1]] = lkg[0]
    end
    
    
    def HybridCross.chi#executing ChiSquare test to check linkage between genes
      @@crosses.each {|p|
        total = p.f2_wild + p.f2_p1 + p.f2_p2 + p.f2_p1p2
        w = total*9/16
        par1 = total*3/16
        par2 = total*3/16
        p1p2 = total*1/16
        chisq = ((p.f2_wild - w)**2/w) + ((p.f2_p1 - par1)**2/par1) + ((p.f2_p2 - par2)**2/par2) + ((p.f2_p1p2 - p1p2)**2/p1p2)
        if chisq > 7.81 #given that there are 3 degrees of freedom and 0,05 cutoff this represents both genes are linked
          lkg = Array.new#creating array to store gene_ids (from SeedStock class) in order to later be transformed into gene_names (from Gene class)
          SeedStock.ret_stock.each do |g|
              if p.p1 == g.seed_stock_id
                lkg.push(g.mutant_gene_id)
              end
              if p.p2 == g.seed_stock_id
                lkg.push(g.mutant_gene_id)
              end
          end
        flkg = Array.new#creating array to store gene_names (from Gene) corresponding to gene_ids(from SeedStock class)
        Gene.ret_genes.each {|g|
            if lkg[0] == g.gene_id
                flkg.push(g.gene_name)
            end
            if lkg[1] == g.gene_id
                flkg.push(g.gene_name)
            end}
        p.store_linked_genes = flkg#storing gene_names into hash
        puts "\nRecording: #{flkg[0]} is genetically linked to #{flkg[1]} with ChiSquareValue #{chisq}"
        end}
    end
    
  
    def self.ret_crosses#in order to retrieve information in other ruby scripts
        return @@crosses    
    end
    
end