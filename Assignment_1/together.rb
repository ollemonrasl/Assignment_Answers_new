class Gene
    attr_accessor :gene_id
    attr_accessor :gene_name
    attr_accessor :mutant_phen
    attr_accessor :linked_to
    @@genes = []
    
    def initialize(params = {})
      @gene_id = params.fetch(:gene_id, "AT0G00000")
 #     abort "Gene identificator #{self.gene_id} not valid" unless @gene_id.match(/A[Tt]\d[Gg]\d\d\d\d\d/)
      @gene_name = params.fetch(:gene_name, "XXX")
      @mutant_phen = params.fetch(:mutant_phen, "xxxxxxxxx")
      @linked_to = params.fetch(:linked_to, "XXXX")
      @@genes << self
    end
    
    def Gene.get_genes(pathG)
      File.readlines(pathG)[1..-1].each do |line|
        gi, gn, mp = line.strip.split("\t")
        Gene.new(gene_id => gi, :gene_name => gn, :mutant_phen => mp)
      end
      if @linked_genes.key?(gi)
        @linked_to = @linked_genes[gi]
      end
      
    end
end

class SeedStock
  attr_accessor :seed_stock_id
  attr_accessor :mutant_gene_id
  attr_accessor :last_planted
  attr_accessor :storage
  attr_accessor :gr_remaining
  @@stocks = []
  
  def initialize(params = {})
    @seed_stock_id = params.fetch(:seed_stock_id, "XOO")
    @mutant_gene_id = params.fetch(:mutant_gene_id, "AT0G00000")
    @last_planted = params.fetch(:last_planted, "DD/MM/YYYY")
    @storage = params.fetch(:storage, "cama00")
    @gr_remaining = params.fetch(:gr_remaining, 00)
    @@stocks << self
  end
  
  def SeedStock.get_seeds(pathS)#reading seed_stock_data.tsv and putting each variable in its corresponding attribute
    File.readlines(pathS)[1..-1].each do |line|
      ssi, mgi, lp, s, gr = line.strip.split("\t")
      SeedStock.new(:seed_stock_id => ssi, :mutant_gene_id => mgi, :last_planted => lp, :storage => s, :gr_remaining => gr)
    end
  end
  
  def SeedStock.planting_seeds(planted_seeds)#if we can, plant the seeds
    if @gr_remaining > planted_seeds
        @gr_remaining-=planted_seeds
        @last_planted = Time.new.strftime("%d/%m/%Y")
    else
      @gr_remaining = 0
      puts "WARNING: We have run out of Seed Stock #{@seed_stock_id}"
    end
  end
  
  def SeedStock.create_database(pathD)
    updb = File.open(pathD)
    updb << "Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining\n"
    @@stocks.each {|p|
      updb.write("#{p.seed_stock_id}\t#{p.mutant_gene_id}\t#{p.last_planted}\t#{p.storage}\t#{p.gr_remaining}")}
  end
  
end

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
      @f2_wild = params.fetch(:f2_wild, "000")
      @f2_p1 = params.fetch(:f2_p1, "000")
      @f2_p2 = params.fetch(:f2_p2, "000")
      @f2_p1p2 = params.fetch(:f2_p1p2, "000")
      @linked_genes = Hash.new
      @@crosses << self
    end
    
    def HybridCross.get_crosses(pathC)
      File.readlines(pathC)[1..-1].each do |line|
        p1, p2, f2w, f2p1, f2p2, f2p1p2 = line.strip.split("\t")
        HybridCross.new(:p1 => p1, :p2 => p2, :f2_wild => f2w, :f2_p1 => f2p1, :f2_p2 => f2p2, :f2_p1p2 => f2p1p2)
      end
    end
    
    def HybridCross.store_linked_genes(g1,g2)
      @linked_genes[g1] = g2
      @linked_genes[g2] = g1
    end
      
    def HybridCross.chi
      @@crosses.each {|p|
        total = p.f2_wild + p.f2_p1 + p.f2_p2 + p.f2_p1p2
        w = total*9/16
        p1 = total*3/16
        p2 = total*3/16
        p1p2 = total*1/16
        chisq = ((p.f2_wild - w)**2/w) + ((p.f2_p1 - p1)**2/p1) + ((p.f2_p2 - p2)**2/p2) + ((p.f2_p1p2 - p1p2)**2/p1p2)
        if chisq > 7.81 #given that there are 3 degrees of freedom and 0,05 cutoff this represents both genes are linked
          puts "#{p1} is genetically linked to #{p2} with chisquarescore #{chisq}"
          store_linked_genes(p1,p2)
        end}
    en    
end

Gene.get_genes(ARGV[0])
SeedStock.get_seeds(ARGV[1])
HybridCross.get_crosses(ARGV[2])

SeedStock.planting_seeds(7)
SeedStock.create_database(ARGV[3])
HybridCross.chi