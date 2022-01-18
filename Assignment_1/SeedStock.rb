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
    @gr_remaining = params.fetch(:gr_remaining, 00).to_i
    @@stocks << self
  end
  
  def SeedStock.get_seeds(pathS)#reading seed_stock_data.tsv and putting each variable in its corresponding attribute.
    stock_file = open(pathS,"r")
    File.readlines(stock_file)[1..-1].each do |line|
      ssi, mgi, lp, s, gr = line.strip.split("\t")
      SeedStock.new(seed_stock_id: ssi, mutant_gene_id: mgi, last_planted: lp, storage: s, gr_remaining: gr)
    end
    stock_file.close
  end
  
  def SeedStock.planting_seeds(planted_seeds)#planting the seeds.
    planting = planted_seeds.to_i
    puts "\n"
    #for each Seed Stock, check whether we can or cannot plant ALL the seeds we want (7). If not, return a warning message. Aslo update the date of planting the seeds.
    @@stocks.each {|e|
      if e.gr_remaining > planting
        e.gr_remaining = e.gr_remaining - planting
        e.last_planted = Time.new.strftime("%d/%m/%Y")
      else
        e.gr_remaining = 0
        e.last_planted = Time.new.strftime("%d/%m/%Y")
        puts "WARNING: We have run out of Seed Stock #{e.seed_stock_id} on date #{e.last_planted}"
      end}
  end
  
  def SeedStock.create_database(pathD)#creating an output database with updated data.
    updb = File.open(pathD,"w")
    updb << "Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining\n"#introducing header.
    @@stocks.each {|p|
      updb.write("#{p.seed_stock_id}\t#{p.mutant_gene_id}\t#{p.last_planted}\t#{p.storage}\t#{p.gr_remaining}\n")}#appending updated attributes.
  end
  
  def self.ret_stock#in order to retrieve information in other ruby scripts
        return @@stocks    
  end
end