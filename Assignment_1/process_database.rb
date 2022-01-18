require "./SeedStock.rb"
require "./Gene.rb"
require "./HybridCross.rb"

Gene.get_genes(ARGV[0])
SeedStock.get_seeds(ARGV[1])
HybridCross.get_crosses(ARGV[2])

SeedStock.planting_seeds(7)
SeedStock.create_database(ARGV[3])
HybridCross.chi
Gene.linkages