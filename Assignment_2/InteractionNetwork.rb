require "rest-client"
require "./Gene.rb"
require "./Network.rb"

file_names = ARGV[0]

gene_file = File.open(file_names)
gene_lines = gene_file.readlines()

gene_lines.each {|gene|
  puts gene
  gene_id = gene.split("\n")[0]
  Gene.gene_info(gene_id)
  }
  
puts "INFO RETRIEVED! STARTING INTERACTOME GENERATION"
Network.interactome