require "./gffscript.rb"
require "rest-client"
require "bio"

gene_file = ARGV[0]

gene_names = File.open(gene_file)
gene_lines = gene_names.readlines()

g_l = "AT5g19120"

