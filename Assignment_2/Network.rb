require "rest-client"
require "./Gene.rb"

class Network
    attr_accessor :net_id
    attr_accessor :net_node_id
    attr_accessor :interactd
    @@genenet = []
    
    def initialize(params = {})
        @net_id = params.fetch(:net_id, 0)
        @net_node_id = params.fetch(:net_node_id, "AT0G00000")
        @interact = params.fetch(:interact, Array.new)
        @@genenet << self
    end
#this function retrieves the information from the gene object and creates different dictionaries to store the information
#in the interaction dictionary, it looks for genes that act as intermediate interactors between genes of the list
#by iterating over the elements of the dictionary
   def self.interactome
       interactiondict = {}
       Gene.get_net.each {|inter|
           interactiondict[inter.gene_id]=inter.network unless inter.gene_id=="AT0G00000"}
       interactiondict.each {|queryid,queryint|
                interactiondict.each {|id2,inter2|
                    unless queryid == id2 or queryint.empty? == true or inter2.empty? == true
                        indirectinteract = []
                        queryint.each {|prot1|
                            inter2.each {|prot2|
                                if prot1 == prot2 and indirectinteract.include?(id2) == false
                                    indirectinteract.push(id2)
                                end
                                }}
                        queryint.push(indirectinteract) if indirectinteract.empty? == false
                    end
                    }
                }
       keggdict = {}
       godict = {}
       Gene.get_net.each {|elem|
           keggdict[elem.gene_id]=elem.kegg
           godict[elem.gene_id]=elem.go}
#this part of the function writes the report by searching for the id in the created dictionaries
       reportfile = File.open("./network_report.txt","a+")
       reportfile.puts "INTERACTION NETWORK REPORT\n"
       reportfile.puts "This file contains, for each gene of the file:\n"
       reportfile.puts "1. A list of the genes with direct interaction\n"
       reportfile.puts "2. A list of the genes OF THE FILE with indirect interactions (1 interactor inbetween)\n"
       reportfile.puts "3. The KEGG and GO annotations for the main gene\n"
       Gene.get_net.each {|g|
            gene = g.gene_id
            reportfile.puts "\n\n\nFOR GENE #{gene}:"
            reportfile.puts "\n1. Direct interactions:"
            interactiondict.fetch(gene).each {|inter|
                reportfile.puts inter unless inter.class == Array}
            reportfile.puts "\n2. Indirect interactions (GENES FROM LIST):"
            reportfile.puts interactiondict.fetch(gene).last if interactiondict.fetch(gene).last.kind_of?(Array) == true
            reportfile.puts "\n3.1. KEGG annotations:"
            keggdict.fetch(gene).each {|k,v|
                reportfile.puts "KEGG id: #{k}, KEGG pathway: #{v}"}
            reportfile.puts "\n3.2. GO annotations:"
            godict.fetch(gene).each {|k2,v2|
                reportfile.puts "GO id: #{k2}, #{v2}"}
       }
   end
    
end