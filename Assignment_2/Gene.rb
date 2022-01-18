require "rest-client"
require "json"

class Gene
  attr_accessor :gene_id
  attr_accessor :network
  attr_accessor :kegg
  attr_accessor :go
  @@netw = []
  
  def initialize(params = {})
    @gene_id = params.fetch(:gene_id,"AT0G00000")
    @network = params.fetch(:network,Array.new)
    @kegg = params.fetch(:kegg,Hash.new)
    @go = params.fetch(:go,Hash.new)
    @@netw << self
  end

#this function retrieves the information for the given gene id from the three different repositories
#initially the parameters are stored in objects that at the end of the function are assigned to the attributes
  def self.gene_info(geneid)
    def self.fetch(url, headers = {accept: "*/*"}, user = "", pass="")
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
        return response
      rescue RestClient::Exception => e
        $stderr.puts e.inspect
        response = false
        return response
      rescue Exception => e
        $stderr.puts e.inspect
        response = false
        return response
    end
    puts "RETRIEVING INFO..."
    address = "http://bar.utoronto.ca:9090/psicquic/webservices/current/search/query/#{geneid}?format=tab25"
    res = fetch(address)
    if res
      net = Array.new
      body = res.body
      bodylines = body.split("\n")
      for line in bodylines do
        threshold = 0.6
        miscore = line.split("\t")[14].split(":")[1]
        if miscore.to_f >= threshold
          geneA = line.split("\t")[2].split(":")[1] if line.split("\t")[2].split(":")[1].match(/[Aa][Tt]\w[Gg]\d\d\d\d\d/)
          geneB = line.split("\t")[3].split(":")[1] if line.split("\t")[3].split(":")[1].match(/[Aa][Tt]\w[Gg]\d\d\d\d\d/)
          if geneA.to_s.downcase == geneid.downcase and not net.include?(geneB)
            net.push(geneB)
          elsif geneB.to_s.downcase == geneid.downcase and not net.include?(geneA)
            net.push(geneA)
          end
        end
      end  
      keggdict = {}
      kegg = "http://togows.org/entry/kegg-genes/ath:#{geneid}/pathways.json"
      reskegg = RestClient::Request.execute(method: :get, url: kegg)
      datakegg = JSON.parse(reskegg.body)
      if datakegg[0]
       datakegg[0].each do |keggid, path_name|  #Extracting kegg info
         keggdict[keggid] = path_name
       end
      end
      godict = {}
      go = "http://togows.org/entry/ebi-uniprot/#{geneid}/dr.json"
      resgo = RestClient::Request.execute(method: :get, url: go)
      datago = JSON.parse(resgo.body)
      if datago[0]["GO"]
       datago[0]["GO"].each do |annotation|  #Extracting go info
        if annotation[1].match(/^P:/)#Extracting only biological processes
          godict[annotation[0]] = annotation[1]
        end
       end
      end
      Gene.new(gene_id: geneid, network: net, kegg: keggdict, go: godict)
    end
  end
  
  def self.get_net
    return @@netw
  end
  
end