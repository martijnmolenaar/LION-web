### R-script, network visualization
### LION/web: www.lipidontology.com
  
### please cite: 
###   LION/web: a web-based ontology enrichment tool for lipidomic data analysis. 
###   Gigascience. 2019 Jun 1;8(6). pii: giz061. doi: 10.1093/gigascience/giz061. 
###   Martijn R. Molenaar, Aike Jeucken, Tsjerk A. Wassenaar, Chris H. A. van de Lest, Jos F. Brouwers, J. Bernd Helms. 
  
  
load(file = 'network.Rdata')
  
library(igraph)
  
plot(net, layout = layout,
     vertex.label.cex=.5, 
     vertex.label.degree = 1.5*pi, edge.arrow.size = .40,
     vertex.label.color='black')
  
library(visNetwork)
  
interactive_network
