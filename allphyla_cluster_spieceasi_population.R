setwd("/project/noujdine_61/mgosborn/Gametophytes/NetworkTopology/allphyla/PopulationComparison")
library(phyloseq)
library(dplyr)
library(tidyr)
library(tidyverse) 
library(microbiome) 
library(vegan)
library(SpiecEasi)
library(igraph)

## This script runs through network and zeta analysis. 

## Remember: Grab abundance table from 
## Command: metaxa2_dc -o AbundanceTable.txt -r "metaxa2_allreads" -p "^[^.]+" *level_7.txt

## Remember: Create updated phenotype data and plug in to line 40 below
## Use R script: CompileMetadataFromFarm.R


###################
###
### Loading abundance table and phenotype data. 
### Preparing phyloseq object.
###
###################

otumat <- read.delim("051721_AbundanceTable.txt", sep = "\t", row.names = 1, check.names = FALSE)
taxmat <- matrix(nrow = nrow(otumat), ncol = 0)
rownames(taxmat) <- rownames(otumat)
taxmat <- as.data.frame(cbind(Taxon = rownames(taxmat), taxmat))
levels <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
taxmat <- taxmat %>% separate(Taxon, levels, sep = ";")
taxmat <- as.matrix(taxmat)

#create phyloseq objects
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)

#load metadata csv file that denotes population
metadata = read.csv("070721_metadata.csv", row.names = 1)
metadata_merged = dplyr::select(metadata, -c("SampleName","Sequencer","Lane","Index"))
metadata_merged = unique(metadata_merged)
rownames(metadata_merged) <- metadata_merged$SampleID

#add population info to physeq object
physeq@sam_data = sample_data(metadata)

#### Remove singletons and doubletons
physeq <- prune_taxa(taxa_sums(physeq) > 2, physeq) 

#normalize abundance by sequencer
sequencers <- c(as.character(unique(metadata$Sequencer)))
for(sequencer in sequencers){
  temp = subset_samples(physeq, Sequencer==sequencer)
  total = median(sample_sums(temp))
  standf = function(x, t=total) round(t * (x / sum(x)))
  temp = transform_sample_counts(temp, standf)
  assign(sequencer,temp)
}

#list = cat(paste(shQuote(sequencers, type="cmd2"), collapse=", "))
physeq_merged = merge_phyloseq(V300059128, V300062015, FP100000947BR, FP100001023TR, V300059139, V300059435, FP100000946BR, V300061963, V300059334, V300060570, V300062007, FP100000945BR, FP100000948BR, FP100001025TR, FP100001024TR, FP100001026TR, FP100000944BR, V300061032, V300059379, DP8400010332BR, V300059000, V300057490, V300058991, V300043045, V300043035, V300059290, V300058990, V300048900, DP8400010343BR, V300049674, V300058998, FP100001022TR)
rm(V300059128, V300062015, FP100000947BR, FP100001023TR, V300059139, V300059435, FP100000946BR, V300061963, V300059334, V300060570, V300062007, FP100000945BR, FP100000948BR, FP100001025TR, FP100001024TR, FP100001026TR, FP100000944BR, V300061032, V300059379, DP8400010332BR, V300059000, V300057490, V300058991, V300043045, V300043035, V300059290, V300058990, V300048900, DP8400010343BR, V300049674, V300058998, FP100001022TR)

#Merge samples with multiple runs (warnings will be produced, don't worry- it is regarding the metadata. will be fixed in next lines)
physeq_merged = merge_samples(physeq_merged, "SampleID",fun = mean)
# Since mean function is additive for abundance counts, divide by number of times sample occurs in dataset
temp_otu <- t(as.data.frame(physeq_merged@otu_table))
temp2 <- as.data.frame(table(metadata$SampleID))
temp2 <- subset(temp2, temp2$Freq>1)
row.names(temp2) <- temp2$Var1
for (z in temp2$Var1) {
  temp_otu[,z] <- round((temp_otu[,z])/temp2[z,"Freq"])
}
temp_otu <- t(temp_otu)

#### Remove singletons and doubletons
physeq <- prune_taxa(taxa_sums(physeq) > 2, physeq)

#update OTU table and metadata in physeq object
OTU = otu_table(temp_otu, taxa_are_rows = FALSE)
physeq_merged@otu_table = OTU
sample_data(physeq_merged) <- metadata_merged

#Only keep  bacteria
physeq_merged = subset_taxa(physeq_merged, Domain=="Bacteria")

# #Only keep top 10 most abundant Phyla. 
# topph = sort(tapply(taxa_sums(physeq_merged), tax_table(physeq_merged)[, "Phylum"], sum), TRUE)[1:10]
# physeq_merged = subset_taxa(physeq_merged, Phylum %in% names(topph))

###################
###
### Creating separate phyloseq objects for all populations and all quartiles. 
### Total of 20 phyloseq objects. 
###
###################

#pops = c("All", "LC", "CB", "CI", "AQ")
pops = c("LC", "CB", "CI", "AQ")
#pops = c("LC")
#levels <- c("Order","Family","Genus","Species")
levels <- c("Genus")

fraction_unclassified = data.frame()
pos_and_neg_edges = data.frame()
hubscores = data.frame()

for (p in pops){
  if (p == "All"){
    physeq_subpop = physeq_merged
  } else
    physeq_subpop = subset_samples(physeq_merged, Population == p)
  
  for (level in levels){
    #Conglomerate to specified taxa level
    physeq_subpop_glom = tax_glom(physeq_subpop, level)
    #Get rid of taxa that aren't present
    physeq_subpop_glom = prune_taxa(taxa_sums(physeq_subpop_glom) > 0, physeq_subpop_glom)
    #Get rid of taxa unclassified at that level
    unclassified = str_subset(as.data.frame(physeq_subpop_glom@tax_table)[[level]], "Unclassified.*")
    physeq_subpop_glom = prune_taxa(!tax_table(physeq_subpop_glom)[,level] %in% unclassified, physeq_subpop_glom)
    
    temp = c(p, level, ntaxa(physeq_subpop_glom), (length(unclassified)/(ntaxa(physeq_subpop_glom)+length(unclassified))))
    fraction_unclassified <- rbind(temp, fraction_unclassified)
    
    # Create separate physeq objects for each quartile (based on total biomass)
    quantiles = quantile(physeq_subpop_glom@sam_data[["average_total_biomass"]], na.rm = TRUE)
    physeq_subpop_glom_1 <- subset_samples(physeq_subpop_glom, average_total_biomass <= quantiles[2])
    physeq_subpop_glom_2 <- subset_samples(physeq_subpop_glom, average_total_biomass > quantiles[2] & average_total_biomass <= quantiles[3])
    physeq_subpop_glom_3 <- subset_samples(physeq_subpop_glom, average_total_biomass > quantiles[3] & average_total_biomass <= quantiles[4])
    physeq_subpop_glom_4 <- subset_samples(physeq_subpop_glom, average_total_biomass > quantiles[4])
    
    #remove zero-count taxa
    physeq_subpop_glom_1 = prune_taxa(taxa_sums(physeq_subpop_glom_1) > 0, physeq_subpop_glom_1)
    physeq_subpop_glom_2 = prune_taxa(taxa_sums(physeq_subpop_glom_2) > 0, physeq_subpop_glom_2)
    physeq_subpop_glom_3 = prune_taxa(taxa_sums(physeq_subpop_glom_3) > 0, physeq_subpop_glom_3)
    physeq_subpop_glom_4 = prune_taxa(taxa_sums(physeq_subpop_glom_4) > 0, physeq_subpop_glom_4)
    
    #Add metadata
    physeq_subpop_glom_1@sam_data$Quantile = paste("Quantile 1 ( <", round(quantiles[2],2), "g)")
    physeq_subpop_glom_2@sam_data$Quantile = paste("Quantile 2 (", round(quantiles[2],2), "g to ",round(quantiles[3], 2), "g)")
    physeq_subpop_glom_3@sam_data$Quantile = paste("Quantile 3 (", round(quantiles[3],2), "g to ",round(quantiles[4], 2), "g)")
    physeq_subpop_glom_4@sam_data$Quantile = paste("Quantile 4 ( >", round(quantiles[4],2), "g)")
    
    quantiles = c("physeq_subpop_glom_1","physeq_subpop_glom_2","physeq_subpop_glom_3","physeq_subpop_glom_4")
    
    #Combine into single physeq object
    physeq_subpop_glom <- merge_phyloseq(physeq_subpop_glom_1, physeq_subpop_glom_2, physeq_subpop_glom_3, physeq_subpop_glom_4)
    
    ###################
    ###
    ### SpiecEasi
    ###
    ###################
    
    ###HUB SCORES
    
    se.mb <- spiec.easi(physeq_subpop_glom, method='mb')#, lambda.min.ratio=1e-2, nlambda =10)#, nlambda = 100)#, nlambda = 20)#, nlambda = 10)#,lambda.min.ratio=1e-2, nlambda = 10)
    print(paste(p, level, "All Quantiles", "All Samples", sep = " "))
    print(getStability(se.mb))
    optbeta <- as.matrix(symBeta(getOptBeta(se.mb)))
    edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]
    
    optbeta <- symBeta(getOptBeta(se.mb))
    edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
    edge_weights <- edge_weights/max(edge_weights)
    ig2.mb <- adj2igraph(getRefit(se.mb),  rmEmptyNodes=TRUE,
                         vertex.attr=list(name=taxa_names(physeq_subpop_glom)),
                         edge.attr=list(color=edge_cols, curved = 0, weight = edge_weights))
    
    if(level != "Class"){
      temp <- as.data.frame(head(sort(hub_score(ig2.mb)$vector, decreasing = TRUE), 10))
      colnames(temp)[1] <-  "Hub Score"
      temp <- tibble::rownames_to_column(temp, "OTU")
      temp$Population = p
      temp$Quantile = "All Quantiles"
      temp$Subsample = "All Samples"
      temp$Level = level
      
      hubscores <- rbind(hubscores, temp)
      
    }
    
    ### Other network features with 30 subsamples
    
    if(level != "Class"){
      
      for(i in 1:100){
        #randomly select 30 samples (100 times)
        temp = physeq_subpop_glom
        temp <- prune_samples(sample(sample_names(temp), 30), temp)
      
        se.mb <- spiec.easi(temp, method='mb')#, lambda.min.ratio=1e-2, nlambda =10)#, nlambda = 100)#, nlambda = 20)#, nlambda = 10)#,lambda.min.ratio=1e-2, nlambda = 10)
        print(paste(p, level, "All Quantiles", "30 subsampled", sep = " "))
        print(getStability(se.mb))
        optbeta <- as.matrix(symBeta(getOptBeta(se.mb)))
        edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]
        
        optbeta <- symBeta(getOptBeta(se.mb))
        edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
        edge_weights <- edge_weights/max(edge_weights)
        ig2.mb <- adj2igraph(getRefit(se.mb),  rmEmptyNodes=TRUE,
                             vertex.attr=list(name=taxa_names(physeq_subpop_glom)),
                             edge.attr=list(color=edge_cols, curved = 0, weight = edge_weights))
        
        TotalNodes = gorder(ig2.mb)
        TotalEdges = gsize(ig2.mb)
        PositiveEdges = sum(optbeta@x>0)
        NegativeEdges = sum(optbeta@x<0)
        PropPosNeg = PositiveEdges/NegativeEdges
        PosTotal = PositiveEdges/TotalEdges
        NegTotal = NegativeEdges/TotalEdges
        AvgPathLength = average.path.length(ig2.mb) #Shortest paths between vertices
        Cluster = cluster_louvain(ig2.mb) #Louvain method is an algorithm to detect communities in large networks. It maximizes a modularity score for each community, where the modularity quantifies the quality of an assignment of nodes to communities.This means evaluating how much more densely connected the nodes within a community are, compared to how connected they would be in a random network.
        Modularity = modularity(Cluster) #Strength of division of network into modules
        DegreeDist = degree(ig2.mb)
        AvgDegree = mean(DegreeDist)
        ClusterCoeff = transitivity(ig2.mb) #Transitivity ("clustering coefficient") measures the probability that the adjacent vertices of a vertex are connected.
        
        #degree heterogeneity
        temp2 = 0
        for(i in 1:length(table(DegreeDist))){
          temp = sum((1-(table(DegreeDist)[i]/sum(table(DegreeDist))))^2)
          temp2 = temp2+temp
        }
        HetIndex = sqrt((1/TotalNodes)*temp2)
        MaxHetInt = sqrt(1 - (3/TotalNodes) + ((TotalNodes+2)/(TotalNodes^3)))
        HetMeasure = HetIndex/MaxHetInt
        
        temp <- c(p, level, "All Quantiles", "30", TotalNodes, TotalEdges, PositiveEdges, NegativeEdges, PropPosNeg, PosTotal, NegTotal, AvgPathLength, Modularity, AvgDegree, HetMeasure, ClusterCoeff)
        pos_and_neg_edges <- rbind(pos_and_neg_edges, temp)
        colnames(pos_and_neg_edges) <- c("Population", "Taxa Level", "Biomass Quantile", "# Subsampled","Total Nodes","Total Edges", "Positive Edges", "Negative Edges", "Pos/Neg","Pos/Total","Neg/Total", "Avg Path Length", "Modularity", "Avg Degree", "Heterogeneity", "Clustering Coefficient")
      }
      
    }
    
    # if(level != "Class"){
    #   for(quantile in quantiles){
    #     temp = get(quantile)
    #     
    #     se.mb <- spiec.easi(temp, method='mb')#, lambda.min.ratio=1e-2, nlambda =10)#, nlambda = 100)#, nlambda = 20)#, nlambda = 10)#,lambda.min.ratio=1e-2, nlambda = 10)
    #     print(paste(p, level, quantile, sep = " "))
    #     print(getStability(se.mb))
    #     optbeta <- as.matrix(symBeta(getOptBeta(se.mb)))
    #     edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]
    #     
    #     optbeta <- symBeta(getOptBeta(se.mb))
    #     edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
    #     edge_weights <- edge_weights/max(edge_weights)
    #     ig2.mb <- adj2igraph(getRefit(se.mb),  rmEmptyNodes=TRUE,
    #                          vertex.attr=list(name=taxa_names(temp)),
    #                          edge.attr=list(color=edge_cols, curved = 0, weight = edge_weights))
    #     
    #     ###HUB SCORES
    #     if(level != "Class"){
    #       temp <- as.data.frame(head(sort(hub_score(ig2.mb)$vector, decreasing = TRUE), 10))
    #       colnames(temp)[1] <-  "Hub Score"
    #       temp <- tibble::rownames_to_column(temp, "OTU")
    #       temp$Population = p
    #       temp$Quantile = quantile
    #       temp$Level = level
    #       
    #       hubscores <- rbind(hubscores, temp)
    #       
    #     }
    #     
    #     # for(i in 1:100){
    #     #   #randomly select 50 samples (100 times)
    #     #   temp = get(quantile)
    #     #   temp <- prune_samples(sample(sample_names(temp), 50), temp)
    #     #   
    #     #   se.mb <- spiec.easi(temp, method='mb')#, lambda.min.ratio=1e-2, nlambda =10)#, nlambda = 100)#, nlambda = 20)#, nlambda = 10)#,lambda.min.ratio=1e-2, nlambda = 10)
    #     #   print(paste(p, level, quantile, sep = " "))
    #     #   print(getStability(se.mb))
    #     #   optbeta <- as.matrix(symBeta(getOptBeta(se.mb)))
    #     #   edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]
    #     #   
    #     #   optbeta <- symBeta(getOptBeta(se.mb))
    #     #   edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
    #     #   edge_weights <- edge_weights/max(edge_weights)
    #     #   ig2.mb <- adj2igraph(getRefit(se.mb),  rmEmptyNodes=TRUE,
    #     #                        vertex.attr=list(name=taxa_names(temp)),
    #     #                        edge.attr=list(color=edge_cols, curved = 0, weight = edge_weights))
    #     #   
    #     #   TotalNodes = gorder(ig2.mb)
    #     #   TotalEdges = gsize(ig2.mb)
    #     #   PositiveEdges = sum(optbeta@x>0)
    #     #   NegativeEdges = sum(optbeta@x<0)
    #     #   PropPosNeg = PositiveEdges/NegativeEdges
    #     #   PosTotal = PositiveEdges/TotalEdges
    #     #   NegTotal = NegativeEdges/TotalEdges
    #     #   AvgPathLength = average.path.length(ig2.mb) #Shortest paths between vertices
    #     #   Cluster = cluster_louvain(ig2.mb) #Louvain method is an algorithm to detect communities in large networks. It maximizes a modularity score for each community, where the modularity quantifies the quality of an assignment of nodes to communities.This means evaluating how much more densely connected the nodes within a community are, compared to how connected they would be in a random network.
    #     #   Modularity = modularity(Cluster) #Strength of division of network into modules
    #     #   DegreeDist = degree(ig2.mb)
    #     #   AvgDegree = mean(DegreeDist)
    #     #   ClusterCoeff = transitivity(ig2.mb) #Transitivity ("clustering coefficient") measures the probability that the adjacent vertices of a vertex are connected.
    #     #   
    #     #   #degree heterogeneity
    #     #   temp2 = 0
    #     #   for(i in 1:length(table(DegreeDist))){
    #     #     temp = sum((1-(table(DegreeDist)[i]/sum(table(DegreeDist))))^2)
    #     #     temp2 = temp2+temp
    #     #   }
    #     #   HetIndex = sqrt((1/TotalNodes)*temp2)
    #     #   MaxHetInt = sqrt(1 - (3/TotalNodes) + ((TotalNodes+2)/(TotalNodes^3)))
    #     #   HetMeasure = HetIndex/MaxHetInt
    #     #   
    #     #   temp <- c(p, level, quantile, TotalNodes, TotalEdges, PositiveEdges, NegativeEdges, PropPosNeg, PosTotal, NegTotal, AvgPathLength, Modularity, AvgDegree, HetMeasure, ClusterCoeff)
    #     #   pos_and_neg_edges <- rbind(pos_and_neg_edges, temp)
    #     #   colnames(pos_and_neg_edges) <- c("Population", "Taxa Level", "Biomass Quantile", "Total Nodes","Total Edges", "Positive Edges", "Negative Edges", "Pos/Neg","Pos/Total","Neg/Total", "Avg Path Length", "Modularity", "Avg Degree", "Heterogeneity", "Clustering Coefficient")
    #     #   
    #     #   
    #     # }
    #     
    #   }
    #   
    # }
  }
}
colnames(fraction_unclassified) <- c("Population","Taxa Level","Number of Taxa Across Samples","Fraction Unclassified")
write.csv(fraction_unclassified, "/project/noujdine_61/mgosborn/Gametophytes/NetworkTopology/allphyla/PopulationComparison/fraction_unclassified.csv")
write.csv(pos_and_neg_edges, "/project/noujdine_61/mgosborn/Gametophytes/NetworkTopology/allphyla/PopulationComparison/network_topology.csv")
write.csv(hubscores, "/project/noujdine_61/mgosborn/Gametophytes/NetworkTopology/allphyla/PopulationComparison/network_hubscores.csv")
