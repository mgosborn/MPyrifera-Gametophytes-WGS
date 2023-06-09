setwd("/Users/Mel/Desktop/Artemis/Kelp Microbiome")
library(phyloseq)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(microbiome)
library(vegan)
library(ALDEx2)

## This script runs through alpha diversity, beta diversity, differential abundance, and random forest analysis. 

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
metadata$Population <- gsub("CB", "CP", metadata$Population) #rename CB as CP
metadata$GametophyteCode <- gsub("CB", "CP", metadata$GametophyteCode)
metadata$Sporophyte <- gsub("CB", "CP", metadata$Sporophyte)
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

#Only keep top 10 most abundant Phyla. 
#topph = sort(tapply(taxa_sums(physeq_merged), tax_table(physeq_merged)[, "Phylum"], sum), TRUE)[1:10]
#physeq_merged = subset_taxa(physeq_merged, Phylum %in% names(topph))

###################
###
### Creating separate phyloseq objects for all populations and all quartiles. 
### Total of 20 phyloseq objects. 
###
###################

pops = c("All")#, "LC", "CB", "CI", "AQ")
#levels <- c("Order","Family","Genus")#,"Species")
levels <- "Order"
#p = "All"

fraction_unclassified = data.frame()
sig_aldex2_corr = data.frame()
sig_aldex2_da = data.frame()
rf_top10_pred = data.frame()

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
    
    # Create separate physeq objects for each quantile (based on total biomass)
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
    #physeq_subpop_glom <- merge_phyloseq(physeq_subpop_glom_1, physeq_subpop_glom_4)
    
    # # Only include taxa >1% relative abundance
    # rel_glom <- transform_sample_counts(physeq_subpop_glom, function(x) x / sum(x) )
    # rel_glom = filter_taxa(rel_glom, function(x) sum(x) > 1, TRUE)
    # physeq_subpop_glom = prune_taxa(taxa_names(physeq_subpop_glom) %in% taxa_names(rel_glom), physeq_subpop_glom)

    ###################
    ###
    ### Alpha Diversity (Taxonomic Richness -- Box plots)
    ###
    ###################
    
    file_location = "/Users/Mel/Desktop/Alpha Div/"

    if (p == pops[1]){
      adiv <- data.frame(
        "Observed Richness" = phyloseq::estimate_richness(physeq_subpop_glom, measures = "Observed"),
        "Shannon Diversity" = phyloseq::estimate_richness(physeq_subpop_glom, measures = "Shannon"),
        "Population" = phyloseq::sample_data(physeq_subpop_glom)$Population)
      a <- adiv %>%
        gather(key = metric, value = value, c("Observed", "Shannon")) %>%
        mutate(metric = factor(metric, levels = c("Observed","Shannon"))) %>%
        ggplot(aes(x = Population, y = value, fill = Population)) +
        theme_bw() +
        geom_boxplot() +
        scale_fill_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
        #geom_jitter(aes(color = Population), height = 0, width = .2) +
        labs(x = "", y = "", title = level) + #title = paste("Population: ", p, ", Taxa Level: ", level, sep = "")) +
        facet_wrap(~ metric, scales = "free") +
        theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position="none")

      #a <- a + stat_compare_means(label.y = min(adiv[1]) - 10)
      #a <- a + stat_compare_means(label.y = 1, label.x=1.3)
      a <- a + stat_compare_means(comparisons = list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)),
                                  symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                     symbols = c("****", "***", "**", "*", "ns")))
      #ns: p > 0.05
      #*: p <= 0.05
      #**: p <= 0.01
      #***: p <= 0.001
      #****: p <= 0.0001

      if(level == "Class"){
        pop_class <- a
      }
      if(level == "Order"){
        pop_order <- a
      }
      if(level == "Family"){
        pop_family <- a
      }
      if(level == "Genus"){
        pop_genus <- a
      }

      # if(level == "Species"){
      #   class_box <- a
      # }

      # print(a)
      # # Save alpha div plot
      # plot_name = paste("Alpha Div",p, "by Pop", level, sep = " ")
      # jpeg(paste(file_location, plot_name, ".jpg", sep = ""))
      # plot(a)
      # dev.off()

      # ggarrange(class_box, unclass_box,
      #           labels = c("A", "B"),
      #           ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")

      # ggarrange(pop_class, pop_order, pop_family, pop_genus,
      #           labels = c("A", "B","C","D"),
      #           ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")

    }

    # adiv <- data.frame(
    #   "Observed Richness" = phyloseq::estimate_richness(physeq_subpop_glom, measures = "Observed"),
    #   "Shannon Diversity" = phyloseq::estimate_richness(physeq_subpop_glom, measures = "Shannon"),
    #   "Quantile" = phyloseq::sample_data(physeq_subpop_glom)$Quantile)
    # a <- adiv %>%
    #   gather(key = metric, value = value, c("Observed", "Shannon")) %>%
    #   mutate(metric = factor(metric, levels = c("Observed","Shannon"))) %>%
    #   ggplot(aes(x = Quantile, y = value, fill = Quantile)) +
    #   theme_bw() +
    #   geom_boxplot() +
    #   scale_fill_manual(values=c("#92E4FF", "#4FB6D8", "#1A6E8A", "#003A4D")) +
    #   #geom_jitter(aes(color = Quantile), height = 0, width = .2) +
    #   labs(x = "Average Biomass Quantiles", y = "", title = level) + #paste("Population: ", p, ", Taxa Level: ", level,sep = "")) +
    #   facet_wrap(~ metric, scales = "free") +
    #   theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    #
    # #a <- a + stat_compare_means(label.y = 1)
    # a <- a + stat_compare_means(comparisons = list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)),
    #                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    #                                                symbols = c("****", "***", "**", "*", "ns")))
    # #ns: p > 0.05
    # #*: p <= 0.05
    # #**: p <= 0.01
    # #***: p <= 0.001
    # #****: p <= 0.0001
    #
    # #print(a)
    # # Save alpha div plot
    # # plot_name = paste("Alpha Div",p, level, sep = " ")
    # # jpeg(paste(file_location, plot_name, ".jpg", sep = ""))
    # # plot(a)
    # # dev.off()
    #
    # if(level == "Class"){
    #   bio_class <- a
    # }
    # if(level == "Order"){
    #   bio_order <- a
    # }
    # if(level == "Family"){
    #   bio_family <- a
    # }
    # if(level == "Genus"){
    #   bio_genus <- a
    # }
    #
    # ggarrange(bio_class, bio_order, bio_family, bio_genus,
    #           labels = c("A", "B","C","D"),
    #           ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
    
    ###################
    ###
    ### Beta Diversity (PCA)
    ###
    ###################
    
    # file_location = "/Users/Mel/Desktop/PCA/"
    # 
    # ## CLR Transform and Ordinate
    # temp <- microbiome::transform(physeq_subpop_glom, 'clr')
    # ord_clr <- phyloseq::ordinate(temp, "RDA", distance = "euclidean")
    # 
    # # Generate distance matrix
    # temp_clr_dist_matrix <- phyloseq::distance(temp, method = "euclidean")
    # 
    # # if (p == pops[1]){
    # #   ## Check if samples cluster beyond that expected by sampling variability w/ PERMANOVA (via vegan adonis)
    # #   ## Difference between centroids
    # #   #ADONIS test
    # #   print(paste("perMANOVA (distance) by Population: ", p, ", Taxa Level: ", level,sep = ""))
    # #   ado <- vegan::adonis(temp_clr_dist_matrix ~ phyloseq::sample_data(temp)$Population)
    # #   print(ado)
    # # 
    # #   # Check for differences in dispersion
    # #   print(paste("Differences in dispersion by Population: ", p, ", Taxa Level: ", level,sep = ""))
    # #   dispr <- vegan::betadisper(temp_clr_dist_matrix, phyloseq::sample_data(temp)$Population)
    # #   print(dispr)
    # #   print(permutest(dispr))
    # # 
    # #   ## Differences in dispersion
    # #   ## ANOSIM
    # #   print(paste("ANOSIM (dispersion) by Population: ", p, ", Taxa Level: ", level,sep = ""))
    # #   anosim <- anosim(temp_clr_dist_matrix, phyloseq::sample_data(temp)$Population)
    # #   print(anosim)
    # # 
    # #   #Scale axes and plot ordination
    # #   clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
    # #   clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
    # #   temp2 <- phyloseq::plot_ordination(physeq_subpop_glom, ord_clr, type="samples", color="Population", label = "", title = level) + #title = paste("Population: ", p, ", Taxa Level: ", level, " Pr(>F): ", ado$aov.tab$`Pr(>F)`[1], sep = "")) +
    # #     theme_bw() +
    # #     geom_point(size = .5) +
    # #     coord_fixed(clr2 / clr1) +
    # #     stat_ellipse(aes(group = Population), linetype = 2) +
    # #     scale_color_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
    # #     theme(plot.title = element_text(hjust = 0.5))
    # #   #dev.new()
    # #   #print(temp2)
    # #   # Save PCA plot
    # #   plot_name = paste("PCA by Pop",p, level, sep = " ")
    # #   jpeg(paste(file_location, plot_name, ".jpg", sep = ""))
    # #   plot(temp2)
    # #   dev.off()
    # # 
    # #   if(level == "Class"){
    # #     poppca_class <- temp2
    # #   }
    # #   if(level == "Order"){
    # #     poppca_order <- temp2
    # #   }
    # #   if(level == "Family"){
    # #     poppca_family <- temp2
    # #   }
    # #   if(level == "Genus"){
    # #     poppca_genus <- temp2
    # #   }
    # #   if(level == "Species"){
    # #     poppca_species <- temp2
    # #   }
    # # 
    # #   ggarrange(poppca_class, poppca_order, poppca_family, poppca_genus,
    # #             labels = c("A", "B", "C","D"), common.legend = TRUE, legend = "bottom",
    # #             ncol = 2, nrow = 2)
    # # 
    # #   }
    # 
    # ## Check if samples cluster beyond that expected by sampling variability w/ PERMANOVA (via vegan adonis)
    # #ADONIS test
    # print(paste("perMANOVA by Quantile for Population: ", p, ", Taxa Level: ", level,sep = ""))
    # ado <- vegan::adonis(temp_clr_dist_matrix ~ phyloseq::sample_data(temp)$Quantile)
    # print(ado)
    # 
    # # Check if ADONIS confounded by differences in dispersion
    # dispr <- vegan::betadisper(temp_clr_dist_matrix, phyloseq::sample_data(temp)$Quantile)
    # print(dispr)
    # print(permutest(dispr))
    # 
    # #Scale axes and plot ordination
    # clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
    # clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
    # temp2 <- phyloseq::plot_ordination(physeq_subpop_glom, ord_clr, type="samples", color="Quantile", label = "") + #, title = level) + #title = paste("Population: ", p, ", Taxa Level: ", level, " Pr(>F): ", ado$aov.tab$`Pr(>F)`[1], sep = "")) +
    #   theme_bw() +
    #   geom_point(size = 2) +
    #   coord_fixed(clr2 / clr1) +
    #   stat_ellipse(aes(group = Quantile), linetype = 2) +
    #   scale_color_manual(values=c("#92E4FF", "#003A4D")) +
    #   theme(plot.title = element_text(hjust = 0.5))
    # 
    # # "#92E4FF", "#4FB6D8", "#1A6E8A", "#003A4D"
    # 
    # #print(temp2)
    # # Save PCA plot
    # plot_name = paste("PCA Quantile",p, level, sep = " ")
    # jpeg(paste(file_location, plot_name, ".jpg", sep = ""))
    # plot(temp2)
    # dev.off()
    # 
    # if(level == "Class"){
    #   biopca_class <- temp2
    # }
    # if(level == "Order"){
    #   biopca_order <- temp2
    # }
    # if(level == "Family"){
    #   biopca_family <- temp2
    # }
    # if(level == "Genus"){
    #   biopca_genus <- temp2
    # }
    # if(level == "Species"){
    #   biopca_species <- temp2
    # }
    # 
    # # ggarrange(biopca_order, biopca_family, biopca_genus,
    # #           labels = c("A", "B", "C"), common.legend = TRUE, legend = "bottom",
    # #           ncol = 1, nrow = 3)

    ###################
    ###
    ### ALDEx2 (Correlation & Differential Abundance)
    ###
    ###################
    
    # ## Correlation against continuous average total biomass (keep only Spearman)
    # aldex2_corr <- ALDEx2::aldex(t(data.frame(phyloseq::otu_table(physeq_subpop_glom))), phyloseq::sample_data(physeq_subpop_glom)$average_total_biomass, test="corr", effect = TRUE, denom="all", cont.var = phyloseq::sample_data(physeq_subpop_glom)$average_total_biomass)
    # #Clean up presentation
    # temp <- aldex2_corr %>%
    #   rownames_to_column(var = "OTU") %>%
    #   filter(spearman.ep < 0.05 | kendall.ep < 0.05) %>%
    #   arrange(spearman.erho, spearman.ep) %>%
    #   dplyr::select(OTU, spearman.erho, spearman.ep, spearman.eBH, kendall.etau, kendall.ep, kendall.eBH)
    # 
    # if(nrow(temp)>0){
    #   temp <- cbind(p, level, temp)
    #   sig_aldex2_corr <- rbind(temp, sig_aldex2_corr)
    # }
    # 
    # ## Differential Abundance (Wilcoxon bottom vs top quantile)
    # aldex_physeq <- subset_samples(physeq_subpop_glom, Quantile != unique(physeq_subpop_glom@sam_data$Quantile)[2])
    # aldex_physeq <- subset_samples(aldex_physeq, Quantile != unique(aldex_physeq@sam_data$Quantile)[2])
    # 
    # aldex2_da <- ALDEx2::aldex(t(data.frame(phyloseq::otu_table(aldex_physeq))), phyloseq::sample_data(aldex_physeq)$Quantile, test="t", effect = TRUE, denom="all")
    # 
    # #Clean up presentation
    # temp <- aldex2_da %>%
    #   rownames_to_column(var = "OTU") %>%
    #   filter(wi.ep < 0.05) %>%
    #   arrange(effect, wi.ep) %>%
    #   dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)
    # #print(sig_aldex2)
    # 
    # if(nrow(temp)>0){
    #   temp <- cbind(p, level, temp)
    #   sig_aldex2_da <- rbind(temp, sig_aldex2_da)
    # }

    ###################
    ###
    ### Multiple Regression (glm with significant taxa)
    ###
    ################### 
    
    # Class: Clostridia, Cytophagia
    # Order: Thermoanaerobacterales, Cytophagales
    # Family: Family III, Flammeovirgaceae, Rhodobacteraceae, Rhodobiaceae, Pseudomonadaceae, Rickettsiaceae
    # Genus: Mesorhizobium, Caldicellulosiruptor, Shinella, Acanthopleuribacter, Parvibaculum, Kangiella, Labrenzia, Pseudomonas, Orientia
    # Species: Mesorhizobium sp., Aquamicrobium sp., Mesorhizobium genosp., Zoogloea ramigera, Nitratireductor sp., Parvibaculum lavamentivorans, Labrenzia sp., Kangiella sp.

    # if(level == "Class"){
    #   glm_physeq <- microbiome::transform(physeq_subpop_glom, 'clr')
    #   glm_physeq = subset_taxa(glm_physeq, Class == "Cytophagia" | Class == "Verrucomicrobiae")
    #   taxa_names(glm_physeq) <- glm_physeq@tax_table[,3]
    #   temp <- as.data.frame(glm_physeq@otu_table)
    #   temp2 <- as.data.frame(glm_physeq@sam_data)
    #   temp2 <- temp2[,"average_total_biomass"]
    #   glm_data <- cbind(temp,temp2)
    #   class_model <- glm(average_total_biomass~.,family = "gaussian", data=glm_data)
    # }
    # 
    # if(level == "Order"){
    #   glm_physeq <- microbiome::transform(physeq_subpop_glom, 'clr')
    #   glm_physeq = subset_taxa(glm_physeq, Order == "Sphingomonadales" | Order=="Thermoanaerobacterales")
    #   taxa_names(glm_physeq) <- glm_physeq@tax_table[,4]
    #   temp <- as.data.frame(glm_physeq@otu_table)
    #   temp2 <- as.data.frame(glm_physeq@sam_data)
    #   temp2 <- temp2[,"average_total_biomass"]
    #   glm_data <- cbind(temp,temp2)
    #   order_model <- glm(average_total_biomass~.,family = "gaussian", data=glm_data)
    # }
    # 
    # if(level == "Family"){
    #   glm_physeq <- microbiome::transform(physeq_subpop_glom, 'clr')
    #   glm_physeq = subset_taxa(glm_physeq, Family=="Family III" | Family== "Brucellaceae")
    #   taxa_names(glm_physeq) <- glm_physeq@tax_table[,5]
    #   temp <- as.data.frame(glm_physeq@otu_table)
    #   temp2 <- as.data.frame(glm_physeq@sam_data)
    #   temp2 <- temp2[,"average_total_biomass"]
    #   glm_data <- cbind(temp,temp2)
    #   family_model <- glm(average_total_biomass~.,family = "gaussian", data=glm_data)
    # }
    # 
    # if(level == "Genus"){
    #   glm_physeq <- microbiome::transform(physeq_subpop_glom, 'clr')
    #   glm_physeq = subset_taxa(glm_physeq, Genus=="Mesorhizobium")
    #   taxa_names(glm_physeq) <- glm_physeq@tax_table[,6]
    #   temp <- as.data.frame(glm_physeq@otu_table)
    #   temp2 <- as.data.frame(glm_physeq@sam_data)
    #   temp2 <- temp2[,"average_total_biomass"]
    #   glm_data <- cbind(temp,temp2)
    #   genus_model <- glm(average_total_biomass~.,family = "gaussian", data=glm_data)
    # }
    # 
    # if(level == "Species"){
    #   glm_physeq <- microbiome::transform(physeq_subpop_glom, 'clr')
    #   glm_physeq = subset_taxa(glm_physeq, Species== "Mesorhizobium sp" & Genus == "Mesorhizobium" | Species==  "Labrenzia sp")
    #   taxa_names(glm_physeq) <- glm_physeq@tax_table[,7]
    #   temp <- as.data.frame(glm_physeq@otu_table)
    #   temp2 <- as.data.frame(glm_physeq@sam_data)
    #   temp2 <- temp2[,"average_total_biomass"]
    #   glm_data <- cbind(temp,temp2)
    #   species_model <- glm(average_total_biomass~.,family = "gaussian", data=glm_data)
    # }

    # summary(class_model)
    # summary(order_model)
    # summary(family_model)
    # summary(genus_model)
    # summary(species_model)
  }
}

###################
###
### Barplots (Taxonomic Richness)
###
###################

### Class (color options: “magma”, “plasma”, “inferno”, “civids”, “mako”, “rocket”, "turbo")
### color options: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

glom <- tax_glom(physeq_subpop_glom, taxrank = 'Class')
# create dataframe from phyloseq object
dat <- data.table(psmelt(glom))
# convert Phylum to a character vector from a factor because R
dat$Class <- as.character(dat$Class)
# group dataframe by Phylum, calculate median rel. abundance
dat[, median := median(Abundance, na.rm = TRUE),
    by = "Class"]
# Change name to remainder of Phylum less than 1%
dat[(median <= 0.01), Class := "< 1% Abundance"]

classplot <- ggplot(data = dat, mapping = aes_string(x = "Population",y = "Abundance")) +
  geom_bar(aes(fill=Class), stat="identity", position="fill") +
  theme_bw() +
  ylab("Relative Abundance") + xlab("") +
  scale_fill_viridis(discrete=TRUE)

# ## If we want to filter out <1% abun. taxa:
# physeq_subpop_glom
# rel_glom <- tax_glom(physeq_subpop_glom, taxrank = 'Class')
# rel_glom <- transform_sample_counts(rel_glom, function(x) x / sum(x) )
# rel_glom = filter_taxa(rel_glom, function(x) sum(x) > 1, TRUE)
#
# classplot <- ggplot(data = psmelt(rel_glom), mapping = aes_string(x = "Population",y = "Abundance")) +
#   geom_bar(aes(fill=Class), stat="identity", position="stack") +
#   ylab("Relative Abundance") +
#   scale_fill_viridis(discrete=TRUE)

### Order

glom <- tax_glom(physeq_subpop_glom, taxrank = 'Order')
# create dataframe from phyloseq object
dat <- data.table(psmelt(glom))
# convert Phylum to a character vector from a factor because R
dat$Order <- as.character(dat$Order)
# group dataframe by Phylum, calculate median rel. abundance
dat[, median := median(Abundance, na.rm = TRUE),
    by = "Order"]
# Change name to remainder of Phylum less than 1%
dat[(median <= 0.01), Order := "< 1% Abundance"]

ordplot <- ggplot(data = dat, mapping = aes_string(x = "Population",y = "Abundance")) +
  geom_bar(aes(fill=Order), stat="identity", position="fill") +
  theme_bw() +
  ylab("") + xlab("") +
  scale_fill_viridis(discrete=TRUE)

### Family

glom <- tax_glom(physeq_subpop_glom, taxrank = 'Family')
# create dataframe from phyloseq object
dat <- data.table(psmelt(glom))
# convert Phylum to a character vector from a factor because R
dat$Family <- as.character(dat$Family)
# group dataframe by Phylum, calculate median rel. abundance
dat[, median := median(Abundance, na.rm = TRUE),
    by = "Family"]
# Change name to remainder of Phylum less than 1%
dat[(median <= 0.01), Family := "< 1% Abundance"]

famplot <- ggplot(data = dat, mapping = aes_string(x = "Population",y = "Abundance")) +
  geom_bar(aes(fill=Family), stat="identity", position="fill") +
  theme_bw() +
  ylab("") + xlab("") +
  scale_fill_viridis(discrete=TRUE)

ggarrange(classplot, ordplot, famplot + rremove("y.text"),
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

### Most abundant overall and by pop

# physeq_subpop_glom_AQ = subset_samples(physeq_subpop_glom, Population == "AQ")
# physeq_subpop_glom_CI = subset_samples(physeq_subpop_glom, Population == "CI")
# physeq_subpop_glom_CP = subset_samples(physeq_subpop_glom, Population == "CP")
# physeq_subpop_glom_LC = subset_samples(physeq_subpop_glom, Population == "LC")
# 
# top_overall <- as.data.frame(sort(tapply(taxa_sums(physeq_subpop_glom), tax_table(physeq_subpop_glom)[, "Species"], sum), TRUE)/sum(taxa_sums(physeq_subpop_glom))*100) %>% filter_at(vars(1), any_vars(. > .1))
# top_AQ <- as.data.frame(sort(tapply(taxa_sums(physeq_subpop_glom_AQ), tax_table(physeq_subpop_glom_AQ)[, "Species"], sum), TRUE)/sum(taxa_sums(physeq_subpop_glom_AQ))*100) %>% filter_at(vars(1), any_vars(. > .1))
# top_CI <- as.data.frame(sort(tapply(taxa_sums(physeq_subpop_glom_CI), tax_table(physeq_subpop_glom_CI)[, "Species"], sum), TRUE)/sum(taxa_sums(physeq_subpop_glom_CI))*100) %>% filter_at(vars(1), any_vars(. > .1))
# top_CP <- as.data.frame(sort(tapply(taxa_sums(physeq_subpop_glom_CP), tax_table(physeq_subpop_glom_CP)[, "Species"], sum), TRUE)/sum(taxa_sums(physeq_subpop_glom_CP))*100) %>% filter_at(vars(1), any_vars(. > .1))
# top_LC <- as.data.frame(sort(tapply(taxa_sums(physeq_subpop_glom_LC), tax_table(physeq_subpop_glom_LC)[, "Species"], sum), TRUE)/sum(taxa_sums(physeq_subpop_glom_LC))*100) %>% filter_at(vars(1), any_vars(. > .1))

# ###################
# ###
# ### Outputs
# ###
# ###################
# 
# colnames(fraction_unclassified) <- c("Population","Taxa Level","Number of Taxa Across Samples","Fraction Unclassified")
# write.csv(fraction_unclassified, "/Users/Mel/Desktop/fraction_unclassified.csv")
#write.csv(sig_aldex2_corr, "/Users/Mel/Desktop/sig_aldex2_corr.csv")
#write.csv(sig_aldex2_da, "/Users/Mel/Desktop/sig_aldex2_da.csv")
#write.csv(rf_top10_pred, "/Users/Mel/Desktop/rf_top10_pred.csv")
