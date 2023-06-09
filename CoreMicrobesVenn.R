#https://microbiome.github.io/tutorials/core_venn.html

#install.packages("eulerr") # If not installed
library(eulerr)
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)

pseq <- physeq_subpop_glom

# simple way to count number of samples in each group
table(meta(pseq)$Population, useNA = "always")
table(meta(pseq)$Quantile, useNA = "always")

pseq.rel <- microbiome::transform(pseq, "compositional")

disease_states <- unique(as.character(meta(pseq.rel)$Population))
print(disease_states)

list_core <- c() # an empty object to store information

for (n in disease_states){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Population == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 75% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

# Specify colors and plot venn
# supplying colors in the order they appear in list_core
mycols <- c(LC="#004d40", CI="#1e88e5", CP="#ffc107", AQ = "#d81b60" )  
unclass_venn <- plot(venn(list_core), fills = mycols, main = "Classified + Unclassified")

list_core<- Reduce(intersect,list_core)
print(list_core)


ggarrange(class_venn, unclass_venn,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
