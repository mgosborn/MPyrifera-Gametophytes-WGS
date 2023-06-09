setwd("/Users/Mel/Desktop")
library(dplyr)
library("ggpubr")

ordfam <- read.csv("network_topology.csv")
levels <- c("Order","Family")
factors <- c("Total.Nodes","Total.Edges","Pos.Neg","Avg.Path.Length","Modularity","Avg.Degree","Heterogeneity","Clustering.Coefficient")

ordfam$Population[ordfam$Population=="CB"]<-"CP"
ordfam$Population <- ordered(ordfam$Population,
                             levels = c("AQ", "CI", "CP","LC"))

levels <- c("Family")


for (level in levels){
  levelsub <- subset(ordfam, Taxa.Level == level) 
  
  for (factor in factors){
    group_by(levelsub, Population) %>%
      summarise(
        count = n(),
        mean = mean(get(factor), na.rm = TRUE),
        sd = sd(get(factor), na.rm = TRUE),
        median = median(get(factor), na.rm = TRUE),
        IQR = IQR(get(factor), na.rm = TRUE)
      )
    
    # a <- ggboxplot(levelsub, x = "Population", y = factor, 
    #           color = "Population", palette = c("#D81B60", "#1E88E5", "#FFC107", "#004D40"),
    #           order = c("AQ", "CI", "CP","LC"),
    #           ylab = factor, xlab = "Population", legend = "right")
    
    
    if(factor == "Total.Nodes"){
      y_lab <- "Total Nodes"
    }

    if(factor == "Total.Edges"){
      y_lab <- "Total Edges"
    }

    if(factor == "Pos.Neg"){
      y_lab <- "Positive to Negative Edge Ratio"
    }

    if(factor == "Avg.Path.Length"){
      y_lab <- "Average Path Length"
    }

    if(factor == "Modularity"){
      y_lab <- "Modularity"
    }


    if(factor == "Avg.Degree"){
      y_lab <- "Average Degree"
    }

    if(factor == "Heterogeneity"){
      y_lab <- "Heterogeneity"
    }

    if(factor == "Clustering.Coefficient"){
      y_lab <- "Clustering Coefficient"
    }
    
    a <- ggboxplot(levelsub, x = "Population", y = factor, 
                   fill = "Population", palette = c("#D81B60", "#1E88E5", "#FFC107", "#004D40"),
                   order = c("AQ", "CI", "CP","LC"),
                   ylab = y_lab, xlab = "Population") + grids() + theme(legend.position = "none")
    
    
    
    print(paste("Kruskal test comparing all populations for", factor, "at the", level, "level",sep = " "))
    print(kruskal.test(get(factor) ~ Population, data = levelsub))
    #pairwise.wilcox.test(levelsub[[factor]], levelsub$Population,
    #                     p.adjust.method = "BH")
    
    a <- a + stat_compare_means(comparisons = list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)),
                                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                                                  symbols = c("****", "***", "**", "*", "ns")))
    print(a)
    
    if(factor == "Total.Nodes"){
      Total.Nodes <- a
    }
    
    if(factor == "Total.Edges"){
      Total.Edges <- a
    }

    if(factor == "Pos.Neg"){
      Pos.Neg <- a
    }

    if(factor == "Avg.Path.Length"){
      Avg.Path.Length <- a
    }

    if(factor == "Modularity"){
      Modularity <- a
    }


    if(factor == "Avg.Degree"){
      Avg.Degree <- a
    }

    if(factor == "Heterogeneity"){
      Heterogeneity <- a
    }

    if(factor == "Clustering.Coefficient"){
      Clustering.Coefficient <- a
    }
    # 
    # ggarrange(Total.Nodes, Total.Edges, Pos.Neg, Avg.Path.Length, Modularity, Avg.Degree, Heterogeneity, Clustering.Coefficient,
    #           labels = c("A","B","C","D","E","F","G","H"),
    #           ncol = 2, nrow = 4)
    # 
    # grid.arrange(p1, p2, nrow = 1)
    # 
    # 
    # ggarrange(Total.Nodes,
    #           labels = c("A")),
    #           ncol = 1, nrow = 1)
    
  }
  
  # ggarrange(Total.Nodes, Total.Edges, Pos.Neg, Avg.Path.Length, Modularity, Avg.Degree, Heterogeneity, Clustering.Coefficient,
  #           labels = c("A","B","C","D","E","F","G","H"),
  #           ncol = 2, nrow = 4, common.legend = TRUE, legend = "right")
  
}



tiff("/Users/Mel/Desktop/figS3_5.tiff", units="in", width=9, height=16, res=300)
ggarrange(Total.Nodes, Total.Edges, Pos.Neg, Avg.Path.Length, Modularity, Avg.Degree, Heterogeneity, Clustering.Coefficient,
          labels = c("A","B","C","D","E","F","G","H"),
          ncol = 2, nrow = 4, common.legend = TRUE, legend = "right")
dev.off()





  


