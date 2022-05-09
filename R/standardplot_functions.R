#Common plots - standardized 

# Load packages -----------------------------------------------------------
library(tidyverse)
library(NbClust)
library(rgl)
library(car)
library(pca3d)
library("RColorBrewer")
library(fpc)


# Unit test: Scaled dataframe? -----------------------------------------------------
ut.scaled.sd <- function(df){
  #Unit test: 
  if(df %>% select_if(negate(is.numeric)) %>% ncol != 0){
    print("Error: non-numeric column in the input dataframe.")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  #nrow SD != 0 or mean != 1
  qc.sd <- apply(df, 2 , sd) %>% data.frame %>% rename(qc = ".") %>% mutate(qc = ifelse(round(qc, digits = 0) == 1, TRUE, FALSE)) %>% filter(qc == FALSE) %>% nrow()
  qc.mean <- apply(df, 2 , mean) %>% data.frame %>% rename(qc = ".") %>% mutate(qc = ifelse(round(qc, digits = 0) == 0, TRUE, FALSE)) %>% filter(qc == FALSE) %>% nrow
  
  qc.out <- ifelse(qc.sd == 0 & qc.mean == 0, TRUE, FALSE)
  
  return(qc.out)
}




# PCA values --------------------------------------------------------------------
#generates PCA values with hierarchical clustering for plotting


pca.hclust.fun <- function(k.in, #k number of clusters to generate. Max number = 8 for this function (limited by Dark2 palette). 
                           df, #dataframe. Conains variables for PCA. Must be scaled and NA removed. 
                           label #each sample group labels in the same order as df observations order. (i.e df$group)
)
  #Requires a scaled dataframe without NA
{
  
  #Is the input a dataframe? 
  if(class(df) == FALSE){
    print("Error: input df is not dataframe class")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  #Unit test: scaled df? 
  if(ut.scaled.sd(df) == FALSE){
    print("Error: dataframe is not scaled")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  #Principal Components
  
  PCA.scaled <<- prcomp(df, scale. = FALSE)
  
  summary(PCA.scaled)
  PCA.scaled
  PCA.scaled.values <- PCA.scaled$x %>% data.frame()#Store PCA values for plotting in scatter plot
  PC1.scaled <- PCA.scaled.values$PC1
  PC2.scaled <- PCA.scaled.values$PC2
  PC3.scaled <- PCA.scaled.values$PC3
  
  
  pca.plot <- data.frame(cbind(PC1.scaled, PC2.scaled, PC3.scaled))
  str(pca.plot)
  
  colors <- brewer.pal(n=k.in, name="Dark2")
  
  #Hclust
  dist <- dist(PCA.scaled.values, method = 'euclidian')
  hclust <- hclust(dist, method = 'ward.D2')
  clust.cut <<- as.factor(cutree(hclust, k = k.in))
  
  pca.plot <<- data.frame(cbind(pca.plot, label, clust.cut))
  
  #Cluster stability
  print(1:k.in)
  clust.stability <<- clusterboot(df, 
                                  B=1000, 
                                  clustermethod=hclustCBI,
                                  method="ward.D", 
                                  k=k.in, 
                                  count=FALSE) 
  clust.stability.df <<- data.frame(clust.stability$bootmean)
  clust.stability.df <<- cbind(1:k.in, clust.stability.df)
  colnames(clust.stability.df) <<- c("Cluster", "Jaccard Stability")
  
} #requires a scaled df



# PCA plotting ------------------------------------------------------------

#2D plot

pca.ggplot <- function(pca.plot.df = pca.plot, #contains PC1.scaled, PC2.scaled and label as a single dataframe
                       prcomp.in = PCA.scaled, #output R PCA function. Used to get variance for axis labels. 
                       ellipses = TRUE,
                       group = "label"
)
{
  #Component variance for axis labels
  extract.variance <- summary(prcomp.in)
  extract.variance <- extract.variance$importance %>% data.frame %>% select(PC1, PC2) %>% rownames_to_column(var = "stat") %>% filter(stat == "Proportion of Variance")
  pc1.variance <- extract.variance %>% pull(PC1) 
  pc2.variance <- extract.variance %>% pull(PC2) 
  
  x.axis.name <- paste0("PC1 (", round(pc1.variance * 100, digits = 1), "%)")
  y.axis.name <- paste0("PC2 (", round(pc2.variance * 100, digits = 1), "%)")
  
  #Plot 
  plot <- ggplot(pca.plot.df, aes(x = PC1.scaled, y = PC2.scaled, color = eval(parse(text = paste(group))))) + geom_point() + theme_classic() +
    xlab(x.axis.name) + ylab(y.axis.name) + theme(legend.title=element_blank()) + 
    {if(ellipses == TRUE)
      stat_ellipse(geom = "polygon", alpha = 0.05, aes(fill =eval(parse(text = paste(group)))))
    } 
  
  return(plot)}



