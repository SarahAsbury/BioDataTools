#Functions for Mixture Model CLustering using RMixMod

# Load libraries --------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(Rmixmod)
library(logr)
library(foreach)
library(parallel)

# RMixModel - Model Optimization - Data Extraction and Plotting  ---------------------------------------------------
mmCluster.OptRes <- function(mod #MixmodCluster object
)
  #Extract model optimization results to dataframe
{
  model.name <- deparse(substitute(mod))
  model.results <- eval(parse(text = paste0(model.name, "@results")), envir = parent.frame())
  n.models <- length(model.results) #total number of models generated

  df <- data.frame()
  for(i in 1:n.models){
    save.row <- cbind(i,
                      model.results[[i]]@model,
                      model.results[[i]]@nbCluster,
                      model.results[[i]]@criterion,
                      model.results[[i]]@criterionValue,
                      model.results[[i]]@likelihood
    )
    df <- rbind(df, save.row)

  }
  colnames(df) <- c("Model.Index", "Model", "k", "Criterion", "CriterionScore", "Likelihood")
  df <- df %>% mutate(CriterionScore = as.numeric(CriterionScore),
                      k = as.factor(as.numeric(k)),
                      Likelihood = as.numeric(Likelihood))
  return(df %>% arrange(CriterionScore))
}


plot.OptRes <- function(OptRes.df #output from mmCluster.OptRes
)
  #Plot model optimization results using the Criterion Score used in model optimization (i.e. mixmodCluster(criterion = ...))
{
  shape.var <- rep(c(16, 15, 17, 18),
                   ceiling(length(OptRes.df$Model %>% unique)/4))

  plot1 <- ggplot(data = OptRes.df, aes(x = k, y = CriterionScore, group = Model, color = Model, fill = Model, shape = Model)) +
    geom_line() +
    geom_point(alpha = 0.8, size = 3) + theme_classic() +
    scale_shape_manual(values = shape.var)

  plot2 <- ggplot(data = OptRes.df, aes(x = k, y = Likelihood, group = Model, color = Model, fill = Model)) +
    geom_line() +
    geom_point(alpha = 0.75, size = 2.5) + theme_classic()

  legend <- get_legend(plot1)
  cowplot::plot_grid(plot1 + theme(legend.position = "none"),
                     plot2 + theme(legend.position = "none"),
                     legend, ncol = 3) #deprecated; plot of both criterion (bic or icl) and likielhood
  return(plot1)
}



# R mix mod - Clustering stability  ---------------------------------------
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


#Functtion: Extract optimization results from rmixmod
extract.rmixmod.clusters <- function(rmixmod.object, #rmixmod clustering object
                                     model.index #index of which models/runs to extract from rmixmod object
)

  #input: rmixmod object and index of which models to extract information from
  #outputs: model cluster assignment, probabilities, and metadata (including criterion score - i.e ICL/BIC)
{

  #Initalize loop
  cluster.assignment <- list()
  cluster.probability <- list()
  meta <- data.frame()
  counter <- 0

  #Loop - extract relevant models
  for (i in model.index){
    counter <- counter + 1
    mod <- rmixmod.object@results[[i]]
    samples <- rmixmod.object@data %>% rownames

    cluster.assignment[[counter]] <- data.frame(samples = samples, cluster = mod@partition)
    cluster.probability[[counter]] <- mod@proba

    metadata <- c(counter, #og.cluster index
                  i, #mod index
                  mod@model, #model at indexed value
                  mod@nbCluster, #k at indexed value
                  mod@criterionValue)
    meta <- rbind(meta, metadata)
  }
  colnames(meta) <- c("cluster.index", "mod.index", "model", "k", "CriterionScore")
  print(meta)

  cluster.results <- list(meta, cluster.assignment,cluster.probability) #combine results to list
  names(cluster.results) <- c("meta", "cluster", "probability")
  str(cluster.results)

  #return
  return(cluster.results)
}




# Logr workaround ---------------------------------------------------------
#Originally had logr output print statements
#No longer a log

bootstrap.rmixmod.nolog <- function(df,
                              nbCluster.input = 2:5,
                              dataType.input = "quantitative",
                              criterion.input = "ICL",
                              models.input = mixmodGaussianModel(),
                              nboot = 1000,
                              bestmods.df, #if only testing specific k and model combinations, provide df
                              directory, #full path of directory to export results
                              numCores = 1, #setup number of cares to use for parallel computing
                              rmixmod.seed = 1, #clustering seed
                              log = FALSE
)

{
  #Set wd
  setwd(directory)

  #Open log
  if(log == TRUE){
    sink(file = paste0("bootstrap_rmixmod_log_", Sys.time(),".txt"))
  }

  #Declare wd
  print("Set directory")
  print(paste("Export results to:", directory))

  #Setup parallel computing
  doParallel::registerDoParallel(numCores)
  print("Register parallel cores")


# === Original Model Clustering ===
#Cluster with full dataset
print("Generate original clusters")
og.mod <- mixmodCluster(df,
                        nbCluster = nbCluster.input,
                        dataType = dataType.input,
                        criterion = criterion.input,
                        models = models.input, #only uses specified models (should be pre-optimized to select for models to test stability on)
                        seed = rmixmod.seed
)

#Index results to be extracted
#Only extract k cluster models in specific.model.df
og.meta <- mmCluster.OptRes(og.mod) %>% mutate(model_k = paste(Model, k, sep = "_")) %>% filter(model_k %in% bestmods.df$model_k)
og.index <- og.meta %>% pull(Model.Index) %>% as.numeric()

#Extract cluster assignment results
og.clusters <- extract.rmixmod.clusters(og.mod, og.index)
print("Original cluster assignment results:")
print(og.clusters)

# === Bootstrapped Model Clustering ===
#Initialize outside loop
boot.cluster.parent <- data.frame()

#Bootstrap
x <- foreach (i = 1:nboot, .combine=rbind) %dopar% {
  print(paste("Resample set", i))

  #initialize inside loop
  cluster.similarity.df <- data.frame()


  #Resample with replacement
  set.seed(i)
  resampled <- df[sample(nrow(df), replace = TRUE), ]
  print(paste("Samples in bootstrapped set:", i))
  print(rownames(resampled))

  #Cluster
  boot.mod <- mixmodCluster(resampled,
                            nbCluster = nbCluster.input,
                            dataType = dataType.input,
                            criterion = criterion.input,
                            models = models.input, #only uses specified models (should be pre-optimized to select for models to test stability on)
                            seed = rmixmod.seed
  )

  #Index results for extraction
  boot.meta <- mmCluster.OptRes(boot.mod) %>% mutate(model_k = paste(Model, k, sep = "_")) %>% filter(model_k %in% bestmods.df$model_k)
  boot.index <- boot.meta %>% pull(Model.Index) %>% as.numeric()

  #Extract cluster assignment results
  boot.clusters <- extract.rmixmod.clusters(boot.mod, boot.index)
  print(paste("Bootstrapped sample", i, "results"))
  print(boot.clusters)

  #Save boostraped cluster assignment
  assign(boot.cluster.parent, rbind(boot.cluster.parent, boot.clusters), envir = parent.env(environment()))


  # == Identify which bootstrapped cluster is most similar to original cluster ==
  #Calculate Jaccard index for each cluster
  print("original cluster metadata:")
  print(head(og.clusters[["meta"]]))
  print("bootstrapped cluster metadata")
  print(head(boot.clusters[["meta"]]))



  # === cluster similarity ===
  for(j in 1:(og.clusters[["meta"]] %>% nrow())) #repeat for each selected model
  {

    #Get i-th model metadata
    model.og.meta <- og.clusters[["meta"]] %>% filter(cluster.index == j) #get information for i-th model

    #Original and bootsrapped i-th model index
    model.og.index <- j #extract i-th model from original clustering
    model.boot.index <- boot.clusters[["meta"]] %>% filter(model == model.og.meta$model & k == model.og.meta$k) %>% pull(cluster.index) %>% as.numeric() #extract i-th model for bootstrapped clustering

    #Extract sample cluster assignment for i-th model
    og.cluster.samples <- og.clusters[["cluster"]][[model.og.index]]
    boot.cluster.samples <- boot.clusters[["cluster"]][[model.boot.index]]


    #Number of clusters in i-th model
    k.clusters <- model.og.meta %>% pull(k)


    # == cluster similarity ==

    # = qc =
    #Can we remove duplicates samples from bootstrapped samples? (i.e. do both occurrences of the sample get assigned to the same cluster?)
    qc.remove.duplicated.1 <- boot.cluster.samples %>%
      rename(bootstrapped.sample.names = samples) %>% #save duplicated samples
      mutate(samples = gsub(pattern = "\\..*", replacement = "", x = .$bootstrapped.sample.names))
    qc.remove.duplicated.2 <- cbind(qc.remove.duplicated.1, duplicated.sample = duplicated(qc.remove.duplicated.1$samples))
    qc.remove.duplicated.3 <- qc.remove.duplicated.2 %>% filter(duplicated.sample == TRUE)
    qc.remove.duplicated.4 <- qc.remove.duplicated.2 %>% filter(samples %in% qc.remove.duplicated.3$samples)

    qc.remove.duplicated.5 <<- qc.remove.duplicated.4 %>%
      mutate(sample.num = (ifelse(duplicated.sample == TRUE, gsub(pattern = "^.*\\.", replacement = "", x = bootstrapped.sample.names), 0)) %>% as.numeric() + 1) %>%
      select(-c(bootstrapped.sample.names, duplicated.sample))

    qc.remove.duplicated.6 <- qc.remove.duplicated.5 %>% pivot_wider(names_from = sample.num, values_from = cluster)


    print("QC results")
    print(qc.remove.duplicated.6)

    qc.remove.duplicated <- qc.remove.duplicated.6 %>% select(-samples)
    qc.ncol <- ncol(qc.remove.duplicated)
    qc.df <- c()
    for(p in 1:(qc.ncol-1)){
      print(p)
      col1 <- qc.remove.duplicated[,p]
      col2 <- qc.remove.duplicated[,(p +1)]
      qc <- data.frame(col1 = col1, col2 = col2) %>%
        mutate(qc = col1 == col2) %>%
        mutate(qc = ifelse(is.na(qc), TRUE, qc))
      print("qc")
      print(qc %>% filter(qc != TRUE))
      qc.df <- append(qc.df, qc$qc %>% all)
    }

    print(qc.df)

    if(all(qc.df) != TRUE){
      stop("Duplicated samples not assigned to the same cluster")
      print("Error: Duplicated samples not assigned to the same cluster. Running UNDECLARED To kill function.")
    }






    # = calculate cluster similarity =
    for(h in 1:k.clusters){
      #Extract j-th cluster in original clustering model
      cluster.h <- og.cluster.samples %>% filter(cluster == h)
      #Compare each of the bootstrapped clusters to original cluster
      #
      for (m in 1:k.clusters){
        print(paste("Original cluster:", h, "Bootstrapped cluster:", m, "Model:", model.og.meta$model, "Total clusters (k):", k.clusters))
        cluster.m <- boot.cluster.samples %>% filter(cluster == m) %>%
          rename(bootstrapped.sample.names = samples) %>% #save duplicated samples
          mutate(samples = gsub(pattern = "\\..*", replacement = "", x = .$bootstrapped.sample.names))



        #Jaccard Similarity
        og.samp.jac <<- cluster.h %>%
          filter(samples %in% rownames(resampled)) %>% #removes samples not present in bootstrapped sample
          pull(samples)

        boot.samp.jac <<- cluster.m %>%
          group_by(samples) %>% slice(1) %>% ungroup() %>% #extract first occurrence (i.e. remove duplicates)
          pull(samples)

        cluster.similarity <- jaccard(og.samp.jac, boot.samp.jac)

        #Save similarity information and metadata
        cluster.similarity.data <- cbind(model.og.meta, data.frame(i,h, m, cluster.similarity))
        cluster.similarity.df <- rbind(cluster.similarity.df, cluster.similarity.data)
      }





    }
  }
  cluster.similarity.df <- cluster.similarity.df %>% rename(nboot = i, original.cluster = h, bootstrap.cluster = m) %>%
    relocate(nboot,.before = cluster.index)
  print("Cluster similarity df:")
  print(cluster.similarity.df)

}

#Find best matching orginal - bootstrap pair of clusters
best.cluster.df <- x %>%
  mutate(model_k = paste(model, k, sep = "_")) %>%
  group_by(nboot, model_k, original.cluster) %>%
  slice(which.max(cluster.similarity))

print("Best matching clusters:")
print(best.cluster.df)

if(log == TRUE){
sink()
}
return(list(best = best.cluster.df, all = x, original = og.clusters, boot = boot.cluster.parent))
}

