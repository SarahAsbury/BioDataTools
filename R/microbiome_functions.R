
# Functions: Alpha Diversity Pipeline -------------------------------
# === Step 1: calculate individual alpha diversity === 
alpha.diversity.df <- function(physeq, 
                               alpha.metrics = c("Shannon", "Chao1", "InvSimpson", "Observed") #any metrics accepted by phyloseq::estimate_richness are accepted here
)
  #Input phyloseq object
  #Returns dataframe with alpha diversity metrics 
{
  alpha <- estimate_richness(physeq, measures= alpha.metrics)
  alpha_meta <- left_join(physeq %>% sample_data %>% data.frame %>% rownames_to_column("join_ID"),
                          alpha %>% rownames_to_column("join_ID"), by = "join_ID") %>% select(-join_ID)
  return(alpha_meta)
}


# === Step 2: calculate alpha diversity group differences
alpha.diversity.pairwise <- function(alpha.diversity.data, #dataframe output from alpha.diversity.df 
                                     alpha.metrics = c("Shannon", "Chao1", "InvSimpson", "Observed"), #should be the same metrics input to alpha.diversity.df 
                                     group #character string; column for alpha diversity comparison; column should only have 2 levels/unique character values
)
  #pairwise analysis for 2 groups only 
{
  #unit test: 2-group comparison 
  unique.group <- eval(parse(text = (paste0("alpha.diversity.data$", group)))) %>% unique
  
  if(length(unique.group) > 2){
    stop("More than 2 categories in group variable. This function is only for pairwise commparison of 2 groups. Please ensure NAs are also removed before running function")
  }
  
  #alpha diversity loop
  kw.res.out <- data.frame()
  for(i in alpha.metrics){
    
    kw <- kruskal.test(eval(parse(text = i))~ eval(parse(text = group)), alpha.diversity.data)
    kw.res <- c(i, unique.group[1], unique.group[2], kw$p.value)
    kw.res.out <- rbind(kw.res.out, kw.res)
  }
  colnames(kw.res.out) <- c("metric", "group1", "group2", "kw.p")
  return(kw.res.out)
}



# === Step 3: plot alpha diversity group differences === 
alpha.diversity.pairwise.plots <- function(alpha.diversity.data, 
                                           alpha.diversity.pairwise.data,
                                           alpha.metrics = c("Shannon", "Chao1", "InvSimpson", "Observed"), #should be the same metrics input to alpha.diversity.df 
                                           group #should be the same group used to generate alpha.diversity.pairwise
){
  
  plot.list <- list()
  for(i in alpha.metrics){
    p.position <- (eval(parse(text = (paste0("alpha.diversity.data$", i)))) %>% max)* 1.1
    
    
    p <- ggboxplot(alpha.diversity.data, 
                   x = group, y = i, 
                   color = group,
                   add = "jitter") + 
      stat_pvalue_manual(alpha.diversity.pairwise.data %>% filter(metric == i) %>% select(-metric) %>% rename(p = kw.p) %>% mutate(p = ifelse(as.numeric(p) > 0.0001, round(as.numeric(p), digits = 4), ">0.0001")), 
                         label = "p", 
                         y.position = p.position) + 
      theme(legend.position = "none")
    plot.list[[i]] <- p
    
  }
  out <- cowplot::plot_grid(plotlist = plot.list, ncol = 2, scale = 0.95)
  return(out)
}



# === Final: full pipeline === 
alpha.diversity.pipeline <- function(physeq, 
                                     alpha.metrics =  c("Shannon", "Chao1", "InvSimpson", "Observed"), 
                                     group #character string; column for alpha diversity comparison; column should only have 2 levels/unique character values
){
  alpha.diversity.data <- alpha.diversity.df(physeq = physeq, 
                                             alpha.metrics = alpha.metrics)
  
  
  alpha.diversity.pairwise.data <- alpha.diversity.pairwise(alpha.diversity.data = alpha.diversity.data, 
                                                            alpha.metrics = alpha.metrics, 
                                                            group = group)
  
  plots <- alpha.diversity.pairwise.plots(alpha.diversity.data = alpha.diversity.data, 
                                          alpha.metrics = alpha.metrics, 
                                          group = group,
                                          alpha.diversity.pairwise.data = alpha.diversity.pairwise.data)
  
  out <- list(alpha.diversity.data, alpha.diversity.pairwise.data, plots)
  names(out) <- c("alpha.df", "alpha.pairwise", "plots")
  return(out)
}




