
###################################################################################################
########################## S1 CODE. R FUNCTIONS FOR NULL MODEL ANALYSES ###########################
###################################################################################################

## PAPER: Elevational gradients in beta-diversity reflect variation in the strength of local 
##  community assembly mechanisms across spatial scales
## JOURNAL: PlosONE
##
## AUTHORS: J. Sebastián Tello, Jonathan A. Myers, Manuel J. Macía, Alfredo F. Fuentes, 
##	Leslie Cayola, Gabriel Arellano, M. Isabel Loza, Vania Torrez, Maritza Cornejo, 
##	Tatiana B. Miranda and Peter M. Jørgensen
##
## CONTACT INFO.: J. Sebastián Tello - jsebastiantello@gmail.com or sebastian.tello@mobot.org 
## DATE: 18-June-2014
##
## DESCRIPTION: The following two functions were used for null-model analyses. 
## 	1. 'assemblages.from.pool.randA' was used to produce null local assemblages expected by random
## 	sampling from observed species pools while eliminating the effect of mechanisms controlling
## 	(1) the distribution of species across local assemblages and (2) the regional species abundances. 
##
## 	2. 'div.partition.nullM' was used to calculate beta-diversity from empirical and null 
## 	assemblages, and to calculate beta-deviations.
##
## These functions are also accompanied by two supplementary functions not meant to be called 
## directly by the user.




####################################################################################################
#############################  FUNCTION 1 - assemblages.from.pool.randA  ###########################
####################################################################################################

# SUMMARY 
# A randomization algorithm that produces matrices of species composition expected by randomly 
# re-assigning individuals from the species pool into local assemblages.

# USAGE
# assemblages.from.pool.randA(compo, rand.N=999, fix.local.abund=TRUE, fix.rSAD=TRUE, 
#  save.output=FALSE, save.format="matrices", path.to.save, show.progress=FALSE)

# DESCRIPTION 
# This function implements a randomization algorithm that creates null "local" assemblages as random 
# samples from species pools. It is designed to study the structure of communities while removing 
# the influence of sampling effects and regional species pools. This algorithm was recently used by
# Kraft et al. (2011), De Cáceres et al. (2012), Mori et al. (2013), and Tello et al. (submitted).
# In the randomization algorithm, all individuals have the same probability of becoming part of any 
# local assemblage. The regional abundances of species can be fixed (constrained) to be the same as 
# in the empirical data (when 'fix.rSAD' is TRUE), or regional abundances can be randomized by also 
# randomly re-assigning individuals into species (when 'fix.rSAD' is FALSE). Here, it is important 
# to note that the function makes sure that every species that is present in the region gets at 
# least one individual. In this way, no species can "disappear" from the null composition matrices 
# if they had at least 1 individual in the empirical data. Additionally, the total number of 
# individuals at each local assemblage can be fixed to be identical to the empirical data (when 
# 'fix.local.abund' is TRUE), or it can be allowed to vary as a result from the random distribution 
# of individuals (when 'fix.local.abund' is FALSE). The output of the function is a series of 
# matrices of species composition produced by the randomization algorithm. 

# ARGUMENTS
# compo - a matrix or data frame where rows are sites and columns species. Values within the 
#   matrix represent abundances as integers. The function will not work with presence/absence 
#   matrices or measures of abundance other than number of individuals.
#
# rand.N - an integer determining the number of null matrices to be produced. The default is ‘999’.
#
# fix.local.abund - TRUE/FALSE argument that determines whether the abundances in local assemblages 
#   (row sums) will be constrained to be the same as in the empirical data. The default is 
#   TRUE.
#
# fix.rSAD - TRUE/FALSE argument that determines whether the regional species abundance distribution 
#   (SAD; column sums) will be constrained to be the same as in the empirical data. If FALSE, 
#   individuals are randomly re-assigned to species before they are re-distributed among local 
#   assemblages. The default is TRUE.
#
# save.output - TRUE/FALSE argument that defines whether results should be saved into files 
#   (when TRUE), or should be returned into the console (when FALSE). The default is FALSE.
#
# save.format - character argument that can take one of two values: "matrices" or "list". This 
#   argument applies only if 'save.output' is TRUE. When the "matrices" option is requested, each
#   of the 'rand.N' null matrices will be written into a file as the function finishes each 
#   randomization. If the "list" option is requested, a single list containing all 'rand.N' null 
#   matrices will be saved after all randomizations have been completed. The default is "matrices".
#
# path.to.save - character argument that gives the address of the folder in the local computer 
#   where output files should be saved. This argument applies only if 'save.output' is TRUE. 
#   The default is the working directory. The path should end with the name of the target folder, 
#   not with '\\' or '/'.
#
# show.progress - TRUE/FALSE argument that defines whether the progress of the function across 
#   randomizations should be shown in the console. The default is FALSE.
  
# OUTPUT
# A list with the following 2 elements:
#
# rand.parameters - a vector containing information on the parameters used for the randomization of 
#   the data.
# 
# rand.datasets - When 'save.output' is FALSE, this is a list of length 'rand.N'. Each element of 
#   the list contains a null matrix of species composition that is of identical dimensions as the 
#   empirical matrix ('compo'). When 'save.output' is TRUE, 'rand.datasets' is simply the path given 
#   by the argument 'path.to.save', and  the composition matrices are saved into either a single 
#   list file (when 'save.format' is "list"), or to individual files for each matrix (when 
#   'save.format' is "matrices").

# REFERENCES
# De Cáceres, M., Legendre, P., Valencia, R., Cao, M., Chang, L.-W., Chuyong, G., et al. (2012). 
#   The variation of tree beta diversity across a global network of forest plots. Glob. Ecol. 
#   Biogeogr., 21, 1191–1202.
#
# Kraft, N.J.B., Comita, L.S., Chase, J.M., Sanders, N.J., Swenson, N.G., Crist, T.O., et al. 
#   (2011). Disentangling the Drivers of beta Diversity Along Latitudinal and Elevational 
#   Gradients. Science, 333, 1755–1758.
#
# Mori, A.S., Shiono, T., Koide, D., Kitagawa, R., Ota, A.T. & Mizumachi, E. (2013) Community 
#   assembly processes shape an altitudinal gradient of forest biodiversity. Global Ecology and 
#   Biogeography, 22, 878–888.

# EXAMPLES
if(FALSE) # Prevents examples from running when loading the functions
{
  library(vegan)
  data(mite)

  # 1. Fixed local abundances and regional species abundance distribution  
  null.res.1 <- assemblages.from.pool.randA(compo=mite, rand.N=3, fix.local.abund=TRUE, 
    fix.rSAD=TRUE, save.output=FALSE) 
  par(mfrow=c(3,2), mar=c(4.5,4.5,1,1))
  for(i in 1:length(null.res.1$rand.datasets))
  {
    print(i)
    plot(colSums(mite)~colSums(null.res.1$rand.datasets[[i]]), ylab="Obs. Abundance", 
      xlab="Null Abundance")
    abline(c(0,1))
    plot(rowSums(mite)~rowSums(null.res.1$rand.datasets[[i]]), ylab="Obs. Local Density", 
      xlab="Null Local Density")
    abline(c(0,1))
  }  

  
  # 2. Random local abundances and fixed regional species abundance distribution  
  null.res.2 <- assemblages.from.pool.randA(compo=mite, rand.N=3, fix.local.abund=FALSE, 
    fix.rSAD=TRUE, save.output=FALSE) 
  par(mfrow=c(3,2), mar=c(4.5,4.5,1,1))
  for(i in 1:length(null.res.2$rand.datasets))
  {
    print(i)
    plot(colSums(mite)~colSums(null.res.2$rand.datasets[[i]]), ylab="Obs. Abundance", 
      xlab="Null Abundance")
    abline(c(0,1))
    plot(rowSums(mite)~rowSums(null.res.2$rand.datasets[[i]]), ylab="Obs. Local Density", 
      xlab="Null Local Density")
    abline(c(0,1))
  }  

  
  # 3. Random local abundances and regional species abundance distribution  
  null.res.3 <- assemblages.from.pool.randA(compo=mite, rand.N=3, fix.local.abund=FALSE, 
    fix.rSAD=FALSE, save.output=FALSE) 
  par(mfrow=c(3,2), mar=c(4.5,4.5,1,1))
  for(i in 1:length(null.res.3$rand.datasets))
  {
    print(i)
    plot(colSums(mite)~colSums(null.res.3$rand.datasets[[i]]), ylab="Obs. Abundance", 
      xlab="Null Abundance")
    abline(c(0,1))
    plot(rowSums(mite)~rowSums(null.res.3$rand.datasets[[i]]), ylab="Obs. Local Density", 
      xlab="Null Local Density")
    abline(c(0,1))
  }
}




### FUNCTION CODE ##################################################################################
assemblages.from.pool.randA <- function(compo, rand.N=999, fix.local.abund=TRUE, fix.rSAD=TRUE, 
  save.output=FALSE, save.format="matrices", path.to.save, show.progress=FALSE) 
{
  ## Makes sure there is a path to save files if needed
  if(save.output==TRUE & missing(path.to.save)) path.to.save <- getwd()

  ## Makes sure 'compo' is in the right format
  compo <- as.matrix(compo, ncol=ncol(compo))
  if(is.null(rownames(compo))) rownames(compo) <- paste("site_", 1:nrow(compoT), sep="")
  if(is.null(colnames(compo))) colnames(compo) <- paste("sp_", 1:nrow(compoT), sep="")
  
  ## Calculates density of individuals per site
  site.densities <- rowSums(compo)
  if(min(site.densities)<=0) warning("Some sites have no species")

  ## Calculates total regional abundance
  regional.abundance <- sum(site.densities)
  
  ## Calculates regional abundances per species - the regional SAD
  spp.abund <- colSums(compo)
  
  ## Checks for a couple of potential problems
  if(length(which((compo - round(compo, 0)) != 0)) > 0) 
    stop("This function can not randomize a matrix with species abundances that are not integers")
  if(min(spp.abund)<=0) 
    warning("Some species have abundances of less than 1")
  
  ## Finds the species with at least one individual
  spp.more.than.0.indices <- which(spp.abund>0)
  spp.more.than.0.names <- colnames(compo)[spp.more.than.0.indices]
  
  ## Produces a vector of species identities for each individual
  individual.id <- as.vector(rep(colnames(compo), spp.abund))

  ## Produces an empty list to save null composition matrices
  rand.datasets <- sapply(rep(NA, rand.N), list)
  names(rand.datasets) <- paste("RandDataset_", 1:rand.N, sep="")
  
  for(i in 1:rand.N)
  {
    if(show.progress==TRUE) 
      print(i)
          
    ## When the regional SAD is NOT fixed, produces a null SAD to be used in analyses by 
    ## randomly assigning individuals to species 
    if(fix.rSAD==FALSE)
    {
      # Assigns individuals to species at random
      null.sp.to.ind.assignation.1 <- spp.more.than.0.names # IMPORTANT - This ensures that all 
                                                            # species in the region receive at least 
                                                            # one individual in the null SAD
      null.sp.to.ind.assignation.2 <- character()
          
      if((regional.abundance-length(spp.more.than.0.names))>0) 
        null.sp.to.ind.assignation.2 <- sample(spp.more.than.0.names, 
          regional.abundance-length(spp.more.than.0.names), replace=TRUE)
          
      null.sp.to.ind.assignation <- c(null.sp.to.ind.assignation.1, null.sp.to.ind.assignation.2)
        
      # Calculates null abundances for species present with at least one individual
      pre.null.spp.abund <- as.numeric(table(null.sp.to.ind.assignation))
      pre.null.spp.abund <- pre.null.spp.abund[sample(1:length(pre.null.spp.abund))]
        
      # Inserts the null abundances into a vector that might contain some zeroes (i.e., if 
      # there are species present in the empirical composition table that do not have any 
      # individuals)
      null.spp.abund <- spp.abund
      null.spp.abund[spp.more.than.0.indices] <- pre.null.spp.abund
      
      # Produces a vector of null species identities for each individual
      null.individual.id <- as.vector(rep(colnames(compo), null.spp.abund))
    }
      
    ## When the regional SAD is fixed, makes the null SAD identical to the empirical SAD  
    if(fix.rSAD==TRUE)
    {
      null.spp.abund <- spp.abund
      null.individual.id <- individual.id
    }
  
    ## Assigns individuals to local sites at random
    if(fix.local.abund==TRUE) 
      null.site.assignation <- sample(rep(rownames(compo), times=site.densities), 
        regional.abundance, replace=FALSE)
    
    if(fix.local.abund==FALSE) 
      null.site.assignation <- sample(rownames(compo), regional.abundance, replace=TRUE)

    ## Creates a null species composition matrix using the null assignation of individuals 
    ## to sites
    null.compo <- tapply(rep(1, regional.abundance), list(null.site.assignation, 
      null.individual.id), sum)
    null.compo[is.na(null.compo)] <- 0

    ## Reshapes the null composition matrix in case there is only 1 species
    if(ncol(compo)==1) 
    {
      null.compo <- matrix(null.compo, ncol=1)
      colnames(null.compo) <- colnames(compo)
      rownames(null.compo) <- unique(null.site.assignation)
    }
    
    ## Inserts empty plots into the null composition matrix in case there are any
    if(nrow(null.compo)<nrow(compo))
    {
      missing.plots <- setdiff(rownames(compo), rownames(null.compo))

      empty.plots <- matrix(0, ncol=ncol(null.compo), nrow=length(missing.plots))
      colnames(empty.plots) <- colnames(null.compo)
      rownames(empty.plots) <- missing.plots

      null.compo <- rbind(null.compo, empty.plots)

      if(ncol(compo)==1) 
      {
        null.compo <- matrix(null.compo, ncol=1)
        colnames(null.compo) <- colnames(compo)
        rownames(null.compo) <- c(unique(null.site.assignation), missing.plots)
      }
    }
    
    ## Matches the column and row names in the null and empirical composition matrices  
    if(ncol(compo)>1) null.compo <- null.compo[,match(colnames(compo), colnames(null.compo))]
    null.compo <- null.compo[match(rownames(compo), rownames(null.compo)),]

    ## Reshapes the null composition matrix in case there is only 1 species
    if(ncol(compo)==1)
    {
      null.compo <- matrix(null.compo, ncol=1)
      colnames(null.compo) <- colnames(compo)
      rownames(null.compo) <- rownames(compo)
    }

    ## Adds names to the species with zero abundances (i.e., empty columns)
    if(min(null.spp.abund)<=0)
    {
      colnames(null.compo)[null.spp.abund<=0] <- colnames(compo)[null.spp.abund<=0]
      null.compo[is.na(null.compo)] <- 0
    }

    ## Adds the null composition matrix to the list that compiles the results
    if(save.output==FALSE | save.format=="list") 
      rand.datasets[[i]] <- null.compo
    
    ## If requested, the null composition matrix is saved as a file
    if(save.output==TRUE & save.format=="matrices") 
      write.table(null.compo, file=paste(path.to.save, "\\", "RandDataset_", i, ".txt", sep=""), 
        quote=FALSE, sep="\t", na="NA", dec=".", row.names=TRUE, col.names=TRUE)
  }

  ## If saving a list is requested, the full list of null composition matrices is saved as a file 
  if(save.format=="list")
    save(rand.datasets, file=paste(path.to.save, "\\", "RandDatasets", sep="")) 
  
  ## Makes a list of the parameters used in the randomization
  rand.parameters <- c(fix.local.abund, fix.rSAD, rand.N)
  names(rand.parameters) <- c("fix.local.abund", "fix.rSAD", "rand.N")

  ## Creates the output to return, depending on whether the main output was saved or not
  if(save.output==TRUE) 
    output <- list(rand.parameters, path.to.save)
  else 
    output <- list(rand.parameters, rand.datasets)
    
  names(output) <- c("rand.parameters", "rand.datasets")
  
  output  
}




####################################################################################################
#################################  FUNCTION 2 - div.partition.nullM  ###############################
####################################################################################################

# SUMMARY 
# A function to partition diversity into alpha, gamma and beta-components from an empirical matrix 
# of species composition and from a series of null matrices produced by a randomization algorithm.

# USAGE
# div.partition.nullM(compo, null.matrices, dist.method="bray", binary.data=FALSE) 

# DESCRIPTION 
# This function takes as its main inputs (1) an empirical matrix of species composition where 
# rows are sites and columns are species, and (2) a list where each element is a null matrix 
# of identical dimensions as the empirical one. These null matrices must result from a null 
# model that randomizes the empirical matrix. The function takes these inputs and partitions 
# diversity into alpha, gamma and beta components, and also calculates mean pair-wise 
# dissimilarities in species composition. From these data, it calculates summary statistics 
# comparing empirical to null partitions of diversity, including standardized effect sizes.

# ARGUMENTS
# compo - a matrix or data frame where rows are sites and columns species. Values within the 
#   matrix represent abundances as integers. The function will not work with presence/absence 
#   matrices or measures of abundance other than number of individuals.
#
# null.matrices - a list where each element contains a null matrix of species composition that 
#   is of identical dimensions as the empirical composition matrix (argument 'compo'). 
#
# dist.method - a method to calculate dissimilarities in species composition among pairs of sites. 
#   It can be "hellinger" or any of the options available for the function 'vegdist' in the package
#   'vegan'. 
#
# binary.data - TRUE/FALSE argument that indicates whether the data are to be treated as 
#   presence/absence only (when TRUE). This argument is passed to the argument 'binary' of the 
#   function 'vegdist' in the package 'vegan'.

# OUTPUT
# A list with the following 4 elements:
#
# empirical.values - a vector with the empirical values of diversity. "q0" and "q1" refer to 
#   diversities of order 0 and 1 respectively following Jost 2007 approach to partition diversity.
#   "prop.spp.turnover" is defined as 1 - mean.local.richness/regional.richness. "mean.compo.dist"
#   is the average pair-wise dissimilarity in species composition among all pairs of sites. 
#
# summary.table - a data frame with summary statistics for each of the empirical values of diversity. 
#   It includes the means and standard deviations of null values, as well as effects sizes 
#   (empirical - mean null) and standardized effect sizes (effect size / standard deviation of null). 
#
# rand.details - a data frame with the values of diversity indices for each repetition 
#   of the null model. 

# REFERENCES
# Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 
#   2427–2439.
#
# Tuomisto, H. (2010). A diversity of beta diversities: straightening up a concept gone awry. 
#   Part 1. Defining beta diversity as a function of alpha and gamma diversity. Ecography, 33, 2–22.

# EXAMPLES
if(FALSE) # Prevents examples from running when loading the functions
{
  library(vegan)
  library(vegetarian)
  
  data(mite)

  # Creates random matrices according to a randomization algorithm
  rand.results <- assemblages.from.pool.randA(compo=mite, fix.local.abund=TRUE, fix.rSAD=TRUE, 
    rand.N=999, save.output=FALSE)
  null.compos <- rand.results$rand.datasets

  # Partitions variation for the empirical and null matrices
  res.nullM <- div.partition.nullM(compo=mite, null.matrices=null.compos, dist.method="bray", 
    binary.data=FALSE)

  # Plots histograms of null values vs. empirical values for various measures of beta-diversity
  par(mfrow=c(1,3))
    hist(res.nullM$rand.details$q1.beta, 20, main="", xlab="Beta-diversity of order 1", 
      col="grey50", border="grey50", xlim=range(c(res.nullM$rand.details$q1.beta, 
      res.nullM$summary.table["q1.beta","empirical"])))
    lines(x=rep(res.nullM$summary.table["q1.beta","empirical"], 2), y=c(0, 1000), 
      col="red", lwd=2)

    hist(res.nullM$rand.details$q0.beta, 20, main="", xlab="Beta-diversity of order 0", 
      col="grey50", border="grey50", xlim=range(c(res.nullM$rand.details$q0.beta, 
      res.nullM$summary.table["q0.beta","empirical"])))
    lines(x=rep(res.nullM$summary.table["q0.beta","empirical"], 2), y=c(0, 1000), 
      col="red", lwd=2)

    hist(res.nullM$rand.details$mean.compo.dist, 20, main="", xlab="Mean compo. distance", 
      col="grey50", border="grey50", xlim=range(c(res.nullM$rand.details$mean.compo.dist, 
      res.nullM$summary.table["mean.compo.dist","empirical"])))
    lines(x=rep(res.nullM$summary.table["mean.compo.dist","empirical"], 2), y=c(0, 1000), 
      col="red", lwd=2)

  # Plots SESs for various measures of beta-diversity
  par(mfrow=c(1,1))
    barplot(res.nullM$summary.table[c("q0.beta", "q1.beta", 
      "prop.spp.turnover"),"SES"], names.arg=c("q0.beta", "q1.beta", 
      "prop.spp.turnover"), ylab="Standardized Effect Size")
    
}



### FUNCTION CODE ##################################################################################
div.partition.nullM <- function(compo, null.matrices, dist.method="bray", binary.data=FALSE, 
  print.progress=FALSE)
{
  ## Opens required packages
  require(vegetarian)
  require(vegan)

  ## Makes sure 'compo' is in the right format
  dim.compo <- length(dim(compo))
  if(dim.compo==0) 
    stop("'compo' does not have more than 1 dimension")
  if(dim.compo>2) 
    stop("This function does not know what to do when 'compo' has more than 2 dimensions")
  if(dim.compo==2) 
  {
    if(is.null(rownames(compo))) 
      rownames(compo) <- paste("site_", 1:nrow(compo), sep="")
    if(is.null(colnames(compo))) 
      colnames(compo) <- paste("sp_", 1:nrow(compo), sep="")
  } 
  
  ## Checks that values in the cells are all positive integers
  if(sum((compo - round(compo, 0)) != 0) > 0) 
    stop("This function cannot use species abundances that are not integers")
  if(sum(compo<0) > 0)    
    stop("This function cannot use negative species abundances")
  if(sum(is.na(compo)) > 0)    
    stop("This function cannot use NAs")

    
  ## Calculates density of individuals per site
  site.densities <- rowSums(compo)
  if(min(site.densities)==0) 
    warning("Some empirical sites seem to be empty of individuals")
  
  ## Calculates regional abundances per species - the regional SAD
  spp.abund <- colSums(compo)
  if(min(spp.abund)<=0) 
    warning("Some species have abundances of less than or equal to zero")
    
  ## Calculates species richness across all sites
  regional.richness <- length(which(spp.abund>0))
  
  ## Calculates species richness at each site
  site.richness <- rowSums(compo>0)
  mean.site.richness <- mean(site.richness, na.rm=TRUE)
  
  
  ## Calculates empirical composition distances among all possible pairs of sites
  if(dist.method=="hellinger") 
    compo.dist.dist.by.site <- vegdist(decostand(compo, method="hellinger"), method="euclidean", 
    binary=binary.data)
  else 
    compo.dist.dist.by.site <- vegdist(compo, method=dist.method, binary=binary.data)
  
  emp.mean.compo.dist <- mean(compo.dist.dist.by.site)  
  
  ## Partitions diversity in the empirical composition data 
  emp.q0.alpha <- d(compo, lev="alpha", q=0)
  emp.q0.beta <- d(compo, lev="beta", q=0)
  emp.q0.gamma <- d(compo, lev="gamma", q=0)

  emp.q1.alpha <- d(compo, lev="alpha", q=1)
  emp.q1.beta <- d(compo, lev="beta", q=1)
  emp.q1.gamma <- d(compo, lev="gamma", q=1)

  emp.prop.spp.turnover <- (1-(mean.site.richness/regional.richness))

  ## Pastes together empirical values 
  empirical.values <- c(mean.site.richness, regional.richness, emp.q0.alpha, emp.q0.beta, 
    emp.q0.gamma, emp.q1.alpha, emp.q1.beta, emp.q1.gamma, emp.prop.spp.turnover, 
    emp.mean.compo.dist)
  names(empirical.values) <- c("mean.site.richness", "regional.richness", "q0.alpha", "q0.beta", 
    "q0.gamma", "q1.alpha", "q1.beta", "q1.gamma", "prop.spp.turnover", "mean.compo.dist")

  
  ## Calculates the number of null matrices
  null.N <- length(null.matrices)

  ## Creates empty objects to hold results from the null matrices
  rand.details <- as.data.frame(matrix(NA, ncol=length(empirical.values), nrow=null.N))
  colnames(rand.details) <- names(empirical.values)
  rownames(rand.details) <- paste("null_", 1:null.N, sep="")
  
  for (i in 1:null.N)
  {
    if(print.progress==TRUE) 
      print(i)

    ## Defines the focal null composition matrix
    null.compo <- as.matrix(null.matrices[[i]])

    ## Calculates null regional abundances per species - the regional SAD
    null.spp.abund <- colSums(null.compo)
    if(identical(null.spp.abund, spp.abund)==FALSE) 
      warning("Simulated and observed species total abundances are not identical")
    
    ## Calculates species richness across all sites
    null.regional.richness <- length(which(null.spp.abund>0))

    ## Calculates species richness at each site
    null.site.richness <- rowSums(null.compo>0)
    null.mean.site.richness <- mean(null.site.richness)

    ## Calculates null composition distances among all possible pairs of sites
    if(dist.method=="hellinger") 
      null.dist.details <- vegdist(decostand(null.compo, method="hellinger"), 
        method="euclidean", binary=binary.data)
    else 
      null.dist.details <- vegdist(null.compo, method=dist.method, binary=binary.data)    
    
    ## Partitions diversity in the null composition data 
    null.q0.alpha <- d(null.compo, lev="alpha", q=0)
    null.q0.beta <- d(null.compo, lev="beta", q=0)
    null.q0.gamma <- d(null.compo, lev="gamma", q=0)

    null.q1.alpha <- d(null.compo, lev="alpha", q=1)
    null.q1.beta <- d(null.compo, lev="beta", q=1)
    null.q1.gamma <- d(null.compo, lev="gamma", q=1)

    null.prop.spp.turnover <- (1-(null.mean.site.richness/null.regional.richness))
      
    rand.details[i,] <- c(null.mean.site.richness, null.regional.richness, null.q0.alpha, 
      null.q0.beta, null.q0.gamma, null.q1.alpha, null.q1.beta, null.q1.gamma, 
      null.prop.spp.turnover, mean(null.dist.details))
  }

  ## Creates the summary table and outputs results
  summary.table <- summary.rand(empirical=empirical.values, rand.details=rand.details, 
    critial.quantiles=c(0.025, 0.975))
  
  output <- list(summary.table, rand.details)
  names(output) <- c("summary.table", "rand.details")

  class(output) <- "rand.test.results"
  output
}




####################################################################################################
######################################  SUPPLEMENTARY FUNCTIONS ####################################
####################################################################################################

### SUPPLEMENTARY FUNCTION 1: summary.rand #########################################################
summary.rand <- function(empirical, rand.details, critial.quantiles=c(0.025, 0.975), boot.details, 
  boot.CI.quantiles=c(0.025, 0.975))
{
	emp.names <- rownames(as.matrix(empirical))
	rand.names <- colnames(rand.details)

	if(length(emp.names)==0 | length(rand.names)==0) 
		stop("Both 'empirical' and 'rand.details' must have 'names' or 'colnames'")
	if(identical(emp.names, rand.names)==FALSE) 
		stop("Names in 'empirical' do not match names in 'rand.details'")

	if(missing(boot.details)==FALSE)
	{
		boot.names <- colnames(boot.details)
		if(length(boot.names)==0) 
			stop("Both 'empirical' and 'rand.details' must have 'names' and 'colnames'")
		if(identical(emp.names, boot.names)==FALSE) 
			stop("Names in 'empirical' do not match names in 'boot.details'")
	}
	
	empirical <- as.numeric(empirical)
	rand.details <- as.matrix(rand.details)

	if(length(empirical)!=ncol(rand.details)) 
		stop("Number of empirical values and of columns in randomization matrix do not match")

	rand.N <- apply(!is.na(rand.details), 2, sum)	
	NA.rand.N <- apply(is.na(rand.details), 2, sum)
		
	mean.rand <- colMeans(rand.details, na.rm=TRUE)
	median.rand <- apply(rand.details, 2, median, na.rm=TRUE)
	sd.rand <- apply(rand.details, 2, sd, na.rm=TRUE)

	critical.values <- as.matrix(t(apply(rand.details, 2, quantile, probs=critial.quantiles, 
    na.rm=TRUE)))
	ps <- t(sapply(1:length(empirical), rand.p, empirical=empirical, rand.details=rand.details))

	ES <- empirical - mean.rand
	SES <- ES/sd.rand
	
	which.inf <- which(sd.rand==0) 
	if(length(which.inf)>0)
	{
		warning("At least some statistics showed no variation expected by the randomization 
			alrorithm - standard deviations of randomized values equal to zero. Corresponding SES 
			values were set to NA", call.=FALSE)
		SES[which.inf] <- NA
	}
	

	summary.results <- data.frame(empirical, rand.N, NA.rand.N, mean.rand, median.rand, sd.rand, ES, 
    SES, critical.values, ps)
	colnames(summary.results) <- c("empirical", "rand.N", "NA.rand.N", "mean.rand", "median.rand", 
    "sd.rand", "ES", "SES", "lower.critical.value", "upper.critical.value", "p.less.than", 
    "p.greater.than")
	rownames(summary.results) <- emp.names
	
	if(missing(boot.details)==FALSE)
	{
		mean.boot <- colMeans(boot.details, na.rm=TRUE)
		CI.boot <- as.matrix(t(apply(boot.details, 2, quantile, probs=boot.CI.quantiles, na.rm=TRUE)))
		colnames(CI.boot) <- c("lower.CI.boot", "upper.CI.boot")
		
		summary.results <- data.frame(summary.results, mean.boot, CI.boot)
	}

	summary.results
}



### SUPPLEMENTARY FUNCTION 2: rand.p ###############################################################
rand.p <- function(i, empirical, rand.details)
{
	emp.i <- empirical[i]
	rands.i <- rand.details[,i]
	rands.i <- rands.i[which(is.na(rands.i)==FALSE)]
	
	p.less.than <- length(which(rands.i <= emp.i)) / length(rands.i)
	p.greater.than <- length(which(rands.i >= emp.i)) / length(rands.i)
		
	rand.ps <- cbind(p.less.than, p.greater.than)
	names(rand.ps) <- c("p.less.than", "p.greater.than")
		
	rand.ps
}

  
