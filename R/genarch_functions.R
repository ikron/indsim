#This file contains scripts required to run the genetic architecture simulations


#Calculate Hardy-Weinberg frequencies, p = A and q = a
HW.genot.freq <- function(p) {
    q <- 1 - p
    AA <- p^2
    Aa <- 2*p*q
    aa <- q^2
    return(c(AA, Aa, aa))
}
###################################

#Initialise the genotypes of individuals
initialise.ind.mat <- function(ind.mat, genotype.mat, N, nloc) {
    
    #Go from genotype frequencies to N genotypes for each locus
    Ngeno <- function(x,N) { round(x*N) } #genotype frequency times number of individuals, rounded
    N.genot <- apply(genotype.mat, MARGIN = 2, Ngeno, N = N)

    #Over loci
    for(j in 1:nloc) {

        #Sample genotypes
        genot <- c(rep("AA",N.genot[1,j]), rep("Aa",N.genot[2,j]), rep("aa",N.genot[3,j]))
        sampled.genot <- sample(genot, size = N, replace = F)
    
        #Over individuals
        for(i in 1:N) {
            if(sampled.genot[i] == "AA") {ind.mat[[i]][,j] <- c(1,1) }
            if(sampled.genot[i] == "Aa") {ind.mat[[i]][,j] <- c(1,0) }
            if(sampled.genot[i] == "aa") {ind.mat[[i]][,j] <- c(0,0) }
        }
    }
    return(ind.mat)
}
##############################################################################################

#Calculate genotypic effects
calc.genot.effects <- function(ind.mat, ae, type, qtl.ind, nloc) {
    if(type == "biallelic") {
    geno.calc <- function(x, ae) {sum(rbind(x[1,]*ae, x[2,]*ae))} #This function sums genotypic effects
    ind.G <- as.vector(unlist(lapply(ind.mat, geno.calc, ae = ae)))
    }

    if(type == "infinite") {
        #First put deleterious recessives to 0
        my.replace <- function(x, qtl.ind, nloc) {x[1:2, !(1:nloc %in% qtl.ind)] <- 0; return(x)} #This function prevents spurious phenotypic effects
        ind.mat <- lapply(ind.mat, my.replace, qtl.ind = qtl.ind, nloc = nloc)
        ind.G <- as.vector(unlist(lapply(lapply(ind.mat, Re), sum)))
    }
    
    return(ind.G)
}

#Calculate genotypic effects with limits (i.e genotypic effects are bounded between limits[1] and limits[2]
calc.genot.effects.lim <- function(ind.mat, ae, type, qtl.ind, nloc, limits) {
    if(type == "biallelic") {
    geno.calc <- function(x, ae) {sum(rbind(x[1,]*ae, x[2,]*ae))} #This function sums genotypic effects
    ind.G <- as.vector(unlist(lapply(ind.mat, geno.calc, ae = ae)))
    }

    if(type == "infinite") {
        #First put deleterious recessives to 0
        my.replace <- function(x, qtl.ind, nloc) {x[1:2, !(1:nloc %in% qtl.ind)] <- 0; return(x)} #This function prevents spurious phenotypic effects
        ind.mat <- lapply(ind.mat, my.replace, qtl.ind = qtl.ind, nloc = nloc)
        ind.G <- as.vector(unlist(lapply(lapply(ind.mat, Re), sum)))
    }

#Set too large or small genotypic values to their limits
    ind.G <- replace(ind.G, ind.G < limits[1], limits[1])
    ind.G <- replace(ind.G, ind.G > limits[2], limits[2])
    
    return(ind.G)
}

#This function calculates effects of epigenetically modified loci
#Does not work at the moment if allelic architecture is biallelic for qtl
calc.genot.effects.epi <- function(ind.mat, qtl.ind, nloc) {
    #Check which loci are epigenetically modified

    
    #Set all other values first to 0 to avoid spurious effects
    my.replace <- function(x, qtl.ind, nloc) {x[1:2, !(1:nloc %in% qtl.ind)] <- 0; return(x)}
    ind.mat <- lapply(ind.mat, my.replace, qtl.ind = qtl.ind, nloc = nloc)
    
    #Calculate epigenetic effects for all inds
    ind.epi <- lapply(ind.mat, Im)
    if(any(unlist(ind.epi) != 0)) {
        epi.cue <- unlist(ind.epi)[which(unlist(ind.epi) != 0)[1]] } else { epi.cue <- 0 }
    #Epigenetic effect is the slope effect of epigenetically modified slope loci * previous parental cue
    epi.check <- lapply(ind.epi, function(x) {x > 0})
    reals <- lapply(ind.mat, Re)
    ind.Gb <- mapply('*', reals, epi.check, SIMPLIFY = FALSE) #
    #Calculate slope effects for all inds
    ind.Gb <- as.vector(unlist(lapply(ind.Gb, sum))) #Calculate slopes
    ind.Gepimod <- ind.Gb*epi.cue
    return(ind.Gepimod)
}

#This functions resets all epigenetic modifications
reset.epigenetics <- function(ind.mat) {
    ind.mat <- lapply(ind.mat, Re)
    return(ind.mat)
}

##This function inserts epigenetic modifications to certain loci
        epigenetic.modification <- function(ind.mat, qtl.ind, cue, ind.Gepi) {
            modify <- function(x, loc, cue) { x[,loc] <- x[,loc] + cue*1i; return(x) }
            N.ind <- length(ind.Gepi)
            epi.mod.ind <- rbinom(n = N.ind, size = 1, prob = ind.Gepi) #Does epigenetic mod happen?
            select.inds <- ind.mat[as.logical(epi.mod.ind)] #
            select.inds <- lapply(select.inds, modify, qtl.ind, cue) #Epigenetically modify the qtl.slope loci
            ind.mat[as.logical(epi.mod.ind)] <- select.inds
            return(ind.mat)           
        }


##########################################################


#Use truncation selection to select top individuals. Returns the index of selected individuals
truncation.selection <- function(ind.mat, ind.P, sel.intensity = 1) {
    #Selecting individuals that have phenotype of i standard deviations larger than the mean
    #sel.intensity <- 1 #implies i = 1 -> S = \sigma_P
    t.point <- mean(ind.P) + sel.intensity*sd(ind.P) #Truncation point
    t.index <- ind.P >= t.point #Indexes of selected individuals
    #sel.mat <- ind.mat[t.index] #Select individuals,( its better to return just the index)
    return(t.index)
}
#########################################################


#Production and random union of gametes
produce.gametes <- function(sel.mat, loc.mat, N) {
    
    Nsel <- length(sel.mat)
    ind.mat <- rep(list("genotype" = loc.mat), N) #List for individuals, initialize again

    #Need to produce N individuals
    for(i in 1:N) {
        #Pick two individuals at random
        ind.index <- sample(1:Nsel, size = 2)
        pind.mat <- sel.mat[ind.index]

        #Produce gametes
        #Gamete 1
        ### This is for free recombination, modify this if linkage exists! ###
        gamete1 <- apply(pind.mat[[1]], MARGIN = 2, sample, size = 1)
        #Gamete 2
        gamete2 <- apply(pind.mat[[2]], MARGIN = 2, sample, size = 1)

        #Union of gametes, aka fertilization <3
        ind.mat[[i]][1,] <- gamete1
        ind.mat[[i]][2,] <- gamete2
    }
    return(ind.mat)
}

###########################################################
#Implement linkage into produce.gametes.W function, 
#Production and random union of gametes, individuals contribute to the next generation weighted by their fitness
    produce.gametes.W <- function(ind.mat, loc.mat, No, ind.W, link.map = F, linkage, nloc.chr, n.chr) {
        
    Nind <- length(ind.mat)
    ind.mat.new <- rep(list("genotype" = loc.mat), No) #List of the next generation of individuals, initialize
    #4-strand index mat
    strand.mat <- matrix(c(1,3,1,4,2,3,2,4), ncol = 4)

    #Need to produce No individuals
    for(i in 1:No) {
        #Pick two individuals at random, weighted by their fitness
        ind.index <- sample(1:Nind, size = 2, prob = ind.W)
        pind.mat <- ind.mat[ind.index]

        #Free recombination
        if(link.map == F) {
        #produce gamates
        #Gamete 1
        ### This is for free recombination, modify this if linkage exists! ###
        gamete1 <- apply(pind.mat[[1]], MARGIN = 2, sample, size = 1)
        #Gamete 2
        gamete2 <- apply(pind.mat[[2]], MARGIN = 2, sample, size = 1)
        }

        #Recombination using a linkage map
        if(link.map == T) {
        #Produce gametes
        ### Gamete 1 ###
        #Make a linkage format data from parent 1, in the four strand stage (i.e. pachytene)
        linkage.gamete1 <- convert.indmat.2.linkage(linkage, pind.mat[[1]], nloc.chr, n.chr, fourstrand = T)
        n.crossover <- gen.n.crossover(n.chr, lambda = 0.56) #This yields 1.56 crossovers per chromosome on avg.
        #Generate crossover position and pick strands that participate
        foo <- list(0)
        xo.positions.list <- rep(foo, n.chr)
        xo.strands.list <- rep(foo, n.chr)
        for(j in 1:n.chr) {
            xo.positions.list[[j]] <- sort(runif(n.crossover[j], min = 0, max = 1))
            xo.strands.list[[j]] <- sample(c(1,2,3,4), size = n.crossover[j], replace = T)
        }
        #Resolving crossovers
        linkage.gamete1 <- resolve.crossovers(meiocyte = linkage.gamete1, xo.positions.list, xo.strands.list, strand.mat, n.chr, n.crossover)
        #Select chromosomes for gametes (i.e. meiotic divisions)
        gamete1 <- make.gamete(linkage.gamete1)

        ### Gamete2 ###
        linkage.gamete2 <- convert.indmat.2.linkage(linkage, pind.mat[[2]], nloc.chr, n.chr, fourstrand = T)
        n.crossover <- gen.n.crossover(n.chr, lambda = 0.56)
        xo.positions.list <- rep(foo, n.chr)
        xo.strands.list <- rep(foo, n.chr)
        for(j in 1:n.chr) {
            xo.positions.list[[j]] <- sort(runif(n.crossover[j], min = 0, max = 1))
            xo.strands.list[[j]] <- sample(c(1,2,3,4), size = n.crossover[j], replace = T)
        }
        #Resolving crossovers and perform meiotic divisions
        linkage.gamete2 <- resolve.crossovers(meiocyte = linkage.gamete2, xo.positions.list, xo.strands.list, strand.mat, n.chr, n.crossover)
        gamete2 <- make.gamete(linkage.gamete2)
        } #Done making gametes

        #Union of gametes, aka fertilization <3
        ind.mat.new[[i]][1,] <- gamete1
        ind.mat.new[[i]][2,] <- gamete2
    } #No individuals produced
        
    return(ind.mat.new)
}
#####################################################################################33


###########################################################
#This function produces gametes on individual basis, 
#Production and random union of gametes, individuals contribute to the next generation weighted by their fitness
    produce.gametes.W.No <- function(ind.mat, loc.mat, No, ind.W, link.map = F, linkage, nloc.chr, n.chr) {
        
    Nind <- length(ind.mat) #Number of parents
    sum.No <- sum(No) #Total number of individuals to be produced
    ind.mat.new <- rep(list("genotype" = loc.mat), sum.No) #List of the next generation of individuals, initialize
    #4-strand index mat
    strand.mat <- matrix(c(1,3,1,4,2,3,2,4), ncol = 4)
  
    n <- 1 #Initialize counter

    #Need to produce No individuals
    for(i in 1:Nind) { #Loop over all individuals (parents)
        if(No[i] > 0) { #Check that parent had more than zero offspring
            for(k in 1:No[i]) { #Loop over all offspring of focal individual
                #Pick the focal individual and another at random weighted by their fitness
                ind.index <- c(i, sample((1:Nind)[-i], size = 1, prob = ind.W[-i])) #Selfing not possible
                pind.mat <- ind.mat[ind.index]

                #Free recombination
                if(link.map == F) {
                #produce gamates
                #Gamete 1
                ### This is for free recombination, modify this if linkage exists! ###
                gamete1 <- apply(pind.mat[[1]], MARGIN = 2, sample, size = 1)
                #Gamete 2
                gamete2 <- apply(pind.mat[[2]], MARGIN = 2, sample, size = 1)
                }

                #Recombination using a linkage map
                if(link.map == T) {
                #Produce gametes
                ### Gamete 1 ###
                #Make a linkage format data from parent 1, in the four strand stage (i.e. pachytene)
                linkage.gamete1 <- convert.indmat.2.linkage(linkage, pind.mat[[1]], nloc.chr, n.chr, fourstrand = T)
                n.crossover <- gen.n.crossover(n.chr, lambda = 0.56) #This yields 1.56 crossovers per chromosome on avg.
                #Generate crossover position and pick strands that participate
                foo <- list(0)
                xo.positions.list <- rep(foo, n.chr)
                xo.strands.list <- rep(foo, n.chr)
                for(j in 1:n.chr) {
                    xo.positions.list[[j]] <- sort(runif(n.crossover[j], min = 0, max = 1))
                    xo.strands.list[[j]] <- sample(c(1,2,3,4), size = n.crossover[j], replace = T)
                }
                #Resolving crossovers
                linkage.gamete1 <- resolve.crossovers(meiocyte = linkage.gamete1, xo.positions.list, xo.strands.list, strand.mat, n.chr, n.crossover)
                #Select chromosomes for gametes (i.e. meiotic divisions)
                gamete1 <- make.gamete(linkage.gamete1)

                ### Gamete2 ###
                linkage.gamete2 <- convert.indmat.2.linkage(linkage, pind.mat[[2]], nloc.chr, n.chr, fourstrand = T)
                n.crossover <- gen.n.crossover(n.chr, lambda = 0.56)
                xo.positions.list <- rep(foo, n.chr)
                xo.strands.list <- rep(foo, n.chr)
                for(j in 1:n.chr) {
                    xo.positions.list[[j]] <- sort(runif(n.crossover[j], min = 0, max = 1))
                    xo.strands.list[[j]] <- sample(c(1,2,3,4), size = n.crossover[j], replace = T)
                }
                #Resolving crossovers and perform meiotic divisions
                linkage.gamete2 <- resolve.crossovers(meiocyte = linkage.gamete2, xo.positions.list, xo.strands.list, strand.mat, n.chr, n.crossover)
                gamete2 <- make.gamete(linkage.gamete2)
                } #Done making gametes

                #Union of gametes, aka fertilization <3
                ind.mat.new[[n]][1,] <- gamete1
                ind.mat.new[[n]][2,] <- gamete2
                n <- n + 1 #update counter        
            } #Close looping over all offspring of one parent
        } #Close if statement
    } #Close looping over all parents, #No individuals produced
        
    return(ind.mat.new)
}
#####################################################################################33



#####################################################################################################
### Functions for recombination #####################################################################
#Function for going from ind.mat format to linkage format, input: linkage, ind.mat and vector of loci numbers per chromosome, and number of chromosomes
#Assuming that loci are in linear order in ind.mat so that consecutive numbering from 1 to n
convert.indmat.2.linkage <- function(linkage, ind.mat, nloc.chr, n.chr, fourstrand = F) {

    #Loop over linkage groups and copy genotypes from ind.mat
    pos1 <- 1 #Initial starting position
    for(i in 1:n.chr) {
        #pos1 <- 1
        pos2 <- pos1 + nloc.chr[i] - 1
        linkage[[i]][1:2,] <- ind.mat[,pos1:pos2]
        pos1 <- pos2 + 1 #update new starting position
    }

    #Do we want four strand stage? If so, duplicate first and second rows
    if(fourstrand == T) {
        for(i in 1:n.chr) {
            linkage[[i]] <- rbind(linkage[[i]][1,], linkage[[i]][1,], linkage[[i]][2,], linkage[[i]][2,], linkage[[i]][3,]) }
    }
          
    return(linkage)
}

#Function for going to from linkage format to ind.mat format
#This is an elegant solution, fully vectorized!
convert.linkage.2.indmat <- function(linkage) {

    mydrop <- function(x) {x[-3,]} #Function for dropping the third row
    ind.mat <- matrix(unlist(lapply(linkage, mydrop)), nrow = 2) #unlist and rebuild matrix
    return(ind.mat)
}

#Generate crossovers
#Function that simulates crossover numbers per chromosome here
gen.n.crossover <- function(n.chr, lambda, lambda.chr = FALSE) {
    #Number of crossovers per chromosome is poisson distributed (minimum is one xo)
    n.xo <- rpois(n.chr, lambda)+1
    return(n.xo)
}

resolve.crossovers <- function(meiocyte, xo.positions.list, xo.strands.list, strand.mat, n.chr, n.crossover) {
    
    #Resolving crossovers
    #Loop over chromosomes
    for(i in 1:n.chr) {
        #Loop over crossover events
        for(j in 1:n.crossover[i]) {
            #Select segments
            cur.pos <- xo.positions.list[[i]][j] #Current position
            cur.strands <- xo.strands.list[[i]][j]
            sel.strands <- meiocyte[[i]][strand.mat[,cur.strands],] #Select strands
            #Crossover indices
            before.xo <- Re(meiocyte[[i]][5,]) < cur.pos
            after.xo <- !before.xo
            #Resolve crossover
            meiocyte[[i]][(strand.mat[,cur.strands][1]),] <- c(sel.strands[1,before.xo], sel.strands[2,after.xo])
            meiocyte[[i]][(strand.mat[,cur.strands][2]),] <- c(sel.strands[2,before.xo], sel.strands[1,after.xo])
        }
    }
    return(meiocyte)
}

#This function converts the four strand format to a gamete, by randomly sampling one chromatid
#returns a gamete in ind.mat format
make.gamete <- function(linkage.gamete) {
    mydrop <- function(x) {x[-5,]} #For dropping the fifth row
    mysample <- function(x) {x[sample(c(1,2,3,4),1),] } #Sampling one the four chromatids
    gamete <- lapply(lapply(linkage.gamete, mydrop), mysample) #Using list apply
    return(unlist(gamete))
}

##################################################################################################


### Functions for generating mutations ###########################################################
#Do mutations occur or not?
gen.mutations <- function(x, N, mu) {
    x <- rbinom(N, size = 1, prob = mu) }

#This function makes mutations with the K-alleles model
        mutate.K.alleles <- function(x, l, alleles) {
            mut.al.ind <- sample(c(1,2),1)
            cur.allele <- x[mut.al.ind,l]
            new.allele <- sample(alleles[[l]][alleles[[l]] !=cur.allele],1) #Select one allele
            x[mut.al.ind,l] <- new.allele
            return(x)
        }

#This function makes mutations with the infinite alleles model
mutate.inf.alleles <- function(x, l) {
    mut.al.ind <- sample(c(1,2),1)
    x[mut.al.ind,l] <- rnorm(1, mean = 0, sd = 1) #This may need to changed accordingly
    return(x)
}

#################################################################################################


#Some additional functions
#Calculate allele frequencies
#For infinite allele model, frequencies of qtls are calculated such that ancestral allele (i.e. 0) is q
calc.al.freq <- function(data, N, nloc, allele.model, loc.attributes) {
    
    allele.freq <- rep(0, nloc) #initialise allele frequencies

    #First the case of all biallelic calculation
    if(allele.model == "biallelic") {

        for(j in 1:nloc) {
            allele1 <- unname(sapply(data, '[[', 2*j-1)) #Get element [1,2*j-1] from each individual
            allele2 <- unname(sapply(data, '[[', 2*j)) #Get element [2,2*j] from each individual

             allele.freq[j] <- (sum(allele1+allele2))/(2*N)
        }
    }

    #Then in the case of infinite allele calculations
    if(allele.model == "infinite") {
            
        for(j in 1:nloc) {
            if(loc.attributes[j] %in% c("qtl", "qtl.slope", "qtl.adj", "qtl.hed", "qtl.epi") == FALSE) { #Frequencies for deleterious recessives and neutral loci
                allele1 <- unname(sapply(data, '[[', 2*j-1)) #Get element [1,2*j-1] from each individual
                allele2 <- unname(sapply(data, '[[', 2*j)) #Get element [2,2*j] from each individual
                
                allele.freq[j] <- (sum(allele1+allele2))/(2*N)
            } else {
                allele1 <- unname(sapply(data, '[[', 2*j-1)) #Get element [1,2*j-1] from each individual
                allele2 <- unname(sapply(data, '[[', 2*j)) #Get element [2,2*j] from each individual
                allele1 <- allele1 != 0 #All other values (alleles) than 0 get a value of 1
                allele2 <- allele2 != 0

                allele.freq[j] <- (sum(allele1+allele2))/(2*N)
            }
        }
    }
    return(allele.freq)
}
########################################################################################

#########################################################
#Calculate effective population size from neutral markers
#Assumes that generations are sampled t generations apart
calc.Ne <- function(af1, af2, N1, N2, nloc, t) {
    ##Using the method of Nei and Tajima 1981 (methoc C: eqs 15 and 16)
    Fc <- rep(0, nloc)
    #Calculating Fc (equation 15) (A measure of inbreeding F-statistics)
    for(i in 1:nloc) {
        Fc[i] <- (1/2) * ((af1[i] - af2[i])^2 / ( (af1[i]+af2[i])/2 - af1[i]*af2[i]) + ((1-af1[i]) - (1-af2[i]))^2 / ( ((1-af1[i]) + (1-af2[i]) )/2 - (1-af1[i])*(1-af2[i]) ))
    }
    #Mean Fc across all loci (since loci are biallelic no weighting is done for number of alleles)
    meanFc <- sum(Fc)/nloc
    #Estimate of Ne (equation 16)
    Ne <- (t-2) / (2*(meanFc - (1/(2*N1) + 1/(2*N2))))
    return(Ne)
}


#########################################################


#Calculate fitness using stabilizing selection 
    calc.fitness <- function(ind.P, opt.pheno, sel.intensity) {
        exp(-((ind.P-opt.pheno)^2)/(2*sel.intensity^2))
    }

#Calculate fitness effects of deleterious mutations
delres.fitness <- function(ind.mat, delres.ef, delres.ind, nloc) {
    mydselect <- function(x, delres.ind, nloc) {x[1:2, (1:nloc %in% delres.ind)]}
    ind.mat <- lapply(ind.mat, mydselect, delres.ind, nloc) #Select only deleterious recessives

    #Need to define some helper function to look if both alleles contain deleterious recessive
    delres.check <- function(ind.mat) {all(ind.mat == 1)}
    delres.check2 <- function(ind.mat) { apply(ind.mat, MARGIN = 2, delres.check) }
    delr.mat <- lapply(ind.mat, delres.check2) #Check all genotypes
    delr.vec <- as.vector(unlist(lapply(delr.mat, sum)))
    return((1-delres.ef)^delr.vec) #Calculate effects on fitness
}
 
#################################################################

#The string #' indicates roxygen parsing

#### Wrapper function for the whole simulation
#' Perform individual based simulations
#'
#' This function performs individual based simulations according to given parameters
#'
#' @param N Starting population size
#' @param generations Number of generations to run the simulations
#' @param sel.intensity Selection intensity for stabilising selection
#' @param init.f Initial allele frequencies for QTL
#' @param init.n Initial allele frequencies for neutral loci
#' @param sigma.e Environmental standard deviations of the phenotypic trait
#' @param K Carrying capacity of the population
#' @param opt.pheno Vector of phenotypic optimums for each generation
#' @param density.reg How population density is regulated. Currently there are two possible values:  "logistic" which implements logistic population density regulation (see manual) and "sugar" which constrains population to carrying capacity by randomly removing individuals when population size is higher than K
#' @param r If density regulation is "logistic", this is the population growth rate.
#' @param B If density regulation is "sugar", this it the mean number of offspring each perfectly adapted individual has
#' @param allele.model Which allele model to use for qtl. Currently the options are: "infinite" which uses an infinite alleles model for qtl, and "biallelic" which uses a biallelic model for qtl
#' @param nqtl Number of loci affecting the phenotype
#' @param mu Mutation rate, which currently is the same for all loci
#' @param n.chr Number of chromosomes (linkage groups)
#' @param nloc.chr Vector of the length of n.chr giving number of loci for each chromosome
#' @param n.delres Number of deleterious recessives, defaults to 0
#' @param delres.ef Effects of deleterious recessives on fitness (these are multiplicative), for example: setting this parameter to 0.1 means that fitness of an individuals that is homozygous for the allele is 1 - 0.1 = 0.9
#' @param n.neutral Number of neutral loci, defaults to 0
#' @param linkage.map If map distances are determined randomly, set this to "random", if fixed linkage map is provided can set this to "user"
#' @param linkage Matrix of 3 rows and nloc columns, first two rows can be zeros and third row determines map distances relative to chromosome coordinates from 0 to 1
#' @export
indsim.simulate <- function(N, generations, sel.intensity, init.f, init.n, a, sigma.e, K, r, B, opt.pheno, density.reg, allele.model, mu, n.chr, nloc.chr, nqtl, n.delres = 0, delres.ef = NULL, n.neutral = 0, linkage.map = "random", linkage = NULL) {

    #Prepare linkage groups
    foo <- list(0)
    nloc <- nqtl + n.delres + n.neutral

    #Perform some checks that initial parameters are sensible
    if(nloc != sum(nloc.chr)) { stop("Number of loci and loci per chromosome don't match!") }

    loc.attributes <- sample(c(rep("qtl", nqtl), rep("neutral", n.neutral), rep("delres", n.delres)), size = nloc, replace = FALSE) #Types for all loci
    delres.ind <- (1:nloc)[loc.attributes == "delres"] #Indices of deleterious recessives
    qtl.ind <- (1:nloc)[loc.attributes == "qtl"] #Indices of loci that affect the phenotype

    alleles <- rep(list(c(0,1)), nloc) #Alleles for all except qtl
    
    if(allele.model == "biallelic") { #Setup biallelic effects
        #alleles <- rep(list(c(0,1)), nloc) #All loci biallelic
        q <- (1:(nqtl+1)/(nqtl+1))[-(nqtl+1)]
        ae <- qnorm(q, mean = 0, sd = 1) #Allelic effects are distributed normally
        #ae <- dnorm(seq(from = -5, to = -1, length.out = nloc), mean = 0, sd = 1)*10
        ae <- sample(ae, length(ae)) #Randomize allelic effects among loci
        temp <- rep(0,nloc)
        temp[qtl.ind] <- ae
        ae <- temp }

    #Setup linkage map
    if(linkage.map == "random") {
        linkage <- rep(foo, n.chr)
        for(i in 1:n.chr) {
            linkage[[i]] <- matrix(rep(0,3*nloc.chr[i]), ncol = nloc.chr[i])
            linkage[[i]][3,] <- sort(runif(nloc.chr[i])) #Map-positions for each locus
        }
    }

    #Simulation with logistic population regulation and natural selection

    #Initialize the results matrix
    #We want to monitor phenotypes, allele freqs, heritabilities, variances etc. etc.
    results.mat.pheno <- matrix(rep(0, generations*7), ncol = 7)
    colnames(results.mat.pheno) <- c("generation", "pop.mean", "sel.mean", "var.a", "h2", "W.mean", "N")
    results.mat.pheno[,1] <- 1:generations #Store generation numbers

    results.mat.alleles <- matrix(rep(0, generations*(nloc+1)), ncol = nloc+1)
    results.mat.alleles[,1] <- 1:generations
    colnames(results.mat.alleles) <- c("generation", c(paste(rep("locus", nloc), 1:nloc, sep = "")))
    #If there are neutral loci in the simulation add Ne column to calculate variance effective size
    if(n.neutral > 0) { results.mat.alleles <- cbind(results.mat.alleles, c(rep(NA,10), rep(0, generations-10)))
                        colnames(results.mat.alleles)[nloc+2] <- "Ne"
                        neutral.index <- (1:nloc)[loc.attributes == "neutral"] + 1 }

    #Initilize the simulation
    loc.mat <- matrix(rep(0,2*nloc), ncol = nloc)
    ind.mat <- rep(list("genotype" = loc.mat), N) #List for individuals

    #Allele frequencies
    init.freq <- rep(init.f,nloc) #for p
    #Here can also make explicit initial frequencies for qtls, delres and neutral markers
    if(n.neutral > 0) { init.freq[loc.attributes == "neutral"] <- init.n } #Neutral markers start at these frequencies

    #Genotype frequencies
    genotype.mat <- matrix(rep(0,3*nloc), ncol = nloc)
    rownames(genotype.mat) <- c("AA", "Aa", "aa")

    #Initialize genotype frequencies
    for(i in 1:nloc) {
        genotype.mat[,i] <- HW.genot.freq(init.freq[i])
    }
    #################################

    ind.mat <- initialise.ind.mat(ind.mat, genotype.mat, N, nloc) #Sample starting genotypes


    #Loop over generations
    for(g in 1:generations) {

        #Check for extinction
        if(N < 2) {stop("Population went extinct! : (")}

        #Calculate allele frequencies
        results.mat.alleles[g,2:(nloc+1)] <- calc.al.freq(ind.mat, length(ind.mat), nloc, allele.model, loc.attributes)

        

        #Argument type needs to have value of either "biallelic" or "infinite"
        ind.G <- calc.genot.effects(ind.mat, ae, type = allele.model, qtl.ind, nloc) #Calculate genotype effects
        #sigma.e <- (abs(ind.G) + 1)*coef.e #Environmental variation scales with genotypic effects
        
        
        #Calculate environmental and phenotypic effects
        ind.E <- rnorm(n = N, mean = 0, sd = sigma.e)
        ind.P <- ind.G + ind.E
        ###############################################

        #Store phenotypes    
        results.mat.pheno[g,2] <- mean(ind.P) #Trait mean
        results.mat.pheno[g,4] <- var(ind.G) #Additive genetic variance
        results.mat.pheno[g,5] <- var(ind.G)/var(ind.P) #Heritability
        
        #Calculate fitnesses
        ind.W <- calc.fitness(ind.P, opt.pheno[g], sel.intensity)
        if(n.delres > 0) {
            #Calculate fitness effects of deleterious recessives
            ind.W <- ind.W*delres.fitness(ind.mat, delres.ef, delres.ind, nloc)
        }
        

        results.mat.pheno[g,6] <- mean(ind.W) #Population mean fitness
        results.mat.pheno[g,7] <- N #Population size

        #Calculate effective population size (if g > 10)
        if(g > 10) { results.mat.alleles[g,(nloc+2)] <- calc.Ne(af1 = results.mat.alleles[(g-10), neutral.index], af2 = results.mat.alleles[g, neutral.index], N1 = results.mat.pheno[g-10,7], N2 = results.mat.pheno[g,7], nloc = n.neutral, t = 10) }
        
        ### Reproduction ###
        #Using logistic density regulation
        if(density.reg == "logistic") {
            #Calculate the number of offspring the population produces
            #Logistic population regulation
            No <- round(N*K*(1+r) / (N*(1+r) - N + K)) #Logistic population growth
            #Scale number of offspring by population mean fitness
            No <- round(mean(ind.W)*No); if(No < 1) {stop("Population went extinct! : (")}
            ind.mat <- produce.gametes.W(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
        }

        #Population regulation via Mikael's method
        if(density.reg == "sugar") {
            #Calculate the number of offspring that each individual produces
            No <- rpois(n = N, lambda = B*ind.W)
            #If number of offspring larger than K, remove some offspring randomly
            if(sum(No) > K) {
                while(sum(No) > K) { #Remove individuals until only K are left, seems slow...
                    no.ind <- (1:length(No))[No > 0] #Indexes of cases where No > 0
                    rem.ind <- sample(no.ind, 1)
                    No[rem.ind] <- No[rem.ind] - 1
                }
            }
            ind.mat <- produce.gametes.W.No(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
        }
        
        #Next generation
        N <- length(ind.mat) #Store new population size

        ### Mutation ###
        #Produce mutations
        mutations <- rep(foo, nloc)
        #Has to be 2*mu since mutations happen on individual basis and each ind has 2 copies
        mutations <- lapply(mutations, gen.mutations, N = N, mu = 2*mu) 
        for(l in 1:nloc) { #Loop over all loci
            if(any(mutations[[l]] == 1) == TRUE) { #Did any mutations happen?
                mut.ind <- ind.mat[mutations[[l]] == 1] #Select individuals to be mutated
                if(allele.model == "biallelic") { #Which alleles model is used?
                    mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) #Mutate alleles
                }
                if(allele.model == "infinite") {
                    if(loc.attributes[l] == "qtl") {
                        mut.ind <- lapply(mut.ind, mutate.inf.alleles, l = l) } else {
                            mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) }
                }
                ind.mat[mutations[[l]] == 1] <- mut.ind #Store results
            }
        }
    ########################    
        
        
    } #Done looping over generations


    #Modify results matrices if necessary

    #Return results
    return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes))
}
    

###This is a function that allows the evolution of phenotypic plasticity
#Need to add roxygen stuff to this

indsim.plasticity.simulate <- function(N, generations, sel.intensity, init.f, init.n, a, sigma.e, K, r, B, opt.pheno, env.cue, density.reg, allele.model, mu, n.chr, nloc.chr, nqtl, nqtl.slope, n.delres = 0, delres.ef = NULL, n.neutral = 0, linkage.map = "random", linkage = NULL) {

    #Prepare linkage groups
    foo <- list(0)
    nloc <- nqtl + n.delres + n.neutral + nqtl.slope

    #Perform some checks that initial parameters are sensible
    if(nloc != sum(nloc.chr)) { stop("Number of loci and loci per chromosome don't match!") }

    loc.attributes <- sample(c(rep("qtl", nqtl), rep("neutral", n.neutral), rep("delres", n.delres), rep("qtl.slope", nqtl.slope)), size = nloc, replace = FALSE) #Types for all loci
    delres.ind <- (1:nloc)[loc.attributes == "delres"] #Indices of deleterious recessives
    qtl.ind <- (1:nloc)[loc.attributes == "qtl"] #Indices of loci that affect the phenotype (via intercept)
    qtl.slope.ind <- (1:nloc)[loc.attributes == "qtl.slope"] #Indices for loci that affect the phenotype (via reaction norm slope)

    alleles <- rep(list(c(0,1)), nloc) #Alleles for all except qtl
    
    if(allele.model == "biallelic") { #Setup biallelic effects
        #alleles <- rep(list(c(0,1)), nloc) #All loci biallelic
        q <- (1:(nqtl+1)/(nqtl+1))[-(nqtl+1)]
        ae <- qnorm(q, mean = 0, sd = 1) #Allelic effects are distributed normally
        #ae <- dnorm(seq(from = -5, to = -1, length.out = nloc), mean = 0, sd = 1)*10
        ae <- sample(ae, length(ae)) #Randomize allelic effects among loci
        temp <- rep(0,nloc)
        temp[qtl.ind] <- ae
        ae <- temp }

    #Setup linkage map
    if(linkage.map == "random") {
        linkage <- rep(foo, n.chr)
        for(i in 1:n.chr) {
            linkage[[i]] <- matrix(rep(0,3*nloc.chr[i]), ncol = nloc.chr[i])
            linkage[[i]][3,] <- sort(runif(nloc.chr[i])) #Map-positions for each locus
        }
    }

    #Simulation with logistic population regulation and natural selection

    #Initialize the results matrix
    #We want to monitor phenotypes, allele freqs, heritabilities, variances etc. etc.
    #Monitoring also population mean intercept and slope (at genotypic values)
    results.mat.pheno <- matrix(rep(0, generations*9), ncol = 9)
    colnames(results.mat.pheno) <- c("generation", "pop.mean", "sel.mean", "var.a", "h2", "W.mean", "N", "G.int", "G.slope")
    results.mat.pheno[,1] <- 1:generations #Store generation numbers

    results.mat.alleles <- matrix(rep(0, generations*(nloc+1)), ncol = nloc+1)
    results.mat.alleles[,1] <- 1:generations
    colnames(results.mat.alleles) <- c("generation", c(paste(rep("locus", nloc), 1:nloc, sep = "")))

    #Initilize the simulation
    loc.mat <- matrix(rep(0,2*nloc), ncol = nloc)
    ind.mat <- rep(list("genotype" = loc.mat), N) #List for individuals

    #Allele frequencies
    init.freq <- rep(init.f,nloc) #for p
    #Here can also make explicit initial frequencies for qtls, delres and neutral markers
    if(n.neutral > 0) { init.freq[loc.attributes == "neutral"] <- init.n } #Neutral markers start at these frequencies

    #Genotype frequencies
    genotype.mat <- matrix(rep(0,3*nloc), ncol = nloc)
    rownames(genotype.mat) <- c("AA", "Aa", "aa")

    #Initialize genotype frequencies
    for(i in 1:nloc) {
        genotype.mat[,i] <- HW.genot.freq(init.freq[i])
    }
    #################################

    ind.mat <- initialise.ind.mat(ind.mat, genotype.mat, N, nloc) #Sample starting genotypes


    #Loop over generations
    for(g in 1:generations) {

        #Check for extinction
        if(N < 2) {stop("Population went extinct! : (")}

        #Calculate allele frequencies
        results.mat.alleles[g,2:(nloc+1)] <- calc.al.freq(ind.mat, length(ind.mat), nloc, allele.model, loc.attributes)
        

        #Argument type needs to have value of either "biallelic" or "infinite"
        #Calculate genotypic effects for intercept effects
        ind.Ga <- calc.genot.effects(ind.mat, ae, type = allele.model, qtl.ind, nloc) #Calculate genotype effects
        #Calculate genotypic effects for slope effects
        ind.Gb <- calc.genot.effects(ind.mat, ae, type = allele.model, qtl.slope.ind, nloc)
        #sigma.e <- (abs(ind.G) + 1)*coef.e #Environmental variation scales with genotypic effects
        #Calculate total genotypic value in the current environment (based on env.cue)
        ind.G <- ind.Ga + ind.Gb*env.cue[g]
        
        #Calculate environmental and phenotypic effects
        ind.E <- rnorm(n = N, mean = 0, sd = sigma.e)
        ind.P <- ind.G + ind.E
        ###############################################

        #Store phenotypes    
        results.mat.pheno[g,2] <- mean(ind.P) #Trait mean
        results.mat.pheno[g,4] <- var(ind.G) #Additive genetic variance
        results.mat.pheno[g,5] <- var(ind.G)/var(ind.P) #Heritability
        results.mat.pheno[g,8] <- mean(ind.Ga) #Mean intercept genotypic value
        results.mat.pheno[g,9] <- mean(ind.Gb) #Mean slope genotypic value
        
        #Calculate fitnesses
        ind.W <- calc.fitness(ind.P, opt.pheno[g], sel.intensity)
        if(n.delres > 0) {
            #Calculate fitness effects of deleterious recessives
            ind.W <- ind.W*delres.fitness(ind.mat, delres.ef, delres.ind, nloc)
        }
        

        results.mat.pheno[g,6] <- mean(ind.W) #Population mean fitness
        results.mat.pheno[g,7] <- N #Population size
        
        ### Reproduction ###
        #Using logistic density regulation
        if(density.reg == "logistic") {
            #Calculate the number of offspring the population produces
            #Logistic population regulation
            No <- round(N*K*(1+r) / (N*(1+r) - N + K)) #Logistic population growth
            #Scale number of offspring by population mean fitness
            No <- round(mean(ind.W)*No); if(No < 1) {stop("Population went extinct! : (")}
            ind.mat <- produce.gametes.W(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
        }

        #Population regulation via Mikael's method
        if(density.reg == "sugar") {
            #Calculate the number of offspring that each individual produces
            No <- rpois(n = N, lambda = B*ind.W)
            #If number of offspring larger than K, remove some offspring randomly
            if(sum(No) > K) {
                while(sum(No) > K) { #Remove individuals until only K are left, seems slow...
                    no.ind <- (1:length(No))[No > 0] #Indexes of cases where No > 0
                    rem.ind <- sample(no.ind, 1)
                    No[rem.ind] <- No[rem.ind] - 1
                }
            }
            ind.mat <- produce.gametes.W.No(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
        }
        
        #Next generation
        N <- length(ind.mat) #Store new population size

        ### Mutation ###
        #Produce mutations
        mutations <- rep(foo, nloc)
        #Has to be 2*mu since mutations happen on individual basis and each ind has 2 copies
        mutations <- lapply(mutations, gen.mutations, N = N, mu = 2*mu) 
        for(l in 1:nloc) { #Loop over all loci
            if(any(mutations[[l]] == 1) == TRUE) { #Did any mutations happen?
                mut.ind <- ind.mat[mutations[[l]] == 1] #Select individuals to be mutated
                if(allele.model == "biallelic") { #Which alleles model is used?
                    mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) #Mutate alleles
                }
                if(allele.model == "infinite") {
                    if(loc.attributes[l] == "qtl" | loc.attributes[l] == "qtl.slope") {
                        mut.ind <- lapply(mut.ind, mutate.inf.alleles, l = l) } else {
                            mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) }
                }
                ind.mat[mutations[[l]] == 1] <- mut.ind #Store results
            }
        }
    ########################    
        
        
    } #Done looping over generations


    #Modify results matrices if necessary

    #Return results
    return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes))
}

#####################################################################################################
##################################################################
##Plasticity simulations using a model similar to Botero et al.###
##################################################################

#Need to incorporate between generation effects, mediated by epigenetics...
#roxygen stuff
indsim.plasticity2.simulate <- function(N, generations, L, sel.intensity, init.f, init.n, a, sigma.e, K, r, B, R, P, density.reg, allele.model, mu, n.chr, nloc.chr, nqtl, nqtl.slope, nqtl.adj, nqtl.hed, nqtl.epi, n.neutral = 0, linkage.map = "random", linkage = NULL, kd, ka, ke) {

    #Initialize time
    timesteps <- generations*L #There are ngenerations*L time steps in the model
    t <- 1:timesteps

    E <- sin((2*pi*t)/(L*R)) #Initialize the environment
    cues <- rnorm(n = timesteps, mean = E*P, sd = (1-P)/3) #Environmental cues
    
    #Prepare linkage groups
    foo <- list(0)
    nloc <- nqtl + n.neutral + nqtl.slope + nqtl.adj + nqtl.hed + nqtl.epi

    #Perform some checks that initial parameters are sensible
    if(nloc != sum(nloc.chr)) { stop("Number of loci and loci per chromosome don't match!") }

    loc.attributes <- sample(c(rep("qtl", nqtl), rep("neutral", n.neutral), rep("qtl.slope", nqtl.slope), rep("qtl.adj", nqtl.adj), rep("qtl.hed", nqtl.hed), rep("qtl.epi", nqtl.epi)), size = nloc, replace = FALSE) #Types for all loci
    qtl.ind <- (1:nloc)[loc.attributes == "qtl"] #Indices of loci that affect the phenotype (via intercept)
    qtl.slope.ind <- (1:nloc)[loc.attributes == "qtl.slope"] #Indices for loci that affect the phenotype (via reaction norm slope)
    qtl.adj.ind <- (1:nloc)[loc.attributes == "qtl.adj"] #Indices for loci for plasticity adjustment
    qtl.hed.ind <- (1:nloc)[loc.attributes == "qtl.hed"] #Indices for loci for canalization / hedging
    qtl.epi.ind <- (1:nloc)[loc.attributes == "qtl.epi"] #Indices for loci for epigenetic modification

    alleles <- rep(list(c(0,1)), nloc) #Alleles for all except qtl
    
    if(allele.model == "biallelic") { #Setup biallelic effects
        #alleles <- rep(list(c(0,1)), nloc) #All loci biallelic
        q <- (1:(nqtl+1)/(nqtl+1))[-(nqtl+1)]
        ae <- qnorm(q, mean = 0, sd = 1) #Allelic effects are distributed normally
        #ae <- dnorm(seq(from = -5, to = -1, length.out = nloc), mean = 0, sd = 1)*10
        ae <- sample(ae, length(ae)) #Randomize allelic effects among loci
        temp <- rep(0,nloc)
        temp[qtl.ind] <- ae
        ae <- temp }

    #Setup linkage map
    if(linkage.map == "random") {
        linkage <- rep(foo, n.chr)
        for(i in 1:n.chr) {
            linkage[[i]] <- matrix(rep(0,3*nloc.chr[i]), ncol = nloc.chr[i])
            linkage[[i]][3,] <- sort(runif(nloc.chr[i])) #Map-positions for each locus
        }
    }

    #Simulation with logistic population regulation and natural selection

    #Initialize the results matrix
    #We want to monitor phenotypes, allele freqs, heritabilities, variances etc. etc.
    #Monitoring also population mean intercept and slope (at genotypic values)
    results.mat.pheno <- matrix(rep(0, generations*10), ncol = 10)
    colnames(results.mat.pheno) <- c("generation", "pop.mean", "var.a", "W.mean", "N", "G.int", "G.slope", "G.adj", "G.hed", "G.epi")
    results.mat.pheno[,1] <- 1:generations #Store generation numbers

    results.mat.alleles <- matrix(rep(0, generations*(nloc+1)), ncol = nloc+1)
    results.mat.alleles[,1] <- 1:generations
    colnames(results.mat.alleles) <- c("generation", c(paste(rep("locus", nloc), 1:nloc, sep = "")))

    #Initilize the simulation
    loc.mat <- matrix(rep(0,2*nloc), ncol = nloc)
    ind.mat <- rep(list("genotype" = loc.mat), N) #List for individuals

    #Allele frequencies
    init.freq <- rep(init.f,nloc) #for p
    #Here can also make explicit initial frequencies for qtls, delres and neutral markers
    if(n.neutral > 0) { init.freq[loc.attributes == "neutral"] <- init.n } #Neutral markers start at these frequencies

    #Genotype frequencies
    genotype.mat <- matrix(rep(0,3*nloc), ncol = nloc)
    rownames(genotype.mat) <- c("AA", "Aa", "aa")

    #Initialize genotype frequencies
    for(i in 1:nloc) {
        genotype.mat[,i] <- HW.genot.freq(init.freq[i])
    }
    #################################

    ind.mat <- initialise.ind.mat(ind.mat, genotype.mat, N, nloc) #Sample starting genotypes


    
    #Loop over generations
    for(g in 1:generations) {

        #Check for extinction
        if(N < 2) {stop("Population went extinct! : (")}

        #Calculate allele frequencies
        #Removing epigenetic marks for allelic effect calculation
        results.mat.alleles[g,2:(nloc+1)] <- calc.al.freq(lapply(ind.mat, Re), length(ind.mat), nloc, allele.model, loc.attributes)

        #Argument type needs to have value of either "biallelic" or "infinite"
        #Calculate genotypic effects for intercept effects
        ind.Ga <- calc.genot.effects(ind.mat, ae, type = allele.model, qtl.ind, nloc) #Calculate genotype effects
        #Calculate genotypic effects for slope effects
        ind.Gb <- calc.genot.effects(ind.mat, ae, type = allele.model, qtl.slope.ind, nloc)

        ind.Gadj <- calc.genot.effects.lim(ind.mat, ae, type = allele.model, qtl.adj.ind, nloc, limits = c(0,1)) #Calculate genotype effects for plasticity adjustment
        ind.Ghed <- calc.genot.effects.lim(ind.mat, ae, type = allele.model, qtl.hed.ind, nloc, limits = c(0,3)) #Calculate genotype effects for canalization / bet-hedging adjustment
        ind.Gepi <- calc.genot.effects.lim(ind.mat, ae, type = allele.model, qtl.epi.ind, nloc, limits = c(0,1)) #Calculate genotype effects for epigenetic modification probability
        
        
        #Calculate environmental effects
        ind.E <- rnorm(n = N, mean = 0, sd = sigma.e + ind.Ghed) #Add sigma.e and bet-hedging effects

        #Calculate phenotypes for the juvenile stage, using epigenetically marked slope loci
        ind.Gepimod <- calc.genot.effects.epi(ind.mat, qtl.slope.ind, nloc)
        ind.P <- ind.Ga + ind.Gepimod + ind.E

        ##Environmental mismatches
        ind.M <- matrix(rep(0, N*L), ncol = L) #Initialize mismatch matrix
        ind.M[,1] <- abs(E[1 + L*(g-1)] - ind.P) #First env. mismatch

        ##Epigenetic effects are reset
        ind.mat <- reset.epigenetics(ind.mat)

        ##Costs of plasticity
        #ind.costs <- rep(0, N) #Initialize costs vector
        ind.costs <- ifelse(ind.Gepimod > 0, ke, 0)
        
        #Calculate phenotypes for first adult stage (developmental plasticity happens)
        ind.P <- ind.Ga + ind.Gb*cues[2 + L*(g-1)] + ind.E #First adult stage (development)

        ind.M[,2] <- abs(E[2 + L*(g-1)] - ind.P) #Second env. mismatch
        ind.costs <- ind.costs + ifelse(ind.Gb > 0, kd, 0)
        
        #Loop over the rest life stages (phenotypic adjustment can happen)
        for(j in 3:L) {
            ind.adjustment <- rbinom(n = N, size = 1, prob = ind.Gadj) #Did adjustment happen
            new.P <- ind.Ga + ind.Gb*cues[j + L*(g-1)] + ind.E #Calculate new phenotypes
            ind.P[ind.adjustment == 1] <- new.P[ind.adjustment == 1] #Adjust phenotypes
            ind.M[,j] <- abs(E[j + L*(g-1)] - ind.P) #Calculate second mismatch
            ind.costs <- ind.costs + ifelse(ind.Gb > 0, ka*ind.adjustment, 0) #Calculate costs
        }

        ### Modify loci epigenetically using the cue in the last life stage = L
        ind.mat <- epigenetic.modification(ind.mat, qtl.slope.ind, cues[5+L*(g-1)], ind.Gepi)

        #Store phenotypes
        results.mat.pheno[g,2] <- mean(ind.P) #Population mean of trait in final life stage
        results.mat.pheno[g,3] <- var(ind.Ga + ind.Gb) #Genetic variance
        results.mat.pheno[g,6] <- mean(ind.Ga) #Mean intercept genotypic value
        results.mat.pheno[g,7] <- mean(ind.Gb) #Mean slope genotypic value
        results.mat.pheno[g,8] <- mean(ind.Gadj) #Mean plasticity adjustment genotypic value
        results.mat.pheno[g,9] <- mean(ind.Ghed) #Mean bet-hedging genotypic value
        results.mat.pheno[g,10] <- mean(ind.Gepi) #Mean epigenetic modification probability (gen.val.)
        
        
        ##Calculate fitness over all life-stages, subtract costs of plasticity
        ind.W <- exp(-sel.intensity * apply(ind.M, MARGIN = 1, sum)) - ind.costs
        #Because costs can bring fitness below 0, all negative values are set to zero
        ind.W <- replace(ind.W, ind.W < 0, 0) #Set negative values to 0

        results.mat.pheno[g,4] <- mean(ind.W) #Population mean fitness
        results.mat.pheno[g,5] <- N #Population size

        ### Reproduction ###
        #Using logistic density regulation
        if(density.reg == "logistic") {
            #Calculate the number of offspring the population produces
            #Logistic population regulation
            No <- round(N*K*(1+r) / (N*(1+r) - N + K)) #Logistic population growth
            #Scale number of offspring by population mean fitness
            No <- round(mean(ind.W)*No); if(No < 1) {stop("Population went extinct! : (")}
            ind.mat <- produce.gametes.W(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
        }

        #Population regulation via Mikael's method
        if(density.reg == "sugar") {
            #Calculate the number of offspring that each individual produces
            No <- rpois(n = N, lambda = B*ind.W)
            #If number of offspring larger than K, remove some offspring randomly
            if(sum(No) > K) {
                while(sum(No) > K) { #Remove individuals until only K are left, seems slow...
                    no.ind <- (1:length(No))[No > 0] #Indexes of cases where No > 0
                    rem.ind <- sample(no.ind, 1)
                    No[rem.ind] <- No[rem.ind] - 1
                }
            }
            ind.mat <- produce.gametes.W.No(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
        }
        
        #Next generation
        N <- length(ind.mat) #Store new population size

         ### Mutation ###
        #Produce mutations
        mutations <- rep(foo, nloc)
        #Has to be 2*mu since mutations happen on individual basis and each ind has 2 copies
        mutations <- lapply(mutations, gen.mutations, N = N, mu = 2*mu) 
        for(l in 1:nloc) { #Loop over all loci
            if(any(mutations[[l]] == 1) == TRUE) { #Did any mutations happen?
                mut.ind <- ind.mat[mutations[[l]] == 1] #Select individuals to be mutated
                if(allele.model == "biallelic") { #Which alleles model is used?
                    mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) #Mutate alleles
                }
                if(allele.model == "infinite") {
                    if(loc.attributes[l] == "qtl" | loc.attributes[l] == "qtl.slope" | loc.attributes[l] == "qtl.adj" | loc.attributes[l] == "qtl.hed" | loc.attributes[l] == "qtl.epi") {
                        mut.ind <- lapply(mut.ind, mutate.inf.alleles, l = l) } else {
                            mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) }
                }
                ind.mat[mutations[[l]] == 1] <- mut.ind #Store results
            }
        }
    ########################    
        
        
    } #Done looping over generations

    #Modify results matrices if necessary

    #Return results
    return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes))
    
}
################################################################################################
###                      End of plasticity simulations                                       ###
################################################################################################


#############################################
###Following are some plotting functions ####
#############################################

#This function plots an overview of the simulation results
#' Plot simulation results
#'
#' This function plots population size, trait mean, and population mean fitness of the simulation
#'
#' @param data Phenotypic result matrix given by a another function (link to be completed)
#' @export
plot_sim.results <- function(data) {

    #First panel is population size
    A <- ggplot2::ggplot(data.frame(data), aes(x = generation, y = N)) +
           ggplot2::geom_line() 
       

    #Second panel is trait mean
    B <- ggplot2::ggplot(data.frame(data), aes(x = generation, y = pop.mean)) +
           ggplot2::geom_line() +
           ggplot2::ylab(label = "Trait mean")

    #Third panel is mean population fitness
    C <- ggplot2::ggplot(data.frame(data), aes(x = generation, y = W.mean)) +
        ggplot2::geom_line() +
        ggplot2::ylab(label = "W") +
        ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0))

    cowplot::plot_grid(A, B, C, align = "h", ncol = 3)

}

#This function plots an overview of allele frequency changes
#' Plot allele frequency results
#'
#' This function plots allele frequencies over time simulation. Note that this only works in the biallelic allele model
#'
#' @param data Allele frequency result matrix given by another function (link to be added)
#' @param nloc Number of loci in the simulation
#' @param loc.attr Types for the different loci, output by indsim.simulate (link to be added)
#' @export
plot_allele.results <- function(data, nloc, loc.attr) {

    #First need to transform data into long format for ggplot
    data <- data.frame(data) #Need to convert to data frame before melt
    data.long <- reshape2::melt(data, id.vars = c("generation"), measure.vars = c(paste(rep("locus", nloc), 1:nloc, sep = "")), variable.name = "locus", value.name = "freq")
    type.indices <- as.numeric(data.long$locus) #Since loci are in order this is OK
    data.long$loctype <- loc.attr[type.indices]

    ggplot2::ggplot(data.long, aes(x = generation, y = freq, group = locus, colour = loctype)) +
        ggplot2::geom_line() +
        ggplot2::xlab("Generation") +
        ggplot2::ylab("Allele frequency") +
        ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0))

}
            
#This function plots an overview of alleles present in the population
#Only works for the infinite alleles model
#' Plot allele effect distribution with the infinite alleles model
#'
#' This function plots the distribution of all allelic effects and their frequencies. Note that this only works in the infinite allele model
#'
#' @param data List of individual genotypes, produced by function (link to be addded)
#' @param nloc Number of loci used in the simulation
#' @param loc.attr Types for the different loci, output by indsim.simulate (link to be added)
#' @export
plot_inf.alleles <- function(data, nloc, loc.attr, plasticity = FALSE) {

    #First transforming the data
    #Select each locus and create a matrix of 2N alleles
    N <- length(data)
    allele.mat <- matrix(rep(0, 2*N*nloc), ncol = nloc)
    myselect <- function(data, col) { data[,col] }
    for(i in 1:nloc) {
        allele.mat[,i] <- as.vector(unlist(lapply(data, myselect, col = i)))
    }

    allele.mat <- data.frame(allele.mat) #Convert to dataframe and give colnames
    colnames(allele.mat) <- c(paste(rep("locus", nloc), 1:nloc, sep = ""))
    #Reshape with melt()
    allele.long <- reshape2::melt(allele.mat, measure.vars = c(paste(rep("locus", nloc), 1:nloc, sep = "")), variable.name = "locus", value.name = "effect")
    type.indices <- as.numeric(allele.long$locus) #Since loci are in order this is OK
    allele.long$loctype <- loc.attr[type.indices]
    if(plasticity == FALSE) {
    allele.long <- allele.long[allele.long$loctype == "qtl",] #Drop all loci that are not qtls
    
    #Plot with ggplot2
    ggplot2::ggplot(allele.long, aes(x = effect)) +
      ggplot2::geom_histogram(aes(y = ..count../(2*N)), colour = "black", fill = "white", binwidth = 0.1) +
      ggplot2::geom_vline(aes(xintercept =0), linetype = "dashed") +
      ggplot2::xlab("Allelic effect") +
      ggplot2::ylab("Allele frequency") +
      ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      ggplot2::facet_wrap(~ locus) }
    ##Plot when plasticity is used
    if(plasticity == TRUE) {
        allele.long <- allele.long[allele.long$loctype == "qtl" | allele.long$loctype == "qtl.slope" | allele.long$loctype == "qtl.adj" | allele.long$loctype == "qtl.hed" | allele.long$loctype == "qtl.epi",]

        #ggplot2
        ggplot2::ggplot(allele.long, aes(x = effect)) +
            ggplot2::geom_histogram(aes(y = ..count../(2*N)), colour = "black", fill = "white", binwidth = 0.1) +
            ggplot2::geom_vline(aes(xintercept =0), linetype = "dashed") +
            ggplot2::xlab("Allelic effect") +
            ggplot2::ylab("Allele frequency") +
            ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
            ggplot2::facet_wrap(loctype ~ locus) }    
    
}

##This function plots the mean reaction norm of the population
##Only works with the plasticity model
#'Plot mean reaction norm of the population in the final generation
#'
#' Line colour is determined by type of plasticity. Red is reversible plasticity, Blue is irreversable plasticity, green is bet-hedging
#' 
#' @param data genotypic values of intercepts and slopes etc.
#' @export
plot_reaction.norm <- function(data) {
    finalgen <- nrow(data)
    data <- as.data.frame(data)
    a <- data$G.int[finalgen]
    b <- data$G.slope[finalgen]
    cue <- seq(-1, 1, by = 0.1)
    phenotype <- a + b*(cue)
    ctype <- "black"
    if(b > 0.1) {
        if(data$G.adj[finalgen] > 0.1) { ctype <- "red" } else { ctype <- "blue" }
    }
    if(data$G.hed[finalgen] > 0.1) { cytpe <- "green" }
    
    plotdata <- data.frame(cue = cue, phenotype = phenotype)
    
    ggplot2::ggplot(plotdata, aes(x = cue, y = phenotype)) +
        ggplot2::geom_line(color = ctype) +
        ggplot2::coord_cartesian(ylim = c(-1,1)) +
        ggplot2::xlab("Environmental cue") +
        ggplot2::ylab("Phenotype")    
}

################################################################

###########################################################
# Generating different colours of environmental variation #
###########################################################

### Autoregressive (AR) process ###
#E_{t+1} = kappa*E_t + omega_t*sqrt(1-kappa^2)
#When kappa < 0 -> blue noise, kappa > 0 -> red noise, and kappa = 0 -> white noise
#' Generate coloured noise using autoregressive process
#'
#' This function generates coloured noise for generating environmental variation. Noise is generated using the equation E_{t+1} = kappa \times E_t + omega_t \times sqrt{1-kappa^2}
#'
#' @param kappa When kappa < 0, blue noise is generated, when kappa > 0 function generates red noise and when kappa = 0 white noise is generated
#' @param mean Mean of the environment
#' @param sd Standard deviation
#' @param generations How many generations (time steps) to generate noise
#' @export
coloured.noise <- function(kappa = 0, mean = 0, sd = 1, generations = 100) {
    E <- rep(0, generations) #Initialize the environmental variable
    E[1] <- rnorm(1, mean = mean, sd = sd) #Initialize first step
    for(i in 2:generations) {
        E[i] <- kappa*E[i-1] + rnorm(1, mean = mean, sd = sd)*sqrt(1-kappa^2) }
    return(E)
}

### 1/f noise ###
#Generating noise using sinusoidal processes, Following Cohen et al. 2009 
f.noise <- function(beta = 0, generations = 100, final.sd = 1) {
   #Need to generate generations/2 waves, where generations is even
   all.waves <- matrix(rep(0, generations*generations/2), ncol = generations)

   #Loop over all waves
   for(i in 1:(generations/2)) {

   f <- i #waves have different periods
   
   #Generate one wave
   #2*pi*(1:generations)/generations standardises 2pi over the entire lenght of the time series
   #f is the period of the time series, for the first wave it is 1
   #add a random phase for each wave + runif(1, min = 0. max = 2*pi)
   #The term 1/(f^(beta/2) is the spectral exponent
   #Note (Ruokolainen et al. 2009, has -beta/2 as exponent, but they probably actually used beta/2
   #without negative sign, beta = -1 gives blue noise, beta = 1 gives pink noise
   all.waves[i,] <- 1/(f^(beta/2)) * sin(2*pi*(1:generations)/generations * f + runif(1, min = 0, max = 2*pi))
   }

   #Sum all waves together
   summed.wave <- apply(all.waves, MARGIN = 2, sum)

   #Scale variance to wanted values
   #Variance scaling using Wichmann et al. 2005
   final.wave <- final.sd/sd(summed.wave)*(summed.wave - mean(summed.wave))

   return(final.wave)
}
#################################################################################################
