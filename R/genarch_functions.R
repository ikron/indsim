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
        epigenetic.modification <- function(ind.mat, qtl.ind, cue, epi.mod.ind) {
            modify <- function(x, loc, cue) { x[,loc] <- x[,loc] + cue*1i; return(x) }
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

###########################################################
#This function produces gametes on individual basis, 
#Production and random union of gametes, individuals contribute to the next generation weighted by their fitness
#This is a faster version of produce.gametes.W.No
    produce.gametes.faster <- function(ind.mat, loc.mat, No, ind.W, link.map = F, linkage, nloc.chr, n.chr) {
        
    Nind <- length(ind.mat) #Number of parents
    sum.No <- sum(No) #Total number of individuals to be produced
    ind.mat.new <- rep(list("genotype" = loc.mat), sum.No) #List of the next generation of individuals, initialize
    #4-strand index mat
    strand.mat <- matrix(c(1,3,1,4,2,3,2,4), ncol = 4)
  
    n <- 1 #Initialize counter

    No.ind <- No > 0 #Those individuals that have more than zero offspring
    Nind2 <- (1:Nind)[No.ind]

    #Pick the focal individual and another at random weighted by their fitness
    #Selfing not possible
    ind.index <- lapply(Nind2, my.sample.parents, Nind = Nind, No = No, ind.W = ind.W)  

    #Need to produce No individuals
    for(i in 1:length(ind.index)) { #Loop over all individuals (parents)
        
            
           # ind.index <- c(i, sample((1:Nind)[-i], size = No[i], prob = ind.W[-i])) #Selfing not possible
            for(k in 2:length(ind.index[[i]])) { #Loop over all offspring of focal individual
                 
                pind.mat <- ind.mat[c(ind.index[[i]][1], ind.index[[i]][k])]

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
                xo.positions.list <- lapply(n.crossover, my.xo.positions)
                xo.strands.list <- lapply(n.crossover, my.xo.strands)
                #Resolving crossovers
                linkage.gamete1 <- resolve.crossovers(meiocyte = linkage.gamete1, xo.positions.list, xo.strands.list, strand.mat, n.chr, n.crossover)
                #Select chromosomes for gametes (i.e. meiotic divisions)
                gamete1 <- make.gamete(linkage.gamete1)

                ### Gamete2 ###
                linkage.gamete2 <- convert.indmat.2.linkage(linkage, pind.mat[[2]], nloc.chr, n.chr, fourstrand = T)
                n.crossover <- gen.n.crossover(n.chr, lambda = 0.56)
                xo.positions.list <- rep(foo, n.chr)
                xo.strands.list <- rep(foo, n.chr)
                xo.positions.list <- lapply(n.crossover, my.xo.positions)
                xo.strands.list <- lapply(n.crossover, my.xo.strands)
                #Resolving crossovers and perform meiotic divisions
                linkage.gamete2 <- resolve.crossovers(meiocyte = linkage.gamete2, xo.positions.list, xo.strands.list, strand.mat, n.chr, n.crossover)
                gamete2 <- make.gamete(linkage.gamete2)
                } #Done making gametes

                #Union of gametes, aka fertilization <3
                ind.mat.new[[n]][1,] <- gamete1
                ind.mat.new[[n]][2,] <- gamete2
                n <- n + 1 #update counter        
            } #Close looping over all offspring of one parent
    } #Close looping over all parents, #No individuals produced
        
    return(ind.mat.new)
}
####################################################################


##Note the function 'sample' has some undesired behaviour, when x is numeric and length(x) = 1, then function actually samples from 1:x .... WHY R WHY!!!!?
##Use a safer version of the function to avoid this problem
resample <- function(x, ...) x[sample.int(length(x), ...)] ##Safer version of sample function
#Helper function, pick parents
my.sample.parents <- function(x, Nind, No, ind.W) { c(x, resample((1:Nind)[-x], size = No[x], replace = T, prob = ind.W[-x])) }



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

#Functions for picking crossover positions and strands
my.xo.positions <- function(n.crossover) { sort(runif(n.crossover, min = 0, max = 1)) }
my.xo.strands <- function(n.crossover) { sample(c(1,2,3,4), size = n.crossover, replace = T) }


##################################################################################################


### Functions for generating mutations ###########################################################
#Do mutations occur or not?
gen.mutations <- function(x, N, mu) {
    x <- rbinom(N, size = 1, prob = mu) }

#Do mutations occur or not. Version for unique mutation rates per locus
gen.mutations.loc <- function(x, N) {
    rbinom(N, size = 1, prob = x) }

#This function makes mutations with the K-alleles model
        mutate.K.alleles <- function(x, l, alleles) {
            mut.al.ind <- sample(c(1,2),1)
            cur.allele <- x[mut.al.ind,l]
            new.allele <- sample(alleles[[l]][alleles[[l]] !=cur.allele],1) #Select one allele
            x[mut.al.ind,l] <- new.allele
            return(x)
        }

#This function makes mutations with the infinite alleles model
mutate.inf.alleles <- function(x, l, sigma.a) {
    mut.al.ind <- sample(c(1,2),1)
    x[mut.al.ind,l] <- rnorm(1, mean = 0, sd = sigma.a) #This may need to changed accordingly
    return(x)
}

#This function makes deleterious mutations according to the empirical distribution
mutate.delres.alleles <- function(x, l) {
    mut.al.ind <- sample(c(1,2),1)
    x[mut.al.ind,l] <- ifelse(rbinom(1,1,0.3), 1, rbeta(1, 1, 6)) #Empirical distribution
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
delres.fitness <- function(ind.mat, delres.ind, nloc, delres.h) {
    mydselect <- function(x, delres.ind, nloc) {x[1:2, (1:nloc %in% delres.ind)]}
    ind.mat <- lapply(ind.mat, mydselect, delres.ind, nloc) #Select only deleterious recessives
    ##Check that both alleles are different from zero, return the smallest value if true, 0 if false
    delres.check <- function(ind.mat) { ifelse(all(ind.mat == 0), 0, ifelse(all(ind.mat != 0), max(ind.mat), delres.h*max(ind.mat))) }
    delres.check2 <- function(ind.mat) { apply(ind.mat, MARGIN = 2, delres.check) }
    delr.mat <- lapply(ind.mat, delres.check2) #Check all genotypes
    ##Calculate effect on fitness for each ind, effects are multiplicative
    delres.ef <- function(delr.mat) { prod(1-delr.mat) }
    delr.vec <- as.vector(unlist(lapply(delr.mat, delres.ef))) #Calculate effects on fitness
    return(delr.vec) 
}



#Old delres.check function, did not take possible partial recessivity into account
#New function does this with the delres.h parameter
#delres.check <- function(ind.mat) { ifelse(all(ind.mat != 0), min(ind.mat), 0) }

#################################################################


##Function to save and load genotypes to a text file
save.genotypes <- function(ind.mat, nloci, parname, generation, ID ) {

    dirname <- paste("par_", parname, "/ID_", ID, sep = "") #Name of directory
    #Then create a directory if it does not exist
    if(!dir.exists(dirname)) { dir.create(dirname, recursive = TRUE) }
    filename <- paste("gen_", generation, sep = "")
    finalname <- paste(dirname, "/", filename, ".csv", sep = "")
    lapply(ind.mat, write, finalname, append = TRUE, ncolumns = 2*nloci)

}

##Function to load genotypes from text file, that was generated using "save.genotypes"
load.genotypes <- function(filename, nloci) {
    conn <- file(filename, open = "r")
    linn <- readLines(conn)
    #Need to close connection
    close(conn)
    
    #Converting from character format to numeric
    tokens <- strsplit(linn, split = " ")
    myconvert <- function(x, nloci) { matrix(as.numeric(x), ncol = nloci) }
    genotypes <- lapply(tokens, myconvert, nloci = nloci)

    return(genotypes)
}

##Function to load loci attributes
load.loc.attributes <- function(file) {

    loc.attributes <- scan(file, what = "char", nlines = 1, quiet = T) #Read loci types
    n.chr <- scan(file, nlines = 1, skip = 1, quiet = T) #Read number of chromosomes
    nloc.chr <- scan(file, nlines = n.chr, skip = 2, quiet = T) #Read n.loci per chromosome
    linkage <- rep(list(0), n.chr) #initialize linkage map
    #Then read map positions chr by chr
    for(i in 1:n.chr) { linkage[[i]] <- scan(file, nlines = 1, skip = (2 + n.chr + i -1), quiet = T) }

    return(list("loci" = loc.attributes, "n.chr" = n.chr, "nloc.chr" = nloc.chr, "linkage" = linkage))
}

###Functions to calculate epistatic effects
check.interaction <- function(x, genot) {
    check <- apply(genot[,x], 2, function(x) {x != 0})
    dose <- apply(check, 2, sum)
    if(any(dose == 0)) {epistat.g <- 0} else {
        if(sum(dose) == 4) {epistat.g <- 2 }
        if(sum(dose) < 4) {epistat.g <- 1 }
    }
    return(epistat.g)
}

calc.epistat.effects <- function(genot, interactions) {
#Check all interactions for a single genotype
epistat.g <- apply(interactions, 1, check.interaction, genot = genot)
return(sum(epistat.g*epi.effects))
}
##################################################
    

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
#' @param mu.qtl Mutation rate for qtl loci
#' @param mu.delres Mutation rate for deleterious recessives
#' @param mu.neutral Mutation rate for neutral loci
#' @param sigma.a Standard deviation of mutational effects
#' @param n.chr Number of chromosomes (linkage groups)
#' @param nloc.chr Vector of the length of n.chr giving number of loci for each chromosome
#' @param n.delres Number of deleterious recessives, defaults to 0
#' @param n.neutral Number of neutral loci, defaults to 0
#' @param delres.h Dominance coefficient of deleterious mutations
#' @param linkage.map If map distances are determined randomly, set this to "random", if fixed linkage map is provided can set this to "user"
#' @param linkage Matrix of 3 rows and nloc columns, first two rows can be zeros and third row determines map distances relative to chromosome coordinates from 0 to 1
#' @param epistasis Whether epistasis is present in the simulation
#' @param Eprob Probability that two loci exhibit epistasis
#' @param Emean Mean of epistasis coefficients
#' @param sigma.epi Standard deviation of epistasis coefficients
#' @param save.all whether genotypes should be saved in each generation?
#' @param parname name of the parameter set used in saving the genotype files
#' @export
indsim.simulate <- function(N, generations, sel.intensity, init.f, init.n, a, sigma.e, K, r, B, opt.pheno, density.reg, allele.model, mu.qtl, mu.delres = 0, mu.neutral = 0, sigma.a = 1, n.chr, nloc.chr, nqtl, n.delres = 0, n.neutral = 0, delres.h = 0.25, linkage.map = "random", epistasis = FALSE, Eprob = 0, Emean = 0, sigma.epi = 1, linkage = NULL, save.all = FALSE, parname = NULL) {

    #Prepare linkage groups
    foo <- list(0)
    nloc <- nqtl + n.delres + n.neutral

    #Perform some checks that initial parameters are sensible
    if(nloc != sum(nloc.chr)) { stop("Number of loci and loci per chromosome don't match!") }
    if(n.delres > 0 & mu.delres == 0) {stop("Mutation rate of deleterious recessives is zero!") }

    loc.attributes <- sample(c(rep("qtl", nqtl), rep("neutral", n.neutral), rep("delres", n.delres)), size = nloc, replace = FALSE) #Types for all loci
    delres.ind <- (1:nloc)[loc.attributes == "delres"] #Indices of deleterious recessives
    qtl.ind <- (1:nloc)[loc.attributes == "qtl"] #Indices of loci that affect the phenotype
    #Locus specific mutation rates
    loc.mu <- rep(0, nloc)
    loc.mu[delres.ind] <- mu.delres
    loc.mu[qtl.ind] <- mu.qtl
    loc.mu[(1:nloc)[loc.attributes == "neutral"]] <- mu.neutral

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

    ##Setup epistatic interactions if they exist in the simulation
    if(epistasis == TRUE) {
        #probability of epistasis = Eprob
        #choose among all of the possible combinations those mutations that are going to have an interaction
        n.comb <- choose(nqtl, 2) #Number of combinations
        loc.comb <- t(combn(nloc, 2)) #Generate all combinations and transpose
        n.epistat <- rbinom(1, n.comb, Eprob) #Number of interactions
        interactions <- loc.comb[sample(1:n.comb, size = n.epistat, replace = FALSE),]
        #interactions <- matrix(sample(1:nloc, size = 2*n.epistat), ncol = 2) #matrix of interacting loci
        epi.effects <- rnorm(n.epistat, mean = Emean, sd = sigma.epi)
    } #Done generating epistatic interactions

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
    ##Note that + 1 in neutral index is because some following calculations use allele frequencies and column numbers need to be adjusted by one because first column is generation number

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

    ### If all genotypes need to be saved, initialize ID and save loci attributes and linkage map
    #When running multiple parallel simulations this is needed so different runs are saved separately
    #
    if(save.all == TRUE) {
        ##Generate a random run ID
        ID <- substr(as.character(runif(1,1,9)), 3, 10)
        #Saving the loc.attributes information
        dirname <- paste("par_", parname, "/ID_", ID, sep = "") #Directory
        if(!dir.exists(dirname)) { dir.create(dirname, recursive = TRUE) }
        filename <- "locattr"
        finalname <- paste(dirname, "/", filename, ".txt", sep = "")
        write(loc.attributes, file = finalname, append = TRUE, ncolumns = nloc) #Save loci attributes
        write(c(n.chr, nloc.chr), file = finalname, append = TRUE, ncolumns = 1) #Save n chr
        for(i in 1:n.chr) { write(linkage[[i]][3,], file = finalname, append = TRUE, ncolumns = nloc.chr[i]) } #save map positions for each locus in each chromosome
    }
    
    #Loop over generations
    for(g in 1:generations) {

        #Check for extinction
        if(N < 2) {stop("Population went extinct! : (")}

        #Calculate allele frequencies
        results.mat.alleles[g,2:(nloc+1)] <- calc.al.freq(ind.mat, length(ind.mat), nloc, allele.model, loc.attributes)

        ### If need to save genotypes do it here ###
        if(save.all == TRUE) {
            save.genotypes(ind.mat, nloci = nloc, parname = parname, generation = g, ID = ID)
        }


        #Argument type needs to have value of either "biallelic" or "infinite"
        ind.G <- calc.genot.effects(ind.mat, ae, type = allele.model, qtl.ind, nloc) #Calculate genotype effects
        #sigma.e <- (abs(ind.G) + 1)*coef.e #Environmental variation scales with genotypic effects
        if(epistasis == TRUE) {
            EPI <- as.numeric(unlist(lapply(ind.mat, calc.epistat.effects, interactions = interactions)))
            ind.G <- ind.G + EPI #Sum genotypic and epistatic effects, EPI part of genotype
        }

        
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
            ind.W <- ind.W*delres.fitness(ind.mat, delres.ind, nloc, delres.h)
        }
        

        results.mat.pheno[g,6] <- mean(ind.W) #Population mean fitness
        results.mat.pheno[g,7] <- N #Population size

        #Calculate effective population size (if g > 10)
        if(g > 10 & n.neutral > 0) { results.mat.alleles[g,(nloc+2)] <- calc.Ne(af1 = results.mat.alleles[(g-10), neutral.index], af2 = results.mat.alleles[g, neutral.index], N1 = results.mat.pheno[g-10,7], N2 = results.mat.pheno[g,7], nloc = n.neutral, t = 10) }
        
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
                while(sum(No) > K) { #Remove individuals until only K are left
                    No.rem <- sum(No) - K
                    no.ind <- (1:length(No))[No > 0] #Indexes of cases where No > 0
                    if(No.rem < length(no.ind)) {   
                        rem.ind <- sample(no.ind, No.rem, replace = FALSE)
                        No[rem.ind] <- No[rem.ind] - 1
                    }
                    if(No.rem >= length(no.ind)) {
                        No[no.ind] <- No[no.ind] - 1
                    }
                }
            }
            ind.mat <- produce.gametes.faster(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
        }
        
        #Next generation
        N <- length(ind.mat) #Store new population size

        ### Mutation ###
        #Produce mutations
        mutations <- rep(foo, nloc)
        #Has to be 2*mu since mutations happen on individual basis and each ind has 2 copies
        mutations <- lapply(2*loc.mu, gen.mutations.loc, N = N) 
        for(l in 1:nloc) { #Loop over all loci
            if(any(mutations[[l]] == 1) == TRUE) { #Did any mutations happen?
                mut.ind <- ind.mat[mutations[[l]] == 1] #Select individuals to be mutated
                if(allele.model == "biallelic") { #Which alleles model is used?
                    mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) #Mutate alleles
                }
                if(allele.model == "infinite") {
                   if(loc.attributes[l] == "qtl") {
                       mut.ind <- lapply(mut.ind, mutate.inf.alleles, l = l, sigma.a = sigma.a) }
                   if(loc.attributes[l] == "delres") {
                       mut.ind <- lapply(mut.ind, mutate.delres.alleles, l = l) }
                   if(loc.attributes[l] == "neutral") {   
                       mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) }
                }
                ind.mat[mutations[[l]] == 1] <- mut.ind #Store results
            }
        }
    ########################    
        
        
    } #Done looping over generations


    #Modify results matrices if necessary

#    if(save.all == TRUE) {
#        return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = #all.genotypes, "loci" = loc.attributes))
#    }

    #Return results
#    if(save.all == FALSE) {
        return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes))
#    }
}

###########################################
### End of individual based simulations ###
####################################################################################################

#####################################################################################################
##################################################################
##Plasticity simulations using a model similar to Botero et al.###
##################################################################

#Need to incorporate between generation effects, mediated by epigenetics...
#roxygen stuff
#This function performs individuals based simulations with model similar to Botero, evolution of plasticity and between generations effects are possible
#' Simulations with plasticity
#'
#' This function performs simulations with plasticity and epigenetic effects in fluctuating environments.
#'
#'
#' @param N Starting population size
#' @param generations Number of generations to run the simulations
#' @param L Number of timesteps in each generation
#' @param sel.intensity Selection intensity for stabilising selection
#' @param init.f Initial allele frequencies for QTL
#' @param init.n Initial allele frequencies for neutral loci
#' @param sigma.e Environmental standard deviations of the phenotypic trait
#' @param K Carrying capacity of the population
#' @param density.reg How population density is regulated. Currently there are two possible values:  "logistic" which implements logistic population density regulation (see manual) and "sugar" which constrains population to carrying capacity by randomly removing individuals when population size is higher than K
#' @param r If density regulation is "logistic", this is the population growth rate.
#' @param B If density regulation is "sugar", this it the mean number of offspring each perfectly adapted individual has
#' @param R Rate of environmental change relative to generation time
#' @param P Predictability of the environment
#' @param allele.model Which allele model to use for qtl. Currently the options are: "infinite" which uses an infinite alleles model for qtl, and "biallelic" which uses a biallelic model for qtl
#' @param nqtl Number of loci affecting the phenotype (intercept)
#' @param nqtl.slope Number of loci affecting phenotype (slope)
#' @param ntql.adj Number of loci affecting phenotype (probability of phenotypic adjustment)
#' @param ntql.hed Number of loci affecting phenotype (bethedging)
#' @param nqtl.epi Number of loci affecting phenotype (probability of epigenetic modification)
#' @param mu Mutation rate, which currently is the same for all loci
#' @param sigma.a Standard deviation of mutational effects
#' @param n.chr Number of chromosomes (linkage groups)
#' @param nloc.chr Vector of the length of n.chr giving number of loci for each chromosome
#' @param n.neutral Number of neutral loci, defaults to 0
#' @param linkage.map If map distances are determined randomly, set this to "random", if fixed linkage map is provided can set this to "user"
#' @param linkage Matrix of 3 rows and nloc columns, first two rows can be zeros and third row determines map distances relative to chromosome coordinates from 0 to 1
#' @param kd Fitness cost of developmental plasticity
#' @param ka Fitness cost of reversible plasticity
#' @param ke Fitness cost of epigenetic adjustment
#' @param sensitivity Whether to draw random parameter values for sensitivity analysis
#' @export
indsim.plasticity2.simulate <- function(N, generations, L, sel.intensity, init.f, init.n, a, sigma.e, K, r, B, R, P, density.reg, allele.model, mu, sigma.a = 1, n.chr, nloc.chr, nqtl, nqtl.slope, nqtl.adj, nqtl.hed, nqtl.epi, n.neutral = 0, linkage.map = "random", linkage = NULL, kd, ka, ke, sensitivity = FALSE) {

    #If using sensitivity analysis draw new random values for some parameters
    if(sensitivity == TRUE) {
        ################################################################
        ### Getting random parameter values for sensitivity analysis ###
        ################################################################
        #sel.intensity <- runif(1, min = 0.01, max = 1)
        sigma.e <- runif(1, min = 0.01, max = 0.5)
        mu <- 10^runif(1, min = -6, max = -3)
        n.chr <- sample(c(1:10), 1)
        nqtl <- sample(c(4:25), 1) #Sampling number of loci from 4 to 25
        nqtl.slope <- sample(c(4:25), 1)
        nqtl.adj <- sample(c(4:25), 1)
        nqtl.hed <- sample(c(4:25), 1)
        nqtl.epi <- sample(c(4:25), 1)
        #Calculating number of loci per chromosome
        perlocus <- (nqtl+nqtl.slope+nqtl.adj+nqtl.hed+nqtl.epi)%/%n.chr #Integer division
        nloc.chr <- rep(perlocus, n.chr)
        nloc.chr[1] <- nloc.chr[1] + (nqtl+nqtl.slope+nqtl.adj+nqtl.hed+nqtl.epi)%%n.chr #Add remainder to chr 1
        sigma.a <- runif(1, min = 0.01, max = 2)
        
        ################################################################

        ##### Inserting some debug stuff ##################################################
        #debugID <- sample(1:30000, 1)
        #debugID <- paste("sensdebug", debugID, sep = "")

        #system(paste("echo 'Parameters in this simulation were:' >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'sigma.e\t", sigma.e, "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'mu\t", mu, "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'n.chr\t", n.chr, "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'nqtl\t", nqtl, "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'nqtl.slope\t", nqtl.slope, "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'nqtl.adj\t", nqtl.adj, "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'nqtl.hed\t", nqtl.hed, "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'nqtl.epi\t", nqtl.epi, "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'nloc.chr\t", paste(nloc.chr, collapse = " "), "'", " >> ", debugID, ".txt", sep = ""))
        #system(paste("echo 'sigma.a\t", sigma.a, "'", " >> ", debugID, ".txt", sep = ""))
        #####################################################################################
        
    }
    
    
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
        if(N < 2) {return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes, "linkagemap" = linkage))}

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
        #Note! Using some thresholds for plasticity costs that is 0.05 below that there are no costs
        ind.costs <- rep(0, N) #Initialize costs vector
        
        
        #Calculate phenotypes for first adult stage (developmental plasticity happens)
        ind.P <- ind.Ga + ind.Gb*cues[2 + L*(g-1)] + ind.E #First adult stage (development)

        ind.M[,2] <- abs(E[2 + L*(g-1)] - ind.P) #Second env. mismatch
        ind.costs <- ind.costs + ifelse(ind.Gb > 0.05 | ind.Gb < -0.05, kd, 0)
        
        #Loop over the rest life stages (phenotypic adjustment can happen)
        for(j in 3:L) {
            ind.adjustment <- rbinom(n = N, size = 1, prob = ind.Gadj) #Did adjustment happen
            new.P <- ind.Ga + ind.Gb*cues[j + L*(g-1)] + ind.E #Calculate new phenotypes
            ind.P[ind.adjustment == 1] <- new.P[ind.adjustment == 1] #Adjust phenotypes
            ind.M[,j] <- abs(E[j + L*(g-1)] - ind.P) #Calculate next mismatch
            ind.costs <- ind.costs + ifelse(ind.Gb > 0.05 | ind.Gb < -0.05, ka*ind.adjustment, 0) #Calculate costs
        }

        ### Modify loci epigenetically using the cue in the last life stage = L
        epi.mod.ind <- rbinom(n = N, size = 1, prob = ind.Gepi) #Does epigenetic modification happen?
        ind.mat <- epigenetic.modification(ind.mat, qtl.slope.ind, cues[5+L*(g-1)], epi.mod.ind)
        #Costs of epigenetic modification
        ind.costs <- ind.costs + ifelse(ind.Gepi > 0.05, ke*epi.mod.ind, 0)

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
            #Check if there are more than 1 offspring, otherwise pop. goes extinct
            if(sum(No) < 2) { return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes, "linkagemap" = linkage)) }
            #If number of offspring larger than K, remove some offspring randomly
            if(sum(No) > K) {
                while(sum(No) > K) { #Remove individuals until only K are left
                    No.rem <- sum(No) - K
                    no.ind <- (1:length(No))[No > 0] #Indexes of cases where No > 0
                    if(No.rem < length(no.ind)) {   
                        rem.ind <- sample(no.ind, No.rem, replace = FALSE)
                        No[rem.ind] <- No[rem.ind] - 1
                    }
                    if(No.rem >= length(no.ind)) {
                        No[no.ind] <- No[no.ind] - 1
                    }
                }
            }
            ind.mat <- produce.gametes.faster(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
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
                        mut.ind <- lapply(mut.ind, mutate.inf.alleles, l = l, sigma.a = sigma.a) } else {
                            mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) }
                }
                ind.mat[mutations[[l]] == 1] <- mut.ind #Store results
            }
        }
    ########################    
        
        
    } #Done looping over generations

    #Modify results matrices if necessary

    ### Sensitivity debug ###
    #if(sensitivity == TRUE) {
    #    system(paste("echo 'Simulation OK'", " >> ", debugID, ".txt", sep = ""))
    #}
    ########################

    #Return results
    return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes, "linkagemap" = linkage))
    
}
################################################################################################
###                      End of plasticity simulations                                       ###
################################################################################################

###Simulations where simulation is started with existing genotypes
#' Simulations with plasticity using starting genotypes
#'
#' This function performs simulations with plasticity and epigenetic effects in fluctuating environments. Using starting genotypes
#'
#'
#' @param N Starting population size
#' @param generations Number of generations to run the simulations
#' @param L Number of timesteps in each generation
#' @param sel.intensity Selection intensity for stabilising selection
#' @param init.f Initial allele frequencies for QTL
#' @param init.n Initial allele frequencies for neutral loci
#' @param sigma.e Environmental standard deviations of the phenotypic trait
#' @param K Carrying capacity of the population
#' @param density.reg How population density is regulated. Currently there are two possible values:  "logistic" which implements logistic population density regulation (see manual) and "sugar" which constrains population to carrying capacity by randomly removing individuals when population size is higher than K
#' @param r If density regulation is "logistic", this is the population growth rate.
#' @param B If density regulation is "sugar", this it the mean number of offspring each perfectly adapted individual has
#' @param R Old rate of environmental change relative to generation time
#' @param P Old predictability of the environment
#' @param allele.model Which allele model to use for qtl. Currently the options are: "infinite" which uses an infinite alleles model for qtl, and "biallelic" which uses a biallelic model for qtl
#' @param mu Mutation rate, which currently is the same for all loci
#' @param sigma.a Standard deviation of mutational effects
#' @param n.neutral Number of neutral loci, defaults to 0
#' @param linkage.map If map distances are determined randomly, set this to "random", if fixed linkage map is provided can set this to "user"
#' @param linkage Matrix of 3 rows and nloc columns, first two rows can be zeros and third row determines map distances relative to chromosome coordinates from 0 to 1
#' @param kd Fitness cost of developmental plasticity
#' @param ka Fitness cost of reversible plasticity
#' @param ke Fitness cost of epigenetic adjustment
#' @param sensitivity Whether to draw random parameter values for sensitivity analysis
#' @param ind.mat List of the genotypes from a previous simulation
#' @param newR New rate of environmental change
#' @param newP New predictability of the environment
#' @param changepoint Time point at which parameters change
#' @param loci Loci attributes of the previous simulation
#' @export
indsim.plasticity.startgenot.simulate <- function(N, generations, L, sel.intensity, init.f, init.n, a, sigma.e, K, r, B, R, P, density.reg, allele.model, mu, sigma.a = 1, n.chr = 0, nloc.chr = 0, nqtl = 0, nqtl.slope = 0, nqtl.adj = 0, nqtl.hed = 0, nqtl.epi = 0, n.neutral = 0, linkage.map = "user", linkage = NULL, kd, ka, ke, sensitivity = FALSE, ind.mat, newR, newP, changepoint, loci) {

    #If using sensitivity analysis draw new random values for some parameters
    if(sensitivity == TRUE) {
        ################################################################
        ### Getting random parameter values for sensitivity analysis ###
        ################################################################
        #sel.intensity <- runif(1, min = 0.01, max = 1)
        sigma.e <- runif(1, min = 0.01, max = 0.5)
        mu <- 10^runif(1, min = -6, max = -3)
        n.chr <- sample(c(1:10), 1)
        nqtl <- sample(c(4:25), 1) #Sampling number of loci from 4 to 25
        nqtl.slope <- sample(c(4:25), 1)
        nqtl.adj <- sample(c(4:25), 1)
        nqtl.hed <- sample(c(4:25), 1)
        nqtl.epi <- sample(c(4:25), 1)
        #Calculating number of loci per chromosome
        perlocus <- (nqtl+nqtl.slope+nqtl.adj+nqtl.hed+nqtl.epi)%/%n.chr #Integer division
        nloc.chr <- rep(perlocus, n.chr)
        nloc.chr[1] <- nloc.chr[1] + (nqtl+nqtl.slope+nqtl.adj+nqtl.hed+nqtl.epi)%%n.chr #Add remainder to chr 1
        sigma.a <- runif(1, min = 0.01, max = 2)
        ################################################################
    }
    
    
    #Initialize time
    timesteps <- generations*L #There are ngenerations*L time steps in the model
    t <- 1:timesteps
    E <- rep(0, timesteps) #Initialize the environment
    cues <- rep(0, timesteps) #Initialize the environmental cues

    E[1:changepoint] <- sin((2*pi*t[1:changepoint])/(L*R)) #Initialize the environment
    E[(changepoint+1):timesteps] <- sin((2*pi*t[(changepoint+1):timesteps])/(L*newR))
    cues[1:changepoint] <- rnorm(n = changepoint, mean = E[1:changepoint]*P, sd = (1-P)/3) #Environmental cues
    cues[(changepoint+1):timesteps] <- rnorm(n = (timesteps - changepoint), mean = E[(changepoint+1):timesteps]*newP, sd = (1-newP)/3)
    
    #Prepare linkage groups and locus information
    loc.attributes <- loci #Types for all loci and their order
    foo <- list(0)
    nloc <- length(loc.attributes)
    nqtl <- sum(loc.attributes == "qtl")
    nqtl.slope <- sum(loc.attributes == "qtl.slope")
    nqtl.adj <- sum(loc.attributes == "qtl.adj")
    nqtl.hed <- sum(loc.attributes == "qtl.hed")
    nqtl.epi <- sum(loc.attributes == "qtl.epi")
    
    #Setup linkage map
    if(linkage.map == "random") {
        linkage <- rep(foo, n.chr)
        for(i in 1:n.chr) {
            linkage[[i]] <- matrix(rep(0,3*nloc.chr[i]), ncol = nloc.chr[i])
            linkage[[i]][3,] <- sort(runif(nloc.chr[i])) #Map-positions for each locus
        }
    }

    n.chr <- length(linkage)
    nloc.chr <- unlist(lapply(linkage, ncol)) #Number of loci per chromosome

    #Perform some checks that initial parameters are sensible
    if(nloc != sum(nloc.chr)) { stop("Number of loci and loci per chromosome don't match!") }

    
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

    #Get the current N
    N <- length(ind.mat)

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
    #ind.mat <- rep(list("genotype" = loc.mat), N) #List for individuals

    #Allele frequencies
    #init.freq <- rep(init.f,nloc) #for p
    #Here can also make explicit initial frequencies for qtls, delres and neutral markers
    #if(n.neutral > 0) { init.freq[loc.attributes == "neutral"] <- init.n } #Neutral markers start at these frequencies

    #Genotype frequencies
    #genotype.mat <- matrix(rep(0,3*nloc), ncol = nloc)
    #rownames(genotype.mat) <- c("AA", "Aa", "aa")

    #Initialize genotype frequencies
    #for(i in 1:nloc) {
    #    genotype.mat[,i] <- HW.genot.freq(init.freq[i])
    #}
    #################################

    #ind.mat <- initialise.ind.mat(ind.mat, genotype.mat, N, nloc) #Sample starting genotypes


    
    #Loop over generations
    for(g in 1:generations) {

        #Check for extinction
        if(N < 2) {return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes, "linkagemap" = linkage))}

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
        #Note! Using some thresholds for plasticity costs that is 0.05 below that there are no costs
        ind.costs <- rep(0, N) #Initialize costs vector
        
        
        #Calculate phenotypes for first adult stage (developmental plasticity happens)
        ind.P <- ind.Ga + ind.Gb*cues[2 + L*(g-1)] + ind.E #First adult stage (development)

        ind.M[,2] <- abs(E[2 + L*(g-1)] - ind.P) #Second env. mismatch
        ind.costs <- ind.costs + ifelse(ind.Gb > 0.05 | ind.Gb < -0.05, kd, 0)
        
        #Loop over the rest life stages (phenotypic adjustment can happen)
        for(j in 3:L) {
            ind.adjustment <- rbinom(n = N, size = 1, prob = ind.Gadj) #Did adjustment happen
            new.P <- ind.Ga + ind.Gb*cues[j + L*(g-1)] + ind.E #Calculate new phenotypes
            ind.P[ind.adjustment == 1] <- new.P[ind.adjustment == 1] #Adjust phenotypes
            ind.M[,j] <- abs(E[j + L*(g-1)] - ind.P) #Calculate next mismatch
            ind.costs <- ind.costs + ifelse(ind.Gb > 0.05 | ind.Gb < -0.05, ka*ind.adjustment, 0) #Calculate costs
        }

        ### Modify loci epigenetically using the cue in the last life stage = L
        epi.mod.ind <- rbinom(n = N, size = 1, prob = ind.Gepi) #Does epigenetic modification happen?
        ind.mat <- epigenetic.modification(ind.mat, qtl.slope.ind, cues[5+L*(g-1)], epi.mod.ind)
        #Costs of epigenetic modification
        ind.costs <- ind.costs + ifelse(ind.Gepi > 0.05, ke*epi.mod.ind, 0)

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
                while(sum(No) > K) { #Remove individuals until only K are left
                    No.rem <- sum(No) - K
                    no.ind <- (1:length(No))[No > 0] #Indexes of cases where No > 0
                    if(No.rem < length(no.ind)) {   
                        rem.ind <- sample(no.ind, No.rem, replace = FALSE)
                        No[rem.ind] <- No[rem.ind] - 1
                    }
                    if(No.rem >= length(no.ind)) {
                        No[no.ind] <- No[no.ind] - 1
                    }
                }
            }
            ind.mat <- produce.gametes.faster(ind.mat, loc.mat, No, ind.W, link.map = T, linkage, nloc.chr, n.chr)
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
                        mut.ind <- lapply(mut.ind, mutate.inf.alleles, l = l, sigma.a = sigma.a) } else {
                            mut.ind <- lapply(mut.ind, mutate.K.alleles, l = l, alleles = alleles) }
                }
                ind.mat[mutations[[l]] == 1] <- mut.ind #Store results
            }
        }
    ########################    
        
        
    } #Done looping over generations

    #Modify results matrices if necessary

    #Return results
    return(list("phenotype" = results.mat.pheno, "alleles" = results.mat.alleles, "genotypes" = ind.mat, "loci" = loc.attributes, "linkagemap" = linkage))
    
}
################################################################################################
###                      End of plasticity simulations with starting genotypes               ###
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


#' Plot phenotype mean
#'
#' This function plots population mean from the simulation
#' @param data Phenotypic result matrix from simulation
#' @export
plot_sim.trait <- function(data) {
    ggplot2::ggplot(data.frame(data), aes(x = generation, y = pop.mean)) +
        ggplot2::geom_line() +
        ggplot2::ylab(label = "Mean phenotype") +
        ggplot2::xlab("Generation")
}

#' Plot mean fitness
#'
#' This function plots population mean fitness from the simulation
#' @param data Phenotypic results matrix from simulation
#' @export
plot_sim.fitness <- function(data) {
    ggplot2::ggplot(data.frame(data), aes(x = generation, y = W.mean)) +
        ggplot2::geom_line() +
        ggplot2::ylab(label = "W") +
        ggplot2::xlab("Generation") +    
        ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0))
}

#' Plot genotypic values
#'
#' This function plots genotypic values of intercept, slope, adjustment, epigenetic mod and bethedging probabilities
#' @param data Phenotypic results matrix from simulation
#' @export
plot_sim.genotypicval <- function(data) {
    data <- data.frame(data[,-c(2:5)]) #Drop columns that are not needed
    data <- tidyr::gather(data, type, genoval, G.int:G.epi, factor_key = TRUE) #Change to long format
    mycolors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")
    
    ggplot2::ggplot(data, aes(x = generation, y = genoval, color = type)) +
        ggplot2::geom_line() +
        ggplot2::ylab("Genotypic value") +
        ggplot2::xlab("Generation") +
        ggplot2::scale_color_manual(values = mycolors, labels = c("Intercept", "Slope", "Plastic\nadjustment", "Developmental\nvariation", "Epigenetic\nmodification")) +
        ggplot2::theme(legend.key.size = unit(2.5, 'lines'))
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
    mycolors <- c("#e41a1c",  "#4daf4a", "#ff7f00", "#984ea3",  "#377eb8")

    ggplot2::ggplot(data.long, aes(x = generation, y = freq, group = locus, colour = loctype)) +
        ggplot2::geom_line() +
        ggplot2::xlab("Generation") +
        ggplot2::ylab("Allele frequency") +
        ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
        ggplot2::scale_color_manual(values = mycolors, labels = c("Intercept", "Plastic\nadjustment", "Epigenetic\nmodification", "Developmental\nvariation",  "Slope"))  
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
#' @param plasticity Whether plasticity simulations were used or not
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
myplot <-  ggplot2::ggplot(allele.long, aes(x = effect)) +
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
myplot <-  ggplot2::ggplot(allele.long, aes(x = effect)) +
            ggplot2::geom_histogram(aes(y = ..count../(2*N)), colour = "black", fill = "white", binwidth = 0.1) +
            ggplot2::geom_vline(aes(xintercept =0), linetype = "dashed") +
            ggplot2::xlab("Allelic effect") +
            ggplot2::ylab("Allele frequency") +
            ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
            ggplot2::facet_wrap(loctype ~ locus) }
    
   print(myplot) 
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
    if((b > 0.1 | b < -0.1) == TRUE) {
        if(data$G.adj[finalgen] > 0.1) { ctype <- "red" } else { ctype <- "blue" }
    }
    if(data$G.hed[finalgen] > 0.1) { ctype <- "green" }
    
    plotdata <- data.frame(cue = cue, phenotype = phenotype)
    
    ggplot2::ggplot(plotdata, aes(x = cue, y = phenotype)) +
        ggplot2::geom_line(color = ctype) +
        ggplot2::coord_cartesian(ylim = c(-1,1)) +
        ggplot2::xlab("Environmental cue") +
        ggplot2::ylab("Phenotype")    
}


################################################################

##################################################################
###         Popgen summary stats calculations                  ###
##################################################################

### Heterozygosity calculations

#Function to check whether a locus is heterozygous (works with infinite alleles model)
is.het <- function(ind.mat) {
   tol <- .Machine$double.eps^0.5
   return(abs(ind.mat[1] - ind.mat[2]) >= tol) }

##This function calculates heterozygosity for all loci (can be selected) and all individuals
## H = 1/NL sum_{i=1}^{N} sum_{j=1}^{L} H_{ij}
## where N = number of inds., L = number of loci, H_ij = whether ind i at locus j is het? 1 or 0
#'Calculate heterozygosity for selected loci and all individuals
#' @param ind.mat a list of individuals (genotypes) used in the calculation
#' @param loc.ind indexes of those loci used in the calculation
#' @param nloc The total number of loci used in the whole simulation
#' @export
het.calc <- function(ind.mat, loc.ind, nloc) {
   mydselect <- function(x, loc.ind, nloc) {x[1:2, (1:nloc %in% loc.ind)]}
   ind.mat <- lapply(ind.mat, mydselect, loc.ind, nloc) #Select only particular loci
####
   current.n <- length(ind.mat) #Number of individuals
   current.nloc <- dim(ind.mat[[1]])[2] #Number of loci in current selection
####
   calc.het.loci <- function(x) { sum(apply(x, 2, is.het)) }
####
   hets <- lapply(ind.mat, calc.het.loci)
####
   het.obs <- (1/(current.n*current.nloc))*sum(unlist(hets))
####
return(het.obs)
##
}

#####################################################################33

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
#' Generate 1/f noise
#' @param beta Noise parameter
#' @param generation Number of generations to generate
#' @param final.sd Final standard deviation
#' @export
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
