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
calc.genot.effects <- function(ind.mat, ae, type) {
    if(type == "biallelic") {
    geno.calc <- function(x, ae) {sum(rbind(x[1,]*ae, x[2,]*ae))} #This function sums genotypic effects
    ind.G <- as.vector(unlist(lapply(ind.mat, geno.calc, ae = ae)))
    }

    if(type == "infinite") {
        ind.G <- as.vector(unlist(lapply(ind.mat, sum)))
    }
    
    return(ind.G)
}


#Calculate genotypic effects, old function
#calc.genot.effects <- function(ind.mat, N, a) {
#    ind.G <- rep(0,N) #init
#    for(i in 1:N) {
#    ind.G[i] <- sum(ind.mat[[i]]*a) #Sum genotypic effects, assuming that all loci have the same allelic effects, for now
#    }
#    return(ind.G)
#}
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
    produce.gametes.W <- function(ind.mat, loc.mat, No, ind.W, link.map = F, linkage) {
        
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
    produce.gametes.W.No <- function(ind.mat, loc.mat, No, ind.W, link.map = F, linkage) {
        
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
            before.xo <- meiocyte[[i]][5,] < cur.pos
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
        mutate.K.alleles <- function(x, l) {
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
calc.al.freq <- function(data, N, nloc) {
    
    allele.freq <- rep(0, nloc) #initialise allele frequencies

    for(j in 1:nloc) {
        allele1 <- unname(sapply(data, '[[', 2*j-1)) #Get element [1,2*j-1] from each individual
        allele2 <- unname(sapply(data, '[[', 2*j)) #Get element [2,2*j] from each individual
        
        allele.freq[j] <- (sum(allele1+allele2))/(2*N)
    }
    return(allele.freq)
}
########################################################################################

#Calculate fitness using stabilizing selection 
    calc.fitness <- function(ind.P, opt.pheno, sel.intensity) {
        exp(-((ind.P-opt.pheno)^2)/(2*sel.intensity^2))
    }
#################################################################

#The string #' indicates roxygen parsing

#This function plots an overview of the simulation results
#' Plot simulation results
#'
#' This function plots population size, trait mean, and population mean fitness of the simulation
#'
#' @param data Result matrix given by a another function (link to be completed)
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
plot_allele.results <- function(data, nloc) {

    #First need to transform data into long format for ggplot
    data <- data.frame(data) #Need to convert to data frame before melt
    data.long <- reshape2::melt(data, id.vars = c("generation"), measure.vars = c(paste(rep("locus", nloc), 1:nloc, sep = "")), variable.name = "locus", value.name = "freq")

    ggplot2::ggplot(data.long, aes(x = generation, y = freq, group = locus)) +
        ggplot2::geom_line() +
        ggplot2::xlab("Generation") +
        ggplot2::ylab("Allele frequency") +
        ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0))

}
            
#This function plots an overview of alleles present in the population
#Only works for the infinite alleles model
plot_inf.alleles <- function(data, nloc) {

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
    
    #Plot with ggplot2
    ggplot2::ggplot(allele.long, aes(x = effect)) +
      ggplot2::geom_histogram(aes(y = (..ncount..)), colour = "black", fill = "white", binwidth = 0.1) +
      ggplot2::geom_vline(aes(xintercept =0), linetype = "dashed") +
      ggplot2::xlab("Allelic effect") +
      ggplot2::ylab("Allele frequency") +
      ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      ggplot2::facet_wrap(~ locus)
    
}

    
################################################################

###########################################################
# Generating different colours of environmental variation #
###########################################################

### Autoregressive (AR) process ###
#E_{t+1} = kappa*E_t + omega_t*sqrt(1-kappa^2)
#When kappa < 0 -> blue noise, kappa > 0 -> red noise, and kappa = 0 -> white noise
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
