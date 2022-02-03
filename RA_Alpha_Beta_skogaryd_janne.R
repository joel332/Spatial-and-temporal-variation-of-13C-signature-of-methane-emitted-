########################################################################################################################
##              Script used to analyse Taxonomy and Test for statistical differnece
##              written by Joel White 2022-02-03
##
########################################################################################################################

#load required libraries
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("readxl")
library("rlang")
library("dplyr")
library("tibble")
library("reshape2")


#Set working directory
setwd()

#Import data in xlsx format. 
otu_mat <- read_excel)
row.names(otu_mat) <- otu_mat$otu
otu_mat <- otu_mat %>% select (-otu)
otu_mat

tax_mat <- read_excel()
row.names(tax_mat) <- tax_mat$otu
tax_mat <- tax_mat %>% select (-otu) 
tax_mat

samples_df<- read_excel()
samples_df
row.names(samples_df) <- samples_df$Sample_ID_plot
samples_df


#OTU's are the absolute abundances of taxa. Columns are counts, Row 1 is OTU1, OTU2, OTU3 .....ect
otu_mat<- as.matrix(otu_mat)
class(otu_mat)


#Taxmat is the imported taxanomic heiarchy save in .txt format
tax_mat<- as.matrix(tax_mat)
class(tax_mat)

#samples_df is the imported taxanomic heiarchy save in .txt format
samples_df <- as.data.frame(samples_df)
head(sam_data)
class(sam_data)

#combine otumat and taxmat into a phyloseq object (Must be as data matrix)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
SAM = sample_data(samples_df)
head(OTU)
head(TAX)
head(SAM)

#plots abundance vs sample
physeq = phyloseq(OTU, TAX, SAM)
physeq

# Saving on object in RData format
save(physeq, file = "data.RData")


########################################################################
###            Data pruning                                          ###
########################################################################

# only look into high-abundance/high-prevelance OTUs over 10
physeq.f <- filter_taxa(physeq,function(x) sum(x >= 10) > (0.01*length(x)), 
                        prune = TRUE)
#oringinal
physeq
#post pruning
physeq.f

# Covert to relative abundance
physeq.f.ra <- transform_sample_counts(physeq.f, function(x) x*100/sum(x))
physeq.f.ra



############################################################################
###                 Plotting                                             ###
###########################################################################


a <- plot_bar(physeq.f.ra, fill="Species", x="Sample_ID") + facet_wrap(~Ecotype, scales = "free", nrow=1)+
      geom_bar(stat="identity", position="stack", colour="black")+
      ggforce::facet_row(vars(Ecotype), scales = 'free', space = 'free')+
      scale_fill_brewer(palette="Dark2") +
theme_bw()
a

#order the stack of species
a$data$Species <- factor(a$data$Species, levels = c("Type II methanotroph", "Type I methanotroph", "Verrucomicrobia", "Acetoclastic methanogen","Methylotrophic methanogen","Hydr/Methyl/Aceto methanogen","Hydrogenotrophic methanogen"))

#re-order the samples
#a$data$Sample_ID_plot <- factor(a$data$Sample_ID_plot, levels = c())

ra <- a + theme(axis.text.x = element_text(angle = 90),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "transparent",colour = NA),
                plot.background = element_rect(fill = "transparent",colour = NA))
ra

prop <- ra + labs(title = NULL) + ylab("Relative Abundance (%)")+ 
guides(fill=guide_legend(title="Functional group"))

prop

#save plot
ggsave("RA_ecotype.tiff", units="in", width=10, height=8, dpi=300, compression = 'lzw')


### ------------ Load data and package -----------------------------------------
### Load data and package for PERMANOVA
require(vegan)
library(dplyr)
library(rstatix)

myckel_sp <- read.csv()
head(myckel_sp)


myckel_env <- read.csv()
head(myckel_env)


### ------------ Distance matrix ----------------------------------------------
### Compute distance matrix using Bray-Curtis 

#use this to ensure transformed data lies close to between 0 - 10. This makes sure highly abundant species dont over influence the results
# otherwise highly abundant genes would dominate the distance measures
#check the range to ensure transformation id appropriate
range(myckel_sp) #original range
range(myckel_sp ^ 0.5) #Root transformation
range(myckel_sp ^ 0.25) # double root transformation


# We use a double root transformation to reduce the range of the data

dist_myckel <- vegdist(myckel_sp ^ 0.25, method = "bray")
dist_myckel


### ------------ PERMANOVA -----------------------------------------------------

#PERMANOVA Ecotype
pmv_ecotype <- adonis(
  myckel_sp ^ 0.25 ~ Ecotype, data = myckel_env,
  permutations = 999,
  method = "bray")
pmv_ecotype



library(RVAideMemoire)
#pairwise permtutation Post hoc test with FDR adjstment
Wilks_pairwise_ecotype <- pairwise.perm.manova(dist(myckel_sp,"euclidean"),myckel_env$Ecotype,nperm=999, test = "Wilks")
Wilks_pairwise_ecotype

#Save the details
results <- capture.output(print(pmv_ecotype), print(Wilks_pairwise_ecotype))
writeLines(results, con = file("output_PERMANOVA_&_pairwise_test_ecotype.txt"))

### ------------ SIMPER --------------------------------------------------------
sim <- simper(myckel_sp, group = myckel_env$Ecotype, permutations = 999)
summary <- summary(sim)
summary
# contr :   contribution to dissimilarity between upstream and downstream
# sd    :   standard deviation of contribution (is the species response consitent?)
# ratio :   ratio between contr and sd (high ratio = high, consisten contribution)
# av.   :   average abundance per groups
# cumsum:   cumulative contribution (rule of thumb : species till 70% are investigated)


#Save the details
taxa_simper_output <- capture.output(print(summary))
writeLines(taxa_simper_output, con = file())
