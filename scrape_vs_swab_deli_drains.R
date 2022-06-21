library(dada2); packageVersion("dada2")
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(reshape2)

#make sure to change the working directory

#The following lines of code are necessary for the pre-processing of the raw reads using 
# DADA2, this was done separately for swab and scrape samples 
#Now that these are done, can start at line 86 with the RDS files.

path <- "./datasets/biofilm_swabs_16S/"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
plotQualityProfile(fnFs[4:5])
plotQualityProfile(fnRs[4:5])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#trim reads so they overlap, and primers are excluded (trimLeft arguement) EE 2,2
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft=20) 

out

errF <- learnErrors(filtFs, multithread=TRUE, verbose = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, verbose = TRUE)

plotErrors(errF, nominalQ=TRUE)

#consider one of the pooling options here!
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo")

dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = "pseudo")

dadaFs

#merging
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample

head(mergers[[1]]) 

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(398,427)]
dim(seqtab2)
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

#The dada2 tutorial uses an old trainingset, Henk created a new trainset:
#path <- "~/Desktop/RDP/RDPClassifier_16S_trainsetNo16_rawtrainingdata"
#dada2:::makeTaxonomyFasta_RDP(file.path(path, "trainset16_022016.fa"), file.path(path, "trainset16_db_taxid.txt"), "~/tax/rdp_train_set_16.fa.gz")
#dada2:::makeSpeciesFasta_RDP("~/Desktop/RDP/current_Bacteria_unaligned.fa", "~/tax/rdp_species_assignment_16.fa.gz")

taxa <- assignTaxonomy(seqtab.nochim, '/Users/hendrikdenbakker/Documents/Faecal_19/16S_all/DADA2_RDP_files/rdp_train_set_18.fa.gz', multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/hendrikdenbakker/Documents/Faecal_19/16S_all/DADA2_RDP_files/rdp_species_assignment_18.fa.gz")
taxa.print <- taxa 
rownames(taxa.print)
head(taxa.print)

# save the frequency table (seqtab.nochim) and taxonomic classifications (taxa) of SVAs
# in R`s serialized format  
#saveRDS(taxa, "taxa_deli_surface_RDPdefaultEE.rds")
#saveRDS(seqtab.nochim, "table_deli_surface_RDPdefaultEE.rds")

############Start here to do analysis.################

#load data
seqtab.swab <- readRDS('table_deli_surface_RDPdefaultEE.rds')
tax.swab <- readRDS('taxa_deli_surface_RDPdefaultEE.rds')

#load metadata ("mapping file")
swab.mapping <- read.table('metadata_swabs.txt', sep = '\t', header = TRUE)

swab.mapping$SampleID <- as.factor(swab.mapping$Sample)
#Creating Phyloseq Object
rownames(swab.mapping) <- as.factor(swab.mapping$SampleID)
swab.mapping$City <- as.factor(swab.mapping$City)
ps.swab = phyloseq(otu_table(seqtab.swab, taxa_are_rows=FALSE), 
                   sample_data(swab.mapping), 
                   tax_table(tax.swab))
ps.swab

#load data
seqtab.biofilm <- readRDS('table_deli_biofilm_RDP_EEdefault.rds')
tax <- readRDS('taxa_deli_biofilm_RDP_EEdefault.rds')

#load metadata ("mapping file")
bio.mapping <- read.table('biofilm_metadata.txt', sep = '\t', header = TRUE)
bio.mapping
rownames(bio.mapping) <- as.factor(bio.mapping$SampleID)
bio.mapping$City <- as.factor(bio.mapping$City)
bio.mapping$SampleID <- as.factor(bio.mapping$SampleID)
ps.biofilm = phyloseq(otu_table(seqtab.biofilm, taxa_are_rows=FALSE), 
                      sample_data(bio.mapping), 
                      tax_table(tax))
ps.biofilm

#merge phyloseq objects 
ps.combined <- merge_phyloseq(ps.swab, ps.biofilm)
#replace mapping file
merged.mapping <- read.table('combined_samples_data.txt', sep = '\t', header = TRUE)
merged.mapping
rownames(merged.mapping) <- as.factor(merged.mapping$SampleID)
sample_data(ps.combined) <- merged.mapping

ps.combined

#Instructions for making bubble plots from https://jkzorz.github.io/2019/06/05/Bubble-plots.html
################bubble plot Fig 1: Shotgun vs 16S amplicon Phyla################
#scrape 16S vs shotgun on phylum level
ps.phylum.biofilm <-tax_glom(ps.biofilm, taxrank = 'Phylum', NArm = TRUE)
ps.phylum.biofilm <- transform_sample_counts(ps.phylum.biofilm, function(OTU) OTU/sum(OTU))
#no need to do the top taxa, there are only a few phyla. Transform counts to get percentages, not fractions
ps.biofilm.phylum_for_bubble <- transform_sample_counts(ps.phylum.biofilm, function(OTU) OTU*100)
#replace ASV names with taxonomic IDs
taxa_names(ps.biofilm.phylum_for_bubble) <- tax_table(ps.biofilm.phylum_for_bubble)[,"Phylum"]
#and write the abundance table to a csv, here commented out to avoid overwriting
#write.csv(otu_table(ps.biofilm.phylum_for_bubble), file = 'scrape_16S_phyla.csv') 

#create scrape_16S_and_shotgun_phyla.csv from scrape_16S_phyla.csv and
# merged_abundance_table_phylum.txt created with metaphlan script
pc = read.csv("scrape_16S_and_shotgun_phyla.csv", header = TRUE)
pc

#convert data frame from a "wide" format to a "long" format
pcm = melt(pc, id = c("sample", "sequence.approach"))

#maybe colorbrewer here?
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")

#keep order of samples
pcm$sample <- factor(pcm$sample,levels=unique(pcm$sample))

pcm



#plot it
xx_16S_shotgun_phyla_scrape = ggplot(pcm, aes(x = sample, y = variable)) + 
  geom_point(aes(size = value, fill = sequence.approach), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "sequnce.approach")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 9), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 11, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +  
  scale_fill_manual(values = c("darkorange", "skyblue"), guide = guide_legend(override.aes = list(size=5))) + 
  scale_y_discrete(limits = rev(levels(pcm$variable))) 

xx_16S_shotgun_phyla_scrape

################bubble plot Fig 2: Shotgun vs 16S amplicon Genera top30#########
#shotgun_vs_16s_forbubble.csv was created from the top 30 genera found in both 
# the 16S amplicon results and the metaphlan (shotgun metagenomic) results 
# (= merged_abundance_table_genus.txt)
pc = read.csv("shotgun_vs_16s_forbubble.csv", header = TRUE)
pc
library(ggplot2)
library(reshape2)

#convert data frame from a "wide" format to a "long" format
pcm = melt(pc, id = c("sample", "type"))

#maybe colorbrewer here?
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")

#keep order of samples
pcm$sample <- factor(pcm$sample,levels=unique(pcm$sample))

pcm

#plot it
xx = ggplot(pcm, aes(x = sample, y = variable)) + 
  geom_point(aes(size = value, fill = type), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "type")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 9), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 11, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +  
  scale_fill_manual(values = c("darkorange", "skyblue"), guide = guide_legend(override.aes = list(size=5))) + 
  scale_y_discrete(limits = rev(levels(pcm$variable))) 

xx


################bubble plot Fig 3: 16S amplicon, swab vs scarape Genera#########
#create 16S-swab vs 16S scrape from top 50
ps.top50.general <- transform_sample_counts(ps.combined.genus, function(OTU) OTU/sum(OTU))
top50.general<- names(sort(taxa_sums(ps.top50.general), decreasing=TRUE))[1:50]
ps.top50.general <- prune_taxa(top50.general, ps.top50.general)
ps.top50.general_bubble <- ps.top50.general
taxa_names(ps.top50.general_bubble) <- tax_table(ps.top50.general_bubble)[,"Genus"]
ps.top50.general_bubble<- transform_sample_counts(ps.top50.general_bubble, function(OTU) OTU*100)
#next commented out to avoid overwriting 
#write.csv(otu_table(ps.top50.general_bubble), file = 'swab_vs_scrape_16S_top50.txt') 

#open swab_vs_scrape_16S_top50.txt to add additional information and sort etc.

pc = read.csv("swab_vs_scrape_16S_top50.csv", header = TRUE)
pc

#convert data frame from a "wide" format to a "long" format
pcm = melt(pc, id = c("sample", "sample_type"))

#maybe colorbrewer here?
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")

#keep order of samples
pcm$sample <- factor(pcm$sample,levels=unique(pcm$sample))

pcm

#plot it
xx_16S_all = ggplot(pcm, aes(x = sample, y = variable)) + 
  geom_point(aes(size = value, fill = sample_type), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Time")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 9), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 11, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +  
  scale_fill_manual(values = c("darkorange", "skyblue"), guide = guide_legend(override.aes = list(size=5))) + 
  scale_y_discrete(limits = rev(levels(pcm$variable))) 

xx_16S_all

################Fig 4: 16S amplicon, alpha diversity############################
ps.combined
############Alpha Diversity############
plot_richness(ps.combined, x = "Type", measures=c("Shannon", "Simpson")) + geom_boxplot()

################Fig 5: 16S amplicon, beta diversity#############################
############Beta Diversity############
wh0.comb = genefilter_sample(ps.combined, filterfun_sample(function(x) x > 5), A=0.5*nsamples(ps.combined))
psc <- prune_taxa(wh0.comb, ps.combined)
psc <- transform_sample_counts(psc, function(x) 1E6 * x/sum(x))
nmds.comb <- ordinate(psc, method="NMDS", distance="bray")
plot_ordination(psc, nmds.comb, color="Type")
p = plot_ordination(psc, nmds.comb, type = "SampleID", color = "Type")
p + geom_point(size = 2) + stat_ellipse(type = "norm")

############STATS############
#Permanova analysis of beta diversity
meta <- as(sample_data(ps.combined), "data.frame")
meta
BC_n.dist <- vegdist(otu_table(ps.combined, taxa_are_rows=FALSE), distance="bray")
permanova <- adonis(BC_n.dist ~ Type, data = meta, permutations=1000, method = "bray")
#permanova <- adonis(t(otu) ~ Type, data = meta, permutations=99, method = "bray")
print(as.data.frame(permanova$aov.tab)["Type", "Pr(>F)"])


