########crAss001 structure study########
#Load required libraries
setwd("~/Downloads/Rsendover")
library(RColorBrewer) #Palettes
library(data.table) #Working with large dataset files
library(ggplot2)
library(gggenes) #ggplot2 extension to plot out genetic maps
library(igraph) #Graphs
library(ape) #Phylogenetic analysis
library(ggforce)
library(tidytext)
library(tidyverse)
library(patchwork)
library(ggdendro)
library(ggnewscale)
library(ggupset)
set.seed(1)
options(stringsAsFactors = FALSE)
########Set Working directory, read in data########
contig_clustering_all <- read.csv("contig_clustering_all.csv", row.names = 1, header = TRUE) #Change path!
contig_clustering <- read.csv("contig_clustering.csv", row.names = 1, header = TRUE) #Change path!
contig_annot_gff <- read.delim("C:/Users/Navi/Downloads/Rsendover/contig_annot_gff.txt")
orthologs <- read.table("orthologs.txt", row.names = 1, header = TRUE) #Change path!
#Convert families and subfamilies to factors with correct order of levels
contig_clustering$nomen_family <- factor(contig_clustering$nomen_family,
                                         levels = c("Intestiviridae", "Crevaviridae", "Suoliviridae", "Steigviridae",
                                                    "Tinaiviridae", "Jelitoviridae", "Marine_1", "Marine_2"))
contig_clustering_all$nomen_family <- factor(contig_clustering_all$nomen_family,
                                         levels = c("Intestiviridae", "Crevaviridae", "Suoliviridae", "Steigviridae",
                                                    "Tinaiviridae", "Jelitoviridae", "Marine_1", "Marine_2"))
contig_clustering$nomen_subfamily <- factor(contig_clustering$nomen_subfamily,
                                         levels = c("Crudevirinae", "Coarsevirinae", "Densevirinae", "Churivirinae",
                                                    "Doltivirinae", "Asinivirinae", "Loutivirinae", "Oafivirinae",
                                                    "Obtuvirinae", "Boorivirinae", "Bearivirinae", "Uncouvirinae",
                                                    "Lumpivirinae", "Grossvirinae", "Marine_1", "Marine_2"))
contig_clustering_all$nomen_subfamily <- factor(contig_clustering_all$nomen_subfamily,
                                            levels = c("Crudevirinae", "Coarsevirinae", "Densevirinae", "Churivirinae",
                                                       "Doltivirinae", "Asinivirinae", "Loutivirinae", "Oafivirinae",
                                                       "Obtuvirinae", "Boorivirinae", "Bearivirinae", "Uncouvirinae",
                                                       "Lumpivirinae", "Grossvirinae", "Marine_1", "Marine_2"))

########Proteins of interest and clusters they belong to########
prots_of_interest <- c("Gp20","Gp21","Gp22","Gp23","Gp24","Gp25","Gp26","Gp29","Gp34","Gp35","Gp40","Gp43","Gp44","Gp38","Gp39")
homologs <- list() #Empty list to store homologs for each protein
genomes <- contig_annot_gff #A copy of 'gff' table
genomes <- genomes[genomes$feature == "CDS",]
#Add in family and subfamily information
genomes$family <- contig_clustering[tolower(genomes$sequence),"nomen_family"]
genomes$subfamily <- contig_clustering[tolower(genomes$sequence),"nomen_subfamily"]
#Filter out unclassified sequences and non-type strains
genomes <- genomes[!is.na(genomes$family),]

#Append homologs to the gff file
genomes$gene_of_interest <- NA
for(i in prots_of_interest) {
  homologs[[i]] <- contig_annot_gff[which(
    contig_annot_gff$orth_group == contig_annot_gff[grep(paste0("Note=",i), contig_annot_gff$attributes),"orth_group"]
  ),"protein_id"]
  genomes[genomes$protein_id %in% homologs[[i]],"gene_of_interest"] <- i
}
genomes <- genomes[genomes$sequence %in% unique(genomes[!is.na(genomes$gene_of_interest),"sequence"]),] #Forgot why I needed this...
length(genomes$gene_of_interest[!is.na(genomes$gene_of_interest)])
prot_all_vs_all_blastp <- fread(file = "all-vs-all.tsv") #Change path!
colnames(prot_all_vs_all_blastp) <-
  c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
prot_all_vs_all_blastp <- prot_all_vs_all_blastp[evalue < 1e-5] #Filter by e-value
prot_all_vs_all_blastp <- prot_all_vs_all_blastp[qseqid != sseqid] #Remove hits on query protein itself
#Select for hits relevant to our proteins of interest
prot_of_interest_blastp <- prot_all_vs_all_blastp[qseqid %in% unique(as.character(unlist(homologs))) | sseqid %in% unique(as.character(unlist(homologs)))]

########Convert BLASTp output into a graph########
prot_of_interest_graph <- graph_from_data_frame(prot_of_interest_blastp[,.(qseqid,sseqid)], directed = FALSE) #Non-directed graph from a dataframe
#Add vertex attributes
V(prot_of_interest_graph)$cluster <-
  orthologs[names(V(prot_of_interest_graph)),"L1"] #Protein clusters
V(prot_of_interest_graph)$subfamily <-
  as.character(contig_clustering_all[orthologs[names(V(prot_of_interest_graph)),"genome"],"nomen_subfamily"]) #Phage subfamilies
V(prot_of_interest_graph)$homolog <- NA #Homologs of proteins of interest
for(i in names(homologs)) {
  V(prot_of_interest_graph)$homolog[names(V(prot_of_interest_graph)) %in% homologs[[i]]] <- i
}
prot_of_interest_df <- igraph::as_data_frame(prot_of_interest_graph, what = "vertices")

#Look for cases when BLAST hits cross from one protein cluster to a different one ("mismatched hits")
prot_of_interest_blastp$mismatched <- FALSE
prot_of_interest_blastp[which(
  prot_of_interest_df[prot_of_interest_blastp$qseqid,"cluster"] != prot_of_interest_df[prot_of_interest_blastp$sseqid,"cluster"]
),"mismatched"] <- TRUE
#Add edge attributes
E(prot_of_interest_graph)$mismatched <- prot_of_interest_blastp$mismatched #"mismatched hits"
E(prot_of_interest_graph)$length <- prot_of_interest_blastp$length #BLAST alignment length
E(prot_of_interest_graph)$pident <- prot_of_interest_blastp$pident #% identity
E(prot_of_interest_graph)$weight <-
  log10(E(prot_of_interest_graph)$length*(E(prot_of_interest_graph)$pident/100)) #Assign edge weight as log-transformed product of length and identity

prot_of_interest_graph <-
  igraph::simplify(prot_of_interest_graph, edge.attr.comb = "sum") #Simplify graph, edge weights and other attributes are summed up
E(prot_of_interest_graph)$mismatched <- as.logical(E(prot_of_interest_graph)$mismatched)
prot_of_interest_layout <- layout_with_kk(prot_of_interest_graph)
#Cluster graph into connected subgraph - "clusters" again, very confusing... Let's call them subgraphs instead
prot_of_interest_clust <- clusters(prot_of_interest_graph)
V(prot_of_interest_graph)$subgraphs <- prot_of_interest_clust$membership[V(prot_of_interest_graph)$name]
genomes$gene_of_interest[genomes$orth_group == "cr63"] = "Gp44"
genomes$gene_of_interest[genomes$orth_group == "cr209"] = "Gp44"
genomes$gene_of_interest[genomes$orth_group == "cr35"] = "Gp44"
genomes$gene_of_interest[genomes$orth_group == "cr165"] = "Gp44"
genomes$gene_of_interest[genomes$orth_group == "cr6"] = "Gp20"



genomes$gene_of_interest[genomes$protein_id %in% V(prot_of_interest_graph)$name[V(prot_of_interest_graph)$subgraphs ==1]] ="Gp20"
genomes$gene_of_interest[genomes$protein_id %in% V(prot_of_interest_graph)$name[V(prot_of_interest_graph)$subgraphs ==2]] ="Gp21-22-23-25-26-29"
genomes$gene_of_interest[genomes$protein_id %in% V(prot_of_interest_graph)$name[V(prot_of_interest_graph)$subgraphs ==3]] ="Gp40"
genomes$gene_of_interest[genomes$protein_id %in% V(prot_of_interest_graph)$name[V(prot_of_interest_graph)$subgraphs ==4]] ="Gp44"
genomes$gene_of_interest[genomes$protein_id %in% V(prot_of_interest_graph)$name[V(prot_of_interest_graph)$subgraphs ==5]] ="Gp43"
genomes$gene_of_interest[genomes$protein_id %in% V(prot_of_interest_graph)$name[V(prot_of_interest_graph)$subgraphs ==6]] ="Gp39"
genomes$gene_of_interest[genomes$protein_id %in% V(prot_of_interest_graph)$name[V(prot_of_interest_graph)$subgraphs ==7]] ="Gp38"
genomes$gene_of_interest[genomes$protein_id %in% V(prot_of_interest_graph)$name[V(prot_of_interest_graph)$subgraphs ==8]] ="Gp35"
unique(genomes$sequence[genomes$gene_of_interest =="Gp20"])
rm(list = c("prot_all_vs_all_blastp", "prot_of_interest_blastp", "prot_of_interest_clust", "prot_of_interest_df", "prot_of_interest_graph", "prot_of_interest_layout"))
gc()
#my laptop has no space so we clean after ourselves in this house
seqlength = data.frame(sequence = rownames(contig_clustering_all), sequence_size= contig_clustering_all$length)

genomes = merge(genomes, seqlength)
for (seq in unique(genomes$sequence)){
  gp20 = genomes[genomes$sequence == seq & genomes$gene_of_interest == "Gp20" & !is.na(genomes$gene_of_interest),]
  if (nrow(gp20) >1){ 
    gp20 = gp20[1,]}
  if (length(gp20[,2]) ==0 ){genomes = genomes[genomes$sequence !=seq,] }
  else{
    if (gp20$strand == -1){
      genomes$dif[genomes$sequence == seq] = sqrt((genomes$end[genomes$sequence == seq] - genomes$start[genomes$sequence == seq])^2)
      genomes$start[genomes$sequence == seq] = genomes$sequence_size[genomes$sequence == seq] - genomes$end[genomes$sequence == seq]
      genomes$end[genomes$sequence == seq] = genomes$dif[genomes$sequence == seq] + genomes$start[genomes$sequence == seq]
      genomes[genomes$sequence  == seq,"strand"] = genomes[genomes$sequence  == seq,"strand"] * -1}
    genomes[genomes$sequence == seq, "start"] = genomes[genomes$sequence == seq, "start"] - gp20$start
    genomes[genomes$sequence == seq, "end"] = genomes[genomes$sequence == seq, "end"] - gp20$start
    genomes[genomes$sequence == seq & genomes$start < 0, "end"] =   genomes[genomes$sequence == seq & genomes$start < 0, "end"] + seqlength$sequence_size[seqlength$sequence==seq]
    genomes[genomes$sequence == seq & genomes$start < 0, "start"] =   genomes[genomes$sequence == seq & genomes$start < 0, "start"] + seqlength$sequence_size[seqlength$sequence==seq]
    genomes[genomes$sequence == seq & genomes$end < 0, "end"] =  genomes[genomes$sequence == seq & genomes$end < 0, "end"]  + seqlength$sequence_size[seqlength$sequence==seq]
    genomes[genomes$sequence  == seq & genomes$strand == -1,"strand"] = genomes[genomes$sequence  == seq & genomes$strand == -1,"strand"] +1
    
  }}
for (seq in unique(genomes$sequence)){
  gp20 = genomes[genomes$sequence == seq & genomes$gene_of_interest == "Gp20" & !is.na(genomes$gene_of_interest),]
  if (nrow(gp20) >1){ 
    gp20 = gp20[1,]}
  if ((gp20$end - gp20$start)> 4000){
    genomes[genomes$sequence == seq & genomes$protein_id == gp20$protein_id & !is.na(genomes$orth_group),"end"] = genomes[genomes$sequence == seq & genomes$protein_id == gp20$protein_id & !is.na(genomes$orth_group),"end"] +genomes[genomes$sequence == seq,"sequence_size"][1]
    tempstart = genomes[genomes$sequence == seq & genomes$orth_group != "cr6" & !is.na(genomes$orth_group),"end"]
    tempend = genomes[genomes$sequence == seq & genomes$orth_group != "cr6" & !is.na(genomes$orth_group),"start"]
    genomes[genomes$sequence == seq & genomes$orth_group != "cr6" & !is.na(genomes$orth_group),"end"] = tempstart
    genomes[genomes$sequence == seq & genomes$orth_group != "cr6" & !is.na(genomes$orth_group),"start"] = tempend
    }
}
#this needs to run twice for some reason, just trust the plan
for (seq in unique(genomes$sequence)){
  gp20 = genomes[genomes$sequence == seq & genomes$gene_of_interest == "Gp20" & !is.na(genomes$gene_of_interest),]
  if (nrow(gp20) >1){ 
    gp20 = gp20[2,]}
  if (length(gp20[,2]) ==0 ){genomes = genomes[genomes$sequence !=seq,] }
  else{
    if (gp20$strand == -1){
      genomes$dif[genomes$sequence == seq] = sqrt((genomes$end[genomes$sequence == seq] - genomes$start[genomes$sequence == seq])^2)
      genomes$start[genomes$sequence == seq] = genomes$sequence_size[genomes$sequence == seq] - genomes$end[genomes$sequence == seq]
      genomes$end[genomes$sequence == seq] = genomes$dif[genomes$sequence == seq] + genomes$start[genomes$sequence == seq]
      genomes[genomes$sequence  == seq,"strand"] = genomes[genomes$sequence  == seq,"strand"] * -1}
    genomes[genomes$sequence == seq, "start"] = genomes[genomes$sequence == seq, "start"] - gp20$start
    genomes[genomes$sequence == seq, "end"] = genomes[genomes$sequence == seq, "end"] - gp20$start
    genomes[genomes$sequence == seq & genomes$start < 0, "end"] =   genomes[genomes$sequence == seq & genomes$start < 0, "end"] + seqlength$sequence_size[seqlength$sequence==seq]
    genomes[genomes$sequence == seq & genomes$start < 0, "start"] =   genomes[genomes$sequence == seq & genomes$start < 0, "start"] + seqlength$sequence_size[seqlength$sequence==seq]
    genomes[genomes$sequence == seq & genomes$end < 0, "end"] =  genomes[genomes$sequence == seq & genomes$end < 0, "end"]  + seqlength$sequence_size[seqlength$sequence==seq]
    genomes[genomes$sequence  == seq & genomes$strand == -1,"strand"] = genomes[genomes$sequence  == seq & genomes$strand == -1,"strand"] +1
    
  }}
for (seq in unique(genomes$sequence)){
  gp20 = genomes[genomes$sequence == seq & genomes$gene_of_interest == "Gp20" & !is.na(genomes$gene_of_interest),]
  if (nrow(gp20) >1){ 
    gp20 = gp20[1,]}
  if ((gp20$end - gp20$start)> 4000){
    genomes[genomes$sequence == seq & genomes$protein_id == gp20$protein_id & !is.na(genomes$orth_group),"end"] = genomes[genomes$sequence == seq & genomes$protein_id == gp20$protein_id & !is.na(genomes$orth_group),"end"] +genomes[genomes$sequence == seq,"sequence_size"][1]
    tempstart = genomes[genomes$sequence == seq & genomes$orth_group != "cr6" & !is.na(genomes$orth_group),"end"]
    tempend = genomes[genomes$sequence == seq & genomes$orth_group != "cr6" & !is.na(genomes$orth_group),"start"]
    genomes[genomes$sequence == seq & genomes$orth_group != "cr6" & !is.na(genomes$orth_group),"end"] = tempstart
    genomes[genomes$sequence == seq & genomes$orth_group != "cr6" & !is.na(genomes$orth_group),"start"] = tempend
  }
}
#some graphical variables
genomes$opac = genomes$gene_of_interest
genomes$opac[is.na(genomes$gene_of_interest)] = 0
genomes$opac[! is.na(genomes$gene_of_interest)] = 1
genomes$opac = as.numeric(genomes$opac)
genomes$Colour = genomes$orth_group
genomes$Colour[! is.na(genomes$gene_of_interest)] = "GOI"
genomesspecies <- genomes # allow more than just type strains
genomes <- genomes[contig_clustering[tolower(genomes$sequence),"type_strain"],]
write.table(genomes, "genomes.tsv", sep = "\t", row.names = FALSE)
write.table(genomesspecies, "genomesspecies.tsv", sep = "\t", row.names = FALSE)
# ░░░░░░░░▄▄▄▀▀▀▄▄███▄░░░░░░░░░░░░░░
# ░░░░░▄▀▀░░░░░░░▐░▀██▌░░░░░░░░░░░░░
# ░░░▄▀░░░░▄▄███░▌▀▀░▀█░░░░░░░░░░░░░
# ░░▄█░░▄▀▀▒▒▒▒▒▄▐░░░░█▌░░░░░░░░░░░░
# ░▐█▀▄▀▄▄▄▄▀▀▀▀▌░░░░░▐█▄░░░░░░░░░░░
# ░▌▄▄▀▀░░░░░░░░▌░░░░▄███████▄░░░░░░
# ░░░░░░░░░░░░░▐░░░░▐███████████▄░░░
# ░░░░░le░░░░░░░▐░░░░▐█████████████▄
# ░░░░toucan░░░░░░▀▄░░░▐█████████████▄ 
# ░░░░░░has░░░░░░░░▀▄▄███████████████ 
# ░░░░░arrived░░░░░░░░░░░░█▀██████░░
#Now that I have your attention 
#You're going to have to do this via python nametry3.py will give you the next file
###
genomesnamed = read.delim("genomesnamed.tsv", sep = "\t")
genomesspeciesnamed = read.delim("genomesspeciesnamed.tsv", sep = "\t")
#creating a copy of the full genomesize so we can cut genomesnamed
genomesfull = genomesnamed
for (seq in genomesfull$sequence){

  gp44 = subset(genomesfull, sequence == seq)
  gp44 = subset(gp44, gene_of_interest == "Gp44")
  #removes all sequences past gp44
  if (dim(gp44)[1] == 1){
  genomesnamed <- subset(genomesnamed, sequence!=seq | (sequence == seq & end <= gp44$end))
  }
}
genomesspeciesfull = genomesspeciesnamed

for (seq in unique(genomesspeciesfull$sequence)){
  
  gp44 = subset(genomesspeciesfull, sequence == seq)
  gp44 = subset(gp44, gene_of_interest == "Gp44")
  #removes all sequences past gp44
  if (dim(gp44)[1] == 1){
  genomesspeciesnamed <- subset(genomesspeciesnamed, sequence!=seq | (sequence == seq & end <= gp44$end))}
  else{
    genomesspeciesnamed <- subset(genomesspeciesnamed, sequence!=seq )
  }
  }
#more cleaning 
rm(list = c("genomes", "genomesspecies", "genomesfull", "genomesspeciesfull"))
gc()
write.table(genomesspeciesnamed, "gp44species.tsv",sep = "\t", row.names = FALSE)
########clustering######## 
#adding the full names here, as they're clustered on this later and this just makes it easier
#this is an inefficent way of doing it no doubt but I'll fix it later
genomesnamed %>% 
  mutate(full_name = paste("Ph:", Genus, "host:", Host)) %>%
  mutate(sequence = paste(sequence, full_name, sep = " ")) %>%
  dplyr::select(c(sequence, orth_group)) %>% 
  group_by(sequence, orth_group) %>% 
  summarise_all(funs(sum))%>% 
  filter(orth_group != "Oneoff") %>%
  filter(!is.na(orth_group)) %>%
  add_column(seen = 1) %>%
  pivot_wider(names_from = orth_group, values_from = seen, values_fill = 0) %>% 
  arrange(sequence) %>%
  column_to_rownames("sequence")-> boega 
dist_boega = dist(boega, method = "binary")
clust_boega = hclust(dist_boega)
#rcreating the counttables
genomesnamed %>% 
  dplyr::select(c(sequence, orth_group)) %>% 
  group_by(sequence, orth_group) %>% 
  summarise_all(funs(sum))%>% 
  filter(orth_group != "Oneoff") %>%
  filter(!is.na(orth_group)) %>%
  add_column(seen = 1) %>%
  pivot_wider(names_from = orth_group, values_from = seen, values_fill = 0) %>% 
  arrange(sequence) %>%
  column_to_rownames("sequence")-> counttable
####
genomesspeciesnamed %>% 
  mutate(full_name = paste("Ph:", Genus, "host:", Host)) %>%
  mutate(sequence = paste(sequence, full_name, sep = " ")) %>%
  dplyr::select(c(sequence, orth_group)) %>% 
  group_by(sequence, orth_group) %>% 
  summarise_all(funs(sum))%>% 
  filter(orth_group != "Oneoff") %>%
  filter(!is.na(orth_group)) %>%
  add_column(seen = 1) %>%
  pivot_wider(names_from = orth_group, values_from = seen, values_fill = 0) %>% 
  arrange(sequence) %>%
  column_to_rownames("sequence")-> speciesboega
species_dist_boega = dist(speciesboega, method = "binary")
species_clust_boega = hclust(species_dist_boega)
genomesspeciesnamed %>% 
  dplyr::select(c(sequence, orth_group)) %>% 
  group_by(sequence, orth_group) %>% 
  summarise_all(funs(sum))%>% 
  filter(orth_group != "Oneoff") %>%
  filter(!is.na(orth_group)) %>%
  add_column(seen = 1) %>%
  pivot_wider(names_from = orth_group, values_from = seen, values_fill = 0) %>% 
  arrange(sequence) %>%
  column_to_rownames("sequence")-> speciescounttable 
###
#removing all singletons as they're just clutter and changing the rownames for ease
counttable = rownames_to_column(counttable[colSums(counttable)>1], "sequence")
counttable <- pivot_longer(counttable, cols = !sequence, names_to = "mcl", values_to = "presence")
counttable <- unique(left_join(genomesnamed[ , c("sequence", "family", "subfamily", "Genus", "Host")], counttable))
counttable <- merge(y = counttable, x = genomesnamed[ , c("sequence", "family", "subfamily", "Genus", "Host")], by = "sequence", all.x=TRUE)
counttable %>%
mutate(full_name = paste("Ph:", Genus.x, "host:", Host.x)) %>%
  mutate(sequence = paste(sequence, full_name, sep = " ")) -> counttable
#ordering the x and y axis based on likeness and accurance respectively 
counttable %>%
  mutate(sequence = factor(counttable$sequence, levels = clust_boega$labels[clust_boega$order])) %>% unique()-> counttable
{data.table(counttable[,c("sequence", "mcl", "presence", "family.x")])} %>%
  pivot_wider( values_from = presence, names_from = mcl) -> familyordered
DT = data.table(counttable[,c("mcl", "presence")])
ordered = DT[ , .(Totalcount = sum(presence)), by = .(mcl)]
ordered = ordered[order(Totalcount, decreasing = TRUE)]
counttable$mcl = factor(counttable$mcl, levels = ordered$mcl)

###
speciescounttable = rownames_to_column(speciescounttable[colSums(speciescounttable)>1], "sequence")
speciescounttable <- pivot_longer(speciescounttable, cols = !sequence, names_to = "mcl", values_to = "presence")
speciescounttable <- left_join(genomesspeciesnamed[ , c("sequence", "family", "subfamily", "Genus", "Host")], y = speciescounttable,)
speciescounttable = unique(speciescounttable)
speciescounttable %>%
  mutate(full_name = paste("Ph:", Genus, "host:", Host)) %>%
  mutate(sequence = paste(sequence, full_name, sep = " ")) -> speciescounttable
speciescounttable %>%
  mutate(sequence = factor(speciescounttable$sequence, levels = species_clust_boega$labels[species_clust_boega$order])) -> speciescounttable
speciesnumtable = speciescounttable[,c("mcl", "presence")]
speciesDT = data.table(speciesnumtable)
speciesordered = speciesDT[ , .(Totalcount = sum(presence)), by = .(mcl)]
speciesordered = speciesordered[order(Totalcount, decreasing = TRUE)]
speciescounttable$mcl = factor(speciescounttable$mcl, levels = speciesordered$mcl)
###
rm(list = c("speciesnumtable", "numtable", "DT", "FT", "speciesDT", "speciesFT", "genomesfull", "genomesspeciesfull"))
gc()
##### plotting the heat map
pdf("Big heat map", width = 22, height = 17)
counttable%>%
  group_by(family.x) %>%
  filter(sum(presence)> 0)%>%

ggplot( aes(y = sequence, x=mcl, fill= cut(presence, breaks=0:1))) +
  geom_tile() +
  scale_fill_brewer(name="Prescence", labels=c("Present", "Absent"), type="qual", palette = "Dark2")  +
  facet_wrap(family.x~., ncol = 1, scales = "free_y")+
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  dev.off()
#Now for the bigger one
  
  pdf("Bigger heat map", width = 22, height = 17)
speciescounttable%>%
  group_by(family) %>%
  filter(sum(presence)> 0)%>%
  
  ggplot( aes(y = sequence, x=mcl, fill= cut(presence, breaks=0:1) )) +
  geom_tile() +
  scale_fill_brewer(name="Prescence", labels=c("Present", "Absent"), type="qual", palette = "Dark2")  +
  facet_wrap(family~., ncol = 1, scales = "free_y")+
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_blank())
dev.off()  
## upsetgraphs

ldatest <- counttable[, c("family.x", "mcl", "presence")]
unique(ldatest) %>%
  pivot_wider(names_from = mcl, values_from = family.x, values_fn = list ) %>%
  pivot_longer(!presence, names_to = "mcl", values_to = "family.x") -> longer
pdf("Big upset plot", width = 22, height =17)
longer[longer$presence == 1, c("mcl", "family.x")] %>%
  ggplot(aes(x=family.x)) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  geom_bar() +
  scale_x_upset(n_intersections = 400)  
dev.off()
##
speciesldatest <- speciescounttable[, c("family", "mcl", "presence")]
unique(speciesldatest) %>%
  pivot_wider(names_from = mcl, values_from = family, values_fn = list ) %>%
  pivot_longer(!presence, names_to = "mcl", values_to = "family") -> specieslonger
pdf("Bigger upset plot", width = 22, height =17)
specieslonger[specieslonger$presence == 1, c("mcl", "family")] %>%
  ggplot(aes(x=family)) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  geom_bar() +
  scale_x_upset(n_intersections = 400)  
dev.off()
### coloured heatmaps
longer%>%
  unnest_wider(family.x) %>%
  filter(presence ==1 & is.na(...2) ) %>%
  rename(colour = ...1) %>%
  select(mcl, colour) %>%
  drop_na()-> merging

merge(y = counttable, x = merging, by=  "mcl", all = TRUE) -> view
unique(view$colour)
view$fill = 0
view[which(view$colour == "Suoliviridae"), "fill"] =1
view[which(view$colour == "Intestiviridae"), "fill"] =2
view[which(view$colour == "Steigviridae"), "fill"] =3
view[which(view$colour == "Crevaviridae"), "fill"] =4
view$mcl = factor(view$mcl, levels = ordered$mcl)

pdf("Big selected heatmap", width= 21, height = 17)
view%>%
  group_by(family.x) %>%
  filter(sum(presence)> 0)%>%
  replace_na(replace = list(none = "colour")) %>%
  ggplot( aes(y = sequence, x=mcl, alpha= presence, fill= colour)) +
  geom_tile() +
  scale_fill_brewer(name="Specific to", na.translate = F, type="qual", palette = "Dark2")  +
  facet_wrap(family.x~., ncol = 1, scales = "free_y")+
  scale_alpha(range = c(0.1, 1), guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_blank())
dev.off()

##
rm(list = c( "view", "speciesldatest"))
specieslonger%>%
  unnest_wider(family) %>%
  filter(presence ==1 & is.na(...2) ) %>%
  rename(colour = ...1) %>%
  select(mcl, colour) %>%
  drop_na()-> speciesmerging

merge(y = speciescounttable, x = speciesmerging, by=  "mcl", all = TRUE) -> speciesview
speciesview = unique(speciesview)
rm(list = c("specieslonger", "speciesmerging"))
gc()
speciesview$mcl = factor(speciesview$mcl, levels = speciesordered$mcl)
write_tsv(speciesview, "speciesview.tsv")
pdf("Bigger selected heatmap", width= 21, height = 17)

speciesview%>%
  group_by(family) %>%
  filter(sum(presence)> 0)%>%
  replace_na(replace = list(none = "colour")) %>%
  ggplot( aes(y = sequence, x=mcl, alpha= presence, fill= colour)) +
  geom_tile() +
  scale_fill_brewer(name="Specific to", na.translate = F, type="qual", palette = "Dark2")  +
  facet_wrap(family~., ncol = 1, scales = "free_y")+
  scale_alpha(range = c(0.1, 1), guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_blank())
dev.off()

##
speciesordered = speciesordered[speciesordered$Totalcount >= 800]
speciescounttable$mcl = factor(speciescounttable$mcl, levels = speciesordered$mcl)
speciescounttable = speciescounttable[speciescounttable$mcl %in% speciesordered$mcl,]
pdf("Bigger heat mapcut off", width = 22, height = 17)
speciescounttable%>%
  group_by(family) %>%
  filter(sum(presence)> 0)%>%
  
  ggplot( aes(y = sequence, x=mcl, fill= cut(presence, breaks=0:1) )) +
  geom_tile() +
  scale_fill_brewer(name="Prescence", labels=c("Present", "Absent"), type="qual", palette = "Dark2")  +
  facet_wrap(family~., ncol = 1, scales = "free_y")+
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_blank())
dev.off()  

#genomesnamed$sequence[clust_boega$order]
#genomessorted = genomesnamed[order("sequence"),]

#p1 =ggplot(segment(ddata1)) +
#  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#  coord_flip() +
#  scale_y_reverse(expand = c(0.2, 0))
#p2 =ggplot(segment(ddata2)) +
#  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#  coord_flip() +
#  scale_y_reverse(expand = c(0.2, 0))
#p3 =ggplot(segment(ddata3)) +
#  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#  coord_flip() +
#  scale_y_reverse(expand = c(0.2, 0))
#p4 =ggplot(segment(ddata4)) +
#  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#  coord_flip() +
#  scale_y_reverse(expand = c(0.2, 0))

  genomesnamed$gene_of_interest[genomesnamed$orth_group == "cr2"] = "Gp32"
  genomesnamed$Colour[genomesnamed$orth_group == "cr2"] = NA
  genomesnamed$opac[!is.na(genomesnamed$gene_of_interest)] =1
  genomesnamed$gene_of_interest[is.na(genomesnamed$gene_of_interest)]= "None"
  prot_all_vs_all_blastp <- fread(file = "all-vs-all.tsv")
  colnames(prot_all_vs_all_blastp) <-
    c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  prot_all_vs_all_blastp <- prot_all_vs_all_blastp[evalue < 1e-5] #Filter by e-value
  prot_all_vs_all_blastp <- prot_all_vs_all_blastp[qseqid != sseqid] #Remove hits on query protein itself
  prot_all_vs_all_blastp <- prot_all_vs_all_blastp[length >= 100]
  genomesnamed %>%
    mutate(full_name = paste("Ph:", Genus, "host:", Host)) %>%
    mutate(sequence = paste(sequence, full_name, sep = " ")) %>%
    mutate(sequence = factor(sequence, levels = clust_boega$labels[clust_boega$order])) ->genomesnamed_1
  genomesnamed_1 %>%
  ggplot(aes(xmin = start, xmax = end, y = sequence, forward = as.logical(strand), fill =Colour, colour =gene_of_interest, size =opac)) +
    scale_size_continuous(range = c(0.1, 1.5))+
    geom_gene_arrow(arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm")) +
    geom_text(inherit.aes = FALSE,aes(x = end - ((end-start)/2), y = sequence, label = orth_group,  size=1)) +
    #geom_gene_label(inherit.aes = FALSE,aes(xmin = start, xmax = end, y = sequence,label=orth_group, color = "#000000"))+
    facet_col(~family, scales = "free_y",space = "free", strip.position = "top") +
    scale_fill_discrete(guide = "none") +
    scale_y_discrete(position = "right") +
    theme_bw() ->plot
  d <- setDT(layer_data(plot))
  d$protein_id = genomesnamed_1$protein_id
  d2 <-d[,c("protein_id", "y", "xmin", "xmax", "PANEL", "forward")]
  colnames(d2) <-c("sprotein_id", "sy", "sxmin", "sxmax", "PANEL", "forward")
  d[d$forward == FALSE, c("xmin", "xmax")] <- d[d$forward == FALSE, c("xmax", "xmin")]
  d2[d2$forward == FALSE, c("sxmin", "sxmax")] <- d2[d2$forward == FALSE, c("sxmax", "sxmin")]
  prot_all_vs_all_blastp[which(qseqid %in% d$protein_id & sseqid %in% d$protein_id),] %>% 
    left_join(d[,c("protein_id", "y", "xmin", "xmax", "PANEL", "forward")], by = c("qseqid" = "protein_id")) %>%
    left_join(d2, by = c("sseqid" = "sprotein_id")) %>% 
    filter(PANEL.x == PANEL.y & (sy-y)^2 ==1) -> selection
  selection$id = as.character( 1:nrow(selection))
  polydf = data.frame("id" = c("-"), "xs" = c("-"), "ys" = c("-"))
  obs ="1"
  for (obs in selection$id) {
    if (selection$forward.x[selection$id == obs] == TRUE){
      polydf = rbind(polydf, c(obs, (selection$qstart[selection$id == obs]+ selection$xmin[selection$id == obs]),selection$y[selection$id == obs]))
      polydf = rbind(polydf, c(obs, (selection$qend[selection$id == obs ]*3+ selection$xmin[selection$id == obs]),selection$y[selection$id == obs]))}
    else {
      polydf = rbind(polydf, c(obs, (selection$xmin[selection$id == obs] - selection$qstart[selection$id == obs]),selection$y[selection$id == obs]))
      polydf = rbind(polydf, c(obs, (selection$xmin[selection$id == obs] - selection$qend[selection$id == obs ]*3),selection$y[selection$id == obs]))}
    if (selection$forward.y[selection$id == obs] == TRUE){
      polydf = rbind(polydf, c(obs, (selection$send[selection$id == obs ]*3+ selection$sxmin[selection$id == obs]),selection$sy[selection$id == obs]))
      polydf = rbind(polydf, c(obs, (selection$sstart[selection$id == obs ]+ selection$sxmin[selection$id == obs]),selection$sy[selection$id == obs]))}
    else {
      polydf = rbind(polydf, c(obs, (selection$sxmin[selection$id == obs] - selection$send[selection$id == obs]*3),selection$sy[selection$id == obs]))
      polydf = rbind(polydf, c(obs, (selection$sxmin[selection$id == obs] - selection$sstart[selection$id == obs ]),selection$sy[selection$id == obs]))}
  }
  datapoly <- merge(selection, polydf, by = c("id"))
  datapoly$xs = as.integer(datapoly$xs)
  datapoly$ys = as.integer(datapoly$ys)
pdf(file="crass_genome_organisationaltregions.pdf", width=22, height=25)
  
genomesnamed_1%>%
    left_join(datapoly, by =c("protein_id" = "qseqid")) %>%
    ggplot( aes(xmin = start, xmax = end, y = sequence, forward = as.logical(strand), fill =Colour, colour =gene_of_interest, size =opac)) +
    
    scale_size_continuous(range = c(0.1, 1.5), guide = "none")+
    geom_gene_arrow(arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm")) +
    #geom_gene_label(inherit.aes = FALSE,aes(xmin = start, xmax = end, y = sequence,label=orth_group, color = "#000000"))+
    facet_col(~family, scales = "free_y",space = "free", strip.position = "top") +
    scale_fill_discrete(guide = "none") +
    scale_color_brewer( type = "qual", palette=3)+
    scale_alpha(guide = "none", range = c(0.1, 0.5)) +
    scale_y_discrete(position = "right") +
    new_scale_colour()+
    new_scale_fill()+
    scale_fill_gradient(low = "#ffffcc", high= "#b10026")+
    scale_color_gradient(low = "#cccca3", high= "#b00002", guide = "none")+
    geom_polygon(aes(x = xs, y=ys, group=id, alpha = 0.5, fill = pident, colour= pident, size=0.5), inherit.aes = F)+
    geom_text(inherit.aes = FALSE,aes(x = end - ((end-start)/2), y = sequence, label = orth_group,  size=1)) +
    
    theme_bw()
dev.off()
  
#   genomesnamed %>%
#     #mutate(name = paste(CRISPR_hit, order(unique(sequence)))) %>%
#     #mutate(name = paste(CRISPR_hit, Habitat))%>% view()
#     #mutate(name = factor(name, levels = sort(unique(name))[clust_boega$order]))%>%
#     #mutate(Genus = factor(Genus, unique(genomesnamed$Genus[match(genomesnamed$sequence , sort(unique(genomesnamed$sequence))[clust_boega$order])]))) %>%
#     mutate(full_name = paste("Ph:", Genus, "host:", Host)) %>%
#     mutate(sequence = paste(sequence, full_name, sep = " ")) -> genomesnames_1
#   write.csv(genomesnames_1, "genomesnames_1")
#   clusttest = data.frame(labels = c(clust_boega$labels), order=  c(clust_boega$order))
#   write.csv(clusttest, "clust_boega")
#   read.csv("genomesnames_1")
#   read.csv("clust_boega")
#   genomesnames_1 %>%
#     mutate(sequence = factor(genomesnames_1$sequence, levels = clust_boega$labels[clust_boega$order])) -> genomesnames_1
# plot <- ggplot(genomesnames_1, aes(xmin = start, xmax = end, y = sequence, forward = as.logical(1), fill =Colour, colour =gene_of_interest, size =opac)) +
#     scale_size_continuous(range = c(0.1, 0.75))+
#     geom_gene_arrow(arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm")) +
#     geom_text(inherit.aes = FALSE,aes(x = end - ((end-start)/2), y = sequence, label = orth_group,  size=1)) +
#     #geom_gene_label(inherit.aes = FALSE,aes(xmin = start, xmax = end, y = sequence,label=orth_group, color = "#000000"))+
#     facet_col(~family, scales = "free_y",space = "free", strip.position = "top") +
#     scale_fill_discrete(guide = "none") +
#     scale_y_discrete(position = "right") +
#     theme_bw()
# 
#   d <- layer_data(plot)
#   setDT(d)
#   length(unique(d$group))
#   length(genomesnames_1$protein_id)
#   d$protein_id = genomesnamed$protein_id
#   d2 <- layer_data(plot)
#   setDT(d2)
#   d2$protein_id = genomesnamed$protein_id
#   prot_all_vs_all_blastp <- fread(file = "all-vs-all.tsv")
#   colnames(prot_all_vs_all_blastp) <-
#     c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
#   prot_all_vs_all_blastp <- prot_all_vs_all_blastp[evalue < 1e-5] #Filter by e-value
#   prot_all_vs_all_blastp <- prot_all_vs_all_blastp[qseqid != sseqid] #Remove hits on query protein itself
#   prot_all_vs_all_blastp <- prot_all_vs_all_blastp[length >= 100]
#   colnames(d2) = c("ssize","sfill","scolour","sy","sxmin","sxmax","sforward","sPANEL","sgroup","salpha","slinetype","sprotein_id")
#   merged = merge(prot_all_vs_all_blastp, d, by.y="protein_id", by.x= "qseqid")
#   merged2 = merge(merged, d2, by.y="sprotein_id", by.x= "sseqid")
#   sortedmerge = merged2[(merged2$y - merged2$sy)^2 == 1 & merged2$PANEL == merged2$sPANEL]
#   fwb = merge(genomesnames_1, sortedmerge, by.x = "protein_id", by.y = "qseqid", all = TRUE)
#   fwb$id = as.character( 1:nrow(fwb))
#   polydf = data.frame("id" = c("-"), "xs" = c("-"), "ys" = c("-"))
#   for (obs in fwb$id) {
#     polydf = rbind(polydf, c(obs, (fwb$qstart[fwb$id == obs]+ fwb$xmin[fwb$id == obs]),fwb$y[fwb$id == obs]))
#     polydf = rbind(polydf, c(obs, (fwb$qend[fwb$id == obs ]*3+ fwb$xmin[fwb$id == obs]),fwb$y[fwb$id == obs]))
#     polydf = rbind(polydf, c(obs, (fwb$send[fwb$id == obs ]*3+ fwb$sxmin[fwb$id == obs]),fwb$sy[fwb$id == obs]))
#     polydf = rbind(polydf, c(obs, (fwb$sstart[fwb$id == obs ]+ fwb$sxmin[fwb$id == obs]),fwb$sy[fwb$id == obs]))
# 
#   }
#   datapoly <- merge(fwb, polydf, by = c("id"))
#   datapoly$xs = as.integer(datapoly$xs)
#   datapoly$ys = as.integer(datapoly$ys)
#   ggplot(datapoly, aes(xmin = start, xmax = end, y = sequence, forward = as.logical(1), fill =Colour, colour =gene_of_interest, size =opac)) +
#     scale_size_continuous(range = c(0.1, 0.75), guide ="none")+
#     geom_gene_arrow(arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm")) +
#     geom_segment(aes(x = xmax - (xmax-xmin)/2, xend = sxmax - (sxmax-sxmin)/2, y = y, yend = sy, alpha= (opac+0.5)/2, size=0.2, colour =gene_of_interest), inherit.aes = FALSE)+
#     #geom_polygon(aes(x = xs, y=ys, group=id, alpha= (opac+0.5)/2), size =0.5)+
#     geom_text(inherit.aes = FALSE,aes(x = end - ((end-start)/2), y = sequence, label = orth_group,  size=1)) +
#     facet_col(~family, scales = "free_y",space = "free", strip.position = "top") +
#     scale_fill_discrete(guide = "none") +
#     scale_y_discrete(position = "right") +
#     #scale_colour_discrete(range(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')))+
#     scale_color_brewer( type = "qual", palette=3)+
#     scale_alpha(guide = "none", range = c(0.1, 0.5)) +
#     xlab("Base Pair") +
#     theme_bw()
#   #(p1/p2/p3/p4)|p5
#       dev.off()
#       pdf(file="crass_genome_organisationaltregions.pdf", width=22, height=17)
#       ggplot(datapoly, aes(xmin = start, xmax = end, y = sequence, forward = as.logical(strand), fill =Colour, colour =gene_of_interest, size =opac)) +
#         scale_size_continuous(range = c(0.1, 0.75), guide ="none")+
#         geom_gene_arrow(arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm")) +
#         #geom_segment(aes(x = xmax - (xmax-xmin)/2, xend = sxmax - (sxmax-sxmin)/2, y = y, yend = sy, alpha= (opac+0.5)/2, size=0.2, colour =gene_of_interest), inherit.aes = FALSE)+
#         geom_polygon(aes(x = xs, y=ys, group=id, alpha= (opac+0.5)/2), size =0.5)+
#         geom_text(inherit.aes = FALSE,aes(x = end - ((end-start)/2), y = sequence, label = orth_group,  size=1)) +
#         facet_col(~family, scales = "free_y",space = "free", strip.position = "top") +
#         scale_fill_discrete(guide = "none") +
#         scale_y_discrete(position = "right") +
#         #scale_colour_discrete(range(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')))+
#         scale_color_brewer( type = "qual", palette=3)+
#         scale_alpha(guide = "none", range = c(0.1, 0.5)) +
#         xlab("Base Pair") +
#         theme_bw()
#       #(p1/p2/p3/p4)|p5
#       dev.off()
#       datapoly$pident
#       pdf(file="crass_genome_organisationaltregionsidentity.pdf", width=22, height=17)
#       ggplot(datapoly, aes(xmin = start, xmax = end, y = sequence, forward = as.logical(1),fill =Colour, colour =gene_of_interest, size =opac)) +
#         scale_size_continuous(range = c(0.1, 0.75), guide ="none")+
#         geom_gene_arrow(arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm")) +
#         #geom_segment(aes(x = xmax - (xmax-xmin)/2, xend = sxmax - (sxmax-sxmin)/2, y = y, yend = sy, alpha= (opac+0.5)/2, size=0.2, colour =gene_of_interest), inherit.aes = FALSE)+
#         scale_fill_discrete(guide = "none") +
#         new_scale_fill()+
#         scale_fill_gradient(low = "#ffffcc", high= "#b10026")+
#         geom_polygon(aes(x = xs, y=ys, group=id, alpha= (opac+0.5)/2, fill=pident, size=0.1))+
#         geom_text(inherit.aes = FALSE,aes(x = end - ((end-start)/2), y = sequence, label = orth_group,  size=1)) +
#         facet_col(~family, scales = "free_y",space = "free", strip.position = "top") +
#         scale_y_discrete(position = "right") +
#         #scale_colour_discrete(range(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')))+
#         scale_color_brewer( type = "qual", palette=3)+
#         scale_alpha(guide = "none", range = c(0.1, 0.5)) +
#         xlab("Base Pair") +
#         theme_bw()
#       #(p1/p2/p3/p4)|p5
#       dev.off()

# selectedmsa2 %>%
#   mutate(sequence = factor(sequence, levels = c("MT774383","MT774386","NG-14253_11xxx_lib223957_5684_NODE_2_length_92241_cov_49.26","EM181_T3_ms_1","MT774377","MN917146","MT774409","MT774379","MT774376","NG-16262_31_lib261986_NODE_2_length_99704_cov_87.7682","Sib2_ms_1","MT774400","MT774404","MT774407","MT774397","MT774398","MT774399"))) -> selectedmsa2
# plot <- ggplot(selectedmsa2, aes(xmin = start, xmax = end, y = sequence, forward = as.logical(1), fill =Colour, colour =gene_of_interest, size =opac)) +
# scale_size_continuous(range = c(0.1, 0.75))+
# geom_gene_arrow(arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm")) +
# geom_text(inherit.aes = FALSE,aes(x = end - ((end-start)/2), y = sequence, label = orth_group,  size=1)) +
# #geom_gene_label(inherit.aes = FALSE,aes(xmin = start, xmax = end, y = sequence,label=orth_group, color = "#000000"))+
# facet_col(~family, scales = "free_y",space = "free", strip.position = "top") +
# scale_fill_discrete(guide = "none") +
# scale_y_discrete(position = "right") +
# theme_bw()
# d <- layer_data(plot)
# setDT(d)
# 
# d$protein_id = selectedmsa2$protein_id
# d2 <- layer_data(plot)
# setDT(d2)
# d2$protein_id = selectedmsa2$protein_id
# prot_all_vs_all_blastp <- fread(file = "all-vs-all.tsv")
# colnames(prot_all_vs_all_blastp) <-
#   c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
# prot_all_vs_all_blastp <- prot_all_vs_all_blastp[evalue < 1e-5] #Filter by e-value
# prot_all_vs_all_blastp <- prot_all_vs_all_blastp[qseqid != sseqid] #Remove hits on query protein itself
# prot_all_vs_all_blastp <- prot_all_vs_all_blastp[length >= 100]
# colnames(d2) = c("ssize","sfill","scolour","sy","sxmin","sxmax","sforward","sPANEL","sgroup","salpha","slinetype","sprotein_id")
# merged = merge(prot_all_vs_all_blastp, d, by.y="protein_id", by.x= "qseqid")
# merged2 = merge(merged, d2, by.y="sprotein_id", by.x= "sseqid")
# sortedmerge = merged2[(merged2$y - merged2$sy)^2 == 1 & merged2$PANEL == merged2$sPANEL]
# fwb = merge(selectedmsa2, sortedmerge, by.x = "protein_id", by.y = "qseqid", all = TRUE)
# fwb$id = as.character( 1:nrow(fwb))
# polydf = data.frame("id" = c("-"), "xs" = c("-"), "ys" = c("-"))
# for (obs in fwb$id) {
#   polydf = rbind(polydf, c(obs, (fwb$qstart[fwb$id == obs]+ fwb$xmin[fwb$id == obs]),fwb$y[fwb$id == obs]))
#   polydf = rbind(polydf, c(obs, (fwb$qend[fwb$id == obs ]+ fwb$xmin[fwb$id == obs]),fwb$y[fwb$id == obs]))
#   polydf = rbind(polydf, c(obs, (fwb$send[fwb$id == obs ]+ fwb$sxmin[fwb$id == obs]),fwb$sy[fwb$id == obs]))
#   polydf = rbind(polydf, c(obs, (fwb$sstart[fwb$id == obs ]+ fwb$sxmin[fwb$id == obs]),fwb$sy[fwb$id == obs]))
#   
# }
# datapoly <- merge(fwb, polydf, by = c("id"))
# datapoly$xs = as.integer(datapoly$xs)
# datapoly$ys = as.integer(datapoly$ys)
# write.table(datapoly, "datapoly.tsv", sep = "\t")
# pdf(file="selectedmsa2.pdf", width=22, height=17)
# ggplot(datapoly, aes(xmin = start, xmax = end, y = sequence, forward = as.logical(1), fill =Colour, colour =gene_of_interest, size =opac)) +
#   scale_size_continuous(range = c(0.1, 0.75), guide ="none")+
#   geom_gene_arrow(arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm")) +
#   geom_text(inherit.aes = FALSE,aes(x = end - ((end-start)/2), y = sequence, label = orth_group,  size=10)) +
#   #geom_segment(aes(x = xmax - (xmax-xmin)/2, xend = sxmax - (sxmax-sxmin)/2, y = y, yend = sy, alpha= (opac+0.5)/2, size=0.2, colour =gene_of_interest), inherit.aes = FALSE)+
#   geom_polygon(aes(x = xs, y=ys, group=id, alpha= (opac+0.5)/2))+
#   facet_col(~family, scales = "free_y",space = "free", strip.position = "top") +
#   scale_fill_discrete(guide = "none") +
#   scale_y_discrete(position = "right") +
#   #scale_colour_discrete(range(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')))+
#   scale_color_brewer( type = "qual", palette=3)+
#   scale_alpha(guide = "none", range = c(0.1, 0.5)) +
#   xlab("Base Pair") +
#   theme_bw()
dev.off()
########All vs all blastp########
prot_all_vs_all_blastp <- fread(file = "all-vs-all.tsv") #Change path!
colnames(prot_all_vs_all_blastp) <-
  c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
prot_all_vs_all_blastp <- prot_all_vs_all_blastp[evalue < 1e-5] #Filter by e-value
prot_all_vs_all_blastp <- prot_all_vs_all_blastp[qseqid != sseqid] #Remove hits on query protein itself
prot_all_vs_all_blastp <- prot_all_vs_all_blastp[length >= 100] #Filter by length of alignment, optionally
prot_all_vs_all_blastp <- prot_all_vs_all_blastp[pident >= 35] #Filter by percent identity of alignment, optionally
#find the real length of the proteins
lengthframe = read.delim("lengthlist.txt", header = TRUE)
#add the size of the query and subject
prot_all_vs_all_blastp <-merge(prot_all_vs_all_blastp, lengthframe, by.x = "qseqid", by.y= "protein_id")
prot_all_vs_all_blastp <-merge(prot_all_vs_all_blastp, lengthframe, by.x = "sseqid", by.y= "protein_id")
#choose the smallest of query and subject size then divide the size of the hit by the length of the chosen protein
prot_all_vs_all_blastp$coverage <- ifelse(prot_all_vs_all_blastp$size.x<prot_all_vs_all_blastp$size.y,(prot_all_vs_all_blastp$qend- prot_all_vs_all_blastp$qstart +1) /prot_all_vs_all_blastp$size.x,(prot_all_vs_all_blastp$send- prot_all_vs_all_blastp$sstart +1) /prot_all_vs_all_blastp$size.y)
#prot_all_vs_all_blastp <- prot_all_vs_all_blastp[coverage >= 0.49] #Filter by percent identity of alignment, optionally

#Select for hits relevant to our proteins of interest
prot_of_interest_blastp <- prot_all_vs_all_blastp[qseqid %in% unique(as.character(unlist(homologs))) | sseqid %in% unique(as.character(unlist(homologs)))]

########Convert BLASTp output into a graph########
prot_of_interest_graph <- graph_from_data_frame(prot_of_interest_blastp[,.(qseqid,sseqid)], directed = FALSE) #Non-directed graph from a dataframe
#Add vertex attributes
V(prot_of_interest_graph)$cluster <-
  orthologs[names(V(prot_of_interest_graph)),"L1"] #Protein clusters
V(prot_of_interest_graph)$subfamily <-
  as.character(contig_clustering_all[orthologs[names(V(prot_of_interest_graph)),"genome"],"nomen_subfamily"]) #Phage subfamilies
V(prot_of_interest_graph)$homolog <- NA #Homologs of proteins of interest
for(i in names(homologs)) {
  V(prot_of_interest_graph)$homolog[names(V(prot_of_interest_graph)) %in% homologs[[i]]] <- i
}
prot_of_interest_df <- igraph::as_data_frame(prot_of_interest_graph, what = "vertices")

#Look for cases when BLAST hits cross from one protein cluster to a different one ("mismatched hits")
prot_of_interest_blastp$mismatched <- FALSE
prot_of_interest_blastp[which(
  prot_of_interest_df[prot_of_interest_blastp$qseqid,"cluster"] != prot_of_interest_df[prot_of_interest_blastp$sseqid,"cluster"]
  ),"mismatched"] <- TRUE
#Add edge attributes
length(prot_of_interest_blastp$mismatched[prot_of_interest_blastp$mismatched])
E(prot_of_interest_graph)$mismatched <- prot_of_interest_blastp$mismatched #"mismatched hits"
E(prot_of_interest_graph)$length <- prot_of_interest_blastp$length #BLAST alignment length
E(prot_of_interest_graph)$pident <- prot_of_interest_blastp$pident #% identity
E(prot_of_interest_graph)$weight <-
  log10(E(prot_of_interest_graph)$length*(E(prot_of_interest_graph)$pident/100)) #Assign edge weight as log-transformed product of length and identity

prot_of_interest_graph <-
  igraph::simplify(prot_of_interest_graph, edge.attr.comb = "sum") #Simplify graph, edge weights and other attributes are summed up
E(prot_of_interest_graph)$mismatched <- as.logical(E(prot_of_interest_graph)$mismatched)
prot_of_interest_layout <- layout_with_kk(prot_of_interest_graph)

#Cluster graph into connected subgraph - "clusters" again, very confusing... Let's call them subgraphs instead
prot_of_interest_clust <- clusters(prot_of_interest_graph)
V(prot_of_interest_graph)$subgraphs <- prot_of_interest_clust$membership[V(prot_of_interest_graph)$name]
V(prot_of_interest_graph)
########Plot graphs########
colors = c(rainbow(85), rep("black", 8))

#Graph colored by protein cluster. Red edges are BLASTp hits crossing between clusters
pdf(file="crass_proteins_graph1.pdf", width=11, height=8.5)
plot(prot_of_interest_graph,
     layout = prot_of_interest_layout,
     vertex.size = 2,
     vertex.color = colors[as.numeric(as.factor(V(prot_of_interest_graph)$cluster))],
     edge.color = c("grey", "red")[as.numeric(E(prot_of_interest_graph)$mismatched)+1],
      vertex.label=NA)
legend(1,1.2, levels(as.factor(V(prot_of_interest_graph)$cluster)),
       fill = colors, cex = 0.4, bty = "n")
dev.off()

#Graph colored by homologs of proteins of interest. Red edges are BLASTp hits crossing between clusters
pdf(file="crass_proteins_graph2.pdf", width=11, height=8.5)
plot(prot_of_interest_graph,
     layout = prot_of_interest_layout,
     vertex.size = 2,
     vertex.color = c(brewer.pal(12, "Paired"),"black")[as.numeric(as.factor(V(prot_of_interest_graph)$homolog))],
     edge.color = c("grey", "red")[as.numeric(E(prot_of_interest_graph)$mismatched)+1],
     vertex.label=NA)
legend(1,-0.3, levels(as.factor(V(prot_of_interest_graph)$homolog)),
       fill = c(brewer.pal(12, "Paired"),"black"), cex = 0.7, bty = "n")
dev.off()

#Graph colored by phage subfamilies of proteins of interest. Red edges are BLASTp hits crossing between clusters
pdf(file="crass_proteins_graph3.pdf", width=11, height=8.5)
plot(prot_of_interest_graph,
     layout = prot_of_interest_layout,
     vertex.size = 2,
     vertex.color = c(brewer.pal(12,"Paired"),"grey","black","darkorchid1","darkmagenta")[as.numeric(
       factor(V(prot_of_interest_graph)$subfamily, levels = levels(contig_clustering$nomen_subfamily))
       )],
     edge.color = c("grey", "red")[as.numeric(E(prot_of_interest_graph)$mismatched)+1],
     vertex.label=NA)
legend(1,-0.3, levels(contig_clustering$nomen_subfamily),
       fill = c(brewer.pal(12,"Paired"),"grey","black","darkorchid1","darkmagenta"), cex = 0.7, bty = "n")
dev.off()

#Graph colored by connected subgraphs.
pdf(file="crass_proteins_graph4.pdf", width=11, height=8.5)
plot(prot_of_interest_graph,
     layout = prot_of_interest_layout,
     vertex.size = 2,
     vertex.color = brewer.pal(12,"Paired")[V(prot_of_interest_graph)$subgraphs],
     edge.color = c("grey", "red")[as.numeric(E(prot_of_interest_graph)$mismatched)+1],
     vertex.label=NA)
legend(1,-0.3, unique(V(prot_of_interest_graph)$subgraphs),
       fill = brewer.pal(12,"Paired"), cex = 0.7, bty = "n")
dev.off()


graphdata = igraph::as_data_frame(prot_of_interest_graph, "vertices")
write.csv(graphdata, "graphdata4.csv")
########Plot a subgraph########
#A zoomed-in view into subgraph #4 - this one contains gp44 and its homologs
prot_of_interest_subgraph <- subgraph(prot_of_interest_graph, which(V(prot_of_interest_graph)$subgraphs == 4))
prot_of_interest_layout2 <- layout_with_fr(prot_of_interest_subgraph)

pdf(file="crass_proteins_graph5.pdf", width=11, height=8.5)
plot(prot_of_interest_subgraph,
     layout = prot_of_interest_layout2,
     vertex.size = 2,
     vertex.color = c(brewer.pal(12, "Paired"),"black")[as.numeric(as.factor(V(prot_of_interest_subgraph)$homolog))],
     edge.color = c("grey", "red")[as.numeric(E(prot_of_interest_subgraph)$mismatched)+1],
     vertex.label=NA)
legend(1,-0.3, levels(as.factor(V(prot_of_interest_subgraph)$homolog)),
       fill = c(brewer.pal(12, "Paired"),"black"), cex = 0.7, bty = "n")
dev.off()

pdf(file="crass_proteins_graph6.pdf", width=11, height=8.5)
plot(prot_of_interest_subgraph,
     layout = prot_of_interest_layout2,
     vertex.size = 2,
     vertex.color = c(brewer.pal(12, "Paired"),"black")[as.numeric(as.factor(V(prot_of_interest_subgraph)$cluster))],
     edge.color = c("grey", "red")[as.numeric(E(prot_of_interest_subgraph)$mismatched)+1],
     vertex.label=NA)
legend(1,-0.3, levels(as.factor(V(prot_of_interest_subgraph)$cluster)),
       fill = c(brewer.pal(12, "Paired"),"black"), cex = 0.7, bty = "n")
dev.off()

save.image()
