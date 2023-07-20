#Code written by Nicolás Mongiardino Koch to support publication of manuscript
#'Embracing the taxonomic and topological stability of phylogenomics' (2023)
#Sci. Rep.

#modify the working directory to the location of this file. All directory
#changes are made relative to this folder, and the structure of the data is
#assumed to be as in the GitHub repository
setwd('')

rm(list=ls())
library(phangorn)
library(tidyverse)
library(combinat)
library(ggplot2)
library(MetBrewer)

#STEP 1: combine sequences into matrices and align------------------------------
setwd('./A_sequences_NCBI')

#list all folders
folders = list.dirs(recursive = F)

#loop through each folder
for(i in 1:length(folders)) {
  #get all file names (== species names) within the folder
  full_files = list.files(folders[i], full.names = T)
  names = list.files(folders[i])
  names = paste(sapply(strsplit(names, '_'), '[', 1), 
                sapply(strsplit(names, '_'), '[', 2), sep = '_')
  
  #loop through the files and combine their content generating a single fasta
  #file
  for(j in 1:length(full_files)) {
    this_file = readLines(full_files[j])
    this_file[1] = paste0('>', names[j])
    
    if(j == 1) {
      all_files = this_file
    } else {
      all_files = c(all_files, this_file)
    }
  }
  
  #create folder where to place
  if(i == 1) dir.create('../B_matrices')
  
  #write the fasta file
  writeLines(all_files, paste0(gsub('./', '../B_matrices/', folders[i]), '_all.fa'))
  
  #...and another one without S. purpuratus, the only sequence added relative to
  #Lee et al. (2023)
  if('Strongylocentrotus_purpuratus' %in% names) {
    start = grep('Strongylocentrotus_purpuratus', all_files)
    if(tail(start == grep('>', all_files), n = 1)) {
      end = length(all_files)
    } else {
      end = grep('>', all_files)[which(grep('>', all_files) == start)+1]-1
    }
    
    all_files = all_files[-(start:end)]
    writeLines(all_files, paste0(gsub('./', '../B_matrices/', folders[i]), '_no_spur.fa'))
  }
}

#print commands used to align loci with MAFFT:
commands = paste('linsi', list.files('../B_matrices'), '>', 
                 gsub('.fa', '_linsi.fa', list.files('../B_matrices')))
commands

#STEP2: convert to phylip, remove duplicated taxon and trim---------------------
rm(list=ls())
setwd('../B_matrices')

#list all aligned fasta files
files = list.files(pattern = 'linsi.fa', full.names = T, recursive = T)

#loop through the files
for(i in 1:length(files)) {
  #load alignment
  aligned_gene = read.phyDat(files[i], format = 'fasta', type = 'DNA')
  gene = as.character(aligned_gene)
  
  #COX1 for Scaphechinus mirabilis came from two non-overlapping sequences, the
  #following code combines them into a single line of the alignment before
  #inference
  if(grepl('cox1', files[i])) {
    Scaph1 = which(rownames(gene) == 'Scaphechinus_mirabilis')
    Scaph2 = which(rownames(gene) == 'Scaphechinus_mirabilis2')
    
    gene[Scaph1, which(gene[Scaph2,] != '-')] = gene[Scaph2, which(gene[Scaph2,] != '-')]
    
    gene = gene[-Scaph2,]
  }
  
  #save as phylip format
  write.phyDat(phyDat(gene, type = 'DNA'), 
               file = gsub('.fa', '_untrimmed.phy', files[i]), 
               format = 'phylip')
  
  #additionally, delete positions with > 50% gaps
  deleted_pos = which(apply(gene, 2, function(x) mean(x == '-')) > 0.5)
  if(length(deleted_pos) > 0) {
    gene = gene[,-deleted_pos]
  }
  
  write.phyDat(phyDat(gene, type = 'DNA'), 
               file = gsub('.fa', '_trimmed.phy', files[i]), 
               format = 'phylip')
}

#STEP 3: concatenate------------------------------------------------------------
rm(list=ls())
#list all phylip alignments
all_files = list.files(pattern = '.phy', full.names = T, recursive = T)

#these correspond to 4 types with/without S. purpuratus + with/without trimming
types = c('all_linsi_untrimmed.phy', 'all_linsi_trimmed.phy', 
          'no_spur_linsi_untrimmed.phy', 'no_spur_linsi_trimmed.phy')

#loop through each file type and concatenate
for(n in 1:length(types)) {
  #list corresponding files
  files = all_files[grep(types[n], all_files)]
  
  #find all taxon names and the length of the loci
  all_taxa = c()
  lengths = c()
  for(i in 1:length(files)) {
    gene = read.phyDat(files[i], format = 'phylip', type = 'DNA')
    gene = as.character(gene)
    all_taxa = c(all_taxa, rownames(gene)[which(!rownames(gene) %in% all_taxa)])
    lengths = c(lengths, dim(gene)[2])
  }
  
  #generate empty alignment with taxon names ordered alphabetically
  all_taxa = all_taxa[order(all_taxa)]
  dataset = matrix(NA, ncol = sum(lengths), nrow = length(all_taxa))
  row.names(dataset) = all_taxa
  
  #loop through each aligned loci and incorporate
  for(i in 1:length(files)) {
    gene = read.phyDat(files[i], format = 'phylip')
    gene = as.character(gene)
    gene = gene[order(rownames(gene)),]
    
    #if a taxon is missing, add '-'
    if(nrow(gene) < nrow(dataset)) {
      missing = rownames(dataset)[which(!rownames(dataset) %in% rownames(gene))]
      for(j in 1:length(missing)) {
        gene = rbind(gene, rep('-', ncol(gene)))
        rownames(gene)[length(rownames(gene))] = missing[j]
      }
      gene = gene[order(rownames(gene)),]
    }
    
    if(i == 1) {
      dataset[,1:lengths[1]] = gene
    } else {
      dataset[,(sum(lengths[1:(i-1)])+1):sum(lengths[1:i])] = gene
    }
  }
  
  #save the supermatrix
  if(n == 1) dir.create('../C_supermatrices')
  write.phyDat(phyDat(dataset, type = 'DNA'), 
               file = paste0('../C_supermatrices/', types[n]), 
               format = 'phylip')
  
  #as well as the associated partition file
  partition = matrix(NA, ncol = 2, nrow = length(files))
  for(i in 1:nrow(partition)) {
    if(i == 1) {
      partition[i,1] = 1
      partition[i,2] = lengths[1]
    } else {
      partition[i,1] = partition[(i-1),2] + 1
      partition[i,2] = partition[i,1] + lengths[i] - 1
    }
  }
  
  partition = apply(partition, 1, paste, collapse = '-')
  names = gsub('.phy', '', sapply(strsplit(files, '/'), '[', 2))
  partitions_tosave = paste0('DNA, ', names, ' = ', partition)
  write(partitions_tosave, 
        file = gsub('.phy', '.txt', paste0('../C_supermatrices/', types[n])))
}

#print commands used to infer trees with iqtree2:
commands = paste('iqtree -s', list.files('../C_supermatrices', pattern = '.phy'), 
                 '-st DNA -spp', list.files('../C_supermatrices', pattern = '.txt'), 
                 '-m MFP+MERGE -bb 1000')
commands

#STEP 4: topology testing-------------------------------------------------------
rm(list=ls())
dir.create('./D_topology_test')
setwd('./D_topology_test')

#copy trimmed alignment with S. purpuratus
file.copy('../C_supermatrices/all_linsi_trimmed.phy', 
          './all_linsi_trimmed.phy')
file.copy('../C_supermatrices/all_linsi_trimmed.txt', 
          './all_linsi_trimmed.txt')

#load alignment and patition fil
data = read.phyDat('all_linsi_trimmed.phy', format = 'phylip', type = 'DNA')
data = as.DNAbin(data)
partitions = read.table('all_linsi_trimmed.txt', sep = ' ')
names = as.character(unlist(enframe(partitions[,(which(partitions[1,] == '=') - 1)], 
                                    name = NULL), use.names = F))
partitions = as_tibble(partitions[,4])
partitions = partitions %>% separate(value, into = c('Start', 'End'), sep = '-') %>%
  mutate_if(is.character, as.numeric)

#extract individual alignments from the supermatrix
for(i in 1:nrow(partitions)) {
  gene = data[,partitions$Start[i]:partitions$End[i]]
  j = 1
  while(j <= nrow(gene)) {
    if(all(as.character(gene[j,]) == '-')) {
      gene = gene[-j,]
    } else {
      j = j + 1
    }
  }
  write.phyDat(phyDat(gene, type = 'DNA'), 
               file = paste(getwd(), '/', gsub('_all', '', names[i]), '.phy', sep = ''), 
               format = 'phylip')
}

#generate constraint topologies
out = '((Colobocentrotus_mertensii,Strongylocentrotus_purpuratus),((Linopneustes_longispinus,Maretia_planulata),('
cassi = '(Rhyncholampas_pacificus,Conolampas_sigsbei,Echinolampas_crassa)'
lagani = '(Echinocyamus_pusillus,Laganum_fudsiyama,Peronella_lesueri,Peronella_japonica)'
scutelli = '(Astriclypeus_mannii,Echinodiscus_bisperforatus,Sculpsitechinus_auritus,Sinaechinocyamus_mai,Dendraster_excentricus,Scaphechinus_mirabilis,Echinarachnius_parma,Mellita_isometra,Mellita_notabilis,Lanthonia_longifissa,Mellitella_stokesii,Leodia_sexiesperforata,Encope_aberrans,Encope_grandis)'
clype = '(Arachnoides_placenta,Clypeaster_reticulatus,Clypeaster_japonicus,Clypeaster_virescens)'
end = ')));'

#list all possible permutations
permn = permn(c('cassi', 'lagani', 'scutelli', 'clype'))
permn_table = matrix(nrow = length(permn), ncol = length(permn[[1]]))
for(i in 1:length(permn)) permn_table[i,] = permn[[i]]

#remove duplicates product of node rotation
#(an error will be reported, this is the expected behavior)
for(i in 1:nrow(permn_table)) {
  if(sum(permn_table[,1] == permn_table[i,1] & permn_table[,2] == permn_table[i,2]) == 2) {
    del = which(permn_table[,1] == permn_table[i,1] & permn_table[,2] == permn_table[i,2])[2]
    permn_table = permn_table[-del,]
  }
}

permn_table = permn_table[order(permn_table[,1]),]

#add three extra slots for non comb-like topologies
permn_table = rbind(permn_table, permn_table[1:3,])

#transform these into topologies
all_trees = c()
for(i in 1:nrow(permn_table)) {
  if(i <= 12) {
    all_trees[i] = paste0(out, get(permn_table[i,1]), ',(', get(permn_table[i,2]), 
                          ',(', get(permn_table[i,3]), ',', get(permn_table[i,4]), 
                          '))', end)
  } else {
    all_trees[i] = paste0(out, '(', get(permn_table[i,1]), ',', get(permn_table[i,2]), 
                          '),(', get(permn_table[i,3]), ',', get(permn_table[i,4]), 
                          ')', end)
  }
  
}

#save to files
writeLines(all_trees, 'all_trees_supermatrix.tre')
all_trees = unroot(read.tree('all_trees_supermatrix.tre'))
write.tree(all_trees, 'all_trees_supermatrix.tre')

#prune to taxa present in each locus
alignments = list.files(pattern = '.phy')
alignments = alignments[-grep('all', alignments)]
all_trees = read.tree('all_trees_supermatrix.tre')

#now prune to the taxa represented in each individual loci
for(i in 1:length(alignments)) {
  this_dataset = read.phyDat(alignments[i], format = 'phylip', type = 'DNA')
  this_dataset = as.DNAbin(this_dataset)
  
  remove = all_trees[[1]]$tip.label[!all_trees[[1]]$tip.label %in% rownames(this_dataset)]
  
  this_trees = drop.tip.multiPhylo(all_trees, tip = remove)
  this_name = paste0(unlist(strsplit(alignments[i], '_'))[1], '_trees.tre')
  
  #remove duplicate trees in case occupancy of one of the clades was zero
  this_trees = unique.multiPhylo(this_trees)
  this_trees = unroot(this_trees)
  write.tree(this_trees, this_name)
  
  #also split into individual files
  this_trees_text = readLines(this_name)
  for(i in 1:length(this_trees_text)) {
    writeLines(this_trees_text[i], paste0(LETTERS[i], '_', this_name))
  }
}

#also split the tree file for the full dataset into individual topologies
all_trees_text = readLines('all_trees_supermatrix.tre')
for(i in 1:length(all_trees_text)) {
  writeLines(all_trees_text[i], paste0(LETTERS[i], '_all_trees_supermatrix.tre'))
}

#generate commands for constrained tree search
alignments = sort(list.files(pattern = 'phy'))
tree_files = sort(list.files(pattern = 'tre'))
tree_files = tree_files[grep(paste0(paste0(LETTERS[1:length(all_trees)], '_'), 
                                    collapse = '|'), tree_files)]

#repeat alignment names to match tree names
alignments = alignments[sapply(sapply(strsplit(tree_files, '_'), '[', 2), 
                               function(x) grep(x, alignments))]


commands = paste('iqtree -s', alignments, '-st DNA -m MFP -g', tree_files, 
                 '-pre', gsub('_trees_supermatrix.tre', '', 
                              gsub('_trees.tre', '', tree_files)))
commands[grep('all', commands)] = gsub('-m MFP', 
                                       '-spp all_linsi_trimmed.txt -m MFP+MERGE', 
                                       commands[grep('all', commands)])

sort(commands)

#after inference, get all optimized trees into single files and perform topology
#testing
alignments = sort(list.files(pattern = 'phy'))
alignments = sapply(strsplit(alignments, '_'), '[', 1)
inferred_trees = list.files(pattern = 'treefile')

for(i in 1:length(alignments)) {
  this_inferred_trees = inferred_trees[grep(alignments[i], inferred_trees)]
  
  for(j in 1:length(this_inferred_trees)) {
    this_tree = read.tree(this_inferred_trees[j])
    
    if(j == 1) {
      all_trees = this_tree
    } else {
      all_trees = c(all_trees, this_tree)
    }
  }
  
  write.tree(all_trees, paste0('inferred_', alignments[i], '_trees.tre'))
  
}

alignments = sort(list.files(pattern = 'phy'))
commands = paste('iqtree -s', alignments, '-m MFP -z', 
                 list.files(pattern = 'inferred'), '-n 0 -zb 10000 -zw -au -pre', 
                 paste0(sapply(strsplit(alignments, '_'), '[', 1), '_test'))
commands[grep('all', commands)] = gsub('-m MFP', '-spp all_linsi_trimmed.txt -m MFP+MERGE', 
                                       commands[grep('all', commands)])

commands

#STEP 5: Estimate logL----------------------------------------------------------
rm(list=ls())
setwd('../E_transcriptomes')

#the directory contains:
# A) A phylogenomic dataset in fasta format ('echinoids3_TS.fa')
# B) A partition file ('echinoids3.txt')
# C) A best-fit partitioned model in nexus format ('echinoids3.txt.best_scheme.nex')
# D) Two constrained topologies in newick format, replicating those of:
#   i. Mongiardino Koch et al. (2022) - ('echinoids3_unconstrained.tre')
#   ii. Lee et al. (2023) - ('echinoids3_constrained.tre')

#first, optimize the branch lengths of these trees, then obtain site-wise
#likelihood under each, using the following commands:

#iqtree -nt AUTO -s echinoids3_TS.fa -st AA -spp echinoids3.txt.best_scheme.nex -te echinoids3_unconstrained.tre -pre unconstrained
#iqtree -nt AUTO -s echinoids3_TS.fa -st AA -spp echinoids3.txt.best_scheme.nex -te echinoids3_constrained.tre -pre constrained

#iqtree -nt AUTO -s echinoids3_TS.fa -st AA -spp unconstrained.best_model.nex -z echinoids3_unconstrained.tre -wsl -n 0 -pre likelihood_unconstrained
#iqtree -nt AUTO -s echinoids3_TS.fa -st AA -spp unconstrained.best_model.nex -z echinoids3_constrained.tre -wsl -n 0 -pre likelihood_constrained
  
#then, analyze results with the code below

#load data, partitions and site-wise logL scores
data = read.phyDat('echinoids3_TS.fa', format = 'fasta', type = 'AA')
data = as.AAbin(data)

partitions = read.table('echinoids3.txt', sep = ' ')
names = as.character(unlist(enframe(partitions[,(which(partitions[1,] == '=') - 1)], 
                                    name = NULL), use.names = F))
partitions = as_tibble(partitions[,4])
partitions = partitions %>% separate(value, into = c('Start', 'End'), sep = '-') %>%
  mutate_if(is.character, as.numeric)

siteL_unconstrained = readLines('likelihood_unconstrained.sitelh')[2]
siteL_unconstrained = unlist(strsplit(siteL_unconstrained, ' '))
siteL_unconstrained = siteL_unconstrained[-which(siteL_unconstrained == '')]
siteL_unconstrained = siteL_unconstrained[-1]
siteL_unconstrained = as.numeric(siteL_unconstrained)

siteL_constrained = readLines('likelihood_constrained.sitelh')[2]
siteL_constrained = unlist(strsplit(siteL_constrained, ' '))
siteL_constrained = siteL_constrained[-which(siteL_constrained == '')]
siteL_constrained = siteL_constrained[-1]
siteL_constrained = as.numeric(siteL_constrained)

#sum the up into gene-wise logL scores
geneL_unconstrained = geneL_constrained = vector(length = nrow(partitions))
for(i in 1:nrow(partitions)) {
  geneL_unconstrained[i] = sum(siteL_unconstrained[partitions$Start[i]:partitions$End[i]])
  geneL_constrained[i] = sum(siteL_constrained[partitions$Start[i]:partitions$End[i]])
}

#compute their difference, i.e., deltaGLS
liks = data.frame(geneL_unconstrained, geneL_constrained)
liks = liks %>% mutate(deltaL = geneL_unconstrained - geneL_constrained, 
                       who_wins = NA, pos = factor(1:nrow(liks)))

#classify loci into those supporting one of the trees (deltaGLS > 2 or < -2) or
#ambiguous (-2 < deltaGLD < 2)
for(i in 1:nrow(liks)) {
  if(liks$deltaL[i] > 2) {
    liks$who_wins[i] = 'A'
  } else {
    if(liks$deltaL[i] < -2) {
      liks$who_wins[i] = 'B'
    } else {
      liks$who_wins[i] = 'C'
    }
  }
}

#plot Fig. 2
ggplot(liks, aes(x = pos, y = deltaL, fill = who_wins)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values=met.brewer("Demuth", n=5)[c(2,3,4)]) + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank()) + theme(legend.position = "none")

#plot pie chart of Fig. 2
data.frame(type = names(table(liks$who_wins)), values = as.numeric(table(liks$who_wins))) %>% 
  ggplot(aes(x = '', y = values, fill = type)) + geom_bar(stat = 'identity', width = 1) + 
  coord_polar('y', start = 0) + scale_fill_manual(values=met.brewer("Demuth", n=5)[c(2,3,4)]) + 
  theme_bw() + theme(legend.position = "none")