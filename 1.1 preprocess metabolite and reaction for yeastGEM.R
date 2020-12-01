# This code is used to prepare the data for the map
# source('model change.R')
# source('transition for cellDesigner.R')
library(FALCONET)
library(stringr)
library(tidyverse)
library(readxl)
library(igraph)
library(networkD3)
library(hongR)

# prepare the reaction format
rxn <- read_excel("data/yeastGEM_latest version1.xls", sheet = "Reaction List")
metabolite <-  read_excel("data/yeastGEM_latest version1.xls", sheet = "Metabolite List")
# Update the metabolite name in rxn sheet
rxn_split <- splitAndCombine0(rxn$Reaction,rxn$Abbreviation,sep=" ")
#using the detailed name to take place of the short name
rxn_split$v3 <- getSingleReactionFormula(metabolite$Description,metabolite$Abbreviation,rxn_split$v1)
for (i in 1:length(rxn_split$v2)){
  if(rxn_split$v3[i]=="NA"){
    rxn_split$v3[i] <- rxn_split$v1[i]
  } else{
    rxn_split$v3[i] <- rxn_split$v3[i]
  }
  
}
rxn_split$v3 <- str_replace_all(rxn_split$v3, " \\[.*?\\]", "")
rxn_split$v3 <- str_trim(rxn_split$v3, side = "both")
rxn_split$compartment <- str_extract(rxn_split$v1, "\\[.*?\\]")
for (i in 1:nrow(rxn_split)){
  if(!is.na(rxn_split$compartment[i])){
    rxn_split$v3[i] <- paste(rxn_split$v3[i],rxn_split$compartment[i], sep = "")
  } else{
    rxn_split$v3[i] <- rxn_split$v3[i]
  }
}
rxn$Description_new <- getMultipleReactionFormula(rxn_split$v3,rxn_split$v2,rxn$Abbreviation)
rxn$Description_new <- str_replace_all(rxn$Description_new,";"," ")
# Update the subsytem in rxn sheet to make sure each reaction belong to one subsystem
subsystem_V3 <- read_excel("data/subsystem for yeast8 map_V3.xlsx")
rxn$Subsystem_new <- getSingleReactionFormula(subsystem_V3$subsystem_map_v3,subsystem_V3$rxnID,rxn$Abbreviation)

# filter rxn based on the fluxes
# rxn <- filter(rxn, Flux != 0) # we can remove the reaction with zero fluxes
# which(rxn$Abbreviation =='r_2034')


# Update the metabolite formula
metabolite$Description <- str_replace_all(metabolite$Description, " \\[.*?\\]", "")
metabolite$Description <- str_trim(metabolite$Description, side = "both")
metabolite$compartment <- str_extract(metabolite$Abbreviation, "\\[.*?\\]")
for (i in 1:nrow(metabolite)){
  if(!is.na(metabolite$compartment[i])){
    metabolite$Description[i] <- paste(metabolite$Description[i],metabolite$compartment[i], sep = "")
  } else{
    metabolite$Description[i] <- metabolite$Description[i]
  }
}




# prepare the rxn format for the cellDesigner
colnames(metabolite) <-  c('Metabolite name','Metabolite description','Metabolite formula','Charge','Compartment','KEGGID','CHEBI')
rxn_split_refine <- splitRxnToMetabolite.Yeast(rxn, metabolite)

# analysis subsystem
analysis_subsystem <- rxn %>%
  count(Subsystem_new) %>%
  arrange(., desc(n)) 

# choose the subsytem
# subsystem1 <- c("glycolysis / gluconeogenesis \\( sce00010 \\)")
subsystem1 <- c("valine, leucine and isoleucine metabolism")
# Define the currency metabolite in each subsystem
currency_metabolites <- DefineCurrencyMet(rxn_split_refine, 
                                          subsystem0=subsystem1,
                                          numberGEM=14,
                                          numberSubsystem=1)

# remove the reactions with only one metabolite
# if we do not remove the currency metabolite in the model then this step is mainly removed exchange reaction
rxn_split_refine <- removeRxnWithSingleMet(rxn_split=rxn_split_refine)


#--------------------------------------------------------------------------------------------
## define the base reactant and product for cellDesigner
#---------------------------------------------------------------------------------------------
rxn_split_refine <- addBaseTypeIntoRxn(rxn_split_refine, metabolite, currency_metabolites)


#---------------------------------------------------
# choose reaction based on the subsystem
#-----------------------------------------------------
rxn_core_carbon <- chooseRxnFromSubsystem_new(rxn_split_refine_inf = rxn_split_refine, subsystem0 = subsystem1)
rxnID_choose <- unique(rxn_core_carbon$v2)


#------------------------------------------------------------------
# produce the met, rxn and gpr used for the map production
#------------------------------------------------------------------
# prepare the metabolites formula
# this funcion is used to prepare the metabolite annotation for cell designer
met_annotation <- prepareMET(rxn_core_carbon, currency_metabolites,rxnID_choose)
# prepare the rxn formula
rxn_core_carbon_cellD0 <- prepareRXN(rxn_core_carbon,met_annotation,currency_metabolites)
# prepare the protein and gene
gpr <- prepareGPR(met_annotation)

#save the exampel data format for cell designer
#write.table(met_annotation,"result/met_annotation for example.txt", row.names = FALSE, sep = "\n")
#write.table(rxn_core_carbon_cellD0,"result/rxn_core_carbon_cellD0 for example.txt", row.names = FALSE, sep = "\n")
#write.table(gpr,"result/gpr for example.txt", row.names = FALSE, sep = "\n")

#------------------------------------------------------------------
# produce the file as the input for the cellDesigner
#------------------------------------------------------------------
produceInputForCellDesigner(met_annotation, 
                            gpr,
                            rxn_core_carbon_cellD0,
                            x_size=1200, 
                            y_size=2000,
                            subsystem0=subsystem1)




# function need to be updated in next time
produceInputForCellDesigner <- function (met_annotation_inf, gpr_inf, rxn_core_carbon_inf, x_size = 2500, 
          y_size = 3600, subsystem0=subsystem1) 
{
  txt1 <- readLines(file("data/template1"))
  txt3 <- readLines(file("data/template3"))
  txt5 <- readLines(file("data/template5"))
  txt7 <- readLines(file("data/template7"))
  txt1 <- getNewSize(x = 2500, y = 3600)
  txt2 <- vector()
  txt3 <- vector()
  txt4 <- vector()
  id <- met_annotation_inf$id
  species <- met_annotation_inf$species
  name <- met_annotation_inf$name
  x <- as.character(met_annotation_inf$x)
  y <- as.character(met_annotation_inf$y)
  for (i in seq(length(id))) {
    txt2 <- c(getPosition(id01 = id[i], species01 = species[i], 
                          x01 = x[i], y01 = y[i]), txt2)
  }
  protein0 <- filter(met_annotation_inf, type == "PROTEIN")
  protein0$type <- "GENERIC"
  protein0$protein_id <- str_replace_all(protein0$rxnID, "r_", 
                                         "pr")
  gene0 <- filter(met_annotation_inf, type == "GENE")
  gene0$gene_id <- str_replace_all(gene0$rxnID, "r_", "gn")
  getProteinGene <- function(id, type, name) {
    txt0 <- vector()
    if (type == "GENERIC") {
      txt0[1] <- paste("<celldesigner:protein id=\"", id, 
                       "\" name=\"", name, "\" type=\"GENERIC\"/>", 
                       sep = "")
    }
    if (type == "GENE") {
      txt0[1] <- paste("<celldesigner:gene id=\"", id, 
                       "\" name=\"", name, "\" type=\"GENE\"/>", sep = "")
    }
    return(txt0)
  }
  txt3_1 <- readLines(file("data/txt3_1.xml"))
  protein_id <- protein0$protein_id
  type_p <- protein0$type
  name_p <- protein0$name
  txt3_2 <- vector()
  for (i in seq(length(protein_id))) {
    txt3_2 <- c(txt3_2, getProteinGene(id = protein_id[i], 
                                       type = type_p[i], name = name_p[i]))
  }
  txt3_3 <- readLines(file("data/txt3_3.xml"))
  gene_id <- gene0$gene_id
  type_g <- gene0$type
  name_g <- gene0$name
  txt3_4 <- vector()
  for (i in seq(length(gene_id))) {
    txt3_4 <- c(txt3_4, getProteinGene(id = gene_id[i], type = type_g[i], 
                                       name = name_g[i]))
  }
  txt3_5 <- readLines(file("data/txt3_5.xml"))
  txt30 <- c(txt3_1, txt3_2, txt3_3, txt3_4, txt3_5)
  met_without_protein0 <- filter(met_annotation_inf, type == 
                                   "SIMPLE_MOLECULE")
  species_unique <- unique(met_without_protein0$species)
  name_unique <- unique(met_without_protein0$name)
  txt4_1 <- vector()
  for (i in seq(length(species_unique))) {
    txt4_1 <- c(txt4_1, getAnnotation(species01 = species_unique[i], 
                                      name01 = name_unique[i]))
  }
  protein_specie <- protein0$species
  protein_id <- protein0$protein_id
  type_p <- protein0$type
  name_p <- protein0$name
  txt4_2 <- vector()
  for (i in seq(length(protein_id))) {
    txt4_2 <- c(txt4_2, getProAnnotation(species01 = protein_specie[i], 
                                         name01 = name_p[i], proteinID01 = protein_id[i]))
  }
  gene_id <- gene0$gene_id
  gene_specie <- gene0$species
  type_g <- gene0$type
  name_g <- gene0$name
  txt4_3 <- vector()
  for (i in seq(length(gene_id))) {
    txt4_3 <- c(txt4_3, getGeneAnnotation(species01 = gene_specie[i], 
                                          name01 = name_g[i], geneID01 = gene_id[i]))
  }
  txt4 <- c(txt4_1, txt4_2, txt4_3)
  txt_final0 <- c(txt1, txt2, txt30, txt4, txt5, txt7)
  R3 <- readLines(file("data/txt6_repeated_part3"))
  R5 <- readLines(file("data/txt6_repeated_part5.xml"))
  R7 <- readLines(file("data/txt6_repeated_part7"))
  rxn <- rxn_core_carbon_inf
  rxn_p0 <- gpr_inf
  rxnID <- unique(rxn$rxnID)
  rxn_annotation <- list()
  rxn_id <- list()
  for (i in 1:length(rxnID)) {
    rxn_id[[i]] <- which(rxn$rxnID %in% rxnID[i] == TRUE)
    rxn_annotation[[i]] <- rxn[rxn_id[[i]], ]
  }
  gpr0 <- list()
  for (i in 1:length(rxn_p0$rxnID)) {
    gpr0[[i]] <- rxn_p0[i, ]
  }
  txt6_new <- vector()
  for (j in 1:length(rxnID)) {
    txt6_new <- c(txt6_new, getRxnInformation(rxn_annotation[[j]], 
                                              rxn_p = gpr0[[j]], R30 = R3, R70 = R7))
  }
  txt_new <- c(txt1, txt2, txt30, txt4, txt5, txt6_new, txt7)
  #writeLines(txt_final0, file("result/model_add_metabolite_protein.xml"))
  outfile_name <- paste("result/",subsystem0, ".xml", sep = "")
  writeLines(txt_new, file(outfile_name))
}









