# An example to show how we can map the predicted gene targets on to the map
# the xml file obtained based on the automatic layout.
# otherwise the followed code can't be used


###### function to splite all the related genes
###### this function will be merged into FALCONET
getSpliteGene <- function(genelist, reactionName){
  ##split of or relation
  #genelist <- yeast_7_7$Gene.reaction.association
  #reactionName <-yeast_7_7$Rxn.name
  
  rxnNum <- length(reactionName)
  
  
  GR <- list() # gene relation
  ss <- list()
  
  
  Add_and <- function(x1){
    tt0 <- vector()
    tt <- x1
    tt_length <- length(tt)
    for (i in 1:tt_length){
      tt0[i] <- paste(i, tt[i], sep = ";")
    }
    return(tt0)
  }
  
  
  Add_or <- function(x1){ #x1 <- GR[[1]]
    tt0 <- vector()
    tt <- x1
    tt_length <- length(tt)
    for (i in 1:tt_length){
      tt0[i] <- paste(i, tt[i], sep = ";")
    }
    return(tt0)
  }
  
  
  ## split the first order of gene relation
  
  splitOrRelation <- function(genelist, reactionName){
    #genelist <- "((YAL023C and YDL095W) or YDL093W or YJR143C or YOR321W)"
    #reactionName <- "r_0016"
    GR <- list() # gene relation
    ss <- list()
    if (!(str_detect(genelist, "\\) or \\(") | str_detect(genelist, "\\) or") | str_detect(genelist, "or \\(") | str_detect(genelist, "\\) and \\(") | str_detect(genelist, "\\) and") | str_detect(genelist, "and \\("))){
      
      GR[1] <- genelist
      ss[[1]] <- paste(reactionName,GR, sep = "@none;" )
      
    } else {
      
      if (str_detect(genelist, "\\) or \\(")){
        GR[1] <- str_split(genelist,"\\) or \\(" )
        
      } else{
        GR[1] <- genelist
      }
      
      if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "\\) or")) {
        GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"\\) or")
      } else{
        GR[1] <- GR[1]
      }
      
      if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "or \\(")){
        GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"or \\(" )
        
      } else{
        GR[1] <- GR[1]
      }
      
      
      
      if (str_detect(genelist, "\\) and \\(")){
        GR[1] <- str_split(genelist,"\\) and \\(" )
        
      } else{
        GR[1] <- GR[1]
      }
      
      if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "\\) and")) {
        GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"\\) and")
      } else{
        GR[1] <- GR[1]
      }
      
      if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "and \\(")){
        GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"and \\(" )
        
      } else{
        GR[1] <- GR[1]
      }
      
      
      
      ##for genelist <- "((YAL023C and YDL095W) or YDL093W or YJR143C or YOR321W)" 
      ##reactionName <- "r_0016"
      ##tt <- GR[1]
      
      sl <-length(GR[[1]])
      for (i in 1:sl){
        if( (str_detect(GR[[1]][i],"\\(\\(")) | (str_detect(GR[[1]][i],"\\)\\)"))){
          GR[[1]][i] <- GR[[1]][i]
        } else if (str_detect(GR[[1]][i], "or")){
          GR[[1]][i]  <- str_split(GR[[1]][i], "or")
        } 
      }
      ##### above is for special condition
      
      
      
      GR[1] <- paste0(unlist(GR[1]),collapse = "@@")
      GR[1] <- str_split(GR[1],"@@")
      
      if (str_detect(genelist, "\\) or \\(") | str_detect(genelist, "\\) or") | str_detect(genelist, "or \\(")){
        ss[[1]] <- paste(reactionName,Add_or(GR[[1]]), sep = "@or" )
      } else {
        ss[[1]] <- paste(reactionName,Add_and(GR[[1]]), sep = "@and" )
      }
    }
    
    return(ss[[1]])
    
  }
  
  for (i in 1:rxnNum){
    
    ss[[i]] <- splitOrRelation(genelist[i],reactionName[i])
    
  }
  
  
  GR1 <- unlist(ss)
  
  GR2 <-  str_replace_all(GR1,"\\(", "")
  GR2 <-  str_replace_all(GR2,"\\)", "")
  GR3 <- str_split(GR2,";")
  tt0 <- length(GR3)
  GR4 <- data.frame( ID= character(tt0),gene=character(tt0), stringsAsFactors = FALSE )
  for (i in 1:tt0){
    GR4$ID[i] <- GR3[[i]][1]
    GR4$gene[i] <- GR3[[i]][2]
    
  }
  
  ## split the first order of gene relation
  rxnNum2 <- length(GR4$ID)
  GR_2 <- list()
  ss2 <- list()
  for (i in 1:rxnNum2){
    if (str_detect(GR4$gene[i], " or ")){
      GR_2[i] <- str_split(GR4$gene[i]," or " )
      ss2[[i]] <- paste(GR4$ID[i],Add_or(GR_2[[i]]), sep = "@or" )
      
    } else if(str_detect(GR4$gene[i], " and ")){
      GR_2[i] <- str_split(GR4$gene[i]," and " )
      ss2[[i]] <- paste(GR4$ID[i],Add_and(GR_2[[i]]), sep = "@and" )
      
    } else{
      GR_2[i] <- GR4$gene[i]
      ss2[[i]] <- paste(GR4$ID[i],GR_2[[i]], sep = "@none;" )
      
    }
  }
  
  
  ## obtain the final formula
  GR_3 <- unlist(ss2)
  GR_4 <- str_split(GR_3,";")
  tt <- length(GR_4)
  GR_5 <- data.frame( ID= character(tt),gene=character(tt), stringsAsFactors = FALSE )
  
  for (i in 1:tt){
    GR_5$ID[i] <- GR_4[[i]][1]
    GR_5$gene[i] <- GR_4[[i]][2]
  }
  
  GR_or_and <- str_split(GR_5$ID,"@")
  
  for (i in 1:tt){
    GR_5$IDnew[i] <- GR_or_and[[i]][1]
    GR_5$R1[i] <- GR_or_and[[i]][2]
    GR_5$R2[i] <- GR_or_and[[i]][3]
  }
  
  GPs_redesign <- select(GR_5,IDnew, R1, R2, gene) ##obtain the new format of GPRs, the can be base to add the new GPRs or correct GPRs
  
  return(GPs_redesign)
  
}

# Function to change the reaction color based on the in silico prediction results
defineRxnColorForDesign <- function(rxn_adjust, over_expression = "OE", down_regulation = "KD", knock_out = "KO") {
  color <- vector()
  for (i in seq_along(rxn_adjust)) {
    #i <- 7
    if (rxn_adjust[i] == over_expression) {
      color[i] <- "ffff0000" # red
    }
    else if (rxn_adjust[i] == down_regulation) {
      color[i] <- "ff0000ff" # blue:ff0000ff  green: ff0eb10e
    }
    else if (rxn_adjust[i] == knock_out) {
      color[i] <- "ff009999" # yellow:ffffff66  grey:ffb3d2ff 
    }
    else {
      color[i] <- "ff000000"
    }
  }
  return(color)
}

# Function to mapping the pathway design onto the map
pathwayDesignMapping <- function(input_map, flux_inf, output_map) {
  yeast_map <- readLines(file(input_map))
  index_rxn <- which(str_detect(yeast_map, "reaction metaid"))
  # extract the rxns
  rxnid <- yeast_map[index_rxn]
  rxnid0 <- str_split(rxnid, " ")
  rxnid1 <- vector()
  for (i in 1:length(rxnid0)) {
    rxnid1[i] <- rxnid0[[i]][3]
  }
  rxnid1 <- str_replace_all(rxnid1, "id=", "") %>%
    str_replace_all(., "\"", "")
  
  flux_data <- data.frame(rxnid = rxnid1, stringsAsFactors = FALSE)
  flux_data$value <- getSingleReactionFormula(flux_inf$actions, flux_inf$Abbreviation, flux_data$rxnid)
  #flux_data$value <- as.numeric(flux_data$value)
  
  # define the line color
  line_color <- defineRxnColorForDesign(rxn_adjust = flux_data$value)
  
  flux_data$line_color <- line_color
  
  
  # from index_rxn choose the two
  # it seems the followed code is not very good as the index_line1 is not consistent with the corresponding reaction ID
  
  
  # condition1--for the simple graph
  #index_line <- which(str_detect(yeast_map, "celldesigner:editPoints"))
  #index_line1 <- index_line + 1
  
  #for (i in 1:length(index_line1)) {
  #  oldline <- "line width=\"1.0\" color=\"ff000000\"/>"
  #  if (line_color[i] != "ff000000") {
  #    newline <- paste("line width=\"3.0\"", " color=\"", line_color[i], "\"/>", sep = "")
  #  } else {
  #    newline <- paste("line width=\"1.0\"", " color=\"", line_color[i], "\"/>", sep = "")
  #   }
  #  yeast_map[index_line1[i]] <- str_replace_all(yeast_map[index_line1[i]], oldline, newline)
  # }
  
  
  # condition2-- for the complex graph
  total_row_num <- length(yeast_map)
  index_rxn_with_last_row <- c(index_rxn, total_row_num)
  # give the line index range for each rxn
  flux_data$line_start_index <- NA
  flux_data$line_end_index <- NA
  for(i in 1:nrow(flux_data)){
    print(i)
    flux_data$line_start_index[i] <- index_rxn_with_last_row[i]+1
    flux_data$line_end_index[i] <- index_rxn_with_last_row[i+1]-1
  }
  
  for(i in 1:nrow(flux_data)){
    print(i)
    #i <- 2
    start0 <-  flux_data$line_start_index[i]
    end0 <- flux_data$line_end_index[i]
    newcolor <- paste("color=\"",line_color[i], "\"", sep = "")
    
    # whether change width of reaction
    if (flux_data$value[i]=="NA"){
      newwidth <- "width=\"1.0\""}
    else{
      newwidth <- "width=\"15.0\""
    }
    
    for(j in start0:end0){
      print(j)
      #j <- 11237
      yeast_map[j]
      print(yeast_map[j])
      if (str_detect(yeast_map[j], "<celldesigner:line width=")){
        ss0 <- unlist(str_split(yeast_map[j], " "))
        ss0[2] <- newwidth
        ss0[3] <- newcolor
        if(length(ss0)==3){
          ss0[3] <- paste(ss0[3],"/>",sep = "" )
        }
        
        yeast_map[j] <- paste0(ss0,collapse = " ")
      }
      
    }
    
  }
  
  
  writeLines(yeast_map, file(output_map))
}




# input the rxn-GPR relations
flux_map <- read_excel("data/flux_map.xlsx")
rxnid_gpr <- flux_map[,c("GPR","Abbreviation")]
rxnid_gpr <- rxnGeneMapping(rxnid_gpr)
colnames(rxnid_gpr) <- c("Gene", "Abbreviation")
# input the gene targets
predicted_targets <- read.table("data/isoleucine_targets.txt", header = TRUE, sep = "\t")
predicted_targets <- predicted_targets[,c("genes","actions")]
# mapping the gene's manipulation onto reactions
rxnid_gpr$actions <- getMultipleReactionFormula(predicted_targets$actions,predicted_targets$genes,rxnid_gpr$Gene)
rxnid_gpr <- rxnid_gpr[!is.na(rxnid_gpr$actions),]
rxnid_gpr$subsystem <- getSingleReactionFormula(flux_map$subsystem,flux_map$Abbreviation,rxnid_gpr$Abbreviation)


# test1
# model_test_check.xml is obtained based on the automatic layout of cellDesigner
pathwayDesignMapping(input_map ="result/model_test_check.xml",
                flux_inf = rxnid_gpr,
                output_map = "result/map_with_ME_strategy.xml")

# test2
# for maps by Zhengming, need some quality improvement
# carbon metabolism.xml is the core metabolite pathway by Zhengming
improveMapConsistency <- function(input_map = "data/carbon metabolism.xml", output_map = "data/carbon metabolism_update.xml") {
  yeast_map <- readLines(file(input_map))
  index_rxn <- which(str_detect(yeast_map, "reaction metaid"))
  # test
  for (i in index_rxn) {
    print(i)
    #i <- 11547
    # extract the rxns
    rxnid <- yeast_map[i]
    # update rxnid
    rxn_inf <- unlist(str_split(rxnid, " "))
    name20 <- rxn_inf[2]
    name2 <- str_replace(rxn_inf[2], "metaid", "id")
    name3 <- rxn_inf[3]
    if (name2 != name3) {
      name3 <- str_replace_all(name20, "metaid", "id")
    }
    rxn_inf[3] <- name3
    rxnid_new <- paste0(rxn_inf, collapse = " ")

    yeast_map[i] <- rxnid_new
  }

  writeLines(yeast_map, file(output_map))
}

improveMapConsistency(input_map = "data/carbon metabolism.xml", output_map = "data/carbon metabolism_update.xml")

pathwayDesignMapping(input_map ="data/carbon metabolism_update.xml",
                     flux_inf = rxnid_gpr,
                     output_map = "result/map_core_metabolic_pathway_ME_strategy.xml")


# test 3 - use a whole map
improveMapConsistency(input_map = "data/Yeast8.xml", output_map = "data/Yeast8_update.xml")

pathwayDesignMapping(input_map ="data/Yeast8_update.xml",
                     flux_inf = rxnid_gpr,
                     output_map = "result/Yeast8_update_ME_strategy.xml")
# note: Yeast8_update2 is a new version of Yeast8 map after refining the layout
pathwayDesignMapping(input_map ="data/Yeast8_v2.xml",
                     flux_inf = rxnid_gpr,
                     output_map = "result/Yeast8_v2_ME_strategy.xml")



#########################################
# generate more maps based on Iven result
#########################################
# data input
all_gene_targets <- read.delim2("/Users/xluhon/Documents/GitHub/CellFactory-yeast-GEM/results/production_targets/targetsMatrix_mech_validated.txt", stringsAsFactors = FALSE)
# input the rxn-GPR relations
flux_map <- read_excel("data/flux_map.xlsx")
rxnid_gpr <- flux_map[,c("GPR","Abbreviation")]
rxnid_gpr <- rxnGeneMapping(rxnid_gpr)
colnames(rxnid_gpr) <- c("Gene", "Abbreviation")
# input the gene targets
product_list <- colnames(all_gene_targets)
product_list0 <- product_list[5:length(product_list)]
# one example
product0 <- product_list0[1]
product0 <- product_list0[45]
predicted_targets <- all_gene_targets[,c("genes", product0)]
predicted_targets <- predicted_targets[predicted_targets[[product0]] !=0,]
predicted_targets$phenylethanol_fam_alc[predicted_targets[[product0]]==3] <- "OE"
predicted_targets$phenylethanol_fam_alc[predicted_targets[[product0]]==2] <- "KD"
predicted_targets$phenylethanol_fam_alc[predicted_targets[[product0]]==1] <- "KO"
colnames(predicted_targets) <- c("genes","shortname", "actions")

# mapping the gene's manipulation onto reactions
rxnid_gpr$actions <- getMultipleReactionFormula(predicted_targets$actions,predicted_targets$genes,rxnid_gpr$Gene)
rxnid_gpr <- rxnid_gpr[!is.na(rxnid_gpr$actions),]
rxnid_gpr$subsystem <- getSingleReactionFormula(flux_map$subsystem,flux_map$Abbreviation,rxnid_gpr$Abbreviation)

# define the output file name
out_dir <- paste("result/Yeast8_v3_ME_strategy_for_",  product0, ".xml", sep = "")
pathwayDesignMapping(input_map ="data/Yeast8_v3.xml",
                     flux_inf = rxnid_gpr,
                     output_map = out_dir)

