# ============================================================================== #
# Script: 19_summary_graphs_methods_markers_popsizes.R
# Title: Multi-Factor Summary Graphs & Comparison Tables
# Description:
#   Aggregates results across methods, marker densities, and population sizes into multi-panel graphical figures and tabular summaries.
#
# Inputs:
#   - results/<Method>_Folds_*_Times_* (All raw cross-validation result text files)
#
# Outputs:
#   - results/Table<namef>.csv (Statistical comparison table)
#   - results/Plot<namef>.png (Accuracy comparison plot)
#   - results/Ability_Plot<namef>.png (Predictive ability plot)
#
# Required R Packages: argparse, ggplot2, reshape2, ggpubr
# ============================================================================== #

# ------------------------------------------------------------------------------
# Directory Setup & Output Path Configuration
# ------------------------------------------------------------------------------
# Automatically create "results" folder for output deliverables
results_dir <- file.path(getwd(), "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}

# Automatically locate shared input data directory
if (dir.exists(file.path(getwd(), "data"))) {
  data_dir <- file.path(getwd(), "data")
} else if (dir.exists(file.path(dirname(getwd()), "data"))) {
  data_dir <- file.path(dirname(getwd()), "data")
} else {
  data_dir <- getwd()
}
data_path <- paste0(data_dir, "/")
# ------------------------------------------------------------------------------


##Once the cross-validation is done, After the results are stored in a text file#####
##########The following code can be used to generate the graphs and tables#####
# ---- Argument Handling (Combined Approach) ----
# setwd() configured via project environment

import_packages<-function(){
  #define packages to install
  packages <- c("argparse","ggplot2","reshape2","ggpubr")
  
  install.packages(setdiff(packages, rownames(installed.packages())))
  
  invisible(lapply(packages, library, character.only=TRUE))
}

import_packages()

# parser <- ArgumentParser(description = "CV graphs")
# 
# # Define Arguments 
# parser$add_argument("--dir", default = "data", required = TRUE,help = "Data directory")
# parser$add_argument("--graph_type", required = TRUE, help = "Graph type")
# parser$add_argument("--trait_name",required = TRUE, help = "Trait name")
# parser$add_argument("--pmethods", nargs = "+", default = c("GBLUP","BA","BB","BC","BL","BRR","EN","RR","RF"), help = "Prediction methods (space-separated)")
# parser$add_argument("--iter", type = "integer", default = 3, help = "Number of iterations")
# parser$add_argument("--fold", type = "integer", default = 3, help = "Number of folds")
# parser$add_argument("--pop_sizes", nargs = "+", default = c("NULL"), help = "Population sizes (space-separated)")
# parser$add_argument("--marker_sizes", nargs = "+", default = c("NULL"), help = "Marker sizes (space-separated)")
# parser$add_argument("--Fixed_effect", nargs = "+", default = c("NULL"), help = "Fixed effects (space-separated)")
# parser$add_argument("--h2", type = "double",required = TRUE, help = "Heritability")
# #parser$add_argument("--marker_start", type = "integer", help = "Starting marker size")
# #parser$add_argument("--marker_end", type = "integer", help = "Ending marker size")


# Parse Arguments
dir<-"data"
graph_type<-"A"
#trait_name<-c("Is_Dead")
#trait_name<-c("DPC_date")
trait_name<-c("DPC_time")
pmethods<-c("PBLUP","GBLUP","BA","BB","BC","BL","BRR","RR","EN","RF")
iter<-10
fold<-5
pop_sizes<-"NULL"
marker_sizes<-c("top_40000")
#marker_sizes<-c("top_40000")
#marker_sizes = c(50,"top_50",100,"top_100",500,"top_500",1000,"top_1000",5000,"top_5000",10000,"top_10000",20000,"top_20000",30000,"top_30000",40000,"top_40000",50000,"top_50000","ALL")
Fixed_effect<-c("NULL")
h2<-1

#sur 0.145
#dpc date 0.186
#dpc time 0.185

# Parse Arguments
dir<-"data"
graph_type<-"C"
#trait_name<-c("Is_Dead")
#trait_name<-c("DPC_date")
trait_name<-c("DPC_time")
pmethods<-c("GBLUP")
iter<-10
fold<-3
#fold<-5
pop_sizes = c("50","100", "150", "200", "250", "300", "350", "400", "474")
marker_sizes = c("top_1000")
#marker_sizes<-c("top_40000")
#marker_sizes = c(50,"top_50",100,"top_100",500,"top_500",1000,"top_1000",5000,"top_5000",10000,"top_10000",20000,"top_20000",30000,"top_30000",40000,"top_40000",50000,"top_50000","ALL")
Fixed_effect<-c("NULL")
h2<-1



# Calculate marker sizes sequence (Added)
#marker_sizes <- unique(c(
#  seq(args$marker_start, min(1000, args$marker_end), by = 100),  # 100 to 1000 or end by 100
#  seq(1000, min(10000, args$marker_end), by = 1000),             # 1000 to 10000 or end by 1000
#  seq(10000, args$marker_end, by = 10000)                         # 10000 to end by 10000
#))

# Check if last marker is below 10000 and add the end marker if needed (Added)
#if (max(marker_sizes) < args$marker_end) {
#  marker_sizes <- c(marker_sizes, args$marker_end) 
#}

data_path = paste0(dir, "/")

if (graph_type=="A"){
  
  trait<-trait_name
  pmethod<-pmethods
  compareM<-matrix(0,iter,length(pmethod))
  for (x in seq_along(pmethod)) {
    datap <- read.table(paste0(data_path,
                               pmethod[x], "_Folds_", fold, "_Times_", iter, "_", 
                               paste(trait, collapse = "_"), "_", paste(Fixed_effect, collapse = "_"), 
                               "_", paste(pop_sizes, collapse = "_"), "_", 
                               paste(marker_sizes, collapse = "_"), ".txt"
    ))
    for (n in 1:iter) {
      compareM[n, x] <- cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
    }
  }
  
  colnames(compareM) <- pmethod
  rownames(compareM) <- paste0("CV", 1:iter)
  
  namef<-paste0("_Folds_", fold, "_Times_", iter, "_", 
                paste(trait, collapse = "_"), "_", paste(Fixed_effect, collapse = "_"), 
                "_", paste(pop_sizes, collapse = "_"), "_", 
                paste(marker_sizes, collapse = ""))
  
  write.table(compareM,paste0(data_path,"Table",namef,".csv"),sep = ",",quote = F)
  data<-reshape2::melt(compareM)
  
  table_graphic <- ggtexttable(round(compareM,4))
  ggsave(paste0(data_path, "Table", namef, ".png"), table_graphic, scale=1, dpi = 300, height = 2, width = 8)
  
  p<-ggplot(data, aes(x=Var2, y=value/sqrt(h2))) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    ylim(-0.1, 0.3)+
    labs(title=paste0("Plot_",trait),x= "Prediction method",fill="Method", y = "Prediction Accuracy")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  ggsave(paste0(data_path,"Plot",namef,".png"),p,scale=1, dpi = 300, height = 4, width = 4)
  
  p<-ggplot(data, aes(x=Var2, y=value)) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    ylim(-0.1, 0.3)+
    labs(title=paste0("Plot_",trait),x= "Prediction method",fill="Method", y = "Prediction Ability")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  ggsave(paste0(data_path,"Ability_Plot",namef,".png"),p,scale=1, dpi = 300, height = 4, width = 4)
}

if (graph_type=="B"){
  
  pmethod<-pmethods
  sname<-trait_name
  marker_sizes<-marker_sizes
  # Use the marker_sizes generated earlier
  compareM <- matrix(0, iter, length(marker_sizes))
  
  for (i in seq_along(marker_sizes)) {  
    filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, "_", sname, "_", paste(Fixed_effect, collapse = "_"), "_", pop_sizes, "_", marker_sizes[i], ".txt")
    datap <- read.table(filename)
    for (n in 1:iter) {
      compareM[n, i] <- cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
    }
  }
  
  # Use marker_sizes for column names
  colnames(compareM) <- paste0("m=", marker_sizes)
  rownames(compareM) <- paste0("CV", 1:iter)
  
  namef<-paste0("_Folds_", fold, "_Times_", iter, "_",sname,"_", paste(Fixed_effect, collapse = "_"), "_", pop_sizes, "_marker")
  
  write.table(compareM,paste0(data_path,"Table_",pmethod,namef,".csv"),sep = ",",quote = F)
  
  data<-reshape2::melt(compareM)
  
  p<-ggplot(data, aes(x=Var2, y=value/sqrt(h2))) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    #ylim(0.37, 0.81)+
    labs(title=paste0("Plot_",sname),x= "Marker No",fill="Markers", y = "Prediction Accuracy")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  ggsave(paste0(data_path,"Plot_",pmethod,namef,".png"),p,scale=1, dpi = 300, height = 4, width = 12)
  
  p<-ggplot(data, aes(x=Var2, y=value)) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    #ylim(0.37, 0.81)+
    labs(title=paste0("Plot_",sname),x= "Marker No",fill="Markers", y = "Prediction Ability")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_text(angle = 45, size = 8),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  ggsave(paste0(data_path,"Ability_Plot_",pmethod,namef,".png"),p,scale=1, dpi = 300, height = 4, width = 12)
  
}
if (graph_type=="C"){
  
  sname<-trait_name
  pop_sizes<-pop_sizes
  pmethod<-pmethods
  compareM<-matrix(0,iter,length(pop_sizes))
  
  for (i in seq_along(pop_sizes)) {  
    filename <- paste0(data_path, pmethod, "_Folds_", fold, "_Times_", iter, "_", sname, "_", Fixed_effect, "_", pop_sizes[i], "_", marker_sizes, ".txt") 
    datap <- read.table(filename)
    for (n in 1:iter) {
      compareM[n, i] <- cor(datap[, paste0("Predicted", n)], datap[, paste0("Actual", n)])
    }
  }
  
  colnames(compareM) <- paste0("s=", pop_sizes)
  rownames(compareM) <- paste0("CV", 1:iter)
  
  namef<-paste0("_Folds_", fold, "_Times_", iter, "_",sname,"_", Fixed_effect, "_pop_", marker_sizes)
  
  write.table(compareM,paste0(data_path,"Table_",pmethod,namef,".csv"),sep = ",",quote = F)
  
  data<-reshape2::melt(compareM)
  
  p<-ggplot(data, aes(x=Var2, y=value/sqrt(h2))) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    #ylim(0.3, 0.8)+
    labs(title=paste0("Plot_",sname),x= "Pop size",fill="Pop size", y = "Prediction Accuracy")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  
  ggsave(paste0(data_path,"Plot_",pmethod,namef,".png"),p,scale=1, dpi = 300, height = 4, width = 12)
  
  p<-ggplot(data, aes(x=Var2, y=value)) + 
    geom_boxplot(aes(fill = Var2),color="black") +
    #ylim(0.3, 0.8)+
    labs(title=paste0("Plot_",sname),x= "Pop size",fill="Pop size", y = "Predictive Ability")+
    stat_summary(fun.y=mean, geom="point", shape=21, size=1, color="black",fill="red") +
    theme(legend.title = element_text(colour="blue", size=12,face="bold"),
          plot.title = element_text(hjust=0.5, size=12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.text.x=element_blank(),
          plot.background = element_rect(fill = "#FFFFFF"),
          legend.background = element_rect(linetype = "dashed"),
          panel.border = element_rect(fill = NA))
  
  ggsave(paste0(data_path,"Ability_Plot_",pmethod,namef,".png"),p,scale=1, dpi = 300, height = 4, width = 12)
}

