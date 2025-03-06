gc()
load_or_install_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% rownames(installed.packages())) {
        warning(pkg, " is not available")
      } else {
        BiocManager::install(pkg, ask = FALSE)
      }
    }
    require(pkg, character.only = TRUE, quietly = TRUE)
  }
}
packages <- c("ggplot2", "dplyr", "tidyverse", "tidyr", "ggupset", "corrplot", "ggvenn", "limma", "graphics", "directPA", "reactome.db", "org.Hs.eg.db", "SummarizedExperiment", "PhosR", "calibrate", "ggpubr", "GGally", "network", "annotate", "stats4", "methods", "datasets", "gridExtra", "venn", "UpSetR", "grid", "reshape2", "extrafont", "PerformanceAnalytics", "xts", "zoo", "ClueR", "parallel", "e1071", "lubridate", "forcats", "stringr", "purrr", "readr", "tibble", "GenomicRanges", "GenomeInfoDb", "MatrixGenerics", "matrixStats", "XML", "AnnotationDbi", "IRanges", "S4Vectors", "Biobase", "BiocGenerics", "stats4", "calibrate", "MASS", "svglite", "pheatmap", "ggplotify","doParallel", "enrichplot" ,"clusterProfiler", "ggrepel", "devtools", "volcano3D","edgeR","RColorBrewer","readxl","ggVennDiagram", "ggbeeswarm", "rlang", "xtable")
load_or_install_packages(packages)

data(PhosphoSitePlus)
Sys.setenv(LANG = "en")


names_map <- c("20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_1_phospho" = "WT_Endogenous_Rep1",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_2_phospho" = "WT_Endogenous_Rep2",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_3_phospho" = "WT_Endogenous_Rep3",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_4_phospho" = "WT_Rescue_Rep1",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_5_phospho" = "WT_Rescue_Rep2",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_6_phospho" = "WT_Rescue_Rep3",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_13_phospho" = "DTU_Endogenous_Rep1",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_14_phospho" = "DTU_Endogenous_Rep2",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_15_phospho" = "DTU_Endogenous_Rep3",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_16_phospho" = "DTU_Rescue_Rep1",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_17_phospho" = "DTU_Rescue_Rep2",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_18_phospho" = "DTU_Rescue_Rep3",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_25_phospho" = "CS_Endogenous_Rep1",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_26_phospho" = "CS_Endogenous_Rep2",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_27_phospho" = "CS_Endogenous_Rep3",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_28_phospho" = "CS_Rescue_Rep1",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_29_phospho" = "CS_Rescue_Rep2",
               "20221127_VC_140min_noFAIMS_1ug_HRMS1_DIA_sample_30_phospho" = "CS_Rescue_Rep3")
old_colnames <- colnames(Keratinocyte_Phosphoproteome_2024)[1:18]
names_map <- c(names_map, setdiff(old_colnames, names(names_map)))
# Use the names_map to rename the columns
colnames(Keratinocyte_Phosphoproteome_2024)[1:18] <- names_map[old_colnames]
renamed_df = Keratinocyte_Phosphoproteome_2024
your_dataframe <- renamed_df %>%
  mutate(Column_3 = str_extract(PTM_collapse_key, "(?<=_).*?(?=_)"),
         Combined_Column = paste(PG.UniProtIds, PG.Genes, Column_3, PEP.StrippedSequence, sep=";"))
your_dataframe <- your_dataframe %>%
  distinct(Combined_Column, .keep_all = TRUE)


ppe <- your_dataframe 
class(ppe)
dim((ppe))
ppe_matrix <- as.matrix(ppe)
grps = gsub("_Rep[0-9]", "", colnames(ppe_matrix[,c(1:18)]))
#grps_list <- lapply(grps, function(grps) {
#  grep(grps, colnames(ppe_matrix[,c(1:36)]), value = TRUE)
#})
#grps= grps_list
# ONLY SELECT NO Serum
#grps = grps[(grep("_NoSerum", (grps)))]

rownames(ppe_matrix) = unique(ppe_matrix[,"Combined_Column"])
#ppe_matrix = ppe_matrix[,(grep("_NoSerum", colnames(ppe_matrix)))]

ppe_matrix = ppe_matrix[,c(1:18)]
ppe_filtered <- selectGrps(ppe_matrix, grps, 0.7, n=1) 
rownames_original <- rownames(ppe_filtered)

ppe_filtered <- apply(ppe_filtered, 2, as.numeric)

# Reassign the original row names to the data frame
# 
rownames(ppe_filtered) <- rownames_original

#dim(ppe_imputed) # log2 before applying
ppe_filter_log = na.omit(((log2(ppe_filtered))))
dim(ppe_filter_log)


FASTA = reformat_fasta_df(uniprot.proteome_Human)

named_FASTA <- fasta_to_named_list(FASTA)

strip_phos_site <- function(GOI) {
  gene_IDs <- gsub("_.*", "", GOI)
  return(gene_IDs)
}
#Assign flanking sequence to phos site instead of tryptic peptide
rownames(ppe_filter_log) = extract_flanking_seq(named_FASTA,rownames(ppe_filter_log))

ppe_filter_log = as.data.frame(ppe_filter_log)

# Create a data.frame with the rownames as a new column
ppe_filter_log_with_id <- ppe_filter_log %>%
  rownames_to_column(var = "ID")

# Modify the ID column to remove the ".1" at the end of the ID if it exists
ppe_filter_log_with_id$ID <- gsub("(\\.1)$", "", ppe_filter_log_with_id$ID)

# Group by the modified ID column and calculate the mean for each group
ppe_filter_log_mean <- ppe_filter_log_with_id %>%
  group_by(ID) %>%
  summarise_all(mean)
# Restore the rownames from the ID column
ppe_filter_log_mean_matrix <- as.matrix(ppe_filter_log_mean)
rownames(ppe_filter_log_mean_matrix) <- ppe_filter_log_mean_matrix[,1]

# Convert columns to numeric using apply()
numeric_columns <- apply(ppe_filter_log_mean_matrix, 2, as.numeric)
# Create a new data frame or matrix with the numeric columns and original rownames
numeric_ppe_filter_log <- as.matrix(numeric_columns)
rownames(numeric_ppe_filter_log) <- rownames(ppe_filter_log_mean_matrix)

ppe_filter_log_mean_matrix = numeric_ppe_filter_log
ppe_filter_log_mean_matrix <- ppe_filter_log_mean_matrix[,-1,drop=FALSE]
# Convert the data frame to a matrix
ppe_filter_log = ppe_filter_log_mean_matrix
rownames(ppe_filter_log) <- gsub(pattern = "\\.", replacement = ";", x = rownames(ppe_filter_log))

#Re iterate to get the flanking sequence of doubles and tryptic peptides
rownames(ppe_filter_log) = extract_flanking_seq(named_FASTA,rownames(ppe_filter_log))

ppe_filter_log = as.data.frame(ppe_filter_log)

# Create a data.frame with the rownames as a new column
ppe_filter_log_with_id <- ppe_filter_log %>%
  rownames_to_column(var = "ID")

# Modify the ID column to remove the ".1" at the end of the ID if it exists
ppe_filter_log_with_id$ID <- gsub("(\\.1)$", "", ppe_filter_log_with_id$ID)

# Group by the modified ID column and calculate the mean for each group
ppe_filter_log_mean <- ppe_filter_log_with_id %>%
  group_by(ID) %>%
  summarise_all(mean)
# Restore the rownames from the ID column
ppe_filter_log_mean_matrix <- as.matrix(ppe_filter_log_mean)
rownames(ppe_filter_log_mean_matrix) <- ppe_filter_log_mean_matrix[,1]

# Convert columns to numeric using apply()
numeric_columns <- apply(ppe_filter_log_mean_matrix, 2, as.numeric)
# Create a new data frame or matrix with the numeric columns and original rownames
numeric_ppe_filter_log <- as.matrix(numeric_columns)
rownames(numeric_ppe_filter_log) <- rownames(ppe_filter_log_mean_matrix)

ppe_filter_log_mean_matrix = numeric_ppe_filter_log
ppe_filter_log_mean_matrix <- ppe_filter_log_mean_matrix[,-1,drop=FALSE]
# Convert the data frame to a matrix
ppe_filter_log = ppe_filter_log_mean_matrix
rownames(ppe_filter_log) <- gsub(pattern = "\\.", replacement = ";", x = rownames(ppe_filter_log))

rownames_split <- strsplit(rownames(ppe_filter_log), ";")
uniprot <- sapply(rownames_split, "[", 1)
gene_symbols <- sapply(rownames_split, "[", 2)
sites <- gsub("[STY]", "", sapply(rownames_split, "[", 3))
sequences <- sapply(rownames_split, "[", 4)
residues <- gsub("[0-9]", "", sapply(rownames_split, "[", 3))
# Create PhosphoExperiment object

ppe_noSerum <- PhosphoExperiment(assays = list(Quantification = ppe_filter_log),UniprotID = uniprot, Site = sites, GeneSymbol = gene_symbols, Residue = residues, Sequence = sequences)


ppe = ppe_noSerum

ppe = ppe_noSerum 
phosphoL6 = (SummarizedExperiment::assay(ppe, "Quantification")) # IS 3997 sites! If not iterate the flanking annotation
# Define the desired order
desired_order <- c("WT_Endogenous_Rep1", "WT_Endogenous_Rep2", "WT_Endogenous_Rep3", 
                   "WT_Rescue_Rep1", "WT_Rescue_Rep2", "WT_Rescue_Rep3",
                   "DTU_Endogenous_Rep1", "DTU_Endogenous_Rep2", "DTU_Endogenous_Rep3",
                   "DTU_Rescue_Rep1", "DTU_Rescue_Rep2", "DTU_Rescue_Rep3",
                   "CS_Endogenous_Rep1", "CS_Endogenous_Rep2", "CS_Endogenous_Rep3",
                   "CS_Rescue_Rep1", "CS_Rescue_Rep2", "CS_Rescue_Rep3")

# Reorder the columns
phosphoL6 <- phosphoL6[, desired_order]
grps = gsub("_Rep[0-9]", "", colnames(ppe))

# Fit the model and contrasts
design = model.matrix(~ grps - 1)
colnames(design) <- c("KO_CS_NoSerum_Control", "KO_CS_NoSerum_Treated", "KO_DTU_NoSerum_Control", "KO_DTU_NoSerum_Treated", "WT_NoSerum_Control", "WT_NoSerum_Treated")
# Update the contrasts for NoSerum data
contrasts <- makeContrasts(
  KO_WT_CS_NoSerum_fold_change = WT_NoSerum_Control- KO_CS_NoSerum_Control,
  KO_WT_DTU_NoSerum_fold_change = WT_NoSerum_Control- KO_DTU_NoSerum_Control,
  KO_CS_NoSerum_fold_change = KO_CS_NoSerum_Treated - KO_CS_NoSerum_Control,
  KO_DTU_NoSerum_fold_change = KO_DTU_NoSerum_Treated - KO_DTU_NoSerum_Control,
  KO_CS_control_WT_control_foldchange = KO_CS_NoSerum_Control - WT_NoSerum_Control,
  KO_DTU_control_WT_control_foldchange = KO_DTU_NoSerum_Control - WT_NoSerum_Control,
  TreatmentEffect_NoSerum_WT = WT_NoSerum_Treated - WT_NoSerum_Control,
  GenotypeTreatmentEffect_NoSerum_DTU = KO_DTU_NoSerum_Treated - WT_NoSerum_Control,
  GenotypeTreatmentEffect_NoSerum_CS = KO_CS_NoSerum_Treated - WT_NoSerum_Control,
  InteractionEffect_CS = (KO_CS_NoSerum_Treated - WT_NoSerum_Treated) - (KO_CS_NoSerum_Control - WT_NoSerum_Control),
  InteractionEffect_DTU = (KO_DTU_NoSerum_Treated - WT_NoSerum_Treated) - (KO_DTU_NoSerum_Control - WT_NoSerum_Control),
  WT_NoSerum_Control = WT_NoSerum_Control - WT_NoSerum_Control,
  CS_DTU_control = KO_CS_NoSerum_Control - KO_DTU_NoSerum_Control,
  levels = design
)


fit <- lmFit(SummarizedExperiment::assay(ppe, "Quantification"), design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Updated tables
table_names <- colnames(coef(fit2))
tables <- lapply(table_names, function(x) topTable(fit2, number=Inf, coef=x))

idx_Up = fit2$coefficients[,c( "KO_WT_CS_NoSerum_fold_change")] > -Inf
phosphoL6.reg_Up <- (SummarizedExperiment::assay(ppe[idx_Up, ], "Quantification"))
desired_order <- c("WT_Endogenous_Rep1", "WT_Endogenous_Rep2", "WT_Endogenous_Rep3", 
                   "WT_Rescue_Rep1", "WT_Rescue_Rep2", "WT_Rescue_Rep3",
                   "DTU_Endogenous_Rep1", "DTU_Endogenous_Rep2", "DTU_Endogenous_Rep3",
                   "DTU_Rescue_Rep1", "DTU_Rescue_Rep2", "DTU_Rescue_Rep3",
                   "CS_Endogenous_Rep1", "CS_Endogenous_Rep2", "CS_Endogenous_Rep3",
                   "CS_Rescue_Rep1", "CS_Rescue_Rep2", "CS_Rescue_Rep3")

# Reorder the columns
phosphoL6.reg_Up <- phosphoL6.reg_Up[, desired_order]

L6.phos.std_Up <- standardise(phosphoL6.reg_Up)

rownames(L6.phos.std_Up) <- paste0(GeneSymbol(ppe), ";", Residue(ppe), Site(ppe), ";")[idx_Up]
L6.phos.seq_Up <- Sequence(ppe)[idx_Up]
empty_or_missing_or_H1_seqs_Up <- sapply(L6.phos.seq_Up, function(x) x == "" || is.na(x) || grepl("S[0-9]", x))
L6.phos.seq_Up <- L6.phos.seq_Up[!empty_or_missing_or_H1_seqs_Up]
L6.phos.std_Up <- L6.phos.std_Up[!empty_or_missing_or_H1_seqs_Up,]
data(PhosphoSitePlus)
L6.matrices_Up_CS <- kinaseSubstrateScore(substrate.list = PhosphoSite.human, 
                                          mat = L6.phos.std_Up, seqs = L6.phos.seq_Up, 
                                          numMotif = 5, numSub = 1, verbose = FALSE, species = "human")

rownames(L6.matrices_Up_CS$combinedScoreMatrix) = gsub(";", " ", rownames(L6.matrices_Up_CS$combinedScoreMatrix))
rownames(L6.matrices_Up_CS$profileScoreMatrix) = gsub(";", " ", rownames(L6.matrices_Up_CS$profileScoreMatrix))
rownames(L6.matrices_Up_CS$motifScoreMatrix) = gsub(";", " ", rownames(L6.matrices_Up_CS$motifScoreMatrix))
rownames(L6.matrices_Up_CS$ksActivityMatrix) = gsub(";", " ", rownames(L6.matrices_Up_CS$ksActivityMatrix))
library(PhosR)
data("KinaseFamily", package = "PhosR")

kinaseHeatmap_Up_CS = (kinaseSubstrateHeatmap_modified(L6.matrices_Up_CS, top =1))



motifScore_Up_CS <- L6.matrices_Up_CS$motifScoreMatrix
profileScore_Up_CS <- L6.matrices_Up_CS$profileScoreMatrix
combinedScore_Up_CS <- L6.matrices_Up_CS$combinedScoreMatrix
kinaseActivityScore_Up_CS <- L6.matrices_Up_CS$ksActivityMatrix

L6.predMat_Up_CS <- kinaseSubstratePred(L6.matrices_Up_CS, cs = 0.5, top = 50, inclusion = 20, iter= 5, verbose = FALSE)
# ksActivityScore violin plots
reshaped_data <- melt(L6.matrices_Up_CS$ksActivityMatrix, variable.name = "Condition", value.name = "Activity")
# Convert Var2 to character type and split into two parts: Genotype and Treatment
reshaped_data$Var2 <- as.character(reshaped_data$Var2)
reshaped_data$Genotype <- sapply(strsplit(reshaped_data$Var2, "_"), function(x) {
  if (x[1] == "WT") {
    return(x[1])
  } else {
    return(x[1])
  }
})
reshaped_data$Treatment <- sapply(strsplit(reshaped_data$Var2, "_"), function(x) {
  if (x[2] %in% c("Endogenous", "Rescue")) {
    return(x[2])
  } else {
    return(x[2])
  }
})
# Replace the values in the Genotype and Treatment columns
reshaped_data$Genotype[reshaped_data$Genotype == "WT"] <- "WT"
reshaped_data$Genotype[reshaped_data$Genotype == "DTU"] <- "DMKN \U03B2\U03B3-/-"
reshaped_data$Genotype[reshaped_data$Genotype == "CS"] <- "DMKN \U03B1\U03B2-/-"

reshaped_data$Treatment[reshaped_data$Treatment == "Control"] <- "Endogenous"
reshaped_data$Treatment[reshaped_data$Treatment == "Rescue"] <- "Supplemented"
reshaped_data$Genotype <- factor(reshaped_data$Genotype, levels = c("WT", "DMKN \U03B2\U03B3-/-", "DMKN \U03B1\U03B2-/-"))
reshaped_data$Treatment <- as.factor(reshaped_data$Treatment)

significant_results <- reshaped_data %>%
  group_by(Var1) %>%
  do({
    data.frame(Genotype = unique(.$Genotype),
               p_value = summary(aov(Activity ~ Genotype, data = .))[[1]]["Pr(>F)"][1,])
  }) %>%
  filter(p_value < 0.05)

# Filter the data based on significant results
significant_genotypes <- significant_results$Genotype
filtered_data <- reshaped_data %>%
  filter(Genotype %in% significant_genotypes & Var1 %in% significant_results$Var1)

significant_results <- significant_results %>%
  mutate(p_value = formatC(p_value, format = "f", digits = 3))

cairo_pdf("O:/MBT-g-Cell-Signaling-Architecture/People Folders/Canbay Vahap/Dermokine_Paper/Supplementary/ksActivityScore_sig_v5.pdf", width = 3.5, height = 8.5, family = "arial")
ksActivityScore <- ggplot(filtered_data, aes(x = Genotype, y = Activity, fill = Treatment)) +
  geom_violin(width = 0.7, alpha = 0.9, color = "black", size = 0.2) +
  ggforce::geom_sina(size = 0.1) +
  facet_wrap(~ Var1, ncol = 3) +  # Set the number of columns to 3
  stat_compare_means(aes(group = Genotype), method = "aov", label.y = 2.0, size = 3) +
  #stat_compare_means(aes(group = Treatment), method = "t.test", label.y = 1.6, size = 4, 
  #label = "p.signif", hide.ns = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Increase font size
        axis.text.y = element_text(size = 10),  # Increase font size
        strip.text = element_text(size = 10, face = "bold"),  # Increase font size
        legend.title = element_text(size = 10, face = "bold"),  # Increase font size
        legend.text = element_text(size = 10),  # Increase font size
        text = element_text(size = 10),  # Increase font size
        legend.position = "bottom") +
  labs(x = "Genotype", y = "Kinase Activity Score", fill = "Treatment")
print(ksActivityScore)
dev.off()



L6.matrices_Up_CS_filter= L6.matrices_Up_CS

# Filter the columns of L6.matrices_Up_CS$combinedScoreMatrix to include only the significant kinases
L6.matrices_Up_CS_filter$combinedScoreMatrix <- L6.matrices_Up_CS_filter$combinedScoreMatrix[, colnames(L6.matrices_Up_CS_filter$combinedScoreMatrix) %in% significant_kinases]
L6.matrices_Up_CS_filter$profileScoreMatrix <- L6.matrices_Up_CS_filter$profileScoreMatrix[, colnames(L6.matrices_Up_CS_filter$profileScoreMatrix) %in% significant_kinases]
L6.matrices_Up_CS_filter$motifScoreMatrix <- L6.matrices_Up_CS_filter$motifScoreMatrix[, colnames(L6.matrices_Up_CS_filter$motifScoreMatrix) %in% significant_kinases]
L6.matrices_Up_CS_filter$ksActivityMatrix <- L6.matrices_Up_CS_filter$ksActivityMatrix[rownames(L6.matrices_Up_CS_filter$ksActivityMatrix) %in% significant_kinases,]

# Filter the columns of L6.predMat_Up_CS to include only the significant kinases
filtered_L6.pred_matrix <- kinaseSubstratePred(L6.matrices_Up_CS_filter, cs = 0.5, top = 50, inclusion = 20, iter= 5, verbose = FALSE)
# Remove rows with missing names from L6.predMat_Up_CS
filtered_L6.pred_matrix <- filtered_L6.pred_matrix[!is.na(rownames(filtered_L6.pred_matrix)), ]
filtered_L6.pred_matrix <- filtered_L6.pred_matrix[!duplicated(rownames(filtered_L6.pred_matrix)), ]

PhosR:::.kinaseNetwork(filtered_L6.pred_matrix, L6.matrices_Up_CS_filter, threskinaseNetwork=0.5, kinase_signalome_color)
PhosR:::.phosphositeClusters(L6.matrices_Up_CS_filter, verbose=T)


Signalomes_results_ENDO_CS_Up_filtered <- Signalomes(KSR=L6.matrices_Up_CS_filter, 
                                                     predMatrix=filtered_L6.pred_matrix, filter =T,
                                                     exprsMat=L6.phos.std_Up, module_res = 8,threskinaseNetwork= 0.95,
                                                     signalomeCutoff = 0.5, KOI = "SRC", verbose = T)

cairo_pdf("C:/Users/vahcan/OneDrive - Danmarks Tekniske Universitet/DMKN Manuscript/Figures/Phospho and KSEA/Signalomes_results_ENDO_CS_Up_filtered.pdf", width = 8, height = 8)
print(Signalomes(KSR=L6.matrices_Up_CS_filter, 
                 predMatrix=filtered_L6.pred_matrix, filter =T,
                 exprsMat=L6.phos.std_Up, module_res = 8,threskinaseNetwork= 0.95,
                 signalomeCutoff = 0.5, KOI = "SRC", verbose = F))
dev.off()

module_inlcuding_proteins = process_signalomes_results_without_converting(Signalomes_results_ENDO_CS_Up_filtered)


moduleUp_entrez <- process_signalomes_results(Signalomes_results_ENDO_CS_Up_filtered)
max_length <- max(sapply(moduleUp_entrez, length))
module_df <- as.data.frame(matrix(NA, nrow = max_length, ncol = length(moduleUp_entrez)))
colnames(module_df) <- names(moduleUp_entrez)

# Fill the data frame with the IDs from each module
for (module_name in names(moduleUp_entrez)) {
  module_df[1:length(moduleUp_entrez[[module_name]]), module_name] <- moduleUp_entrez[[module_name]]
}


# Convert the wide data frame to a long format
long_module_df <- module_df %>%
  gather(module, Entrez, everything(), na.rm = TRUE)
library(clusterProfiler) ; library(enrichplot) ; library(org.Hs.eg.db)

# Use compareCluster with the reshaped data frame
ENDO_CS_UP_cluster_significant_GO <- compareCluster(Entrez ~ module, data = long_module_df, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                                                    ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.01,
                                                    qvalueCutoff = 0.05, minGSSize = 20, readable = TRUE)

# Rename the modules; adjust depending on howm any clusters were generated
new_module_names <- c("1", "2", "3","4", "5", "6", "7", "8")
# Update the Cluster column in the compare_result object
ENDO_CS_UP_cluster_significant_GO@compareClusterResult <- ENDO_CS_UP_cluster_significant_GO@compareClusterResult %>%
  mutate(Cluster = str_replace(Cluster, "module_", ""))

ENDO_CS_UP_cluster_significant_clusterPLOT = dotplot(ENDO_CS_UP_cluster_significant_GO, font.size = 6, size = NULL, x = "GeneRatio",label_format =50, showCategory=3,includeAll = T,split="ONTOLOGY", shape =T)  +
  scale_size_continuous(range = c(1, 4)) + 
  theme( strip.text = element_text(face="bold"), axis.text.x = element_text(face="bold", size =6),
         axis.title.x = element_text(size=6,face="bold" ), 
         legend.key.size = unit(0.5, "cm"),
         legend.text=element_text(size=6),
         legend.position = "bottom",
         legend.box = "vertical",
         text = element_text(size=6),
         strip.background = element_rect()) + labs(x="Gene Ratio") +facet_grid(ONTOLOGY~.,scales = "free")

# Filter rows where 'Description' contains any of the terms
filtered_GO_terms <- ENDO_CS_UP_cluster_significant_GO@compareClusterResult[grep("cadherin binding|adherens junction", ENDO_CS_UP_cluster_significant_GO@compareClusterResult$Description), ]

# Split 'geneID' column into lists of genes
gene_lists_per_term <- strsplit(filtered_GO_terms$geneID, "/")

# Find the intersecting genes
intersecting_genes_among_terms <- Reduce(intersect, gene_lists_per_term)
all_genes <- unlist(gene_lists_per_term)

# Count the frequency of each gene
gene_frequency <- table(all_genes)

# Convert the table to a data frame and sort by frequency
gene_frequency_df <- data.frame(Gene = names(gene_frequency), Frequency = as.integer(gene_frequency))
gene_frequency_df <- gene_frequency_df[order(-gene_frequency_df$Frequency, gene_frequency_df$Gene), ]
# Get the column names from gene_frequency
all_genes <- names(gene_frequency)

# Filter rows of L6.matrices_Up_CS_filter$combinedScoreMatrix based on all_genes
filtered_score_matrix <- L6.matrices_Up_CS_filter$combinedScoreMatrix[grep(paste(all_genes, collapse="|"), rownames(L6.matrices_Up_CS_filter$combinedScoreMatrix)), ]

# Extract the phosphosite part from the row names of filtered_score_matrix
phosphosites <- sapply(strsplit(rownames(filtered_score_matrix), " "), `[`, 2)
# Create a new data frame to store the results
kinase_sites_df <- data.frame()

# Create phosphosite_df from filtered_score_matrix
phosphosite_df <- data.frame(
  Gene = sapply(strsplit(rownames(filtered_score_matrix), " "), `[`, 1),
  Phosphosite = phosphosites,
  
  stringsAsFactors = FALSE
)
# Loop over the kinaseSubstrates list
for (kinase in names(Signalomes_results_ENDO_CS_Up_filtered$kinaseSubstrates)) {
  # Get the substrates for the current kinase
  kinase_substrates <- Signalomes_results_ENDO_CS_Up_filtered$kinaseSubstrates[[kinase]]
  
  # Loop over the rows of phosphosite_df
  for (i in 1:nrow(phosphosite_df)) {
    # Create a string that matches the format in kinaseSubstrates
    phosphosite_string <- paste0(phosphosite_df$Gene[i], " ", phosphosite_df$Phosphosite[i], " ")
    
    # Check if phosphosite_string is in kinase_substrates
    if (phosphosite_string %in% kinase_substrates) {
      # If it is, add a row to kinase_sites_df
      kinase_sites_df <- rbind(kinase_sites_df, data.frame(Kinase=kinase, Gene=phosphosite_df$Gene[i], Phosphosite=phosphosite_df$Phosphosite[i]))
    }
  }
}
# Get the cluster information
clusters <- Signalomes_results_ENDO_CS_Up_filtered$proteinModules

# Split the names of clusters into Gene and Phosphosite
cluster_info <- data.frame(Gene = sapply(strsplit(names(clusters), " "), `[`, 1), 
                           Phosphosite = sapply(strsplit(names(clusters), " "), `[`, 2), 
                           Cluster = as.character(clusters))
# Function to annotate phosphosites with clusters
annotate_with_clusters <- function(row) {
  phosphosites <- strsplit(row$Phosphosite, ", ")[[1]]
  clusters <- sapply(phosphosites, function(phosphosite) {
    cluster <- cluster_info$Cluster[cluster_info$Gene == row$Gene & cluster_info$Phosphosite == phosphosite]
    if (length(cluster) > 0) {
      return(paste0(phosphosite, " (Cluster ", cluster, ")"))
    } else {
      return(phosphosite)
    }
  })
  return(paste(clusters, collapse = ", "))
}

# Merge phosphosite_df with cluster_info
phosphosite_df_with_clusters <- merge(phosphosite_df, cluster_info, by = c("Gene", "Phosphosite"), all.x = TRUE)

# Merge phosphosite_df_with_clusters with kinase_sites_df
phosphosite_df_with_clusters <- merge(phosphosite_df_with_clusters, kinase_sites_df, by = c("Gene", "Phosphosite"), all.x = TRUE)

# Annotate phosphosites with clusters and kinases
phosphosite_df_with_clusters$Phosphosite <- paste(phosphosite_df_with_clusters$Phosphosite, " (", phosphosite_df_with_clusters$Cluster, "|", phosphosite_df_with_clusters$Kinase, ")", sep="")

# Aggregate phosphosites for each gene
aggregated_phosphosites <- aggregate(Phosphosite ~ Gene, phosphosite_df_with_clusters, paste, collapse = ", ")

# Merge gene_frequency_df with aggregated_phosphosites
merged_df_with_clusters <- merge(gene_frequency_df, aggregated_phosphosites, by = "Gene")
merged_df_with_clusters <- rename(merged_df_with_clusters, Protein = Gene, `Phosphosite (Cluster|Kinase)` = Phosphosite)
# Sort merged_df by Frequency and select the top 10 rows
top_10_merged_df <- head(merged_df_with_clusters[order(-merged_df_with_clusters$Frequency), ], 10)

# Select the first 10 phosphosites in the `Phosphosite (Cluster|Kinase)` column
top_10_merged_df <- top_10_merged_df %>%
  mutate(`Phosphosite (Cluster|Kinase)` = sapply(strsplit(`Phosphosite (Cluster|Kinase)`, ", "), function(x) str_c(head(x, 10), collapse = ", ")))

# Convert top_10_merged_df to a LaTeX table
latex_table <- xtable(top_10_merged_df)

# Print the LaTeX table
print(latex_table, type = "latex", include.rownames=FALSE)


# Create a copy of the data frame for export
export <- L6.matrices_Up_CS$combinedScoreMatrix
export <- as.data.frame(export)

# Add row names as a new column
export <- cbind(Phosphosite = rownames(export), export)
reshaped_data_for_export <- rename(reshaped_data, Kinases = Var1, Sample = Var2)
merged_df_with_clusters_export = rename(merged_df_with_clusters, 'Frequency among Clusters and GO Terms'= Frequency)
# Install the writexl package if not already installed
if (!require(writexl)) {
  install.packages("writexl")
}


kinases_dataframe <- aggregated_phosphosites %>%
  mutate(PhosphositeList = strsplit(Phosphosite, ", ")) %>%
  unnest(PhosphositeList) %>%
  separate(PhosphositeList, into = c("Phosphosite", "ClusterKinase"), sep = " \\(") %>%
  separate(ClusterKinase, into = c("Cluster", "Kinase"), sep = "\\|") %>%
  mutate(Kinase = sub("\\)", "", Kinase))


kinases_dataframe <- kinases_dataframe %>%
  group_by(Gene, Kinase) %>%
  mutate(AffectedSites = paste(Phosphosite, collapse = ", ")) %>%
  ungroup()

kinase_counts_per_protein <- kinases_dataframe %>%
  group_by(Gene, Kinase, AffectedSites) %>%
  summarise(KinaseCount = n())

# Sort the dataframe by the KinaseCount column in descending order
sorted_kinase_counts_per_protein <- kinase_counts_per_protein %>%
  arrange(desc(KinaseCount))

# Convert the dataframe to a LaTeX table
latex_table_kinases <- xtable(top_5_kinase_counts_per_protein <- kinase_counts_per_protein %>%
                                arrange(desc(KinaseCount)) %>%
                                head(5))
print(latex_table_kinases, type = "latex", include.rownames=FALSE)
