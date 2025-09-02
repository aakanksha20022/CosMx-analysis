library(DESeq2)     
library(limma)      
library(edgeR)      
library(pheatmap)   
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(sva)

epithelial_mtx= read.csv("Epithelial_pseudobulk.csv", row.names=1, header=T)

#SAMPLE SHEET
col_data <- data.frame(
  sample = colnames(epithelial_mtx),
  sample_group = c("positive", "negative", "positive", "negative", "positive", "negative", "positive",
                   "negative", "positive", "negative", "positive", "negative", "positive", "negative",
                   "positive", "negative", "positive", "negative", "positive"),
  patient_id = c("INC0027", "INC0032", "INC0032", "INC0046", "INC0046", "INC0074", "INC0074", 
                 "INC0223", "INC0223", "INC0315", "INC0315", "INC0510", "INC0510", "INC1662", 
                 "INC1662", "INC2655", "INC2655", "INC2842", "INC2842")
)

# Replace dots with underscores in sample names
col_data$sample <- gsub("\\.", "_", col_data$sample)

------#EPITHELIAL CELLS--------
epithelial_corrected= process_em("Epithelial_pseudobulk.csv")
epithelial_dds = do_de(epithelial_corrected, col_data)
epithelial_res= get_res(epithelial_dds)
epithelial_de= format_de(epithelial_res)
epithelial_norm_em= format_norm_em(epithelial_dds)
write.table(as.data.frame(epithelial_de), "/home/ac715y/scripts/tsv_output/DE_epithelial.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(epithelial_norm_em), "/home/ac715y/scripts/tsv_output/Epithelial_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)


------#INTESTINAL EPITHELIAL CELLS-------
int_ep_corrected= process_em("Intestinal Epithelial_pseudobulk.csv")
int_ep_dds = do_de(int_ep_corrected, col_data)
int_ep_res= get_res(int_ep_dds)
int_ep_de= format_de(int_ep_res)
int_ep_norm_em= format_norm_em(int_ep_dds)
write.table(as.data.frame(int_ep_de), "/home/ac715y/scripts/tsv_output/DE_int_epithelial.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(int_ep_norm_em), "/home/ac715y/scripts/tsv_output/Intestinal_Epithelial_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

-----#IMMUNE T_NK_CELLS---------
immune_t_nk_corrected= process_em("Immune_T_NK_pseudobulk.csv")
immune_t_nk_dds = do_de(immune_t_nk_corrected, col_data)
immune_t_nk_res= get_res(immune_t_nk_dds)
immune_t_nk_de= format_de(immune_t_nk_res)
immune_t_nk_norm_em= format_norm_em(immune_t_nk_dds)
write.table(as.data.frame(immune_t_nk_de), "/home/ac715y/scripts/tsv_output/DE_Immune_T_NK.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(immune_t_nk_norm_em), "/home/ac715y/scripts/tsv_output/Immune_T_NK_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

----#EPITHELIAL INFLAM CELLS-------
ep_inflam_corrected= process_em("Epithelial_Inflam_pseudobulk.csv")
ep_inflam_dds = do_de(ep_inflam_corrected, col_data)
ep_inflam_res= get_res(ep_inflam_dds)
ep_inflam_de= format_de(ep_inflam_res)
ep_inflam_norm_em= format_norm_em(ep_inflam_dds)
write.table(as.data.frame(ep_inflam_de), "/home/ac715y/scripts/tsv_output/DE_epithelial_inflam.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(ep_inflam_norm_em), "/home/ac715y/scripts/tsv_output/Epithelial_Inflam_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

----#PLASMA_B CELLS-------
plasma_b_corrected= process_em("Plasma_B_pseudobulk.csv")
plasma_b_dds = do_de(plasma_b_corrected, col_data)
plasma_b_res= get_res(plasma_b_dds)
plasma_b_de= format_de(plasma_b_res)
plasma_b_norm_em= format_norm_em(plasma_b_dds)
write.table(as.data.frame(plasma_b_de), "/home/ac715y/scripts/tsv_output/DE_plasma_b.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(plasma_b_norm_em), "/home/ac715y/scripts/tsv_output/Plasma_B_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

-----#SMOOTH MUSCLE CELLS-----
smc_fibroblasts_corrected= process_em("smooth_muscle_fibroblasts_pseudobulk.csv")
smc_fibroblasts_dds = do_de(smc_fibroblasts_corrected, col_data)
smc_fibroblasts_res= get_res(smc_fibroblasts_dds)
smc_fibroblasts_de= format_de(smc_fibroblasts_res)
smc_fibroblasts_norm_em= format_norm_em(smc_fibroblasts_dds)
write.table(as.data.frame(smc_fibroblasts_de), "/home/ac715y/scripts/tsv_output/DE_smooth_muscle_fibroblasts.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(smc_fibroblasts_norm_em), "/home/ac715y/scripts/tsv_output/smooth_muscle_fibroblasts_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

------##MYOBLASTS FIBROBLASTS CELLS------
myo_fibro_stromal_corrected= process_em("Myo_fibroblasts & Stromal_pseudobulk.csv")
myo_fibro_stromal_dds = do_de(myo_fibro_stromal_corrected, col_data)
myo_fibro_stromal_res= get_res(myo_fibro_stromal_dds)
myo_fibro_stromal_de= format_de(myo_fibro_stromal_res)
myo_fibro_stromal_norm_em= format_norm_em(myo_fibro_stromal_dds)
write.table(as.data.frame(myo_fibro_stromal_de), "/home/ac715y/scripts/tsv_output/DE_myo_fibro_stromal.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(myo_fibro_stromal_norm_em), "/home/ac715y/scripts/tsv_output/Myo_Fibro_stromal_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

------##MACROPHAGE DENDRITIC MONOCYTES CELLS--------
macro_mono_den_corrected= process_em("Macrophage_Dendritic_Monocytes_pseudobulk.csv")
macro_mono_den_dds = do_de(macro_mono_den_corrected, col_data)
macro_mono_den_res= get_res(macro_mono_den_dds)
macro_mono_den_de= format_de(macro_mono_den_res)
macro_mono_den_norm_em= format_norm_em(macro_mono_den_dds)
write.table(as.data.frame(macro_mono_den_de), "/home/ac715y/scripts/tsv_output/DE_macro_mono_den.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(macro_mono_den_norm_em), "/home/ac715y/scripts/tsv_output/Macrophage_Dendritic_Monocytes_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

------#B CELLS--------
b_corrected= process_em("B_pseudobulk.csv")
b_dds = do_de(b_corrected, col_data)
b_res= get_res(b_dds)
b_de= format_de(b_res)
b_norm_em= format_norm_em(b_dds)
write.table(as.data.frame(b_de), "/home/ac715y/scripts/tsv_output/DE_B.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(b_norm_em), "/home/ac715y/scripts/tsv_output/B_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

------#ENDOTHELIAL CELLS--------
endo_corrected= process_em("Endothelial_pseudobulk.csv")
endo_dds = do_de(endo_corrected, col_data)
endo_res= get_res(endo_dds)
endo_de= format_de(endo_res)
endo_norm_em= format_norm_em(endo_dds)
write.table(as.data.frame(endo_de), "/home/ac715y/scripts/tsv_output/DE_endo.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(endo_norm_em), "/home/ac715y/scripts/tsv_output/Endothelial_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

------#MAST CELLS--------
mast_corrected= process_em("Mast_pseudobulk.csv")
mast_dds = do_de(mast_corrected, col_data)
mast_res= get_res(mast_dds)
mast_de= format_de(mast_res)
mast_norm_em= format_norm_em(mast_dds)
write.table(as.data.frame(mast_de), "/home/ac715y/scripts/tsv_output/DE_mast.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(mast_norm_em), "/home/ac715y/scripts/tsv_output/Mast_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

------#NEURON GLIAL CELLS--------
neuron_glial_corrected= process_em("Neuron_Glial_pseudobulk.csv")
neuron_glial_dds = do_de(neuron_glial_corrected, col_data)
neuron_glial_res= get_res(neuron_glial_dds)
neuron_glial_de= format_de(neuron_glial_res)
neuron_glial_norm_em= format_norm_em(neuron_glial_dds)
write.table(as.data.frame(neuron_glial_de), "/home/ac715y/scripts/tsv_output/DE_neuron_glial.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(neuron_glial_norm_em), "/home/ac715y/scripts/tsv_output/Neuron_Glial_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

------#ALL CELLS--------
all_corrected= process_em("All_pseudobulk.csv")
all_dds = do_de(all_corrected, col_data)
all_res= get_res(all_dds)
all_de= format_de(all_res)
all_norm_em= format_norm_em(all_dds)
write.table(as.data.frame(all_de), "/home/ac715y/scripts/tsv_output/DE_all.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(all_norm_em), "/home/ac715y/scripts/tsv_output/All_pseudobulk.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)


