# Analysis of CnT, ENCODE and TIP-seq using EpiCompare

# CnT peakfiles
active_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peak_files/ActiveMotif_SEACR.bed", as="GRanges")
active_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peak_files/ActiveMotif_MACS2.bed", as="GRanges")
abcam_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peak_files/Abcam-ab4729_SEACR.bed", as="GRanges")
abcam_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peak_files/Abcam-ab4729_MACS2.bed", as="GRanges")
CnT_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peak_files/H3K27ac_SRR8383507_SEACR.bed", as="GRanges")
CnT_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peak_files/H3K27ac_SRR8383507_MACS2.bed", as="GRanges")

# ENCODE peakfiles
H3K27me3_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peak_files/H3K27me3_SEACR.bed", as="GRanges")
H3K27me3_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peak_files/H3K27me3_MACS2.bed", as="GRanges")
encode_H3K27ac_full <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/ENCODE_H3K27ac.bed", as="GRanges")

# TIP-seq peakfiles
tip_R1 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/tip_peakfiles/S_1_R1.peaks.bed.stringent.bed", as="GRanges")
tip_R2 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/tip_peakfiles/S_1_R2.peaks.bed.stringent.bed", as="GRanges")
tip <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/tip_peakfiles/S_4_R1.peaks.bed.stringent.bed", as="GRanges")


# peaklist
peaklist <- list(active_macs, active_seacr, abcam_macs, abcam_seacr, CnT_macs, CnT_seacr, H3K27me3_macs, H3K27me3_seacr, encode_H3K27ac_full, tip_R1, tip_R2, tip)
names(peaklist) <- c("H3K27ac_CnT_ActiveMotif_MACS2", "H3K27ac_CnT_ActiveMotif_SEACR",
                     "H3K27ac_CnT_Abcam-ab4729_MACS2", "H3K27ac_CnT_Abcam-ab4729_SEACR",
                     "H3K27ac_CnT_SRR8383507_MACS2", "H3K27ac_CnT_SRR8383507_SEACR",
                     "H3K27me3_ENCODE_MACS2.bed", "H3K27me3_ENCODE_SEACR.bed",
                     "H3K27ac_ENCODE",
                     "H3K27ac_TIP_Abcam.phase_1_05_jan_2022.S_1_R1", "H3K27ac_TIP_Abcam.phase_1_05_jan_2022.S_1_R2",
                     "H3K27ac_TIP_Abcam.phase_2_03_feb_2022.S_4_R1")

# reference peakfile: ENCODE H3K27ac (data accession: ENCFF044JNJ)
reference <- encode_H3K27ac_full
names(reference) <- "ENCODE_H3K27ac"
# blacklist
data("hg19_blacklist")

# EpiCompare
EpiCompare::EpiCompare(peakfiles = peaklist,
           blacklist = hg19_blacklist,
           reference = reference,
           stat_plot = TRUE,
           save_output = FALSE,
           chrmHMM_plot = TRUE,
           chrmHMM_annotation = "K562",
           chipseeker_plot = TRUE,
           enrichment_plot = TRUE,
           interact = TRUE,
           output_dir = "/Users/serachoi/Documents/EpiCompare")


