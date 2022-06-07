## EpiCompare analysis of CUT&Tag, CUT&Run, ENCODE and TIP-seq

## Load peakfiles in BED format
# CUT&Tag
active_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/ActiveMotif_SEACR.bed", as="GRanges")
active_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/ActiveMotif_MACS2.bed", as="GRanges")
abcam_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/Abcam-ab4729_SEACR.bed", as="GRanges")
abcam_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/Abcam-ab4729_MACS2.bed", as="GRanges")
kayaOkur_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/H3K27ac_SRR8383507_SEACR.bed", as="GRanges")
kayaOkur_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/H3K27ac_SRR8383507_MACS2.bed", as="GRanges")

# CUT&Run
meers_ac_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/H3K27ac_SRR8581604_CnR_MACS2.bed", as="GRanges")
meers_ac_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/H3K27ac_SRR8581604_CnR_SEACR.bed", as="GRanges")
meers_me_macs <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/H3K27me3_SRR9073702_CnR_MACS2.bed", as="GRanges")
meers_me_seacr <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/CnT/H3K27me3_SRR9073702_CnR_SEACR.bed", as="GRanges")

# ENCODE
encode_H3K27me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/ENCODE/ENCODE_H3K27me3.bed", as="GRanges")
encode_H3K27ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/ENCODE/ENCODE_H3K27ac.bed", as="GRanges")

# TIP-seq
tip_1_05_01 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_1_05_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_2_28_01 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_2_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_3_28_01 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_3_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_4_03_02 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_4_03_02_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_4_28_01 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_4_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_5_28_01 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_5_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_6_28_01 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_6_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")

# create peaklist
peaklist <- list(active_seacr, active_macs, abcam_seacr, abcam_macs, kayaOkur_seacr, kayaOkur_macs,
                 meers_ac_seacr, meers_ac_macs, meers_me_seacr, meers_me_macs,
                 encode_H3K27me, encode_H3K27ac,
                 tip_1_05_01, tip_4_03_02, tip_2_28_01, tip_3_28_01, tip_4_28_01, tip_5_28_01, tip_6_28_01)
# set names
my_label <- c("H3K27ac_CnT_ActiveMotif_SEACR", "H3K27ac_CnT_ActiveMotif_MACS2",
              "H3K27ac_CnT_Abcamab4729_SEACR", "H3K27ac_CnT_Abcamab4729_MACS2",
              "H3K27ac_CnT_KayaOkur_SEACR", "H3K27ac_CnT_KayaOkur_MACS2",
              "H3K27ac_CnR_Meers_SEACR", "H3K27ac_CnR_Meers_MACS2",
              "H3K27me3_CnR_Meers_SEACR", "H3K27me3_CnR_Meers_MACS2",
              "H3K27me3_ENCODE", "H3K27ac_ENCODE",
              "H3K27ac_TIP_Abcam.phase_1_05_jan_2022.S_1_R1",
              "H3K27ac_TIP_Abcam.phase_2_03_feb_2022.S_4_R1",
              "H3K27ac_TIP_Diagenode.phase_2_28_jan_2022.S_2_R1",
              "H3K27me3_TIP_Diagenode.phase_2_28_jan_2022.S_3_R1",
              "H3K27ac_TIP_Diagenode.phase_2_28_jan_2022.S_4_R1",
              "H3K27me3_TIP_Diagenode.phase_2_28_jan_2022.S_5_R1",
              "H3K27ac_TIP_Diagenode.phase_2_28_jan_2022.S_6_R1")
names(peaklist) <- my_label

# reference file
reference <- list("ENCODE_H3K27ac" = encode_H3K27ac)
# hg19 blacklist
data("hg19_blacklist")

# Epicompare
EpiCompare(peakfiles = peaklist,
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



