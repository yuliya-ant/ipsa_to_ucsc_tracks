library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(ggrepel)
library(tibble)

setwd("~/Data/Lab2/")

# path1 = 'srsf3_kd_cassete_exons.tsv'
# path2 = 'rbfox2_kd_cassete_exons.tsv'
# path_out <- 'srsf3_rbfox2_kd_cassette_exons.bed'


path1 = 'chx_cassette_exons_mouse.csv'
path_out = 'chx_cassette_exons_mouse.bed'

# path1 = 'chx_cassette_exons_human.csv'
# path_out = 'chx_cassette_exons_human.bed'


# UCSC track columns:
# chrom chromStart chromEnd name score strand 
# thickStart thickEnd itemRgb blockCount blockSizes blockStarts
#
# https://genome.ucsc.edu/goldenPath/help/customTrack.html - here find 'BED custom track with multiple blocks'
# 

# psi = 2e / (2e + i1 + i2)

  
make_ce_bed = function(ce, colorx) {
  ce_bed = ce %>% 
    mutate(
         psi_ctrl =  2*e_ctrl / (2*e_ctrl + i1_ctrl + i2_ctrl), 
         psi_kd =  2*e_kd / (2*e_kd + i1_kd + i2_kd), 
         score = round(psi_ctrl - psi_kd, 2), 
         
         chromStart = A - 2, chromEnd = D + 2, 
         thickStart = A - 2, thickEnd = D + 2, 
         blockCount = 3, itemRgb = colorx,
         blockSizes = paste(2, C - B, 2, sep=','), 
         blockStarts = paste(0, B - A, D - A + 2, sep=','), 
         name = paste0(
           'dpsi=', round(score, 2), ',psi_kd=', round(psi_kd, 2), 
           ',psi_ctrl=', round(psi_ctrl, 2), 
           ',excl_kd=', e_kd, ',excl_ctrl=', e_ctrl, 
           ',i1_kd=', i1_kd, ',i1_ctrl=', i1_ctrl, 
           ',i2_kd=', i2_kd, ',i2_ctrl=', i2_ctrl
         )
         ) %>% 
    dplyr::select(chrom, chromStart, chromEnd, name, score,
                strand, thickStart, thickEnd, itemRgb, 
                blockCount, blockSizes, blockStarts)
  return(ce_bed)
}

color1 = "103,137,33"
color2 = "71,48,120"

ce <- read.table(path1, header=1)
ce_bed1 = make_ce_bed(ce, color1) 
ce_bed1$score[is.na(ce_bed1$score)] = 0

ce <- read.table(path2, header=1)
ce_bed2 = make_ce_bed(ce, color2)

View(ce_bed1)

header <- 'track name="chx_cassette" description="chx_cassette_exons, PSI = 2e / (2e + i1 + i2)" visibility=2 color=103,137,33 itemRgb="on" useScore=1'

header <- 'track name="srsf3_kd_cassette" description="srsf3_kd_cassette_exons, PSI = 2e / (2e + i1 + i2)" visibility=2 color=103,137,33 itemRgb="on" useScore=1'
cat(header, '\n', file=path_out, append = FALSE)
write.table(ce_bed1,
            file = path_out, append = TRUE, sep=' ', 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

header <- 'track name="rbfox2_kd_cassette" description="rbfox2_kd_cassette_exons, PSI = 2e / (2e + i1 + i2)" visibility=2 color=71,48,120 itemRgb="on" useScore=1'
cat(header, '\n', file=path_out, append = TRUE)
write.table(ce_bed2,
            file = path_out, append = TRUE, sep=' ', 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

View(ce_bed1)
