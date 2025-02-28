library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(stringr)
# library(GenomicRanges)
# library(ggrepel)
# library(tibble)

setwd("~/Data/Lab2/")

get_ce_from_ipsa <- function(df){

  t <- data.frame(matrix(nrow=0, ncol=12))
  
  for(i in 1:nrow(df)) {
    row <- df[i,]
    dfc <- df %>% filter(chrom == row$chrom & 
                           start >= row$start & end <= row$end)
    dff <- dfc %>% filter(chrom == row$chrom & 
                           start == row$start & end < row$end) %>% 
      distinct(chrom, start, end, strand, .keep_all=TRUE)
    # dff
    if (dim(dff)[1] > 0) {
      # print(i)
      A = row$start
      D = row$end
      E = row$a_chx
      Ec = row$a_ctrl
      print(paste(i, A, D))
      for (j in dim(dff)[1]) {
        incl1 <- dff[j, ]
        B = incl1$end
        I1 = incl1$a_chx
        I1c = incl1$a_ctrl
        print(paste(i, A, B, D))
        dfff <- dfc %>% filter(chrom == row$chrom & start > B & end == D)%>% 
          distinct(chrom, start, end, strand, .keep_all=TRUE)
        if (dim(dfff)[1] > 0) {
          for (k in dim(dfff)[1]) {
            incl2 <- dfff[k, ]
            C = incl2$start
            I2 = incl2$a_chx
            I2c = incl2$a_ctrl
            print(paste(i, A, B, C,  D))
            new_row <- c(row$chrom, A, B, C, D, 
                         Ec, E, I1c, I1, I2c, I2, row$strand)
            t <- rbind(t, new_row)
          }
        }
      }
    }
  }
  colnames(t) <- c('chrom', 'A', 'B', 'C', 'D',
                   'e_ctrl', 'e_kd', 'i1_ctrl', 'i1_kd', 
                   'i2_ctrl', 'i2_kd', 'strand')
  t <- t %>% distinct(.keep_all=TRUE)
  return(t)
}

load('cycloheximide/df_all2_a.rda')
df <- df_all2 %>% 
  dplyr::rename(a_chx = TOTAL_COUNT_chm, a_ctrl = TOTAL_COUNT_ctrl)
t <- get_ce_from_ipsa(df)
View(df)
save(t, file='t_H_CHX.rda')
write.table(t, 'chx_cassette_exons_human.csv', sep=' ', 
            row.names=FALSE, quote= FALSE) 


clk <- read.table('clk_coordinates_plus_20k.bed')
clk[1, 'V2']
load('cycloheximide/df4_mouse_FULL.rda')

df <- df4_mouse %>% arrange(chrom, start, end) %>% 
  mutate(start = as.integer(start), 
         end = as.integer(end)) %>%
  filter((chrom == clk[1, 'V1'] & start > clk[1, 'V2'] & end < clk[1, 'V3']) |
           (chrom == clk[2, 'V1'] & start > clk[2, 'V2'] & end < clk[2, 'V3']) |
           (chrom == clk[3, 'V1'] & start > clk[3, 'V2'] & end < clk[3, 'V3']) |
           (chrom == clk[4, 'V1'] & start > clk[4, 'V2'] & end < clk[4, 'V3']))
View(df)
rm(df4_mouse)
t <- get_ce_from_ipsa(df)
save(t, file='t_M_CHX.rda')
write.table(t, 'chx_cassette_exons_mouse.csv', sep=' ', 
            row.names=FALSE, quote= FALSE) 
Sys.time()
View(df)



