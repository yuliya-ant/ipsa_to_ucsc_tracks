library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
# library(GenomicRanges)
# library(ggrepel)
# library(tibble)

getwd()
setwd("~/Data/Lab2/")

# ---
# 1. Define paths
# 2. Make metadata table
# 3. Merge IPSA output (junctions and sites) to one table
# ---


# 1. Define paths

path <- '~/Data/Lab2/cycloheximide/pyIPSA_Masha/'


# 2. Make metadata table
#    columns: files, group

# HUMAN CHX METADATA

# files = list.files(paste0(path, 'J4/'))
# chx_files = gsub('.J4.gz', '', files[str_detect(files, "CHX_WT") &
#                     str_detect(files, "J4") ])
# control_files = gsub('.J4.gz', '', files[str_detect(files, "SS_WT") &
#                      str_detect(files, "J4") ])
# meta = rbind(
#   data.frame(files = chx_files, group = 'chx'),
#   data.frame(files = control_files, group = 'ctrl'))
# meta

# MOUSE CHX METADATA

# meta <- data.frame(
#   sample=c("SRX25290400", "SRX25290401", "SRX25290402", 
#            "SRX25290403", "SRX25290404", "SRX25290405"), 
#   group = c("ctrl", "ctrl", 
#             "chx_3h", "chx_3h", 
#             "chx_6h", "chx_6h"), 
#   files=c(
#           "SRX25290400", 
#           "SRX25290401", 
#           "SRX25290402", 
#           "SRX25290403", 
#           "SRX25290404", 
#           "SRX25290405")
#   )

# MOUSE SRSF3 KD METADATA

# meta1 <- read.table('srsf2_kd_mm_ctrl.csv', header = 1)
# meta2 <- read.table('srsf2_kd_mm_ko.csv', header = 1)
# meta1$group = 'ctrl'
# meta2$group = 'kd'
# meta = rbind(meta1, meta2)
# colnames(meta) = c('files', 'group')

# MOUSE RBFOX2 KD METADATA

# meta <- DataFrame(files = paste0('SRR', 1165129:1165136), 
#                   group = c(rep('ctrl', 4), rep('kd', 4)), 
#                   group2 = c(rep('ctrl_luc', 2), 
#                              rep('ctrl_gfp', 2), rep('kd', 4)))
# meta <- as.data.frame(meta)


# 3. Merge IPSA output (junctions and sites) to one table
# 3.1 Junctions: J6 (filtered by entropy), J4 (not filtered)


j4_path <- paste0(path, 'J4/')
j6_path <- paste0(path, 'J6/')

colnames_j <- c('junction_id',	'total_count',	
                'staggered_count',	'entropy',	'annotation_status', 
                'splice_site')


# meta$files[1]
# # [1] "CHX_WT_16_trimmed_sorted_merged"
j4or6 = 'J4'
make_junction_table <- function(meta, path, j4or6){
  
  if (j4or6=='J4') {
    j_path = paste0(path, 'J4/')
  } else {
    j_path = paste0(path, 'J6/')
  }
  
  junctions = data.frame(matrix(nrow=0, ncol = 8))
  colnames(junctions) = colnames_j
  
  for (f0 in meta$files) {
    # HERE UPDATE FILENAMES IF NEEDED, for KDs:
    # f <- paste0(f0, 'Aligned.sortedByCoord.out.J4.gz')
    f <- paste0(f0, '.', j4or6, '.gz')
    ff <- read.csv(paste0(j_path, f), sep='\t', header = FALSE)
    colnames(ff) <- colnames_j
    ff <- ff %>% mutate(file = f0, 
                        chx = meta %>% as.data.frame() %>%
                          dplyr::filter(files == f0) %>% 
                          .$group)
    junctions <- rbind(junctions, ff)
  }
  # View(junctions)

  junctions <- junctions %>% select(!staggered_count & !entropy) %>% 
    group_by(junction_id) %>% 
    summarise(a_chx = sum(total_count[chx %in% c('kd')]), 
              a_ctrl = sum(total_count[chx == 'ctrl']), 
              splice_site = paste(unique(splice_site), collapse = '_'),
              annotation_status = max(annotation_status)
    ) %>% 
    arrange(junction_id) %>% ungroup() %>% 
    separate(junction_id, into=c('chrom', 'start', 'end', 'strand'), 
             remove=FALSE, sep='_') 

# j_test <- junctions %>% filter(junction_id %in% j)
  junctions <- junctions %>% 
    group_by(chrom, start, end) %>% 
    summarise(junction_id = junction_id[1],
              strand = strand[1],
              a_chx = sum(a_chx), 
              a_ctrl = sum(a_ctrl), 
              splice_site = splice_site[1],
              annotation_status = max(annotation_status)
    ) %>% 
    arrange(junction_id) 
  
  return(junctions)
}

junctions <- make_junction_table(meta, 'J4')
junctions <- make_junction_table(meta, 'J6')

View(junctions)
# save(junctions, file='junctions_MOUSE_SRSF2_KD_FULL_J4.rda')
save(junctions, file='junctions_MOUSE_RBFOX2_KD_FULL_J4.rda')

table(junctions$annotation_status)



# s2_path <- '~/Data/Lab2/pyIpsa_srsf2_kd/S2/'
s2_path <- paste0(path, 'S2/')

# files_s <- list.files(s2_path) %>%  as.data.frame()
# colnames(files_s) <- c('files')

colnames_s <- c('site_id',	'total_count',	'staggered_count',	'entropy')
sites = data.frame(matrix(nrow=0, ncol = 6))
colnames(sites) = colnames_s
# save(meta, file='meta_MOUSE_CHX_J4S2.rda')
# load('meta_MOUSE_CHX_J4S2.rda')

for (f0 in meta$files) {
  f <- paste0(f0, 'Aligned.sortedByCoord.out.S2.gz')
  ff <- read.csv(paste0(s2_path, f), sep='\t', header = FALSE)
  colnames(ff) <- colnames_s
  ff <- ff %>% mutate(file = f, 
                      chx = meta %>% filter(files == f0) %>% .$group)
  sites <- rbind(sites, ff)
}

View(sites)
table(sites$chx)


sites <- sites %>% select(!entropy & !staggered_count) %>% 
  group_by(site_id) %>% 
  summarise(c_chx = sum(total_count[chx %in% c('kd')]), 
            c_ctrl = sum(total_count[chx == 'ctrl'])) %>% 
  arrange(site_id) %>% 
  separate(site_id, into=c('chrom', 'site', 'strand'), 
           remove=FALSE, sep='_')

sites <- sites %>%  group_by(chrom, site) %>% 
  summarise(site_id = site_id[1],
            strand = strand[1],
            c_chx = sum(c_chx), 
            c_ctrl = sum(c_ctrl)) %>% 
  arrange(site_id) 

View(sites)

# save(sites, file='sites_MOUSE_SRSF2_KD_FULL_S2.rda')
save(sites, file='sites_MOUSE_RBFOX2_KD_FULL_S2.rda')

df4_hela <- junctions %>% 
  left_join(sites, by=c('chrom', 'start'='site'), 
            relationship = "many-to-many") %>% 
  left_join(sites, by=c('chrom', 'end'='site'), 
            relationship = "many-to-many", 
            suffix = c('_start', '_end')) %>% 
  dplyr::rename(strand = strand.x, 
                strand_start = strand.y, 
                strand_end = strand) %>% 
  rowwise() %>% 
  mutate(sites_counted = 2-sum(is.na(c(site_id_start, site_id_end)))) 

View(df4_hela)

df4_hela <- df4_hela %>% mutate(across(starts_with('c_'), as.numeric)) %>% 
  mutate(across(starts_with('c_'), ~replace_na(., 0))) %>% 
  mutate(c_ctrl = c_ctrl_start + c_ctrl_end, 
         c_chx = c_chx_start + c_chx_end) %>% 
  mutate(cosi_ctrl = round(a_ctrl / (a_ctrl + 0.5*c_ctrl), 2),
         cosi_chx = round(a_chx / (a_chx + 0.5*c_chx), 2))

df4_hela <- df4_hela %>% mutate(cosi_ctrl = replace_na(cosi_ctrl, 1), 
                                cosi_chx = replace_na(cosi_chx, 1), 
                                delta_cosi = cosi_chx - cosi_ctrl, 
                                log2fc_local = log2((a_chx + c_chx) / (a_ctrl + c_ctrl)))

View(df4_mouse)
df4_mouse <- df4_hela
# save(df4_mouse, file="df4_mouse_FULL_SRSF2_KD.rda")
save(df4_mouse, file="df4_mouse_FULL_RBFOX2_KD.rda")
load('df4_mouse_FULL_SRSF2_KD.rda')
View(df4_mouse)



# BED of ALL mouse introns

df4_mouse$ten_reads = ifelse(
  (df4_mouse$a_ctrl + 0.5*df4_mouse$c_ctrl >= 10) | 
    (df4_mouse$a_chx + 0.5*df4_mouse$c_chx >= 10), 
  TRUE, FALSE)

str_which_exp = 'psi_RBFOX2_KD='
str_which_exp = 'psi_SRSF3_KD='

df4_mouse <- df4_mouse[1:100, ]

df4_mouse$name=paste(df4_mouse$strand, ', dpsi=', round(df4_mouse$delta_cosi,2), ', ',
                     'psi_ctrl=', df4_mouse$a_ctrl, '/', 
                     round(df4_mouse$a_ctrl + 0.5*df4_mouse$c_ctrl, 2), '=', 
                     round(df4_mouse$cosi_ctrl,2), ', ',
                     str_which_exp, df4_mouse$a_chx, '/', 
                     round(df4_mouse$a_chx + 0.5*df4_mouse$c_chx, 2), '=', 
                     round(df4_mouse$cosi_chx,2), 
                     sep='')

# chrom chromStart chromEnd name score strand 
# thickStart thickEnd itemRgb blockCount blockSizes blockStarts

df4_mouse_bed <- df4_mouse %>% 
  # rownames_to_column('name') %>% 
  separate(junction_id, 
           into=c('chrom', 'start', 'end', 'strand'), 
           sep='_') %>% 
  mutate(score = round(delta_cosi,2)) %>%
  dplyr::select(chrom:end, name, score, strand, ten_reads) %>% 
  distinct(chrom, start, end, strand, .keep_all = TRUE)
View(df4_mouse_bed)

border_width = 10

df4_mouse_bed2 <- df4_mouse %>% 
  separate(junction_id, 
           into=c('chrom', 'chromStart', 'chromEnd', 'strand'), 
           sep='_') %>% 
  mutate(score = round(delta_cosi,2)) %>%
  dplyr::select(chrom:chromEnd, name, score, strand, ten_reads) %>% 
  mutate(chromStart = as.numeric(chromStart) - border_width, 
         chromEnd = as.numeric(chromEnd) + border_width, 
         thickStart = chromStart, thickEnd = chromEnd, itemRgb = 0, 
         blockSizes = paste(border_width, border_width, sep=','), 
         blockCount = 2) %>% rowwise() %>% 
  mutate(
    blockStarts = paste(0, chromEnd - border_width - chromStart, 
                        sep = ',')) %>%
  distinct(chrom, chromStart, chromEnd, strand, .keep_all = TRUE)
View(df4_mouse_bed2)

# test
file_name <- 'mouse_srsf3_kd_all_introns2_TEST.bed'
header <- 'track name="mouTEST" description="m" visibility=2 color=154,205,50 useScore=1 type=bedDetail'
header <- 'track name="BED track" description="BED format custom track example" visibility=2 color=0,128,0 useScore=1 type=bedDetail \n#chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts'

cat(header, '\n', file=file_name, 
    append = FALSE)
write.table(df4_mouse_bed2 %>% 
              dplyr::filter(abs(score) >= 0.05 & ten_reads) %>% 
              dplyr::select(chrom, chromStart, chromEnd, name, score,
                            strand, thickStart, thickEnd, itemRgb, 
                            blockCount, blockSizes, blockStarts
              ),
            file = file_name, append = TRUE,
            sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)



file_name <- 'mouse_srsf3_kd_all_introns2.bed'
file_name = 'mouse_rbfox2_kd_all_introns.bed'
header <- 'track name="mouse_srsf3_kd_introns" description="mouse_srsf3_kd_introns, |dpsi| >= 0.05 & denominator in kd or ctrl >= 10" visibility=3 color=154,205,50 type=bedDetail'
header <- 'track name="mouse_rbfox2_kd_introns" description="mouse_rbfox2_kd_introns, |dpsi| >= 0.05 & denominator in kd or ctrl >= 10 (PSI = split_reads / (split_reads + 0.5*continuous_reads_in_site1 + 0.5*continuous_reads_in_site2))" visibility=3 color=137,104,205 type=bedDetail'
cat(header, '\n', file=file_name, 
    append = FALSE)
write.table(df4_mouse_bed2 %>% 
              dplyr::filter(abs(score) >= 0.05 & ten_reads) %>% 
              dplyr::select(!ten_reads),
            file = file_name, append = TRUE,
            sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

header2 <- 'track name="mouse_srsf3_kd_rest_introns" description="mouse_srsf3_kd_rest_introns, |dpsi| < 0.05 or denominators in kd and ctrl < 10" visibility=3 color=214,235,174 type=bedDetail'
header2 <- 'track name="mouse_rbfox2_kd_rest_introns" description="mouse_rbfox2_kd_rest_introns, |dpsi| < 0.05 or denominators in kd and ctrl < 10" visibility=3 color=196,179,230 type=bedDetail'
cat(header2, '\n', file=file_name, 
    append = TRUE)
write.table(df4_mouse_bed2 %>% 
              dplyr::filter(!(abs(score) >= 0.05 & ten_reads)) %>% 
              dplyr::select(!ten_reads),
            file = file_name, append = TRUE,
            sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


load('df4_mouse_FULL_SRSF2_KD.rda')

View(df4_mouse)
introns <- df4_mouse %>%  filter(chrom == 'chr15') %>% 
  dplyr::select(chrom, start, end, strand, a_chx) %>% 
  dplyr::rename(reads = a_chx) %>% 
  group_by(chrom, strand) %>%
  arrange(start, end)
# Соединяем интроны по принципу (предыдущий интрон + следующий интрон = экзон между ними)
exon_inclusion <- introns %>% 
  # filter(chrom == 'chr15') %>% 
  # group_by(chrom, strand) %>%
  # arrange(start, end) %>%
  mutate(next_start = ifelse(start != lead(start),
                             lead(start), 
                             next_end = lead(end),
                             I1 = reads,                  # Риды, поддерживающие первый интрон
                             I2 = lead(reads),            # Риды, поддерживающие второй интрон
                             E = reads[ifelse(!is.na(intersect(match(next_end, introns$end), 
                                                               match(start, introns$start))[1]), 
                                              intersect(match(next_end, introns$end), 
                                                        match(start, introns$start))[1], NA
                             )] # Риды, поддерживающие экзонное исключение (перепрыгивающие экзон)
  )) #%>%
# filter(!is.na(next_start) & !is.na(E)) %>%  # Исключаем последние интроны
# mutate(exon_start = end, 
#        exon_end = next_start,
#        PSI = (I1 + I2) / (2 * E + I1 + I2)) %>% # Считаем PSI
# select(chrom, exon_start, exon_end, strand, I1, I2, E, PSI)

# Просмотр результатов

# Загружаем необходимые библиотеки
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

introns <- introns %>% arrange(chrom, start, end)

# Находим потенциальные экзоны и считаем PSI
exon_inclusion0 <- introns %>%
  group_by(chrom, strand, start) %>%
  count() %>% ungroup()

View(exon_inclusion)
View(introns)




# BED --- 


cosi_thr = 0.05
expr_thr = 0.1
reads_thr = 10

df4_mouse_filt <- df4_mouse %>% filter(delta_cosi >= cosi_thr & 
                                         log2fc_local >= expr_thr & 
                                         # type == 'utr3' & 
                                         a_chx + c_chx >= reads_thr)
dim(df4_mouse_filt)
# save(df4_mouse_filt, file="df4_mouse_filt.rda")








## mm10 annotation to regions





load("/home/yulia_ant/Data/Lab2/cycloheximide/qki_ilf2/df_all_filt_1052.rda")
human_df_filt <- as.data.frame(df_all_filt)
rm(df_all_filt)
human_df_filt_bed <- human_df_filt %>% 
  rownames_to_column('name') %>% 
  separate(junction_id, 
           into=c('chrom', 'start', 'end', 'strand'), 
           sep='_') %>% 
  dplyr::select(chrom:end, name, strand) %>% #View()
  mutate(score = '.') %>%
  dplyr::select(chrom:end, name, score, strand) %>% 
  distinct(chrom, start, end, strand, .keep_all = TRUE)



dim(human_df_filt_bed)

View(human_df_filt_bed)
write.table(human_df_filt_bed, 
            "/home/yulia_ant/Data/Lab2/cycloheximide_mouse/human_filt_1052.bed", 
            sep = '\t', quote = FALSE, 
            row.names = FALSE, col.names = FALSE)



# MOUSE junctions (mapped)

mouse_junctions <- read.table(
  "/home/yulia_ant/Data/Lab2/hg_to_mm_filt_1052_10.bed")
colnames(mouse_junctions) <- c('chrom', 'start', 'end', 
                               'name', 'score', 'strand')
View(mouse_junctions)
mouse_junctions <- mouse_junctions %>% 
  unite(col='junction_id', chrom:end, strand, 
        sep='_', remove = FALSE) 
# mouse_junctions %>% 
#   dplyr::select(chrom:end, strand) %>% 
#   write.table(
#     "~/Data/Lab2/cycloheximide_mouse/mouse_liftover_junctions_83.bed", 
#     quote = FALSE, sep='\t', row.names=FALSE)

# REGIONS

regions <- read.csv('~/Data/Lab2/regions_long_mm10.tsv', sep='\t')
View(regions)
dim(regions)
regions <- drop_na(regions, c(chrom, start, end))
regions <- regions  %>%
  mutate(start1 = ifelse(start <= end, start, end), 
         end1 = ifelse(start <= end, end, start)) %>% 
  mutate(start=NULL, end=NULL) %>% 
  dplyr::rename(c(start=start1, end=end1)) 
dim(regions)

df4_gr <- makeGRangesFromDataFrame(df4_mouse, 
                                   keep.extra.columns=TRUE)
df_gr <- makeGRangesFromDataFrame(mouse_junctions, 
                                  keep.extra.columns=TRUE)
regions_gr <- makeGRangesFromDataFrame(regions, 
                                       keep.extra.columns=TRUE)
View(data.frame(regions_gr))

# overlaps_utr <- findOverlaps(df_gr, regions_utr_gr, type='within')
overlaps <- findOverlaps(df_gr, regions_gr, type='within')
overlaps_df4 <- findOverlaps(df4_gr, regions_gr, type='within')
View(data.frame(overlaps))

overlaps_any <- findOverlaps(df_gr, regions_gr, type='any')
View(data.frame(overlaps_any))

# p_utr <- Pairs(df_gr, regions_utr_gr, hits=overlaps_utr)
p <- Pairs(df_gr, regions_gr, hits=overlaps)
p_df4 <- Pairs(df4_gr, regions_gr, hits=overlaps_df4)
View(data.frame(p))

gene_ids_mouse <- p@second$gene_id %>% unique()


# overlaps_df0 <- data.frame(p)
overlaps_df <- data.frame(junction_id = p@first$junction_id,
                          gene_id=p@second$gene_id,
                          type = p@second$type)
overlaps_df4 <- data.frame(junction_id = p_df4@first$junction_id,
                           gene_id=p_df4@second$gene_id,
                           type = p_df4@second$type)


# # overlaps_utr3_df
# # save(overlaps_utr3_df, file='overlaps_utr3_df.rda')
# save(overlaps_df, file='mouse_regions_overlaps_df.rda')
# # pintersect(p)

load('mouse_regions_overlaps_df.rda')




# View(overlaps_utr3_df)
# dim(overlaps_utr3_df)
# overlaps_utr3_df$junc %>% unique() %>% length()
# 
# 
junc_type <- overlaps_df4 %>% group_by(junction_id)  %>%
  summarize(gene_id = paste(sort(unique(gene_id)), collapse = '_'),
            type=paste(sort(unique(type)), collapse = '_'))


junc_type$multiple_type = grepl('_', junc_type$type, fixed = TRUE)
junc_type$multiple_gene = grepl('_', junc_type$gene_id, fixed = TRUE)


junc_type  %>% filter(multiple_gene == TRUE)  %>% head()

View(junc_type)




load("df4_mouse_filt.rda")
load("df4_mouse_FULL.rda")


df <- df4_mouse %>%
  left_join(junc_type, by='junction_id')
df_filt <- df %>%  
  filter(delta_cosi >= cosi_thr & 
           log2fc_local >= expr_thr & 
           type == 'utr3' &
           a_chx + c_chx >= reads_thr)

View(df_filt)
dim(df_filt)
sum(is.na(df$type))
df %>% filter(! is.na(type)) %>% View()
best_j <- df %>% filter(! is.na(type)) %>% .$junction_id

best_j %in% df4_mouse_filt$junction_id

df <- df %>%
  left_join(gtf, by='gene_id') %>%
  mutate(log2fc_local = log2((a_chx + c_chx) / (a_ctrl + c_ctrl)))
# 
# View(df)
# save(df, file='df_A549_2.rda')
# load('df_A549_2.rda')



### APPROX COORDS

mouse_junctions <- mouse_junctions %>%
  left_join(junc_type, by='junction_id')
View(mouse_junctions1)

tbl <- mouse_junctions %>% group_by(type) %>% count() 
tbl$type[is.na(tbl$type)] <- 'NA'
tbl$name <- paste0(tbl$type, '\nN=', tbl$n)
pie(tbl$n, tbl$name)




df_gr <- makeGRangesFromDataFrame(df4_mouse, 
                                  keep.extra.columns=TRUE)
pred_gr <- makeGRangesFromDataFrame(mouse_junctions, 
                                    keep.extra.columns=TRUE)
View(data.frame(pred_gr))

# overlaps_utr <- findOverlaps(df_gr, regions_utr_gr, type='within')
overlaps <- findOverlaps(df_gr, pred_gr, type='within')
View(data.frame(overlaps))

overlaps_any <- findOverlaps(df_gr, pred_gr, type='any')
View(data.frame(overlaps_any))

# p_utr <- Pairs(df_gr, regions_utr_gr, hits=overlaps_utr)
p <- Pairs(df_gr, pred_gr, hits=overlaps)
p_any <- Pairs(df_gr, pred_gr, hits=overlaps_any)
View(data.frame(p))

# gene_ids_mouse <- p_any@second$gene_id %>% unique()
p <- data.frame(p)

p$name=paste('dpsi=', round(p$first.delta_cosi,2), ', ',
             'psi=', p$first.a_ctrl, '/', round(p$first.a_ctrl + 0.5*p$first.c_ctrl, 2), '=', 
             round(p$first.cosi_ctrl,2), ', ',
             'psi_chx=', p$first.a_chx, '/', round(p$first.a_chx + 0.5*p$first.c_chx, 2), '=', 
             round(p$first.cosi_chx,2), 
             sep='')
p_any <- data.frame(p_any)

p_any$name=paste('dpsi=', round(p_any$first.delta_cosi,2), ', ',
                 'psi=', p_any$first.a_ctrl, '/', round(p_any$first.a_ctrl + 0.5*p_any$first.c_ctrl, 2), '=', 
                 round(p_any$first.cosi_ctrl,2), ', ',
                 'psi_chx=', p_any$first.a_chx, '/', round(p_any$first.a_chx + 0.5*p_any$first.c_chx, 2), '=', 
                 round(p_any$first.cosi_chx,2), 
                 sep='')

# p_track <- p %>% 
#   dplyr::select(first.X.seqnames , first.X.start, first.X.end, 
#                 name, first.X.delta_cosi, first.X.strand)

p_filt <- p %>% 
  filter(first.delta_cosi >= cosi_thr & 
           first.log2fc_local >= expr_thr & 
           # second.type == 'utr3' &
           first.a_chx + first.c_chx >= reads_thr)

p_any_filt <- p_any  %>% 
  filter(first.delta_cosi >= cosi_thr & 
           first.log2fc_local >= expr_thr & 
           # second.type == 'utr3' &
           first.a_chx + first.c_chx >= reads_thr)

View(p_any_filt)




p_track_real <- p %>% 
  filter(first.delta_cosi >= cosi_thr & 
           first.log2fc_local >= expr_thr & 
           # second.type == 'utr3' &
           first.a_chx + first.c_chx >= reads_thr) %>% 
  dplyr::select(first.X.seqnames , first.X.start, first.X.end, 
                name, first.X.delta_cosi, first.X.strand)
p_track_pred <- p %>% 
  filter(first.delta_cosi >= cosi_thr & 
           first.log2fc_local >= expr_thr & 
           # second.type == 'utr3' &
           first.a_chx + first.c_chx >= reads_thr) %>% 
  dplyr::select(second.X.seqnames , second.X.start, second.X.end, 
                name, first.X.delta_cosi, first.X.strand)

p_track_real_any <- p_any %>% 
  filter(first.delta_cosi >= cosi_thr & 
           first.log2fc_local >= expr_thr & 
           # second.type == 'utr3' &
           first.a_chx + first.c_chx >= reads_thr) %>% 
  dplyr::select(first.X.seqnames , first.X.start, first.X.end, 
                name, first.X.delta_cosi, first.X.strand)
p_track_pred_any <- p_any %>% 
  filter(first.delta_cosi >= cosi_thr & 
           first.log2fc_local >= expr_thr & 
           # second.type == 'utr3' &
           first.a_chx + first.c_chx >= reads_thr) %>% 
  dplyr::select(second.X.seqnames , second.X.start, second.X.end, 
                name, first.X.delta_cosi, first.X.strand)


p_any$first.delta_cosi %>% hist()


file_name <- 'mouse_chx_liftover_all_regions.bed'
header <- 'track name="mouse_chx_whithin_real" description="mouse CHX data (whithin)" visibility=3 color=154,205,50 type=bedDetail'
cat(header, '\n', file=file_name, 
    append = FALSE)
write.table(p_track_real,
            file = file_name, append = TRUE,
            sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

header <- 'track name="mouse_chx_whithin_liftover" description="liftovered introns (whithin)" visibility=3 color=130,0,160 type=bedDetail'

cat(header, '\n', file=file_name, append = TRUE)
write.table(p_track_pred,
            file = file_name, append = TRUE,
            sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


#any


header <- 'track name="mouse_chx_any_real" description="mouse CHX data (any)" visibility=3 color=34,139,34 type=bedDetail'
cat(header, '\n', file=file_name, append = TRUE)
write.table(p_track_real_any,
            file = file_name, append = TRUE,
            sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

header <- 'track name="mouse_chx_any_liftover" description="liftovered introns (any)" visibility=3 color=75,0,130 type=bedDetail'
cat(header, '\n', file=file_name, append = TRUE)
write.table(p_track_pred_any,
            file = file_name, append = TRUE,
            sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)



