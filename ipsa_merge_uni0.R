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

write.table(df4_mouse, 'df4_mouse.tsv', quote=FALSE, sep='\t')
