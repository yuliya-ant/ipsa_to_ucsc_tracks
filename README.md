# ipsa_to_ucsc_tracks

1. fastqc
2. mapping
3. ipsa (snakemake --cores 12)

scripts:

1. ipsa_merge_uni0.R - merges IPSA results into one data.frame
2. make_ce_list_from_ipsa.R - gets coordinates and counts of cassette exons
3. make_cassette_exon_track.R - makes a cassette exon track
