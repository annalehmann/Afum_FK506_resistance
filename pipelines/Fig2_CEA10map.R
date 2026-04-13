#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("patchwork")
library(patchwork)

setwd("")
bednames_CEA10_map <- list.files(pattern = ".bed")
samples_CEA10_map <- sub(".regions.bed", "", bednames_CEA10_map)
bedfiles_CEA10_map <- list()

# chromosome start/end coords
chr_position_CEA10_map <- c(0, 4910137, 9697127, 13967564, 17787680, 21587582, 25488251, 27373204)
chr_name_CEA10_map <- c("CP097563.1", "CP097564.1", "CP097565.1", "CP097566.1", "CP097567.1", "CP097568.1", "CP097569.1", "CP097570.1")

# read in files
for (f in seq_along(bednames_CEA10_map)) {
  d <- read.delim(bednames_CEA10_map[f], 
                  sep = "\t", 
                  header = FALSE, 
                  col.names = c("Chromosome", "Coordinate", "end", "average_coverage"))
  global_coords <- numeric(nrow(d)) 
  for (i in chr_name_CEA10_map) {
    global_coord_index <- d$Chromosome == i
    global_coords[global_coord_index] <- d$Coordinate[global_coord_index] + chr_position_CEA10_map[match(i, chr_name_CEA10_map)]
  }
  
  d$global_coord <- global_coords
  
  global_avg <- mean(d$average_coverage)
  
  d$Normalized_Coverage <- d$average_coverage / global_avg
  bedfiles_CEA10_map[[f]] <- d
}

names(bedfiles_CEA10_map) <- samples_CEA10_map

# last coord of each chr
breaks_CEA10_map <- c(0, 4910137, 9697127, 13967564, 17787680, 21587582, 25488251, 27373204, 29322347)

colors_CEA10_map <- c()

colors_CEA10_map[grepl("AF100-12_5", names(bedfiles_CEA10_map))] <- "#bcbcbc"
colors_CEA10_map[grepl("AF100-12_5_R1", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-12_5_R2", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-12_5_R3", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-12_5_R4", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-1-3", names(bedfiles_CEA10_map))] <- "#bcbcbc"
colors_CEA10_map[grepl("AF100-1-3_R1", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-1-3_R1REV", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-1-3_R2", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-1-3_R3", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-1-3_R3REV", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("AF100-1-3_R4", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("A_fum_Oryx_1", names(bedfiles_CEA10_map))] <- "#bcbcbc"
colors_CEA10_map[grepl("A_fum_Oryx_1_R2", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("A_fum_Oryx_1_R3", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("A_fum_Oryx_1_R3REV", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("A_fum_Oryx_1_R4", names(bedfiles_CEA10_map))] <- "#012169"
colors_CEA10_map[grepl("A_fum_Oryx_1_R4REV", names(bedfiles_CEA10_map))] <- "#012169"
names(colors_CEA10_map) <- names(bedfiles_CEA10_map)

# make plots
for (f in names(bedfiles_CEA10_map)) {
  p <- ggplot(bedfiles_CEA10_map[[f]], aes(x = global_coord, y = Normalized_Coverage)) + 
    geom_point(aes(fill = Chromosome), color = colors_CEA10_map[f], size = 0.001) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", 
                                      fill=NA, 
                                      linewidth = 0.25),
          axis.ticks.y = element_line(color = "black",
                                      linewidth = 0.25)) +
    geom_vline(xintercept = breaks_CEA10_map,
               linetype = "dashed",
               linewidth = 0.25) +
    xlim(0, 30527940) +
    ylim(0,2.5)
  assign(f, p)
}

figS2_CEA10_map <- `AF100-12_5` /
  `AF100-12_5_R1` /
  `AF100-12_5_R2` /
  `AF100-12_5_R3` /
  `AF100-12_5_R4` /
  `AF100-1-3_R1REV` /
  `AF100-1-3_R2` /
  `AF100-1-3_R3REV` /
`A_fum_Oryx_1_R3REV` /
 `A_fum_Oryx_1_R4REV`

fig2_plot1 <- `AF100-1-3` /
  `AF100-1-3_R1` /
  `AF100-1-3_R3` /
  `AF100-1-3_R4`

fig2_plot3 <- `A_fum_Oryx_1` /
  `A_fum_Oryx_1_R2` /
  `A_fum_Oryx_1_R3` /
  `A_fum_Oryx_1_R4`

ggsave(file="../../plots/fig2_plot1.png", plot=fig2_plot1, width=7, height=4)
ggsave(file="../../plots/fig2_plot3.png", plot=fig2_plot3, width=7, height=4)
ggsave(file="../../plots/figS2_CEA10_map.png", plot=figS2_CEA10_map, width=7, height=10)
