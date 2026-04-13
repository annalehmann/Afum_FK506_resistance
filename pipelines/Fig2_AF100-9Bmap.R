#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("patchwork")
library(patchwork)

setwd("")
bednames_AF1009B_map <- list.files(pattern = ".bed")
samples_AF1009B_map <- sub(".regions.bed", "", bednames_AF1009B_map)
bedfiles_AF1009B_map <- list()

# chromosome start/end coords
chr_position_AF1009B_map <- c(0, 4610729, 9429261, 13623370, 17299250, 21059275, 24865876, 26605102)
chr_name_AF1009B_map <- c("CM000169.1", "CM000170.1", "CM000171.1", "CM000172.1", "CM000173.1", "CM000174.1", "CM000175.1", "CM000176.1")

# read in files
for (f in seq_along(bednames_AF1009B_map)) {
  d <- read.delim(bednames_AF1009B_map[f], 
                  sep = "\t", 
                  header = FALSE, 
                  col.names = c("Chromosome", "Coordinate", "end", "average_coverage"))
  global_coords <- numeric(nrow(d)) 
  for (i in chr_name_AF1009B_map) {
    global_coord_index <- d$Chromosome == i
    global_coords[global_coord_index] <- d$Coordinate[global_coord_index] + chr_position_AF1009B_map[match(i, chr_name_AF1009B_map)]
  }
  
  d$global_coord <- global_coords
  
  global_avg <- mean(d$average_coverage)
  
  d$Normalized_Coverage <- d$average_coverage / global_avg
  bedfiles_AF1009B_map[[f]] <- d
}

names(bedfiles_AF1009B_map) <- samples_AF1009B_map

# last coord of each chr
breaks_AF1009B_map <- c(0, 4610729, 9429261, 13623370, 17299250, 21059275, 24865876, 26605102, 28371557)

colors_AF1009B_map <- c()

colors_AF1009B_map[grepl("AF100-12-8", names(bedfiles_AF1009B_map))] <- "#bcbcbc"
colors_AF1009B_map[grepl("AF100-12-8_R1", names(bedfiles_AF1009B_map))] <- "#012169"
colors_AF1009B_map[grepl("AF100-12-8_R1REV", names(bedfiles_AF1009B_map))] <- "#012169"
colors_AF1009B_map[grepl("AF100-12-8_R2", names(bedfiles_AF1009B_map))] <- "#012169"
colors_AF1009B_map[grepl("AF100-12-8_R2REV", names(bedfiles_AF1009B_map))] <- "#012169"
colors_AF1009B_map[grepl("AF100-12-8_R3", names(bedfiles_AF1009B_map))] <- "#012169"
colors_AF1009B_map[grepl("AF100-12-8_R3REV", names(bedfiles_AF1009B_map))] <- "#012169"
colors_AF1009B_map[grepl("AF100-12-8_R4", names(bedfiles_AF1009B_map))] <- "#012169"
names(colors_AF1009B_map) <- names(bedfiles_AF1009B_map)

# make plots
for (f in names(bedfiles_AF1009B_map)) {
  p <- ggplot(bedfiles_AF1009B_map[[f]], aes(x = global_coord, y = Normalized_Coverage)) + 
    geom_point(aes(fill = Chromosome), color = colors_AF1009B_map[f], size = 0.001) +
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
    geom_vline(xintercept = breaks_AF1009B_map,
               linetype = "dashed",
               linewidth = 0.25) +
    xlim(0, 30527940) +
    ylim(0,2.5)
  assign(f, p)
}

figS2_AF1009B_map <- `AF100-12-8_R1REV` /
  `AF100-12-8_R2REV` /
  `AF100-12-8_R3REV` /
  `AF100-12-8_R4`

fig2_plot2 <- `AF100-12-8` /
  `AF100-12-8_R1` /
  `AF100-12-8_R2` /
  `AF100-12-8_R3`

ggsave(file="../../plots/fig2_plot2.png", plot=fig2_plot2, width=7, height=4)
ggsave(file="../../plots/figS2_AF1009B_map.png", plot=figS2_AF1009B_map, width=7, height=4)


