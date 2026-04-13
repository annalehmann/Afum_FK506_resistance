#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("patchwork")
library(patchwork)

setwd("")
bednames_W72310map <- list.files(pattern = ".bed")
samples_W72310map <- sub(".regions.bed", "", bednames_W72310map)
bedfiles_W72310map <- list()

# chromosome start/end coords
chr_position_W72310map <- c(0, 4929293, 9626092, 13460185, 17062713, 20901619, 24709728, 26402470)
chr_name_W72310map <- c("CM000169.1", "CM000170.1", "CM000171.1", "CM000172.1", "CM000173.1", "CM000174.1", "CM000175.1", "CM000176.1")

# read in files
for (f in seq_along(bednames_W72310map)) {
  d <- read.delim(bednames_W72310map[f], 
                  sep = "\t", 
                  header = FALSE, 
                  col.names = c("Chromosome", "Coordinate", "end", "average_coverage"))
  global_coords <- numeric(nrow(d)) 
  for (i in chr_name_W72310map) {
    global_coord_index <- d$Chromosome == i
    global_coords[global_coord_index] <- d$Coordinate[global_coord_index] + chr_position_W72310map[match(i, chr_name_W72310map)]
  }
  
  d$global_coord <- global_coords
  
  global_avg <- mean(d$average_coverage)
  
  d$Normalized_Coverage <- d$average_coverage / global_avg
  bedfiles_W72310map[[f]] <- d
}

names(bedfiles_W72310map) <- samples_W72310map

# last coord of each chr
breaks_W72310map <- c(0, 4929293, 9626092, 13460185, 17062713, 20901619, 24709728, 26402470, 28150988)

colors_W72310map <- c()

colors_W72310map[grepl("W72310", names(bedfiles_W72310map))] <- "#bcbcbc"
colors_W72310map[grepl("W72310_R1", names(bedfiles_W72310map))] <- "#012169"
colors_W72310map[grepl("W72310_R2", names(bedfiles_W72310map))] <- "#012169"
colors_W72310map[grepl("W72310_R3", names(bedfiles_W72310map))] <- "#012169"
names(colors_W72310map) <- names(bedfiles_W72310map)

# make plots
for (f in names(bedfiles_W72310map)) {
  p <- ggplot(bedfiles_W72310map[[f]], aes(x = global_coord, y = Normalized_Coverage)) + 
    geom_point(aes(fill = Chromosome), color = colors_W72310map[f], size = 0.001) +
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
    geom_vline(xintercept = breaks_W72310map,
               linetype = "dashed",
               linewidth = 0.25) +
    xlim(0, 30527940) +
    ylim(0,2.5)
  assign(f, p)
}

figS2_ATCC46645_map <- `W72310` /
  `W72310_R1` /
  `W72310_R2` /
  `W72310_R3` 

ggsave(file="../plots/figS2_W72310_map_cov.png", plot=figS2_W72310_map, width=7, height=4)
