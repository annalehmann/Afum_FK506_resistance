#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("patchwork")
library(patchwork)

setwd("")
bednames_ATCC46645_map <- list.files(pattern = ".bed")
samples_ATCC46645_map <- sub(".regions.bed", "", bednames_ATCC46645_map)
bedfiles_ATCC46645_map <- list()

# chromosome start/end coords
chr_position_ATCC46645_map <- c(0, 5144364, 10028163, 14408132, 18458229, 22855349, 26866601, 28721294)
chr_name_ATCC46645_map <- c("CM000169.1", "CM000170.1", "CM000171.1", "CM000172.1", "CM000173.1", "CM000174.1", "CM000175.1", "CM000176.1")

# read in files
for (f in seq_along(bednames_ATCC46645_map)) {
  d <- read.delim(bednames_ATCC46645_map[f], 
                  sep = "\t", 
                  header = FALSE, 
                  col.names = c("Chromosome", "Coordinate", "end", "average_coverage"))
  global_coords <- numeric(nrow(d)) 
  for (i in chr_name_ATCC46645_map) {
    global_coord_index <- d$Chromosome == i
    global_coords[global_coord_index] <- d$Coordinate[global_coord_index] + chr_position_ATCC46645_map[match(i, chr_name_ATCC46645_map)]
  }
  
  d$global_coord <- global_coords
  
  global_avg <- mean(d$average_coverage)
  
  d$Normalized_Coverage <- d$average_coverage / global_avg
  bedfiles_ATCC46645_map[[f]] <- d
}

names(bedfiles_ATCC46645_map) <- samples_ATCC46645_map

# last coord of each chr
breaks_ATCC46645_map <- c(0, 5144364, 10028163, 14408132, 18458229, 22855349, 26866601, 28721294, 30527940)

colors_ATCC46645_map <- c()

colors_ATCC46645_map[grepl("08-19-02-10", names(bedfiles_ATCC46645_map))] <- "#bcbcbc"
colors_ATCC46645_map[grepl("08-19-02-10_R1", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("08-19-02-10_R1REV", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("08-19-02-10_R2", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("08-19-02-10_R3", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("08-19-02-10_R3REV", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("08-19-02-10_R4", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("08-19-02-10_R4REV", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("10-01-02-27", names(bedfiles_ATCC46645_map))] <- "#bcbcbc"
colors_ATCC46645_map[grepl("10-01-02-27_R1", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("10-01-02-27_R2", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("10-01-02-27_R3", names(bedfiles_ATCC46645_map))] <- "#012169"
colors_ATCC46645_map[grepl("10-01-02-27_R4", names(bedfiles_ATCC46645_map))] <- "#012169"
names(colors_ATCC46645_map) <- names(bedfiles_ATCC46645_map)

# make plots
for (f in names(bedfiles_ATCC46645_map)) {
  p <- ggplot(bedfiles_ATCC46645_map[[f]], aes(x = global_coord, y = Normalized_Coverage)) + 
    geom_point(aes(fill = Chromosome), color = colors_ATCC46645_map[f], size = 0.001) +
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
    geom_vline(xintercept = breaks_ATCC46645_map,
               linetype = "dashed",
               linewidth = 0.25) +
    xlim(0, 30527940) +
    ylim(0,2.5)
  assign(f, p)
}

figS2_ATCC46645_map <- `08-19-02-10_R1REV` /
  `08-19-02-10_R2` /
  `08-19-02-10_R3REV` /
  `08-19-02-10_R4REV` /
  `10-01-02-27` /
  `10-01-02-27_R1` /
  `10-01-02-27_R2` /
  `10-01-02-27_R3` /
  `10-01-02-27_R4` 

fig2_plot4 <- `08-19-02-10` /
  `08-19-02-10_R1` /
  `08-19-02-10_R3` /
  `08-19-02-10_R4`

ggsave(file="../plots/fig2_plot4_cov.png", plot=fig2_plot4, width=7, height=4)
ggsave(file="../plots/figS2_ATCC46645_map_cov.png", plot=figS2_ATCC46645_map, width=7, height=9)

