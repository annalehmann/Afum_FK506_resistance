#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("patchwork")
library(patchwork)

setwd("")
bednames <- list.files(pattern = ".bed")
samples <- sub(".regions.bed", "", bednames)
bedfiles <- list()

# chromosome start/end coords
chr_position <- c(0, 4910137, 9697127, 13967564, 17787680, 21587582, 25488251, 27373204)
chr_name <- c("CP097563.1", "CP097564.1", "CP097565.1", "CP097566.1", "CP097567.1", "CP097568.1", "CP097569.1", "CP097570.1")

# read in files
for (f in seq_along(bednames)) {
  d <- read.delim(bednames[f], 
                  sep = "\t", 
                  header = FALSE, 
                  col.names = c("Chromosome", "Coordinate", "end", "average_coverage"))
  global_coords <- numeric(nrow(d)) 
  for (i in chr_name) {
    global_coord_index <- d$Chromosome == i
    global_coords[global_coord_index] <- d$Coordinate[global_coord_index] + chr_position[match(i, chr_name)]
  }
  
  d$global_coord <- global_coords
  
  global_avg <- mean(d$average_coverage)
  
  d$Normalized_Coverage <- d$average_coverage / global_avg
  bedfiles[[f]] <- d
}

names(bedfiles) <- samples

# last coord of each chr
breaks <- c(0, 4910137, 9697127, 13967564, 17787680, 21587582, 25488251, 27373204, 29322347)

colors <- c()

colors[grepl("A1163", names(bedfiles))] <- "#a7adba"  
colors[grepl("0AR", names(bedfiles))] <- "#684fa1" 
colors[grepl("0AS", names(bedfiles))] <- "#b4a8d2" 
colors[grepl("6BR", names(bedfiles))] <- "#39773a"
names(colors) <- names(bedfiles)

# make plots
for (f in names(bedfiles)) {
  p <- ggplot(bedfiles[[f]], aes(x = global_coord, y = Normalized_Coverage)) + 
    geom_point(aes(fill = Chromosome), color = colors[f], size = 0.001) +
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
    geom_vline(xintercept = breaks,
               linetype = "dashed",
               linewidth = 0.25) +
    xlim(0, 29322347) +
    ylim(0,2.5)
  assign(f, p)
}

stacked_plot <- A1163 /
  `0AR` /
  `0AS` /
  `6BR` 

ggsave(file="../Chr7_cov.svg", plot=stacked_plot, width=6, height=3)
