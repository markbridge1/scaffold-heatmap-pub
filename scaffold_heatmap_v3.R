# Plot significance heatmap using Scaffold outputs
# Updated 2.4.2021 by Brandon Chen

## Usage ##

# 0. This script requires the 'igraph' and 'pheatmap' packages. If you see
#    an error like "there is no package called `igraph`", install the missing
#    package(s) and re-run the script using the relevant command below:
#   a. install.packages("igraph")
#   b. install.packages("pheatmap")

# 1. Set your working directory. Directory should contain:
#   a. Signif.csv files
#   b. clustered.txt files
#   c. scaffold file
# 2. Run the script.
#   a. Note: Loading a large .scaffold file may take some time.
# 3. In RStudio, save the output using the "Export" button above the plot.

## Dependencies ##

library(igraph)

library(pheatmap)

## Functions ##install

GetNearestLandmarkData <- function(dir) {
  # Get nearest landmark data from scaffold file
  my_load <- function(f_name) {
    con <- file(f_name, "rb")
    retval <- unserialize(con)
    close(con)
    return(retval)
  }
  sc.data <- my_load(dir(dir, pattern = "\\.scaffold")[[1]])
  G <- sc.data$graphs[[1]]
  ee <- get.edgelist(G)
  ee <- ee[V(G)[V(G)$type == 2]$highest_scoring_edge,]
  ee <- ee[, 1]
  return(ee)
}

GetSignifData <- function(dir, log.fold.change = T, mask.nonsig.lfc = T) {
  # Get significance data from Signif.csv files
  filenames <- dir(path = dir, pattern = "Signif.csv")
  m <- c()
  for (i in 1:length(filenames)) {
    csv <- read.csv(filenames[[i]])
    if (i == 1) {
      m <- cbind(m, csv[, 1])
    }
    if (log.fold.change) {
      lfc <- csv[, 4]
      if (mask.nonsig.lfc) {
        non.sig <- csv[, 2] == 0.5
        lfc[non.sig] <- 0
      }
      m <- cbind(m, lfc)
    } else {
      m <- cbind(m, csv[, 2])
    }
  }
  features <- gsub(".csv", "", filenames)
  colnames(m) <- c("cell.type", features)
  return(m)
}

GetCountData <- function(dir) {
  # Get count data from clustered.txt files
  filenames2 <- dir(pattern = "clustered.txt")
  df2 <- c()
  for (f in filenames2) {
    if (grepl(".clustered.txt.scaffold", f) == FALSE) {
      txt <- read.table(f, sep = "\t", header = TRUE)
      df2 <- cbind(df2, txt[, "popsize"])
    }
  }
  popsize <- rowSums(df2)
  return(popsize)
}

PlotSignificanceHeatmap <- function(dir, nearest.landmarks,
                                    use.log.fold.change = T,
                                    mask.nonsig.lfc = T,
                                    max.log.fold.change = NA,
                                    keep.freq.column = T,
                                    cluster.cols = T,
                                    cell.width = NA, cell.height = NA) {

  # Get significance and count data
  signif.matrix <- GetSignifData(dir, use.log.fold.change, mask.nonsig.lfc)
  if (keep.freq.column == FALSE) {
    signif.matrix <- signif.matrix[, colnames(signif.matrix) != "freqSignif"]
  }
  counts <- GetCountData(dir)

  # Aggregate and sort data
  df <- data.frame(signif.matrix, nearest.landmarks, counts)
  names(df) <- gsub("BooleanSignif", "", names(df))
  df <- df[with(df, order(counts, decreasing = T)), ]
  df <- df[with(df, order(nearest.landmarks)), ]
  rownames(df) <- paste0("row", 1:nrow(df))

  # Prep input matrix
  m <- df[!(names(df) %in% c("cell.type", "nearest.landmarks", "counts"))]
  m <- as.matrix(m)

  # Prep row gaps
  row.gaps <- table(df$nearest.landmarks)
  for (i in 2:length(row.gaps)) {
    row.gaps[i] <- row.gaps[i] + row.gaps[i - 1]
  }

  # Prep row annotation
  row.annotation <- data.frame(Count = df$counts,
                               Population = df$nearest.landmarks)
  rownames(row.annotation) <- rownames(df)

  # Prep row annotation colors
  n <- length(unique(df$nearest.landmarks))
  okabeito.pal2 <- c("#424242", "#E69F00", "#56B4E9", "#009E73",
                     "#F0e442", "#0072B2", "#D55E00", "#CC79A7")
  pop.pal <- rep(okabeito.pal2, length.out = n)
  names(pop.pal) <- unique(df$nearest.landmarks)
  ann.colors <- list(Count = c("white", "black"),
                     Population = pop.pal)

  # Prep cluster ids
  cluster.ids <- paste0("c", df$cell.type)


  # Prep breaks (max log fold change)
  if (!is.na(max.log.fold.change)) {
    max.lfc <- max.log.fold.change
  } else {
    max.lfc <- max(abs(m))
  }

  # Plot heatmap
  if (use.log.fold.change) {
    # with log fold change values
    pheatmap(mat = m,
             color = colorRampPalette(c("blue", "#e6e6e6", "red"))(101),
             breaks = seq(-1*max.lfc, max.lfc, length.out = 102),
             border_color = "grey60",
             cellwidth = cell.width,
             cellheight = cell.height,
             cluster_rows = F,
             cluster_cols = cluster.cols,
             annotation_row = row.annotation,
             annotation_colors = ann.colors,
             main = "Scaffold Significance",
             fontsize_row = 5,
             angle_col = 270,
             gaps_row = row.gaps,
             labels_row = cluster.ids)
  } else {
    # without log fold change values
    pheatmap(mat = m,
             color = colorRampPalette(c("blue", "#e6e6e6", "red"))(3),
             breaks = c(0, 0.45, 0.55, 1),
             border_color = "grey60",
             cellwidth = cell.width,
             cellheight = cell.height,
             cluster_rows = F,
             cluster_cols = cluster.cols,
             legend_breaks = c(0, 0.5, 1),
             legend_labels = c("Decreased", "Non-significant", "Increased"),
             annotation_row = row.annotation,
             annotation_colors = ann.colors,
             main = "Scaffold Significance",
             fontsize_row = 5,
             angle_col = 270,
             gaps_row = row.gaps,
             labels_row = cluster.ids)
  }
}

## Script ##

# Run once per dataset
dir <- getwd()
nearest.landmarks <- GetNearestLandmarkData(dir)  # Most time-intensive step

# Re-run just this line to tweak heatmap parameters
PlotSignificanceHeatmap(dir, nearest.landmarks,
                        use.log.fold.change = T,
                        mask.nonsig.lfc = T,
                        max.log.fold.change = 3,
                        keep.freq.column = T,
                        cluster.cols = F,
                        cell.width = NA, cell.height = NA)  # defaults: NA

