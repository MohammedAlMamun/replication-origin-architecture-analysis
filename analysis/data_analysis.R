# Replication Origin Architecture
# Quantitative analyses for Figures 2 and 3
#
# This script reproduces the statistical analyses and figure panels described
# in the manuscript "Replication Origin Architecture: A Multidimensional
# Framework for Genome Duplication".
#
# Run from the repository root with:
#
#   Rscript analysis/data_analysis.R
#
# The yeast input spreadsheet is not redistributed in this repository. Place
# the unchanged Dukaj and Rhind (2021) supplementary file at data/OriDep.xlsx.
# See data/README.md for the complete citation and file details.


#==================================================
# Analysis settings
#==================================================

Nsim <- as.integer(
  Sys.getenv(
    "NSIM",
    unset = "10000"
  )
)

seed <- 123L

data_file <- Sys.getenv(
  "DATA_FILE",
  unset = file.path("data", "OriDep.xlsx")
)

output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  unset = "results"
)

if(is.na(Nsim) || Nsim < 1L)
  stop("NSIM must be a positive integer")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


#==================================================
# Packages
#==================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if(length(missing_packages) > 0){
  stop(
    "Install the required R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)


#==================================================
# Shared functions
#==================================================

# Calculate the coefficient of variation and largest inter-origin gap for a
# circular archaeal chromosome, then compare them with uniformly randomized
# origin arrangements that preserve chromosome length and origin number.

Archaea_R_value <- function(
  GenLength,
  OriPos,
  n_sim = Nsim
){
  OriPos <- sort(OriPos)

  GetDistances <- function(x){
    x <- sort(x)

    if(length(x) < 2)
      stop("Need at least 2 origins")

    # Include the interval that crosses the circular coordinate boundary.
    c(
      diff(x),
      (GenLength - max(x)) + min(x)
    )
  }

  #----------------------------------
  # Observed architecture
  #----------------------------------

  ObsDist <- GetDistances(OriPos)
  ObsCV <- sd(ObsDist) / mean(ObsDist)
  ObsMax <- max(ObsDist)

  #----------------------------------
  # Randomized architectures
  #----------------------------------

  SimCV <- numeric(n_sim)
  SimMax <- numeric(n_sim)

  for(i in seq_len(n_sim)){
    # Hashed sampling avoids constructing a full vector of every base-pair
    # coordinate while retaining uniform sampling without replacement.
    RandPos <- sort(
      sample.int(
        n = GenLength,
        size = length(OriPos),
        replace = FALSE,
        useHash = TRUE
      )
    )

    RandDist <- GetDistances(RandPos)

    SimCV[i] <- sd(RandDist) / mean(RandDist)
    SimMax[i] <- max(RandDist)
  }

  # The one-sided empirical tests ask whether the observed archaeal
  # arrangement is more regular, and has a smaller largest gap, than expected
  # under uniform random placement.

  Summary <- data.frame(
    CoefV = ObsCV,
    RanCoefMean = mean(SimCV),
    RanCoefSd = sd(SimCV),
    MaxDist = ObsMax,
    RanMaxMean = mean(SimMax),
    RanMaxSd = sd(SimMax),
    CV_p = (sum(SimCV <= ObsCV) + 1) / (n_sim + 1),
    MaxGap_p = (sum(SimMax <= ObsMax) + 1) / (n_sim + 1)
  )

  return(
    list(
      Summary = Summary,
      SimCV = SimCV,
      SimMax = SimMax,
      ObsDist = ObsDist
    )
  )
}


# Calculate adjacent inter-origin distances for linear yeast chromosomes.
# Distances are calculated within each chromosome and then pooled across the
# genome. Intervals between chromosome ends and the nearest origins are not
# included, in agreement with the manuscript methods.

DistVarGen <- function(dat){
  dat$Origin_Position_bp <- round(
    (dat$start + dat$end) / 2
  )

  GetDistances <- function(pos){
    pos <- sort(pos)

    if(length(pos) < 2)
      return(NULL)

    diff(pos)
  }

  dist_list <- lapply(
    split(
      dat$Origin_Position_bp,
      dat$chr
    ),
    GetDistances
  )

  dist_all <- unlist(
    dist_list,
    use.names = FALSE
  )

  if(length(dist_all) == 0)
    stop("No within-chromosome inter-origin distances could be calculated")

  Mean_IOD <- mean(dist_all)
  SD_IOD <- sd(dist_all)

  data.frame(
    N_origins = nrow(dat),
    N_distances = length(dist_all),
    Mean_IOD = Mean_IOD,
    SD_IOD = SD_IOD,
    CV_IOD = SD_IOD / Mean_IOD,
    Max_IOD = max(dist_all),
    Min_IOD = min(dist_all),
    Median_IOD = median(dist_all)
  )
}


#==================================================
# Figure 2
# Archaeal origin spacing and cross-system comparison
#==================================================

Figure_2 <- {

  #--------------------------------------------------
  # Archaeal chromosome data
  #--------------------------------------------------

  # Genome lengths and origin coordinates are the values analysed in the
  # manuscript and summarized in Table 1. The order here is retained from the
  # seeded analysis; each list element represents one circular chromosome.

  ArchaeaData <- data.frame(
    Species = c(
      "P.cal",
      "H.med",
      "H.vol",
      "S.isl",
      "S.aci",
      "S.sol",
      "H.nrc1",
      "H.his",
      "A.per"
    ),
    GenLength = c(
      2009313,
      2948884,
      2847757,
      2465177,
      2225959,
      2992245,
      2014239,
      2995271,
      1669695
    )
  )

  ArchaeaData$OriPos <- list(
    c(334, 520242, 891354, 1621662),
    c(576248, 1741600, 2948067),
    c(1594016, 570803, 2847458),
    c(800, 1154508, 1594240),
    c(579741, 1197431, 2225845),
    c(222537, 2008966, 740382),
    c(921939, 1806836),
    c(2994804, 1800696),
    c(445186, 1260265)
  )

  # Keep the simulation order above unchanged because the random seed is set
  # once for the complete analysis. The separate order below matches the
  # revised Figure 2C and 2D panels without altering any randomized result.

  ArchaeaPlotOrder <- c(
    "A.per",
    "H.his",
    "H.med",
    "H.nrc1",
    "H.vol",
    "P.cal",
    "S.aci",
    "S.isl",
    "S.sol"
  )

  # Set the seed once before the complete nine-chromosome analysis. This
  # reproduces the randomization sequence described in the manuscript.

  set.seed(seed)

  ArchaeaSimulations <- lapply(
    seq_len(nrow(ArchaeaData)),
    function(i){
      Archaea_R_value(
        GenLength = ArchaeaData$GenLength[i],
        OriPos = ArchaeaData$OriPos[[i]]
      )
    }
  )

  names(ArchaeaSimulations) <- ArchaeaData$Species

  #--------------------------------------------------
  # Archaeal empirical-test results
  #--------------------------------------------------

  Archaea <- do.call(
    rbind,
    lapply(
      ArchaeaSimulations,
      function(x) x$Summary
    )
  )

  Archaea <- as.data.frame(Archaea)
  Archaea$Species <- rownames(Archaea)
  rownames(Archaea) <- NULL

  Archaea$CV_q_BH <- p.adjust(
    Archaea$CV_p,
    method = "BH"
  )

  Archaea$MaxGap_q_BH <- p.adjust(
    Archaea$MaxGap_p,
    method = "BH"
  )

  Archaea_results <- data.frame(
    Species = Archaea$Species,
    CV_test = "one-sided: lower than randomized",
    Observed_CV = Archaea$CoefV,
    Randomized_mean_CV = Archaea$RanCoefMean,
    CV_p_nominal = Archaea$CV_p,
    CV_q_BH = Archaea$CV_q_BH,
    MaxGap_test = "one-sided: lower than randomized",
    Observed_MaxGap_bp = Archaea$MaxDist,
    Randomized_mean_MaxGap_bp = Archaea$RanMaxMean,
    MaxGap_p_nominal = Archaea$MaxGap_p,
    MaxGap_q_BH = Archaea$MaxGap_q_BH,
    N_simulations = Nsim,
    Seed = seed
  )

  Archaea_results_output <- Archaea_results

  Archaea_results_output[c(
    "Observed_CV",
    "Randomized_mean_CV"
  )] <- lapply(
    Archaea_results_output[c(
      "Observed_CV",
      "Randomized_mean_CV"
    )],
    round,
    digits = 4
  )

  Archaea_results_output[c(
    "Observed_MaxGap_bp",
    "Randomized_mean_MaxGap_bp"
  )] <- lapply(
    Archaea_results_output[c(
      "Observed_MaxGap_bp",
      "Randomized_mean_MaxGap_bp"
    )],
    round,
    digits = 0
  )

  Archaea_results_output[c(
    "CV_p_nominal",
    "CV_q_BH",
    "MaxGap_p_nominal",
    "MaxGap_q_BH"
  )] <- lapply(
    Archaea_results_output[c(
      "CV_p_nominal",
      "CV_q_BH",
      "MaxGap_p_nominal",
      "MaxGap_q_BH"
    )],
    round,
    digits = 6
  )

  write.csv(
    Archaea_results_output,
    file.path(output_dir, "archaeal_empirical_results.csv"),
    row.names = FALSE
  )

  #--------------------------------------------------
  # Figure 2C: archaeal inter-origin distance CV
  #--------------------------------------------------

  CV_df <- do.call(
    rbind,
    lapply(
      names(ArchaeaSimulations),
      function(SpeciesName){
        data.frame(
          Species = SpeciesName,
          SimCV = ArchaeaSimulations[[SpeciesName]]$SimCV
        )
      }
    )
  )

  CV_df$Species <- factor(
    CV_df$Species,
    levels = ArchaeaPlotOrder
  )

  ObsCV <- data.frame(
    Species = factor(
      Archaea$Species,
      levels = ArchaeaPlotOrder
    ),
    Observed = Archaea$CoefV,
    P = Archaea$CV_p
  )

  p2c <- ggplot(
    CV_df,
    aes(
      x = Species,
      y = SimCV
    )
  ) +
    geom_boxplot(
      width = 0.6,
      fill = "grey85",
      colour = "black",
      linewidth = 0.7,
      outlier.shape = NA
    ) +
    geom_point(
      data = ObsCV,
      aes(
        x = Species,
        y = Observed
      ),
      shape = 21,
      size = 3.5,
      fill = "black",
      colour = "black"
    ) +
    labs(
      x = "",
      y = "Inter-origin distance CV"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      )
    )

  #--------------------------------------------------
  # Figure 2D: archaeal largest inter-origin gap
  #--------------------------------------------------

  Max_df <- do.call(
    rbind,
    lapply(
      names(ArchaeaSimulations),
      function(SpeciesName){
        data.frame(
          Species = SpeciesName,
          SimMax = ArchaeaSimulations[[SpeciesName]]$SimMax / 1e6
        )
      }
    )
  )

  Max_df$Species <- factor(
    Max_df$Species,
    levels = ArchaeaPlotOrder
  )

  ObsMax <- data.frame(
    Species = factor(
      Archaea$Species,
      levels = ArchaeaPlotOrder
    ),
    Observed = Archaea$MaxDist / 1e6,
    P = Archaea$MaxGap_p
  )

  p2d <- ggplot(
    Max_df,
    aes(
      x = Species,
      y = SimMax
    )
  ) +
    geom_boxplot(
      width = 0.6,
      fill = "grey85",
      colour = "black",
      linewidth = 0.7,
      outlier.shape = NA
    ) +
    geom_point(
      data = ObsMax,
      aes(
        x = Species,
        y = Observed
      ),
      shape = 21,
      size = 3.5,
      fill = "black",
      colour = "black"
    ) +
    labs(
      x = "",
      y = "Largest inter-origin gap (Mbp)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      )
    )

  #--------------------------------------------------
  # Figure 2B: origin density and spacing variability
  #--------------------------------------------------

  # These 18 values are the cross-system dataset presented in Table 1C. The
  # comparison is descriptive because the underlying studies identify
  # different classes of replication origins.

  DensityCV <- data.frame(
    Species = c(
      "A.per",
      "P.cal",
      "H.med",
      "H.vol",
      "S.aci",
      "S.isl",
      "S.sol",
      "H.nrc1",
      "H.his",
      "K.lac",
      "L.wal",
      "L.klu",
      "S.pom",
      "S.cer",
      "A.tha",
      "D.mel",
      "IMR90",
      "HeLa"
    ),
    Group = c(
      rep("Archaea", 9),
      rep("Yeast", 5),
      "Plant",
      "Drosophila",
      "Human",
      "Human"
    ),
    Origins_per_Mb = c(
      1.20,
      1.99,
      1.02,
      1.05,
      1.35,
      1.22,
      1.00,
      0.99,
      0.67,
      13.83,
      19.02,
      22.30,
      36.51,
      39.93,
      13.38,
      52.76,
      61.29,
      57.16
    ),
    CV = c(
      0.033,
      0.330,
      0.358,
      0.366,
      0.335,
      0.438,
      0.418,
      0.172,
      0.287,
      0.55,
      0.58,
      0.46,
      0.86,
      0.70,
      0.98,
      0.82,
      2.31,
      2.52
    )
  )

  DensityCV$Group <- factor(
    DensityCV$Group,
    levels = c(
      "Archaea",
      "Yeast",
      "Plant",
      "Drosophila",
      "Human"
    )
  )

  DensityCV$LogDensity <- log10(
    DensityCV$Origins_per_Mb
  )

  GroupCols <- c(
    "Archaea" = "grey30",
    "Yeast" = "grey40",
    "Plant" = "black",
    "Drosophila" = "black",
    "Human" = "grey20"
  )

  GroupShapes <- c(
    "Archaea" = 16,
    "Yeast" = 7,
    "Plant" = 3,
    "Drosophila" = 15,
    "Human" = 17
  )

  p2b <- ggplot(
    DensityCV,
    aes(
      x = LogDensity,
      y = CV
    )
  ) +
    geom_point(
      aes(
        colour = Group,
        shape = Group
      ),
      size = 3.5,
      stroke = 1
    ) +
    geom_text(
      data = subset(
        DensityCV,
        Species %in% c("IMR90", "HeLa")
      ),
      aes(label = Species),
      nudge_y = -0.09,
      size = 2.5
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      colour = "grey50",
      alpha = 0.6,
      linewidth = 0.6,
      linetype = 2
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 2.55,
      label = expression(rho == 0.84),
      hjust = 1,
      size = 4
    ) +
    annotate(
      "text",
      x = 0.55,
      y = 2.40,
      label = "P < 0.001",
      hjust = 1,
      size = 4
    ) +
    scale_x_continuous(
      breaks = log10(c(1, 3, 6, 10, 30, 60)),
      labels = c("1", "3", "6", "10", "30", "60"),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    scale_colour_manual(
      values = GroupCols
    ) +
    scale_shape_manual(
      values = GroupShapes
    ) +
    labs(
      x = "Origin density (origins/Mbp)",
      y = "Inter-origin distance CV"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.title = element_blank(),
      legend.position = "right"
    )

  # Preserve the full statistical output used for the manuscript statement.

  density_statistics <- capture.output({
    cat("Spearman correlation: log10 origin density and CV\n\n")
    print(
      cor.test(
        DensityCV$LogDensity,
        DensityCV$CV,
        method = "spearman",
        exact = FALSE
      )
    )

    cat("\nLinear model: all 18 datasets\n\n")
    print(
      summary(
        lm(
          CV ~ LogDensity,
          data = DensityCV
        )
      )
    )

    cat("\nLinear model: human cell lines excluded\n\n")
    print(
      summary(
        lm(
          CV ~ LogDensity,
          data = subset(
            DensityCV,
            Group != "Human"
          )
        )
      )
    )
  })

  writeLines(
    density_statistics,
    file.path(output_dir, "origin_density_statistics.txt")
  )

  #--------------------------------------------------
  # Figure 2 assembly
  #--------------------------------------------------

  # Panel A is the manuscript schematic and is therefore represented by a
  # labelled blank area in the analysis-generated figure assembly.

  BlankA <- ggplot() +
    theme_void() +
    labs(tag = "A")

  p2b <- p2b +
    labs(tag = "B") +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      aspect.ratio = 0.55
    )

  p2c <- p2c +
    labs(tag = "C") +
    theme(legend.position = "none")

  p2d <- p2d +
    labs(tag = "D") +
    theme(legend.position = "none")

  Fig2 <- wrap_plots(
    A = BlankA,
    B = p2b,
    C = p2c,
    D = p2d,
    design = "
      AA
      BB
      CD
    ",
    heights = c(0.85, 1.15, 1)
  ) &
    theme(
      plot.tag = element_text(
        face = "bold",
        size = 18
      )
    )

  ggsave(
    filename = file.path(output_dir, "figure_2.pdf"),
    plot = Fig2,
    width = 10,
    height = 12,
    units = "in"
  )
}


#==================================================
# Figure 3
# Licensing strength and spatial organization in yeast
#==================================================

Figure_3 <- {

  #--------------------------------------------------
  # Read the published yeast dataset
  #--------------------------------------------------

  # Source: manuscript reference 14
  #
  # Dukaj, L. and Rhind, N. (2021) The capacity of origins to load MCM
  # establishes replication timing patterns. PLOS Genetics, 17, e1009467.
  #
  # OriDep.xlsx is the unchanged supplementary spreadsheet downloaded from the
  # published study. It is intentionally not redistributed here. Only the two
  # untreated 0 micromolar auxin replicates are used in the present analysis.

  if(!file.exists(data_file)){
    stop(
      "Yeast input file not found: ",
      data_file,
      "\nSee data/README.md for download and placement instructions."
    )
  }

  dat <- read_excel(data_file)

  if(ncol(dat) != 10L){
    stop(
      "OriDep.xlsx should contain 10 columns; found ",
      ncol(dat),
      "."
    )
  }

  if(nrow(dat) != 421L){
    warning(
      "The manuscript analysis used 421 origins; this file contains ",
      nrow(dat),
      "."
    )
  }

  colnames(dat) <- c(
    "chr",
    "start",
    "end",
    "ARS",
    "0uM_r1",
    "30uM_r1",
    "500uM_r1",
    "0uM_r2",
    "30uM_r2",
    "500uM_r2"
  )

  #--------------------------------------------------
  # Rank origins by mean untreated MCM signal
  #--------------------------------------------------

  # Average the two 0 micromolar auxin controls. This mean is used as an
  # indicator of relative licensing strength throughout Figure 3.

  dat_ranked <- dat %>%
    mutate(
      WT_MCM = (
        `0uM_r1` +
        `0uM_r2`
      ) / 2
    ) %>%
    arrange(desc(WT_MCM)) %>%
    mutate(
      Rank = row_number(),
      Quartile = ntile(-WT_MCM, 4)
    )

  dat_ranked$Quartile <- factor(
    dat_ranked$Quartile,
    levels = c(1, 2, 3, 4),
    labels = c(
      "Strongest",
      "Strong",
      "Weak",
      "Weakest"
    )
  )

  # The four quartiles are displayed in Figure 3B. For the spatial analyses in
  # Figure 3C and 3D, they are combined into the upper and lower halves used in
  # the manuscript. No separate four-quartile randomization is performed.

  dat_ranked$Cohort <- ifelse(
    dat_ranked$Quartile %in% c("Strongest", "Strong"),
    "Upper",
    "Lower"
  )

  dat_ranked$Cohort <- factor(
    dat_ranked$Cohort,
    levels = c("Upper", "Lower")
  )

  quartile_summary <- dat_ranked %>%
    group_by(Quartile) %>%
    summarise(
      n = n(),
      mean_MCM = mean(WT_MCM),
      median_MCM = median(WT_MCM),
      .groups = "drop"
    )

  print(quartile_summary)

  write.csv(
    dat_ranked,
    file.path(output_dir, "DukajRhind2021_ranked_origins.csv"),
    row.names = FALSE
  )

  #--------------------------------------------------
  # Figure 3A: distribution of mean MCM signal
  #--------------------------------------------------

  p3a <- ggplot(
    dat_ranked,
    aes(x = WT_MCM)
  ) +
    geom_histogram(
      bins = 30,
      fill = "grey85",
      colour = "black",
      linewidth = 0.4
    ) +
    labs(
      x = "Mean MCM ChIP signal",
      y = "Number of origins"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      )
    )

  #--------------------------------------------------
  # Figure 3B: licensing quartiles
  #--------------------------------------------------

  quartile_counts <- dat_ranked %>%
    count(Quartile)

  quartile_labels <- c(
    Strongest = paste0(
      "Strongest\n(n=",
      quartile_counts$n[quartile_counts$Quartile == "Strongest"],
      ")"
    ),
    Strong = paste0(
      "Strong\n(n=",
      quartile_counts$n[quartile_counts$Quartile == "Strong"],
      ")"
    ),
    Weak = paste0(
      "Weak\n(n=",
      quartile_counts$n[quartile_counts$Quartile == "Weak"],
      ")"
    ),
    Weakest = paste0(
      "Weakest\n(n=",
      quartile_counts$n[quartile_counts$Quartile == "Weakest"],
      ")"
    )
  )

  p3b <- ggplot(
    dat_ranked,
    aes(
      x = Quartile,
      y = WT_MCM
    )
  ) +
    geom_boxplot(
      width = 0.65,
      fill = "grey85",
      colour = "black",
      linewidth = 0.7,
      outlier.shape = 16,
      outlier.size = 1.2
    ) +
    labs(
      x = "Origin licensing cohort",
      y = "Mean MCM loading"
    ) +
    scale_x_discrete(
      labels = quartile_labels
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(
        size = 11,
        colour = "black"
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      )
    )

  #--------------------------------------------------
  # Upper and lower licensing cohorts
  #--------------------------------------------------

  UpperHalf <- subset(
    dat_ranked,
    Cohort == "Upper"
  )

  LowerHalf <- subset(
    dat_ranked,
    Cohort == "Lower"
  )

  ObsUpper <- DistVarGen(UpperHalf)
  ObsLower <- DistVarGen(LowerHalf)

  print(ObsUpper)
  print(ObsLower)

  #--------------------------------------------------
  # Fixed-location cohort randomization
  #--------------------------------------------------

  # Cohort membership is reassigned among the 421 fixed origin positions.
  # Upper and lower sample sizes are preserved in every randomization.

  set.seed(seed)

  UpperN <- nrow(UpperHalf)
  LowerN <- nrow(LowerHalf)

  SimCV_Half <- data.frame(
    Upper = numeric(Nsim),
    Lower = numeric(Nsim)
  )

  SimMax_Half <- data.frame(
    Upper = numeric(Nsim),
    Lower = numeric(Nsim)
  )

  for(i in seq_len(Nsim)){
    if(i %% 1000 == 0)
      cat("Yeast cohort randomization", i, "of", Nsim, "\n")

    RandOrder <- sample(
      seq_len(nrow(dat_ranked))
    )

    SimData <- dat_ranked[RandOrder, ]

    SimData$SimHalf <- c(
      rep("Upper", UpperN),
      rep("Lower", LowerN)
    )

    ResUpper <- DistVarGen(
      subset(
        SimData,
        SimHalf == "Upper"
      )
    )

    ResLower <- DistVarGen(
      subset(
        SimData,
        SimHalf == "Lower"
      )
    )

    SimCV_Half[i, "Upper"] <- ResUpper$CV_IOD
    SimCV_Half[i, "Lower"] <- ResLower$CV_IOD

    SimMax_Half[i, "Upper"] <- ResUpper$Max_IOD
    SimMax_Half[i, "Lower"] <- ResLower$Max_IOD
  }

  #--------------------------------------------------
  # One-sided empirical tests
  #--------------------------------------------------

  # The upper cohort is tested for lower CV and smaller maximum gap than
  # randomized. The lower cohort is tested in the opposite direction.

  CV_p_upper <- (
    sum(SimCV_Half$Upper <= ObsUpper$CV_IOD) + 1
  ) / (Nsim + 1)

  CV_p_lower <- (
    sum(SimCV_Half$Lower >= ObsLower$CV_IOD) + 1
  ) / (Nsim + 1)

  Max_p_upper <- (
    sum(SimMax_Half$Upper <= ObsUpper$Max_IOD) + 1
  ) / (Nsim + 1)

  Max_p_lower <- (
    sum(SimMax_Half$Lower >= ObsLower$Max_IOD) + 1
  ) / (Nsim + 1)

  Yeast_results <- data.frame(
    Cohort = c("Upper", "Lower"),
    CV_test = c(
      "one-sided: lower than randomized",
      "one-sided: higher than randomized"
    ),
    Observed_CV = c(
      ObsUpper$CV_IOD,
      ObsLower$CV_IOD
    ),
    Randomized_mean_CV = c(
      mean(SimCV_Half$Upper),
      mean(SimCV_Half$Lower)
    ),
    CV_p_nominal = c(
      CV_p_upper,
      CV_p_lower
    ),
    CV_q_BH = p.adjust(
      c(CV_p_upper, CV_p_lower),
      method = "BH"
    ),
    MaxGap_test = c(
      "one-sided: lower than randomized",
      "one-sided: higher than randomized"
    ),
    Observed_MaxGap_bp = c(
      ObsUpper$Max_IOD,
      ObsLower$Max_IOD
    ),
    Randomized_mean_MaxGap_bp = c(
      mean(SimMax_Half$Upper),
      mean(SimMax_Half$Lower)
    ),
    MaxGap_p_nominal = c(
      Max_p_upper,
      Max_p_lower
    ),
    MaxGap_q_BH = p.adjust(
      c(Max_p_upper, Max_p_lower),
      method = "BH"
    ),
    N_simulations = Nsim,
    Seed = seed
  )

  Yeast_results_output <- Yeast_results

  Yeast_results_output[c(
    "Observed_CV",
    "Randomized_mean_CV"
  )] <- lapply(
    Yeast_results_output[c(
      "Observed_CV",
      "Randomized_mean_CV"
    )],
    round,
    digits = 4
  )

  Yeast_results_output[c(
    "Observed_MaxGap_bp",
    "Randomized_mean_MaxGap_bp"
  )] <- lapply(
    Yeast_results_output[c(
      "Observed_MaxGap_bp",
      "Randomized_mean_MaxGap_bp"
    )],
    round,
    digits = 0
  )

  Yeast_results_output[c(
    "CV_p_nominal",
    "CV_q_BH",
    "MaxGap_p_nominal",
    "MaxGap_q_BH"
  )] <- lapply(
    Yeast_results_output[c(
      "CV_p_nominal",
      "CV_q_BH",
      "MaxGap_p_nominal",
      "MaxGap_q_BH"
    )],
    round,
    digits = 6
  )

  print(Yeast_results_output)

  write.csv(
    Yeast_results_output,
    file.path(output_dir, "yeast_licensing_empirical_results.csv"),
    row.names = FALSE
  )

  #--------------------------------------------------
  # Figure 3C: inter-origin distance CV
  #--------------------------------------------------

  CVplot <- rbind(
    data.frame(
      Cohort = "Upper",
      CV = SimCV_Half$Upper
    ),
    data.frame(
      Cohort = "Lower",
      CV = SimCV_Half$Lower
    )
  )

  CVplot$Cohort <- factor(
    CVplot$Cohort,
    levels = c("Upper", "Lower")
  )

  ObservedCV <- data.frame(
    Cohort = factor(
      c("Upper", "Lower"),
      levels = c("Upper", "Lower")
    ),
    CV = c(
      ObsUpper$CV_IOD,
      ObsLower$CV_IOD
    )
  )

  p3c <- ggplot(
    CVplot,
    aes(
      x = Cohort,
      y = CV
    )
  ) +
    geom_boxplot(
      fill = "grey85",
      colour = "black",
      width = 0.60,
      outlier.shape = NA,
      linewidth = 0.8
    ) +
    geom_point(
      data = ObservedCV,
      aes(
        x = Cohort,
        y = CV
      ),
      shape = 21,
      size = 4.5,
      fill = "black",
      colour = "black",
      stroke = 1
    ) +
    labs(
      x = "Origin licensing cohort",
      y = "Inter-origin distance CV"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      )
    )

  #--------------------------------------------------
  # Figure 3D: largest inter-origin gap
  #--------------------------------------------------

  MaxPlot <- rbind(
    data.frame(
      Cohort = "Upper",
      MaxGap = SimMax_Half$Upper
    ),
    data.frame(
      Cohort = "Lower",
      MaxGap = SimMax_Half$Lower
    )
  )

  MaxPlot$Cohort <- factor(
    MaxPlot$Cohort,
    levels = c("Upper", "Lower")
  )

  ObservedMax <- data.frame(
    Cohort = factor(
      c("Upper", "Lower"),
      levels = c("Upper", "Lower")
    ),
    MaxGap = c(
      ObsUpper$Max_IOD,
      ObsLower$Max_IOD
    )
  )

  p3d <- ggplot(
    MaxPlot,
    aes(
      x = Cohort,
      y = MaxGap / 1000
    )
  ) +
    geom_boxplot(
      fill = "grey85",
      colour = "black",
      width = 0.60,
      outlier.shape = NA,
      linewidth = 0.8
    ) +
    geom_point(
      data = ObservedMax,
      aes(
        x = Cohort,
        y = MaxGap / 1000
      ),
      shape = 21,
      size = 4.5,
      fill = "black",
      colour = "black",
      stroke = 1
    ) +
    labs(
      x = "Origin licensing cohort",
      y = "Largest inter-origin gap (kb)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      )
    ) +
    coord_cartesian(
      ylim = c(150, 650)
    )

  #--------------------------------------------------
  # Figure 3 assembly
  #--------------------------------------------------

  Fig3 <- (p3a | p3b) /
    (p3c | p3d) +
    plot_layout(
      heights = c(1, 1.5),
      widths = c(1, 1.5)
    ) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(
          size = 18,
          face = "bold"
        )
      )
    )

  ggsave(
    filename = file.path(output_dir, "figure_3.pdf"),
    plot = Fig3,
    width = 10,
    height = 8,
    units = "in"
  )
}


#==================================================
# Record the analysis environment
#==================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt")
)

cat("Analysis complete. Results written to", output_dir, "\n")
