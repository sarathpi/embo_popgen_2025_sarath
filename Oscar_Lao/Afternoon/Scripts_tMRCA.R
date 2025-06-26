rm(list=ls());

folder.fastSimcoal2 <- "C:\\Users\\u9424\\OneDrive - Universitat Pompeu Fabra\\Grants\\2025\\EMBO Naples June\\Second day\\Afternoon\\Practical\\";

setwd(folder.fastSimcoal2);

# b) assuming a constant population size A with a recent population change
model.single.pop.with.recent.population.change.tMRCA <- function(population_size_in_present, effective_population_size_1, time_change_population_size, number_of_blocks)
{
  ex <- paste(folder.fastSimcoal2,"DemographicModelSplitR\\" ,"DemographicModelSplitR_mrca.txt",sep="");
  if (file.exists(ex)) {
    file.remove(ex)
  }  
  r <- (effective_population_size_1/population_size_in_present);
  lines <- c(
    "//Number of population samples (demes)",
    "1",
    "//Population effective sizes (number of genes)",
    population_size_in_present,
    "//Sample sizes",
    "2",
    "//Growth rates\t: negative growth implies population expansion",
    "0",
    "//Number of migration matrices : 0 implies no migration between demes",
    "0",
    "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
    "1 historical event",
    paste(time_change_population_size, "0 0 1", r,"0 0",sep = " "),
    "//Number of independent loci [chromosome]",
    paste(number_of_blocks,"0",sep=" ")
  )
  
  # Repeating the block 22 times
  block <- c(
    "//Per chromosome: Number of linkage blocks",
    "1",
    "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
    "DNA 1 0 0"
  )
  
  # Append number_of_blocks blocks
  for (i in 1:number_of_blocks) {
    lines <- c(lines, block)
  }

  # Write to file
  writeLines(lines, "DemographicModelSplitR.par")
  
  exe <- ".\\fsc28.exe"
  args <- c("-i", ".\\DemographicModelSplitR.par", "-x", "-s0", "-n", "1", "-q", "--recordMRCA")
  # Execute
  system2(exe, args = args);
  data.t <- read.table(file=paste(folder.fastSimcoal2,"\\DemographicModelSplitR\\DemographicModelSplitR_mrca.txt", sep=""), header = T, skip = 4);
  return(data.t);
}

# Example of execution
data <- model.single.pop.with.recent.population.change.tMRCA(10000, 1000, 500, 1000)
hist(data$Time_.gen., breaks = 50, col = "skyblue", main = "TMRCA Distribution", xlab = "Generations")
abline(v = mean(data$Time_.gen.), col = "red", lwd = 2)

# Exercise 1). Run 100 times the model with Ne current = 5000, Ne ancestral = 10.
# Change at random the time of change of Ne, from 10 to 1000 generations.
# Run 10000 independent fragments. Read the output file and recover the tmrca
# Compute the mean(log(tmrca)). Plot this value against the time of change of Ne

time.change <- rep(NA,100);
mean.log.tmrca <- rep(NA,100);
for(i in 1:100)
{
  time.c <- floor(runif(1,10,1000));
  time.change[i] <- time.c;
  d <- model.single.pop.with.recent.population.change.tMRCA(5000,10,time.c,10000);
  mean.log.tmrca[i] <- mean(log(d$Time_.gen.));
}

# Exercise 2. Make the inference of the inverse coalescence rate, 
# which can be approximated to Ne under the assumption of panmixia and lack of
# population substructure.
# Do some simulations changing the time when there is the reduction in population
# size. Is it working?

iicr <- function(data)
{
  tmrca_values <- data$Time_.gen.
  # Define time bins
  time_bins <- seq(0, 9000, by = 100)
  data$time_bin <- cut(data$Time_.gen., breaks = time_bins, include.lowest = TRUE, right = FALSE)
  
  # PDF of TMRCA (count per bin)
  library(dplyr)
  pdf_tbl <- data %>%
    group_by(time_bin) %>%
    summarise(pdf = n(), .groups = "drop")
  
  # Compute midpoints of time bins
  bin_mids <- sapply(strsplit(as.character(pdf_tbl$time_bin), ","), function(x) {
    as.numeric(gsub("\\[|\\)|\\]", "", x)) |> mean()
  })
  
  pdf_tbl$bin_mid <- bin_mids
  
  #CDF (cumulative count)
  pdf_tbl <- pdf_tbl %>%
    mutate(cdf = cumsum(pdf) - pdf)
  
  # Calculate proportions and IICR
  total <- nrow(data)
  pdf_tbl <- pdf_tbl %>%
    mutate(
      cdf_prop = cdf / total,
      pdf_prop = pdf / total,
      coalescence_rate = ifelse((1 - cdf_prop) > 0, pdf_prop / (1 - cdf_prop), NA),
      IICR = ifelse(coalescence_rate > 0, 1 / coalescence_rate, NA)
    )
  
  # Plot
  library(ggplot2)
  ggplot(pdf_tbl, aes(x = bin_mid, y = IICR)) +
    geom_line(color = "steelblue") +
    labs(x = "Time (generations ago)", y = "IICR", title = "Inverse Instantaneous Coalescence Rate") +
    theme_minimal()  
}

data <- model.single.pop.with.recent.population.change.tMRCA(2000, 500, 2000, 100000)
iicr(data)

# Exercise 3. Goal: Analyze how the strength and timing of bottlenecks impact coalescence times.
# Instructions:
# Run several simulations with effective_population_size_1 = 1,000, 3,000, 5,000, and 10,000.
# Repeat for different time_change_population_size values: 100, 500, 1,000.
# Plot average TMRCA for each combination.
# Questions:
# How does bottleneck timing affect TMRCA more than its severity?
# For which combination does TMRCA reduce the most?


# Exercise 4: Match Simulation to Analytical Expectation
# Goal: Validate coalescence theory under constant size using expected values.
# Instructions:
# Set effective_population_size_1 = population_size_in_present = N
# Under a constant N, the expected TMRCA of two haploid individuals is 
# 2N generations.
# Simulate and compute the empirical mean TMRCA from the file.
# Questions:
# How close is the simulated mean TMRCA to 2N?
# How many blocks are needed for it to converge?
 
model.single.pop.with.bottleneck.tMRCA <- function(population_size_in_present, ne_bottleneck, time_change_population_size, time_duration_bottleneck, number_of_blocks)
{
  ex <- paste(folder.fastSimcoal2,"DemographicModelSplitR\\" ,"DemographicModelSplitR_mrca.txt",sep="");
  if (file.exists(ex)) {
    file.remove(ex)
  }  
  r <- (effective_population_size_1/population_size_in_present);
  lines <- c(
    "//Number of population samples (demes)",
    "1",
    "//Population effective sizes (number of genes)",
    population_size_in_present,
    "//Sample sizes",
    "2",
    "//Growth rates\t: negative growth implies population expansion",
    "0",
    "//Number of migration matrices : 0 implies no migration between demes",
    "0",
    "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
    "2 historical event",
    paste(time_change_population_size, "0 0 1", r,"0 0",sep = " "),
    paste((time_change_population_size + time_duration_bottleneck), "0 0 1", 1/r,"0 0",sep = " "),
    "//Number of independent loci [chromosome]",
    paste(number_of_blocks,"0",sep=" ")
  )
  
  # Repeating the block 22 times
  block <- c(
    "//Per chromosome: Number of linkage blocks",
    "1",
    "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
    "DNA 1 0 0"
  )
  
  # Append number_of_blocks blocks
  for (i in 1:number_of_blocks) {
    lines <- c(lines, block)
  }
  
  # Write to file
  writeLines(lines, "DemographicModelSplitR.par")
  
  exe <- ".\\fsc28.exe"
  args <- c("-i", ".\\DemographicModelSplitR.par", "-x", "-s0", "-n", "1", "-q", "--recordMRCA")
  # Execute
  system2(exe, args = args);
  data.t <- read.table(file=paste(folder.fastSimcoal2,"\\DemographicModelSplitR\\DemographicModelSplitR_mrca.txt", sep=""), header = T, skip = 4);
  return(data.t);
}

# Exercise 5: Bottleneck Detectability Threshold
# Goal: Explore the conditions under which a bottleneck becomes "invisible."
#Instructions:
#  Fix present Ne = 10,000
# Run simulations with:
# Bottleneck Ne = 500
# Bottleneck duration = 50 generations
# Bottleneck timing = vary (recent vs. ancient)
# Use iicr() to compute and plot Ne.
# Find the minimum bottleneck duration or severity that shows a visible dip in IICR.  

# Exercise 6: Explore Forward vs. Backward Thinking
# Goal: Reinforce the idea that the coalescent is a backward-time process.
# Instructions:
#  Simulate a recent population expansion:
#  Ancestral Ne = 1,000
# Present-day Ne = 10,000
# Change at 500 generations ago
# Ask:
#  Is the mean TMRCA affected more by the ancestral or present Ne?
#  What part of the IICR curve reflects present-day size?
#  Discuss the “memory” of TMRCA about the past more than the present.
