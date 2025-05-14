<!--   주석처리  toc:table of content    
author: "Sung Young Kim, MD, PhD^[mailto:palelamp@gmail.com]</center>"
date: "03 June 2020"
output: 
    html_document:
        toc: false  
        toc_float: 
           collapsed: false
        toc_depth: 3
        
output: 
    md_document:
        variant: markdown_github
-->
<!--
This page provides the supplementary R code and data to reproduce the experiments in the following paper : **Accurate prediction of acquired EGFR TKIs resistance using a pathway-based individualized machine learning approach**  
-->

## RClass2025

**Main Dataset**
----------------

-   Preprocessed dataset can be downloaded from [here](https://github.com/centromics-m/RClass/raw/refs/heads/main/data/LC_NT_RClass.rda)
-   The pathways used for model construction can be downloaded from [here](https://github.com/centromics-m/RClass/raw/refs/heads/main/data/pathwayDB_KEGG_202411_RClass.rda)

------------------------------------------------------------------------

**Install Dependencies**
------------------------

-   If don't have the R software installed in our computer, download and install it (check out the [R home page](http://www.r-project.org/))

``` r
if (!require("devtools")) install.packages("devtools")
```

------------------------------------------------------------------------

**Install package via Github**
--------------------------------------------

To install the latest version of package via Github, run following commands in R:

``` r
devtools::install_github("centromics-m/RClass")
```

------------------------------------------------------------------------

## 202505

``` r
options(stringsAsFactors = FALSE) 
getwd(); # setwd()

pkgs = c("BiocManager", "MASS", "ggplot2", "gridExtra", "ggplotify", "pheatmap", "factoextra", "shiny", "glmnet", "reshape" )
bioc_pkgs = c("limma", "clusterProfiler", "enrichplot", "preprocessCore")

library("RClass")

install_if_missing(pkgs, bioc_pkgs)
RClass:::install_if_missing(pkgs, bioc_pkgs)



####################################################### 
# 1-1. Sample, Quantile & Scale
####################################################### 


# Set seed for reproducible results
set.seed(123)

# Simulating large populations based on birth year and sex
# Generate male 1995 population from a normal distribution (mean 1.1, SD 0.1)
population.male.1995 = rnorm(500000/2, 1.1, 0.1) # Generates 250,000 values
# Generate female 1995 population from a normal distribution (mean 0.9, SD 0.1)
population.female.1995 = rnorm(500000/2, 0.9, 0.1) # Generates 250,000 values

# Generate male 1971 population from a normal distribution (mean 1.05, SD 0.2)
population.male.1971 = rnorm(1000000/2, 1.05, 0.2) # Generates 500,000 values
# Generate female 1971 population from a normal distribution (mean 0.85, SD 0.2)
population.female.1971 = rnorm(1000000/2, 0.85, 0.2) # Generates 500,000 values


# Draw 1000 samples from the simulated populations
sample.male.1971 <- sample(population.male.1971, 1000)
sample.female.1971 <- sample(population.female.1971, 1000)



# Task 1: Understanding data characteristics using osteoporosis simulation data 
# 1. Boxplot
boxplot(sample.male.1971, horizontal = T)

# Code to calculate IQR:
q1 <- quantile(sample.male.1971, 0.25) # First Quartile (Q1)
q3 <- quantile(sample.male.1971, 0.75) # Third Quartile (Q3)
iqr <- q3 - q1

# Maximum/minimum boundaries for whiskers (Q1 - 1.5*IQR, Q3 + 1.5*IQR)
lower_limit <- q1 - 1.5 * iqr
upper_limit <- q3 + 1.5 * iqr


# 2. SE = sample standard deviation / sqrt(sample size)
# (Using Z-distribution for CI since sample is large enough - common approximation)
hist(sample.male.1971)
x_bar = mean(sample.male.1971)
abline(v=x_bar, col=2) # Add line for sample mean
se <- sd(sample.male.1971)/ sqrt(length(sample.male.1971)); se # Calculate Standard Error (SE)


# 3. Calculate 95% Confidence Interval for the mean (Using Z-distribution critical value as sample is large enough)
alpha = 0.05; # Significance level (1 - confidence level)
z_critical = qnorm(1 - alpha/2) # Same as qnorm(0.975), 95% confidence level two-tailed critical value (approx. 1.96)

lCI=mean(sample.male.1971)-z_critical * se # Lower CI
uCI=mean(sample.male.1971)+z_critical * se # Upper CI

abline(v=lCI, lty=2, col=3) # Add line for Lower CI bound
abline(v=uCI, lty=2, col=3) # Add line for Upper CI bound




# 2-1  --- Simulation to demonstrate SE ---
# Simulated SE
sample_n_sim <- 100 # Size of each sample for simulation (e.g., drawing 100 people)
num_simulations <- 10000 # Number of times to draw samples and calculate the mean (e.g., repeating 10,000 times)
sample_means <- c()
for (i in 1:num_simulations) {
  sample_data <- sample(population.male.1971, sample_n_sim) # Note: Using sample_n_sim
  sample_means[i] <- mean(sample_data)
}
simulated_se <- sd(sample_means); 
simulated_se 

hist(sample_means)
abline(v = mean(population.male.1971), col = "red", lty = 2) 


# SE value calculated using one sample's standard deviation.
one_sample <- sample(population.male.1971, sample_n_sim) # Note: Using sample_n_sim
estimated_se <- sd(one_sample) / sqrt(sample_n_sim)
estimated_se

# 'Theoretical' SE value calculated using population standard deviation (most accurate)
# Typo fix: sample.male.male.1971 - This seems to be an old comment line
theoretical_se <- sd(population.male.1971) / sqrt(sample_n_sim) 
theoretical_se



# --- Standardization and T-score Calculations ---

# Standardize the male 1995 population (Z-transformation) and show first 10 values
((population.male.1995-mean(population.male.1995))/sd(population.male.1995))[1:10]
scale(population.male.1995)[1:10]

# Calculate T-scores for the male 1971 sample using the 1995 male population as reference
male.1971.tscores = (sample.male.1971-mean(population.male.1995))/sd(population.male.1995)

# Plot histogram of the standardized male 1995 population
hist(scale(population.male.1995))
# Add a line for the mean T-score of the male 1971 sample relative to 1995 male population reference
abline(v=mean(male.1971.tscores), col=2)

# Calculate T-scores for the female 1971 sample using the 1995 female population as reference
female.1971.tscores = (sample.female.1971-mean(population.female.1995))/sd(population.female.1995)

# Add a line for the mean T-score of the female 1971 sample relative to 1995 female population reference
# Note: This line is plotted on the histogram of the *standardized male 1995 population*, which might be confusing visually.
abline(v=mean(female.1971.tscores), col=3)




####################################################### 
# 1-2. Permutation test 
####################################################### 

# Task 2: Performing hypothesis testing using permutation test 
# Generate sample data
# set.seed(123)  # Set seed for reproducibility
# group_A <- rnorm(20, mean = 5, sd = 1)
# group_B <- rnorm(20, mean = 6, sd = 1)
group_A <- male.1971.tscores
group_B <- female.1971.tscores

# Calculate observed mean difference
obs_diff <- mean(group_A) - mean(group_B)

# Set up permutation t-test
n_permutations <- 9999  # Number of permutations
perm_diffs <- numeric(n_permutations)  # Vector to store permutation differences

# Perform permutations
for (i in 1:n_permutations) {
  combined <- sample(c(group_A, group_B))  # Shuffle combined data
  perm_diffs[i] <- mean(combined[1:length(group_A)]) - mean(combined[(length(group_A)+1):length(combined)  ]) # Calculate mean difference
}

hist(perm_diffs)
abline(v=abs(obs_diff))
# Calculate p-value
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q3
p_value





####################################################### 
# 2. Standard Error & Central Limit Theorem Simulation 
####################################################### 

# Task 3: Understanding CLT
# Generate population data for different distributions
n <- 1000
population.list <- list(
  unif = runif(n, 1, 100),
  chisq = rchisq(n, 3),
  lnorm = log(rlnorm(n)),
  nbinom = rnbinom(n, size = 3, mu = 20)
)

# Function to simulate Central Limit Theorem
clt.sim <- function(population.list, iter = 1000, sample.N = 100) {
  par(mfrow = c(length(population.list), 2))


 # population.list=population.list
 # iter = 1000
 # sample.N = 100

  
  for (k in names(population.list)) {
    hist(population.list[[k]], main = k, xlab = "x")
    population <- population.list[[k]]
    
    # Sample means and standard errors
    x_bar <- lapply(1:iter, function(i) sample(population, sample.N))
    x.mean <- sapply(x_bar, mean)
    x.sd <- sapply(x_bar, sd)
    x.se <- x.sd / sqrt(sample.N)
    
    # Confidence interval coverage
    num <- sum((x.mean + 1.96 * x.se > mean(population)) & 
                 (x.mean - 1.96 * x.se < mean(population)))
    cat("1.96 estimated SE:", round(num / iter, 2), "\n")
    
    # Plot sample means
    hist(x.mean, probability = TRUE, main = paste(k, "1.96 estimated SE:", round(num / iter, 2)), xlab = "x_bar.mean")
    abline(v = mean(population) + 1.96 * mean(x.se), col = 2, lty = 2)
    abline(v = mean(population) - 1.96 * mean(x.se), col = 2, lty = 2)
    lines(density(x.mean), col = 2)
  }
  
  par(mfrow = c(1, 1))
}

# Shiny App UI and Server for interactive simulation
library(shiny)
ui <- fluidPage(
  sliderInput("iter", label = "simCLT", min = 0, max = 1000, value = 50),
  plotOutput("plot", height = "750px")
)
server <- function(input, output) {
  v <- reactiveValues()
  
  observeEvent(input$iter, {
    v$data <- input$iter
  })
  
  output$plot <- renderPlot({
    if (is.null(v$data)) return()
    clt.sim(population.list, v$data)
  })
}

# Launch the Shiny App
shinyApp(ui, server)









``` 





