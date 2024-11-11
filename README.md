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

## RClass2024

**Main Dataset**
----------------

-   The main method function R file can be downloaded from [here](http://centromics.org/info/142sup/mainFunctions.R)
-   Preprocessed dataset can be downloaded from [here](http://centromics.org/info/142sup/EGFRTKIs_8set.RData)
-   The pathways used for model construction can be downloaded from [here](http://centromics.org/info/142sup/p.KEGG.PID.BioCarta.RData)

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
devtools::install_github("centromics-m/RClass2024")
```

------------------------------------------------------------------------

## 20241112

``` r
options(stringsAsFactors = FALSE) 
getwd(); setwd()

pkgs = c("BiocManager", "MASS", "preprocessCore", "ggplot2", "gridExtra", "ggplotify", "pheatmap", "factoextra", "shiny", "glmnet" )
bioc_pkgs = c("limma")
install_if_missing(pkgs, bioc_pkgs)

library("RClass")

####################################################### 
# 1. Normal Distribution 
####################################################### 

# Uniform distribution (ontinuous & discrete)
continuous_uniform <- runif(10, min = 0, max = 1) 
discrete_uniform <- sample(1:10, size = 10, replace = T) 

# Normal distribution
# Set parameters for the normal distribution
n=100; mean = 0; sd = 1
n=1000; mean = 100; sd = 10

# Generate data
x <- rnorm(n = n, mean = mean, sd = sd)
hist_data <- hist(x)

# Probability density function (PDF) values
y <- dnorm(x, mean = mean, sd = sd)

# Quantile function values
u.q <- qnorm(0.975, mean = mean, sd = sd) # Upper 2.5% quantile
l.q <- qnorm(0.025, mean = mean, sd = sd) # Lower 2.5% quantile

# Cumulative distribution function (CDF) between quantiles
prob_95 <- pnorm(u.q, mean = mean, sd = sd) - pnorm(l.q, mean = mean, sd = sd)

# Create a data frame for plotting
df <- data.frame(x, y)

# Plot histogram of normal distribution
p1 <- ggplot(df, aes(x)) + 
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Normal Distribution",
       x = "X values",
       y = "Frequency")

# Plot PDF of normal distribution with quantile lines
p2 <- ggplot(df, aes(x, y)) +  
  geom_line() + 
  geom_vline(xintercept = c(mean, l.q, u.q), color = c("red", "blue", "blue"), linetype = "dashed") + 
  labs(title = "Probability Density Function",
       x = "x",
       y = "PDF(x)")

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)





####################################################### 
# 2. Standard Error & Central Limit Theorem Simulation 
####################################################### 

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





####################################################### 
# 3. Chi-square test
####################################################### 

# Observed frequency matrix
observed <- matrix(c(30, 70, 50, 150), 2)
rownames(observed) <- c("Male", "Female")
colnames(observed) <- c("Smoke", "Non-smoke")
observed

# Calculate row sums, column sums, and total sum
row_sums <- rowSums(observed)
col_sums <- colSums(observed)
total_sum <- sum(observed)

# Create expected frequency matrix
expected <- matrix(0, nrow = nrow(observed), ncol = ncol(observed))

# Calculate expected frequencies using a for loop 
expected   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q1

# Calculate Chi-squared statistic
chi_squared <- sum((observed - expected)^2 / expected)
print(paste("Chi-squared statistic:", chi_squared))

# Calculate degrees of freedom
df <- (nrow(observed) - 1) * (ncol(observed) - 1)
print(paste("Degrees of freedom:", df))

# Calculate p-value
p_value <- 1 - pchisq(chi_squared, df)
print(paste("p-value:", p_value))





####################################################### 
# 4. U test
####################################################### 

U_test <- function(x, y, alternative = "two.sided") {
  n_x <- length(x)
  n_y <- length(y)
  
  # Calculate ranks for combined data and U statistic for Group x
  rank_data <- rank(c(x, y))
  U_x <- sum(rank_data[1:n_x]) - (n_x * (n_x + 1)) / 2
  
  # Calculate U statistic for Group y and select the minimum U
  U_y <- n_x * n_y - U_x
  U <- min(U_x, U_y)
  
  # Calculate the mean and standard deviation of U under the null hypothesis
  mean_U    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q2
  sd_U <- sqrt((n_x * n_y * (n_x + n_y + 1)) / 12)
  
  # Calculate p-value based on test alternative
  z <- (U - mean_U) / sd_U
  if (alternative == "two.sided") {
    p_value <- 2 * pnorm(-abs(z))  
  } else if (alternative == "greater") {
    p_value <- pnorm(-z)
  } else if (alternative == "less") {
    p_value <- pnorm(z)
  } else {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  
  return(p_value)
}

# Test the function with example data
x <- c(1, 2, 3, 4)
y <- c(5, 6, 7, 8)
U_test(x, y, alternative = "two.sided")
U_test(x, y, alternative = "greater")
U_test(x, y, alternative = "less")

# Comparison with R's Base Functions
getS3method("wilcox.test", "default")
wilcox.test(x, y, paired=F, correct = F, exact = F, alternative = "two.sided")
wilcox.test(x, y, paired=F, correct = F, exact = F, alternative = "greater")
wilcox.test(x, y, paired=F, correct = F, exact = F, alternative = "less")





####################################################### 
# 5-1. Permutation t-test 
####################################################### 

# Generate sample data
set.seed(123)  # Set seed for reproducibility
group_A <- rnorm(20, mean = 5, sd = 1)
group_B <- rnorm(20, mean = 6, sd = 1)

# Calculate observed mean difference
obs_diff <- mean(group_A) - mean(group_B)

# Set up permutation t-test
n_permutations <- 9999  # Number of permutations
perm_diffs <- numeric(n_permutations)  # Vector to store permutation differences

# Perform permutations
for (i in 1:n_permutations) {
  combined <- sample(c(group_A, group_B))  # Shuffle combined data
  perm_diffs[i] <- mean(combined[1:20]) - mean(combined[21:40])  # Calculate mean difference
}

# Calculate p-value 
p_value   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q3






####################################################### 
# 5-2. geneSetTest example using limma-like logic
####################################################### 

# Testing whether a subset of gene statistics is higher than expected by chance

# Generate example data
statistics <- rnorm(100)  # All gene statistics
index <- 1:10  # Indices for the gene set
statistics[index] <- statistics[index] + 1  # Modify gene set statistics for testing
alternative <- "mixed"  # Testing alternative hypothesis
nsim <- 9999  # Number of simulations

# Extract and filter the subset of statistics
ssel <- statistics[index]
ssel <- ssel[!is.na(ssel)]
nsel <- length(ssel)
if (nsel == 0) return(1)  # If no selected statistics, return p-value of 1

# Set up statistic and function for alternatives
stat <- statistics[!is.na(statistics)]
msel <- mean(ssel)
if (alternative == "either") { 
  posstat <- abs 
} else { 
  posstat <- function(x) x 
}
msel <- posstat(msel)

# Simulation to compute p-value
ntail <- 1
for (i in 1:nsim) {
  if (posstat(mean(sample(stat, nsel))) >= msel) {
    ntail <- ntail + 1
  }
}
p.value <- ntail / (nsim + 1)
p.value





####################################################### 
# 6. PCA Analysis and Visualizatio
####################################################### 

# Basic PCA Example
p <- prcomp(matrix(rnorm(100), 5))  # Sample data for demonstration
plot(p$rotation)  # Plot rotation matrix
plot(p$x)  # Plot principal components
biplot(p, scale = FALSE)  # Biplot
prcomp(matrix(rnorm(100), 20))  # Demonstrates PCA with scale option explained

# Loading necessary library
require("factoextra")

# Custom dataset creation and PCA
x <- t(makeSimData2(25, 50)$exp)  # Transposed data for PCA
group.col.N <- NULL
scale <- FALSE
addEllipses <- FALSE
gradient.cols <- c("#00AFBB", "#E7B800", "#FC4E07")  # Color scheme for plots

# Run PCA
pca <- prcomp(x, scale = scale)  # Set scale = TRUE for standardization if variables have different units
summary(pca)
print(pca$rotation)

# Eigenvalue Visualization
p1 <- fviz_eig(pca)  # Show % of variances explained by each PC
print(p1)

# Define color by groups or quality of representation (cos2)
if (!is.null(group.col.N) && is.numeric(group.col.N)) {
  groups <- as.factor(x[[group.col.N]])
  col.ind <- groups
} else {
  col.ind <- "cos2"
}

# Set label options based on data size
label.ind <- if (dim(x)[1] > 30) "none" else "all"
label.var <- if (dim(x)[2] > 30) "none" else "all"
label.bi <- if (dim(x)[1] > 30 || dim(x)[2] > 30) "none" else "all"

# Sample Plot
p2 <- fviz_pca_ind(pca, col.ind = col.ind, gradient.cols = gradient.cols, 
                   repel = TRUE, addEllipses = addEllipses, label = label.ind)
print(p2)

# Variable Plot
p3 <- fviz_pca_var(pca, col.var = "contrib", gradient.cols = gradient.cols, 
                   repel = TRUE, label = label.var)
print(p3)






####################################################### 
# 7-1. Scale and Quantile Normalization 
#######################################################

# Quantile Normalization

# Sample data matrix
x <- data.frame(matrix(sample(12, replace = TRUE), 4))

# Parameters
tied <- "average"  # Options: "average", "max", etc.
verbose <- TRUE

# Rank the data
rank <- apply(x, 2, rank, ties.method = "min")

# If necessary, apply max method for ties
if (any(tied %in% c("average", "max"))) {
  rank.max <- apply(x, 2, rank, ties.method = "max")
}

# Sort the data
sorted <- apply(x, 2, sort)

# Calculate row-wise means
sorted.row.mean <- apply(sorted, 1, mean)

# Apply the rank-based transformation
x2 <- apply(rank, 2, function(x) sorted.row.mean[x])

# If using max for ties, apply it
if (any(tied %in% c("average", "max"))) {
  x2.max <- apply(rank.max, 2, function(x) sorted.row.mean[x])
}

# Combine results depending on the chosen tie-breaking method
if (tied == "average") {
  x2 <- (x2 + x2.max) / 2
} else if (tied == "max") {
  x2 <- x2.max
}

# Convert back to data frame or matrix if needed
if (any(class(x) == "data.frame")) {
  x2 <- as.data.frame(x2)
  rownames(x2) <- rownames(x)
} else if (any(class(x) == "matrix")) {
  x2 <- as.matrix(x2)
  rownames(x2) <- rownames(x)
}

# If verbose output is enabled, show intermediate steps
if (verbose) {
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), mar = c(3, 3, 1, 1))
  
  cat('Original matrix or data.frame\n')
  peep(x, TRUE)
  cat('Sort the original matrix from lowest to highest\n')
  peep(rank)
  cat('Determine the ranks of original matrix\n')
  peep(sorted)
  cat('\nCalculate the means\n\n')
  peep(sorted.row.mean)
  cat('\n\nFinally substitute the means into our ranked matrix\n')
  peep(x2, TRUE)
  cat(sprintf('If the values were tied, %s is used\n\n', tied))
  
  par(op)
}

# Final output
x2




####################################################### 
# 7-2. Scale and normalization for different techniques
#######################################################

sim <- makeSimData2(25, 50)
sim <- makeSimData2(25, 50, seed = 15, technique = "rna_seq_counts")
sim <- makeSimData2(25, 50, seed = 15, technique = "sc_rna_seq_counts")

# Heatmap of the expression data
pheatmap(sim$exp, cluster_rows = FALSE, cluster_cols = FALSE)

# Extract expression data
x <- sim$exp

# Perform different normalization techniques
q <- normalize.quantiles(as.matrix(x))  # Quantile normalization
xx <- list(
  x = x, 
  l = log2(x),  # Log transformation
  q = q,        # Quantile normalized data
  z = scale(x),  # Z-score normalization
  tzt = t(scale(t(x))),  # Transpose, scale and then transpose back
  qz = scale(q),  # Z-score normalization of quantile normalized data
  gz = (x - mean(x)) / sd(x),  # Alternative Z-score normalization
  ql = log2(q),  # Log transformation after quantile normalization
  qltzt = t(scale(t(log2(q))))  # Transpose, log transform, scale and transpose back
)

# Boxplot visualization of the normalized data
par(mfrow = c(2, 4), mar = c(3, 4, 2, 1))
for (xxx in names(xx)) { 
  boxplot(xx[[xxx]], main = xxx) 
}

# Create heatmaps for each normalization method
a <- list()
for (i in names(xx)) {
  a[[i]] <- as.ggplot(pheatmap(xx[[i]], main = i, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE))
}

# Arrange heatmaps in a grid
gridExtra::grid.arrange(grobs = a)

``` 





