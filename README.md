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

-   Preprocessed dataset can be downloaded from [here](https://github.com/centromics-m/RClass2024/raw/refs/heads/main/data/LC_NT_RClass.rda)
-   The pathways used for model construction can be downloaded from [here](https://github.com/centromics-m/RClass2024/raw/refs/heads/main/data/pathwayDB_KEGG_202411_RClass.rda)

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
getwd(); # setwd()

pkgs = c("BiocManager", "MASS", "preprocessCore", "ggplot2", "gridExtra", "ggplotify", "pheatmap", "factoextra", "shiny", "glmnet", "reshape" )
bioc_pkgs = c("limma", "clusterProfiler", "enrichplot")

library("RClass")

install_if_missing(pkgs, bioc_pkgs)
RClass:::install_if_missing(pkgs, bioc_pkgs)


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
for (i in 1:nrow(observed)) {
  for (j in 1:ncol(observed)) {
    expected[i, j] <- (row_sums[i] * col_sums[j]) / total_sum
  }
}  
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
  mean_U  <- (n_x * n_y)/2   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q2
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
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
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
sorted.row.mean <-     # ----------------------------------------------------------------------------- Q3

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





####################################################### 
# 8. AUC & Confusion_matrix
#######################################################

calculate_auc <- function(true_labels, probabilities) {
  # Calculate the ranks of the probabilities (lower values have lower ranks)
  ranks <- rank(probabilities)
  
  # Get the ranks of the positive class instances
  pos_ranks <- ranks[true_labels == 1]
  
  # Count the number of negative class instances
  neg_count <- sum(true_labels == 0)
  
  # Calculate the sum of ranks for the positive class
  rank_sum_pos <- sum(pos_ranks)
  
  # Calculate the U statistic for rank-sum
  u_statistic <-   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q4
  
  # Calculate AUC
  auc <- u_statistic / (length(pos_ranks) * neg_count)
  return(auc)
}


true_labels <- c(1, 0, 1, 0, 1, 0, 1, 0)
probabilities <- c(0.9, 0.4, 0.8, 0.3, 0.6, 0.2, 0.7, 0.1)
# probabilities <- c(0.51, 0.49, 0.51, 0.49, 0.51, 0.49, 0.51, 0.49)
# probabilities <- c(0.51, 0.49, 0.51, 0.49, 0.51, 0.49, 0.51, 0.6)

# Calculate AUC
auc_value <- calculate_auc(true_labels, probabilities)
cat("AUC:", auc_value, "\n")

# Confusion matrix based on a threshold (0.5 in this case)
threshold <- 0.5
predictions <- ifelse(probabilities >= threshold, 1, 0)

# Create confusion matrix
confusion_matrix <- table(Predicted = predictions, Actual = true_labels)
print("Confusion Matrix:")
print(confusion_matrix)


# Extract values from confusion matrix
TP <- confusion_matrix[2, 2]  # True Positives
TN <- confusion_matrix[1, 1]  # True Negatives
FP <- confusion_matrix[2, 1]  # False Positives
FN <- confusion_matrix[1, 2]  # False Negatives

# Calculate performance metrics
precision <-            # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q5
accuracy <- 
recall <-            
specificity <-      
false_discovery_rate <-   






####################################################### 
# 8. Logistic regression
#######################################################
# Sample data
data <- data.frame(
  smoke = c(0, 1, 0, 1, 0, 1, 0, 0, 1, 1),    # Binary outcome (smoking or not)
  age = c(22, 45, 30, 50, 27, 37, 26, 34, 48, 29), # Predictor variable: age
  income = c(3, 2, 4, 1, 4, 3, 3, 4, 1, 2)   # Predictor variable: income level (1 = low, 5 = high)
); data 

# Fit binomial logistic regression model
model <- glm(smoke ~ age + income, family = binomial(link = "logit"), data = data)

# View model summary
summary(model)

# Predict probabilities of smoking
# This will output the predicted probabilities (from 0 to 1) for each individual in the dataset. These probabilities indicate the likelihood of each individual being a smoker based on their age and income.
predicted_probs <- predict(model, type = "response")
print(predicted_probs)

# Convert probabilities to binary predictions with a 0.5 threshold
predicted_labels <- ifelse(predicted_probs > 0.5, 1, 0)
print(predicted_labels)

# Create the confusion matrix
true_labels <- data$smoke  # Actual values of the response variable
confusion_matrix <- table(Predicted = predicted_labels, Actual = true_labels)
print(confusion_matrix)

auc_value <-   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q6
print(paste("AUC:", auc_value))





####################################################### 
# 8. Logistic regression - 10 fold CV
#######################################################

# Set up random seed for reproducibility
set.seed(123)

# Generate sample data
n <- 200
data <- data.frame(
  x1 = rnorm(n),             # Predictor variable 1 (normally distributed)
  x2 = rnorm(n),             # Predictor variable 2 (normally distributed)
  y = rbinom(n, 1, prob = 0.5)  # Binary outcome variable (0 or 1) with a 50% probability
)

# Split data into training and test sets (70% training, 30% test)
train_index <- sample(1:n, size = 0.7 * n, replace = FALSE)  # Randomly select indices for training set
train_data <- data[train_index, ]  # Training data subset
test_data <- data[-train_index, ]  # Test data subset

# Set up 10-fold cross-validation
k <- 10
folds <- sample(rep(1:k, length.out = nrow(train_data)))  # Randomly assign each row to one of the 10 folds

# Initialize vector to store AUC values from each fold
cv_auc <- c()

# Perform cross-validation to evaluate AUC
for (i in 1:k) {
  # Divide training data into fold-specific training and validation sets
  fold_train <- train_data[folds != i, ]  # Data for training on this fold
  fold_val <- train_data[folds == i, ]    # Data for validation on this fold
  
  # Fit logistic regression model on fold training data
  model <- glm(y ~ x1 + x2, data = fold_train, family = "binomial")
  
  # Predict probabilities on fold validation data
  pred_probs <- predict(model, newdata = fold_val, type = "response")
  
  # Calculate AUC for this fold using U-statistic
  fold_auc <- calculate_auc(fold_val$y, pred_probs)
  cv_auc <- c(cv_auc, fold_auc)  # Store AUC for this fold
}

# Print mean AUC from 10-fold cross-validation
cat("10-Fold Cross-Validation AUC:", mean(cv_auc), "\n")

# Train the final model on the entire training set
final_model <- glm(y ~ x1 + x2, data = train_data, family = "binomial")

# Predict probabilities on the test set
test_probs <- predict(final_model, newdata = test_data, type = "response")

# Calculate and print AUC for test set
test_auc <- calculate_auc(test_data$y, test_probs)
cat("Test Set AUC:", test_auc, "\n")






####################################################### 
# 9. Gene Set Enrichment Analysis
#######################################################

# Access expression data from TCGA liver cancer
data("LC_NT_RClass.rda")
data("pathwayDB_KEGG_202411_RClass.rda")

str(LC_NT_RClass)
str(pathwayDB_KEGG_202411_RClass)

# Create a gene list ordered by expression level in descending order
geneList = LC_NT_RClass$expr[order(LC_NT_RClass$expr[,1], decreasing = T), 1]

# Reshape the pathway database into a long format and keep specific columns
df <- reshape::melt(pathwayDB_KEGG_202411_RClass)
df <- df[, c(2, 1)]

# Run Gene Set Enrichment Analysis (GSEA) with the gene list and pathway data
# and calculate pairwise term similarities
GSEA.results <- pairwise_termsim(GSEA(geneList, TERM2GENE=df))

# Plot GSEA results for the top pathway
gseaplot2(GSEA.results, GSEA.results$ID[1])

# Generate various visualizations of GSEA results
emapplot(GSEA.results)     # Enrichment map plot
cnetplot(GSEA.results)     # Concept network plot
heatplot(GSEA.results)     # Heatmap plot
ridgeplot(GSEA.results)    # Ridge plot
















``` 





