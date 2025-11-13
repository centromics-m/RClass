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

pkgs = c("BiocManager", "MASS", "ggplot2", "dplyr", "gridExtra", "ggplotify", "pheatmap", "factoextra", "shiny", "glmnet", "reshape" )
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

# Maximum/minimum boundaries for whiskers
lower_limit <-     # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q1
upper_limit <- 


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
# 1-2. Permutation test 1-5조
####################################################### 

# Task 2: Performing hypothesis testing using permutation test 

# Task 2-1
set.seed(123)
# Task 2-1
set.seed(123)
group1 <- seq(0.1, 1, by = 0.1); group1
group2 <- seq(0.6, 1.5, by = 0.1); group2


observed_diff <- mean(group2) - mean(group1)
combined <- c(group1, group2)

# permutation test
n_perm <- 10000
perm_diffs <- numeric(n_perm)

for (i in 1:n_perm) {
  permuted <- sample(combined)  
  perm_group1 <- permuted[1:3]
  perm_group2 <- permuted[4:6]
  perm_diffs[i] <- mean(perm_group2) - mean(perm_group1)
}

p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
p_value

hist(perm_diffs, breaks = 30, main = "Permutation Distribution",
     xlab = "Mean Difference", col = "lightblue")
abline(v = observed_diff, col = 2, lwd = 3)




# Task 2-2
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
p_value <-    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q2
p_value





##############################################################
# 2. Standard Error & Central Limit Theorem Simulation  1-5조
##############################################################

# Task 3: Understanding CLT


# Generate population data for different distributions
# runif(10, min = 0, max = 10) # 균등분포: 0 이상 10 이하 구간에서 아무 숫자나 나올 수 있음, 모두 동일한 확률
# rpois(10, lambda = 3) # 포아송분포: 어떤 콜센터에 시간 평균 3통의 전화 
# rnbinom(10, size = 2, mu = 5) # 음이항분포: 두 통의 불만 전화를 받기까지, 평균적으로 5개의 일반 전화가 오는 분포 
# rchisq(10, df = 3) # 카이제곱 분포, 자유도 3

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
# 3. Chi-square test  4조
####################################################### 

# Observed frequency matrix
observed <- matrix(c(70, 200, 80, 600), 2)
rownames(observed) <- c("Osteoporosis", "Non-osteoporosis")
colnames(observed) <- c("Smoke", "Non-smoke")
observed

# Calculate row sums, column sums, and total sum
row_sums <- rowSums(observed)
col_sums <- colSums(observed)
total_sum <- sum(observed)

# Create expected frequency matrix            
expected <- matrix(0, nrow = nrow(observed), ncol = ncol(observed))
for (i in 1:nrow(observed)) {
  for (j in 1:ncol(observed)) {
    expected[i, j] <- (row_sums[i] * col_sums[j]) / total_sum
  }
}  
expected    

# Calculate Chi-squared statistic
chi_squared <- sum((observed - expected)^2 / expected)
chi_squared

# Calculate degrees of freedom
df <- (nrow(observed) - 1) * (ncol(observed) - 1)

# Calculate p-value
p_value <- 1 - pchisq(chi_squared, df)
print(paste("p-value:", p_value))










####################################################### 
# 4. Fisher’s Exact Test
####################################################### 

observed <- matrix(c(70, 200, 80, 600), nrow = 2)
rownames(observed) <- c("Osteoporosis", "Non-osteoporosis")
colnames(observed) <- c("Smoke", "Non-smoke")
row_sums <- rowSums(observed)
col_sums <- colSums(observed)
total_sum <- sum(observed)

# Combination function using product instead of factorial
nCk <- function(n, k) {
  if (k == 0) return(1)
  prod((n - k + 1):n) / prod(1:k)
}

# Observed count of 'a' = Osteoporosis & Smoke
a <- observed[1, 1]
# Fisher’s Exact Test calculates the probabilities for all possible values that ‘a’ can take. Since the total sum and each row and column sum are fixed, knowing only ‘a’ determines all the other cells. Therefore, calculating the probability using just ‘a’ is sufficient. It computes the probabilities for all cases where ‘a’ ranges from its minimum to maximum possible value. Then, it sums the probabilities of the observed ‘a’ value and all more extreme cases to get the p-value. Fisher’s Exact Test calculates the probability of observing the value of ‘a’ (or a more extreme value) by chance. If this probability is very low, we conclude that such a result is unlikely to occur randomly, indicating a statistically significant difference. Since knowing ‘a’ determines ‘b’, ‘c’, and ‘d’, calculating the test based on ‘b’ would yield the same probability distribution.

# Minimum and maximum possible values for 'a'
min_a <- max(0, col_sums[1] - row_sums[2])
max_a <- min(row_sums[1], col_sums[1])

# Calculate the probability of the observed table
observed_prob <- (nCk(col_sums[1], a) * nCk(col_sums[2], row_sums[1] - a)) / nCk(total_sum, row_sums[1])

# Calculate probabilities for all possible values of 'a'
probs <- numeric(max_a - min_a + 1)
idx <- 1
for (x in min_a:max_a) {
  probs[idx] <- (nCk(col_sums[1], x) * nCk(col_sums[2], row_sums[1] - x)) / nCk(total_sum, row_sums[1])
  idx <- idx + 1
}

# Calculate two-sided p-value by summing probabilities less than or equal to the observed probability
p_value <- sum(probs[probs <= observed_prob])







####################################################### 
# 5. U test  2조
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
  
mean_U  <- (n_x * n_y)/2     
sd_U <- sqrt((n_x * n_y * (n_x + n_y + 1)) / 12)

z <- (U - mean_U) / sd_U
z_critical <- qnorm(0.975)   # 1 - 0.05/2 = 0.975

if (abs(z) > z_critical) {
  cat("유의함 (귀무가설 기각)\n")
} else {
  cat("유의하지 않음 (귀무가설 채택)\n")
}
}









####################################################### 
# 6. AUC & Confusion_matrix  2조
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
  #  << Try writing the code yourself >>
  #  u_statistic <-              # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q5
  
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
#  << Try writing the code yourself >>
# precision <-         # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q6
accuracy <- (TP + TN) / (TP + TN + FP + FN)
recall <- TP / (TP + FN)            # Sensitivity or True Positive Rate
specificity <- TN / (TN + FP)       # True Negative Rate
false_discovery_rate <- FP / (TP + FP)

# Print the results
cat("Precision:", precision, "\n")
cat("Accuracy:", accuracy, "\n")
cat("Recall (Sensitivity):", recall, "\n")
cat("Specificity:", specificity, "\n")
cat("False Discovery Rate:", false_discovery_rate, "\n")





####################################################### 
# TCGA data load: Liver cancer  
#######################################################

# setwd("~/MS/edu/RClass/data-raw")
getwd() # Place the downloaded file in the current working directory and load it using the load command.
# load("LC_NT_RClass.rda")
# load("pathwayDB_KEGG_202411_RClass.rda")

str(LC_NT_RClass) # structure
names(LC_NT_RClass)

expr <- LC_NT_RClass$expr
head(expr)
View(expr)

meta <- LC_NT_RClass$meta
View(meta)



####################################################### 
# 7. Logistic regression 4조
#######################################################
# Sample data
set.seed(123)
n <- 50
data <- data.frame(
  osteoporesis = rbinom(n, 1, 0.4),       # 0/1 outcome, 40% 발생
  age = round(rnorm(n, mean=50, sd=10)),  # 나이 30~70 정도
  income = sample(1:5, n, replace=TRUE),  # 소득 1~5
  gender = sample(c("M","F"), n, replace=TRUE)  # 범주형 변수
)
head(data)

# Fit binomial logistic regression model
model <- glm(smoke ~ age + income, family = binomial(link = "logit"), data = data)
# model <- glm(income ~ age + smoke, data = data)  # default: family = gaussian
# plot(data$age, data$income, col=ifelse(data$smoke == 1, 1, 2))
# abline(glm(income ~ age, data = data))
# newdata <- data.frame(age = 35, smoke = 1)
# predict(model, newdata)

# View model summary
summary(model)
exp(coef(model))


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

# AUC?



# Fit binomial logistic regression model using the loaded TCGA liver cancer data.  # <<<<<<<<<<<<<<<<<<<<<<< Q7





####################################################### 
# 8. PCA Analysis and Visualization  1조
####################################################### 


# Basic PCA Example
p <- prcomp(matrix(rnorm(100), 5), scale. = T)  
plot(p$rotation)  # Plot rotation matrix
plot(p$x)  # Plot principal components
biplot(p, scale = FALSE)  # Biplot


# Draw a PCA plot using the loaded TCGA liver cancer data.  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Q8
# Loading necessary library
require("factoextra")

# Custom dataset creation and PCA
x <- t(expr) # Transposed data for PCA
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
# 9. Quantile Normalization  3조
#######################################################

  x= data.frame(matrix(sample(12, replace = T), 4))
  x <- x[rowSums(x)>0, ]

  tied = "average"
  rank <- apply(x, 2, rank,ties.method="min"); 
  if(any(tied %in% c("average", "max"))) rank.max <- apply(x, 2, rank,ties.method="max"); 
  sorted <- apply(x, 2, sort)
  sorted.row.mean <- apply(sorted, 1, mean); 
  x2 <- apply(rank, 2, function(x) sorted.row.mean[x])
  if(any(tied %in% c("average", "max"))) x2.max <- apply(rank.max, 2, function(x) sorted.row.mean[x])
  if(tied=="average") { x2 <- (x2+x2.max)/2 } else if (tied=="max"){x2 <- x2.max } else { }
  
  if( class(x) == "data.frame") { x2 <- as.data.frame(x2); rownames(x2) <- rownames(x) }
  x2










####################################################### 
# 10. Monty Hall problem: 5조
#######################################################











####################################################### 
# 11. geneSetTest example using limma-like logic
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
# 12. Gene Set Enrichment Analysis
#######################################################

if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")


# load("LC_NT_RClass.rda")
# load("pathwayDB_KEGG_202411_RClass.rda")


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





































https://github.com/centromics-m/RClass 

# ---- 1. 간단한 plot 그려보기 ----


# ---- 2. 변수의 선언과 할당 ----


# ---- 3. 주석달아보기 ----

sample_1 <- 10000
# sample_2 <- 20


# ---- 4. 데이터 타입 알아보기 ----


# ---- 5. 연산자 알아보기 ----


# ---- 6.데이터 구조 알아보기 ----

# 6-1.벡터 (Vector)
x <- c(1, 2, 3) # 숫자형 벡터
y <- c(TRUE, FALSE) # 논리형 벡터

# 6-2.팩터 (Factor)
factor_1 <- factor(c("low", "high", "medium"))
factor_1
levels(factor_1)[1]

factor_2 <- factor(c("control", "treat"))
factor_2
levels(factor_2)[1]

factor_3 <- factor(c("male", "female"))
levels(factor_3)
factor_3
levels(factor_3)[1]

# 6-3.리스트 (list)
list_1 <- list(name="yejin", gpa = 2.0, passed = FALSE)
list_1

# 6-4. 매트릭스 (Matrix)
mat_1 <- matrix(1:6, nrow=2)
mat_2 <- matrix(1:24, ncol=3)

# 6-5. 배열 (Array)
array_1 <- array(1:24, dim = c(2, 3, 4))

# 6-6. 데이터 프레임
df <- data.frame(name = c("yejin", "vinnuri"),
                 gpa = c(2.0, 4.5),
                 passed = c(FALSE, TRUE))

df

# ---- 7.제어문 ----

x <- 10
# if문
if (x > 5) {
  print("x는 5보다 크다")
}

# if ... else 문
x <- 3
if (x > 5) {
  print("x는 5보다 크다")
} else {
  print("x는 5 이하이다")
}

# else if 문
x <- 5
if (x > 5) {
  print("x는 5보다 크다")
} else if (x == 5) {
  print("x는 5이다")
} else {
  print("x는 5보다 작다")
}

# for 반복문
for (i in 1:3) {
  print(paste("현재 i =", i))
}

# while 반복문
x <- 1
while (x <= 3) {
  print(paste("x는", x, "입니다"))
  x <- x + 1
}

# ---- 8.함수 ----
# 8-1. 내장함수 사용해보기
# plot()
x <- 1:10
y <- x^2
plot(x, y, main = "squre function", col = "blue", type = "b")

# hist()
data <- rnorm(100)  # 평균 0, 표준편차 1인 정규분포 샘플 100개
hist(data, main = "histogram", col = "skyblue")

# 8-2.함수만들기


# ---- 9. 패키지 설치하기 ----
# 9-1. 패키지 설치
install.packages("ggplot2")

# 9-2. 패키지 불러오기
library(ggplot2)

# ---- 10. 실습 ----








