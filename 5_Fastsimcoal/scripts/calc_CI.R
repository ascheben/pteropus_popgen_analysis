# Load necessary library
library(dplyr)

# Function to calculate 95% confidence interval using mean and standard error
calculate_ci <- function(x) {
  n <- length(x)
  mean_x <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(n)
  error_margin <- qnorm(0.975) * se  # 1.96 * SE for 95% CI
  lower_bound <- mean_x - error_margin
  upper_bound <- mean_x + error_margin
  return(c(mean = mean_x, lower_CI = lower_bound, upper_CI = upper_bound))
}

# Function to calculate 95% quantile-based confidence interval
calculate_quantile_ci <- function(x) {
  lower_bound <- quantile(x, 0.025, na.rm = TRUE)
  upper_bound <- quantile(x, 0.975, na.rm = TRUE)
  return(c(lower_quantile_CI = lower_bound, upper_quantile_CI = upper_bound))
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide an input file name.")
}
file_path <- args[1]

# Generate output file name
output_file <- sub("\\..*$", ".confidence_intervals.txt", file_path)

# Read the input file with a header
data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convert all columns to numeric
data <- data %>% mutate(across(everything(), as.numeric))

# Get column names
param_names <- colnames(data)

# Apply confidence interval calculation to each column
ci_results <- sapply(data, calculate_ci)
quantile_ci_results <- sapply(data, calculate_quantile_ci)

# Convert results to data frame
ci_df <- as.data.frame(t(ci_results))
quantile_ci_df <- as.data.frame(t(quantile_ci_results))

# Add parameter names
final_results <- cbind(Parameter = param_names, ci_df, quantile_ci_df)

# Print or save the results
print(final_results)
write.table(final_results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
