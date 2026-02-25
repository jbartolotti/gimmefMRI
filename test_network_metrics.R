source('R/utils.R')
source('R/process.R')
library(igraph)

results <- computeNetworkMetrics('models/first_model', verbose=TRUE)
print("Network-level metrics (first 3 rows):")
print(head(results$network_metrics, 3))
print("\nNode-level metrics (first 10 rows):")
print(head(results$node_metrics, 10))
