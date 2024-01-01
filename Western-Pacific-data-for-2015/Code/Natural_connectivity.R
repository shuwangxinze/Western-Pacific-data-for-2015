#Adjacency matrix

adj_matrix <- read.delim(file.choose(), row.names = 1, sep = '\t', check.names = FALSE)
#Calculate natural connectivity
natural.connectivity(adj_matrix, eig = NULL, norm = F)

# List of nodes to be removed
node.attri<-read.delim(file.choose(), sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
nodes_to_remove<-as.character(node.attri$nodes_id)


# function that removes the specified node
remove_node <- function(adj_matrix, node) {
  adj_matrix <- adj_matrix[!(rownames(adj_matrix) %in% node), !(colnames(adj_matrix) %in% node), drop = FALSE]
  return(adj_matrix)
}

#Create a data frame to store the results
result_df <- data.frame(Node = character(), Connectivity = numeric(), stringsAsFactors = FALSE)
# Remove the node and calculate the natural connectivity after removal
for (node in nodes_to_remove) {
  adj_matrix <- remove_node(adj_matrix, node)
  connectivity <- natural.connectivity(adj_matrix, eig = NULL, norm = F)
  result_df <- rbind(result_df, data.frame(Node = node, Connectivity = connectivity))
}

write.table(result_df, 'groupX_Natural_connectivity.txt', sep = '\t', row.names = FALSE, quote = FALSE)


