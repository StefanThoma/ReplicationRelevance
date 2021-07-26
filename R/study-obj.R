setClass("study", 
         slots = list(
           name = "character", 
           manyLabs = "character",
           data.revised = "data.frame", 
           data.replication = "data.frame",
           original = "list", 
           variables = "list",
           var.of.interest = "character",
           relevance.threshold = "numeric",
           family = "character",
           mixedModels = "list",
           table = "data.frame",
           difference.table = "list"
         )
)