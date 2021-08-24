#' Define the constant which will be used to round digits.
digitsForRounding <- 3

# study -------------------------------
#' An S4 class to represent a many-labs replication study.
#'
#'
#' @slot name Name of study in question
#' @slot manyLabs From which many-labs project is the data (so far only ML1 & ML5 are supported)
#' @slot data.revised Data generated with revised protocol (ML5)
#' @slot data.replication Data generated with replication protocol (ML1 & ML5)
#' @slot original Parameter of the original study. Best to use output of relevance::inference function and numer ob observations
#' @slot variables X, Y, and type of Measure (OR, SMD or Drop)
#' @slot var.of.interest Depends on the Measure. Usually just one predictorn.
#' @slot relevance.threshold Supply study specific relevance threshold
#' @slot family Model family forwarded to glm and glmer
#' @slot mixedModels Here you save the output of the MixedModels.study function
#' @slot table Here you save the output of the relevance_table.study function
#' @slot difference.table Here you save the output of the effect_differences.study function
#' @export


study <- setClass("study",
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


