#' Revision protocol dataframe of replication project ML5
#' Study: Alb5
#'
#' The experiment 5 of albarracin (2008) is a classic two group experimental setup.
#' The research hypothesis was whether priming a person with *action* vs priming with *inaction*
#' (the experimental `Condition`) would increase the number of correct solutions on a cognitive test
#' featuring 21 questions assessing verbal ability and quantitative ability (`SATTotal`).
#' albarracin (2008) found that participants primed with *action* achieved a significantly
#' higher number of correct results

#' The main difference between the two study protocols was that the revised protocol was *not*
#' conducted online but in person, just like the original study.
#' For more differences, see the [osf page](https://osf.io/x3gce/).
#' This study was replicated within the Many Labs 5 project in eight locations using
#' the revised protocol (with sample sizes $N$ ranging from 81 to 174), and online
#' through Amazon Mechanical Turk (*MTurk*) using the replication protocol, ebersole, 2020.
#' We used this study to exemplify the Relevance procedure for two independent samples, estimating the SMD.
#'
#' @format A data frame with 884 rows and 4 variables:
#' \describe{
#'   \item{ResponseId}{unique ID of participant}
#'   \item{Location}{Grouping variable for where the replication of the experiment took place}
#'   \item{SATTotal}{Sum of all SAT variables, and dependent variable in our experiment}
#'   \item{Condition}{Experimental condition of the experiment}
#'   ...
#' }
#' @source \url{https://osf.io/64njs/}
"alb5_rev"
