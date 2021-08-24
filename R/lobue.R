#' Replication protocol dataframe of replication project ML5
#'
#' Text from my thesis:
#' The study conducted by Schwarz et al. (1985) investigated the effect of the
#' design of questionnaire-scales (scalesgroup) on the information reported by
#' participants (unintuitively labelled scales).
#' Speciﬁcally, they randomly assigned N = 132 participants either into
#' the low category or the high category group.
#' Participants were asked to rate how many hours of tv they watched per
#' day on a questionnaire.
#'
#' Full code book can be found here: https://osf.io/68f7m/
#'
#' @format A data frame with 473 rows and 7 variables:
#' \describe{
#'   \item{ID}{unique ID of participant}
#'   \item{Location}{Grouping variable for where the replication of the experiment took place. In this case just mTurk}
#'   \item{snake_experience}{Has person had previous experience with snake? (e.g. seen before)}
#'   \item{RT.correct}{Reaction time. Dependent variable in the experiment}
#'   \item{protocol}{Either revised protocol (RP) or replication protocol (NP)}
#'   \item{child_parent}{Is Person the child or the parent}
#'   \item{target_stimulus}{Indicates, whether snake (fear relevant) or decoy stimuli was shown}
#'   \item{Condition}{Experimental condition of the experiment}
#'   ...
#' }
#' @source \url{https://osf.io/z2fu9/}
"lobue"
