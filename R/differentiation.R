#-------------------------------------------------------

#  Copyright (C) 2017 José Flávio de Souza Dias Júnior
#  
#  This file is part of "transmod" - <https://github.com/joseflaviojr/transmod>.
#  
#  "transmod" is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  "transmod" is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#  
#  You should have received a copy of the GNU Lesser General Public License
#  along with "transmod". If not, see <http://www.gnu.org/licenses/>.

#-------------------------------------------------------

#' Performs statistical test to determine the differentiation level between two sets of gene expression samples.
#' Performs univariate (rowwise) Welch tests.
#' 
#' @param expressions Gene expression matrix: row = gene, column = sample and element = expression level.
#' @param control_group Column indexes (samples) that are part of the differentiation control group.
#' @param target_group Column indexes (samples) that are part of the differentiation target group.
#' @return List of gene differentiation levels, in the order of the lines of the input matrix.
#' @seealso See example of use in the site \url{https://github.com/joseflaviojr/transmod}
#' @family differentiation
#' @export
differentiate <- function( expressions, control_group, target_group ){

    if( ! requireNamespace("GeneSelector", quietly=TRUE) ){
        stop("Package required: GeneSelector", call.=FALSE)
    }

    expressions <- expressions[,c(control_group,target_group)]

    classes <- c(rep(1,length(control_group)), rep(2,length(target_group)))
    ranking <- GeneSelector::RankingWelchT(as.matrix(expressions), classes, type="unpaired")
    ranking <- GeneSelector::toplist(ranking, top=nrow(expressions))
    
    ranking <- abs(ranking[order(ranking$index),]$statistic)
    ranking[is.na(ranking)] <- 0

    ranking

}

#-------------------------------------------------------