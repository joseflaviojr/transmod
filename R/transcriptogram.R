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

#' Normalizes gene expression profiles to transcriptogram curves.
#' The rows of the matrix (genes) must be in the correct order (seriation).
#' 
#' @param expressions Seriated gene expressions: row = gene, column = sample and element = expression level.
#' @param window Window size for calculating the mean of adjacent expression.
#' @return Transcriptograms matrix.
#' @seealso See example of use in the site \url{https://github.com/joseflaviojr/transmod}
#' @family transcriptogram
#' @export
transcriptogram <- function( expressions, window=251 ){

    genes   <- nrow(expressions)
    samples <- ncol(expressions)
    half    <- floor( window / 2 )
    tgram   <- expressions

    for( gene in 1:genes ){
        wi <- gene - half
        wf <- wi + window - 1
        if( wi < 1 ) wi <- 1
        if( wf > genes ) wf <- genes
        for( sample in 1:samples ){
            tgram[gene,sample] <- mean(expressions[wi:wf,sample], na.rm=TRUE)
        }
    }

    tgram

}

#-------------------------------------------------------

#' Normalizes gene expression profile to transcriptogram curve.
#' The elements of the list (genes) must be in the correct order (seriation).
#' 
#' @param expression Seriated gene expression list: element = expression level.
#' @param window Window size for calculating the mean of adjacent expression.
#' @return Transcriptogram.
#' @seealso See example of use in the site \url{https://github.com/joseflaviojr/transmod}
#' @family transcriptogram
#' @export
transcriptogram_series <- function( expression, window=251 ){

    total <- length(expression)
    half  <- floor( window / 2 )
    tgram <- expression
    
    for( i in 1:total ){
        wi <- i - half
        wf <- wi + window - 1
        if( wi < 1 ) wi <- 1
        if( wf > total ) wf <- total
        tgram[i] <- mean(expression[wi:wf], na.rm=TRUE)
    }

    tgram

}

#-------------------------------------------------------