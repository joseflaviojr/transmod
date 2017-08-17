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

#' Clustering of adjacent values in a series based on peaks and valleys detection.
#' Module: series interval with level elevation.
#' 
#' @param series Series of proportional values/levels that have gradual elevations.
#' @param height_min Minimum peak height (percentage in relation to the highest value).
#' @return List that numerically indicates in which module each element of the series is.
#' @seealso See example of use in the site \url{https://github.com/joseflaviojr/transmod}
#' @family modularity
#' @export
modularize <- function( series, height_min=0.05 ){

    tops <- order(-series)
    total <- length(series)
    modules <- rep(0,total)
    height_min <- series[tops[1]] * height_min

    for( t in tops ){

        if( modules[t] != 0 ) next

        # Begin
        
        begin <- t

        if( t > 1 ){
            
            smaller <- t
            smaller_level <- series[smaller]

            for( i in (t-1):1 ){

                if( modules[i] != 0 ) break

                level <- series[i]

                if( level <= smaller_level ){
                    smaller <- i
                    smaller_level <- series[smaller]
                    begin <- smaller
                }else if( ( level - smaller_level ) >= height_min ){
                    break
                }

            }

        }

        # End
        
        end <- t

        if( t < total ){
            
            smaller <- t
            smaller_level <- series[smaller]

            for( i in (t+1):total ){

                if( modules[i] != 0 ) break

                level <- series[i]

                if( level <= smaller_level ){
                    smaller <- i
                    smaller_level <- series[smaller]
                    end <- smaller
                }else if( ( level - smaller_level ) >= height_min ){
                    break
                }

            }

        }

        # Module

        for( i in begin:end ){
            modules[i] <- t
        }

    }

    # Final adjustment: correct sequence and no gaps

    id <- 1
    current <- modules[1]
    modules[1] <- id

    for( i in 2:total ){
        if( modules[i] == current ){
            modules[i] <- id
        }else{
            id <- id + 1
            current <- modules[i]
            modules[i] <- id
        }
    }

    modules

}

#-------------------------------------------------------

#' Summarizes the main features of each module.
#' 
#' @param series Series of proportional values/levels that have gradual elevations.
#' @param modules List that numerically indicates in which module each element of the series is.
#' @return Summary matrix.
#' @seealso See example of use in the site \url{https://github.com/joseflaviojr/transmod}
#' @family modularity
#' @export
summarize_modules <- function( series, modules ){

    m <- matrix(0,0,7)
    colnames(m) <- c("module","begin","end","size","mean","max","min")

    current <- modules[1]
    begin   <- 1
    end     <- 1

    for( i in 2:length(series) ){
        if( modules[i] == current ){
            end <- i
        }else{
            module <- series[begin:end]
            m <- rbind(
                m,
                c(
                    current,
                    begin,
                    end,
                    length(module),
                    mean(module, na.rm=TRUE),
                    max(module, na.rm=TRUE),
                    min(module, na.rm=TRUE)
                )
            )
            current <- modules[i]
            begin   <- i
            end     <- i
        }
    }

    module <- series[begin:end]
    m <- rbind(
        m,
        c(
            current,
            begin,
            end,
            length(module),
            mean(module, na.rm=TRUE),
            max(module, na.rm=TRUE),
            min(module, na.rm=TRUE)
        )
    )

    as.data.frame(m)

}

#-------------------------------------------------------

#' Selects elements of a series according to modular levels.
#' The modules are classified according to the peak level (the higher the better).
#' 
#' @param series Series of proportional values/levels that have gradual elevations.
#' @param modules List that numerically indicates in which module each element of the series is.
#' @param select Number of elements to select.
#' @param inspect Number of modules to consider. 0 = all modules.
#' @return Index (from 1) of each selected element.
#' @seealso See example of use in the site \url{https://github.com/joseflaviojr/transmod}
#' @family modularity
#' @export
select_from_modules <- function( series, modules, select=100, inspect=0 ){
    
    sel  <- c()

    m <- summarize_modules(series, modules)
    m <- m[order(-m$max),]

    total <- 0
    lines <- 0

    if( inspect == 0 ){
        total <- length(series)
        lines <- nrow(m)
    }else{
        total <- sum(m$size[1:inspect])
        lines <- inspect
    }

    perc <- select / total
    
    for( i in 1:lines ){
        if( select == 0 ) break
        mi <- m[i,]
        ms <- order( - series[mi$begin:mi$end] ) + ( mi$begin - 1 )
        mt <- ceiling( mi$size * perc )
        if( mt > select ) mt <- select
        select <- select - mt
        sel <- c(sel, ms[1:mt])
    }

    sel

}

#-------------------------------------------------------