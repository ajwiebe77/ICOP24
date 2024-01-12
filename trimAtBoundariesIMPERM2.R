### trimAtBoundariesIMPERM2.R
### Andrew J. Wiebe, 14 Dec 2023
#
# Remove points outside aquifer for third side of wedge-shaped/triangular-shaped aquifer that is impermeable.
# 
# Code modified from the script trimAtBoundaries.R at https://github.com/ajwiebe77/theYWVC (Wiebe and McKenzie, 2022)
#
# Revised on 3 Jan 2023 to remove points above y = x * tan(alpha) line for case 2 option 1; revised further on 4 Jan 2023 to restrict the 3 Jan 2023 changes for options 1 to 3
# Revised on 13 Jan 2023 to account for straight line boundary for impermeable boundary case
# Revised on 14 Dec 2023 to account for dimensioned data with x-axis length distance
#
# Parameters:
#     list1 - a list of sets of points
#     case - the case from Nagheli et al. (2020)
#     option - the option from Nagheli et al. (2020) demonstration examples, needed for triangular aquifers
#     alpha - the angle of the triangular aquifer
#     axisLen - the length of the x-axis of the aquifer
#     xy - boolean value denoting whether the x and y vectors in the list are labelled or not
#     mirrorIMPERMEABLE - boolean value denoting whether the third side of the aquifer is impermeable (TRUE) or not (FALSE)
# 
# Return values:
#    list2 <- a list of the sets of points that are within the aquifer
#
# References:
#
#    Nagheli, S., Samani, N., and Barry, D.A., 2020. Capture zone models of a multi-well system in aquifers bounded with regular and irregular inflow boundaries. J. Hydrol. X 7, 100053. doi:10.1016/j.hydroa.2020.100053.
#    Wiebe, A.J., and McKenzie, J.M., 2022. An Open-Source Web Tool for Visualizing Estimates of Well Capture Zones Near Surface Water Features. Poster presentation at: AGU Fall Meeting 2022, Chicago, IL, USA, 12-16 Dec 2022. doi:10.22541/essoar.167267811.10671930/v1.
#

trimAtBoundariesIMPERM2 <- function(list1, case, option, alpha, axisLen, xy, mirrorIMPERMEABLE){
  list2 <- list1 
  
  if(case == 1){
    # ignore this case - use original trimAtBoundaries.R function
    
  }else if(case == 2 && (option == 1 || option == 2 || option == 3 || option == 4)){
    
    for(i in 1:length(list1)){
      if(xy == TRUE){
        for(ii in 1:length(list1[[i]]$x)){ # for each point
          if( ! is.na(list1[[i]]$x[ii])){
            if(list1[[i]]$x[ii] > 0){
              if(list1[[i]]$x[ii] > axisLen || list1[[i]]$y[ii] < 0 || (mirrorIMPERMEABLE == TRUE && list1[[i]]$x[ii] <= axisLen * cos(alpha) && list1[[i]]$y[ii] > list1[[i]]$x[ii] * tan(alpha)) || (mirrorIMPERMEABLE == TRUE && list1[[i]]$x[ii] > axisLen * cos(alpha) && list1[[i]]$y[ii] > (list1[[i]]$x[ii] * tan(pi/2 - alpha/2))) || (mirrorIMPERMEABLE == TRUE && list1[[i]]$x[ii] > axisLen * cos(alpha) && list1[[i]]$x[ii] > (axisLen - list1[[i]]$y[ii] * cos(pi/2 - alpha/2)/sin(pi/2 - alpha/2)))){ # revised 13 Jan 2023
                list2[[i]]$x[ii] <- NaN
                list2[[i]]$y[ii] <- NaN
              }
            }else if(list1[[i]]$x[ii] < 0){ # x < 0
              if(list1[[i]]$x[ii] < axisLen * cos(alpha) || list1[[i]]$y[ii] < list1[[i]]$x[ii] * axisLen * tan(alpha) || list1[[i]]$y[ii] > axisLen * sin(acos(list1[[i]]$x[ii]))){
                list2[[i]]$x[ii] <- NaN
                list2[[i]]$y[ii] <- NaN
              }
            }else if(list1[[i]]$x[ii] == 0){
              if(list1[[i]]$y[ii] > 1 || list1[[i]]$y[ii] < 0 || (option < 3 && list1[[i]]$y[ii] > 0)){ # revised 4 Jan 2023
                list2[[i]]$x[ii] <- NaN
                list2[[i]]$y[ii] <- NaN
              }
            }
          }
        }
      }else{
        # ignore this case - use original trimAtBoundaries.R function
      }
    }
  }
  
  
  return(list2)
}