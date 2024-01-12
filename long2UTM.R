# O'Brien, J., 2012.
# https://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude
#
# Unmodified from version coded at https://github.com/ajwiebe77/theYWVC (Wiebe and McKenzie, 2022)
#
# References:
#    
#    O'Brien, 2012. Response to "Determining UTM zone (to convert) from longitude/latitude." https://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude. Cited 19 Dec 2022.
#    Wiebe, A.J., and McKenzie, J.M., 2022. An Open-Source Web Tool for Visualizing Estimates of Well Capture Zones Near Surface Water Features. Poster presentation at: AGU Fall Meeting 2022, Chicago, IL, USA, 12-16 Dec 2022. http://dx.doi.org/10.22541/essoar.167267811.10671930/v1.
#    
#
long2UTM <- function(long) {
  (floor((long + 180)/6) %% 60) + 1
}
