# default color for displaying XCO2, DW

def.col <- function(){
  return(c('black', '#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8',
           '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))
}

ggdef.col <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

def.colb = function() {
  c("#000000", "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}