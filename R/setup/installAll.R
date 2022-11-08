#============================================================================
# Script to install the following required packages to the default directory
# - scatterplot3d
# - ks
# - matlab
# - numDeriv
# - nlme
# - KernSmooth
# Requires internet connection.
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
# List of packages to install
packages = c("scatterplot3d",
             "ks", 	"matlab", 
             	"numDeriv", "nlme", 
             	"KernSmooth")

# Path to the default personal library
Library = .libPaths()[1]
cat('* Installing all packages to', Library, '*(Internet connection required)\n')
ans <- readline(">>> Please press <Enter> to continue (q to quit): ")
if(ans == 'q') {
  stop(">>> Installation cancelled by user.") 
} else {
  for(i in seq(packages)){
    cat("\n >>>>> Installing package ", packages[i], '\n')
    install.packages(packages[i], lib=Library)
  }
  cat("\n\n >>>>>> Done.\n\n")
}


