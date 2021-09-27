## The packages required to run this script can be installed by uncommenting the appropriate lines of package installation in the code (lines 30 and 33). 
## Packages requiring installation are the following:
## - RColorBrewer
## - cmocean.
## The rest of the script requires only standard R functions [pdf(),smoothScatter(), cor(), text(), read.csv()]. 
## The code uses the data from the folder "survey_model_matrix_comparisons" provided through this repository (https://github.com/mobs-lab/mixing-patterns/tree/main/analysis_results/survey_model_matrix_comparisons). 
## The folder should be added to your working directory, that also can be set modifying the first line after this comment [setwd(".")]. 
## Resulting figure will be created in an automatically generated folder "fig" inside "survey_model_matrix_comparisons" folder.

setwd(".")
dir.create("./survey_model_matrix_comparisons/fig")

all.m.file = list.files("./survey_model_matrix_comparisons", pattern = "([A-Z])([a-z]*)(.csv)", full.names = T)
# check that only necessary files with the following pattern are added to the folder and for each survey a unique model file exists
survey.m.file = all.m.file[grep("(/Survey_)", all.m.file)]
model.m.file = all.m.file[grep("(/Model_)", all.m.file)]

all.files = c(survey.m.file, model.m.file)
print(matrix(all.files, ncol = 2))

m.all = list()
matr.norm = list()
for(i in 1:length(all.files)){
  m.all[[i]] = read.csv(all.files[i], sep = " ", header = F)
  matr.norm[[i]] = m.all[[i]]/sum(m.all[[i]])
}

## Panel A
n=1000
#install.packages("RColorBrewer")
library(RColorBrewer)
k = 100
#install.packages("cmocean")
library(cmocean)
my.cols = cmocean('deep')(k)[k:1]


for(i in 1:length(survey.m.file)){
  countrydata1 = as.numeric(as.matrix(m.all[[i]]))
  idx.0 = which(countrydata1 == 0)
  countrydata2 = as.numeric(as.matrix(m.all[[i+3]]))
  mydata2 = as.data.frame(cbind(log10(countrydata1/sum(countrydata1)), 
                                log10(countrydata2/sum(countrydata2))))
  colnames(mydata2) = c("Survey","Synthetic")
  if(length(idx.0) > 0){
    countrydata1[idx.0] = 0.000001
  }
  ro = cor(log10(countrydata1/sum(countrydata1)), log10(countrydata2/sum(countrydata2)))
  ro0 = cor((countrydata1/sum(countrydata1)), (countrydata2/sum(countrydata2)))
  print(cor.test((countrydata1/sum(countrydata1)), (countrydata2/sum(countrydata2))))
  fullname = unlist(strsplit(all.files[i], "_"))
  pdf(paste0("./survey_model_matrix_comparisons/fig/fig3A_",
             gsub(".csv","",fullname[length(fullname)]),".pdf"),
      width = 6, height = 6)
  par(mar = c(5,4,4,5) + .1)
  #quartz()
  smoothScatter(mydata2, nrpoints=n, colramp=colorRampPalette(my.cols), pch=19, cex=.2,
                nbin = 128, bandwidth = 0.08, col = "grey30")
  text(min(mydata2[which(is.finite(mydata2[,1])),1]),max(mydata2[,2])-0.1, pos=4,
       paste("Pearson's r:", round(ro0,2)),
       col="white")
  dev.off()
}