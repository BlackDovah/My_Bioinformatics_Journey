pollutantmean <- function(directory, pollutant, id = 1:332) {
        mainDir = getwd()
        if (getwd() != "/home/blackdovah/Work and Education/Python39/Bioinformatics/GenomicDataScienceSpecialization/RStuff/specdata") {
                setwd(directory)
        }
        x = c()
        for (i in dir()[id]){
                x =  c(x,read.csv(i)[[pollutant]][!is.na(read.csv(i)[[pollutant]])])
        }
        setwd(mainDir)
        return(mean(x))
}

