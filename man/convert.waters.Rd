\name{convert.waters}
\alias{convert.waters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion of Waters' .raw data.
}
\description{
This functions converts waters' .raw files to mzData. This function is slow and only exists because other conversion software does not currently correctly remove lockmass scans (masswolf) or does not give calibrated data (ProteoWizard).
The function uses masswolf which needs to be in path. MassLynx also needs to be installed.
This function does the conversion and subsequently removes the lockscan masses that are intertwined in the data from masswolf.
}
\usage{
convert.waters(indir,outdir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indir}{Folder with waters .raw files.}
  \item{outdir}{Folder where the mzData files should be written.}
}

\value{None}

\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
\note{
This should probably be changed to accept files instead of folders. But then it would also be needed to be able to specify an output folder for each file.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
## Convert all files in a single folder
convert.waters(indir,outdir)
}




\dontrun{
## Convert files in parallel
library(chemhelper)
library(stringr)
library(parallel)


#list of folders with files to convert
indir_list=c("Standards_May2013/Standards4/Data",
             "Standards_Sept-12/Data"
        )

indir=paste("D:/root_path/raw/",indir_list,"/",sep="")

# Now we make a list of matching output folders
outdir=paste("D:/root_path/Converted/",indir_list,"/",sep="")
outdir=str_replace(outdir,"/Data","")
outdir=str_replace(outdir,"/data","")


# combine indir and outdir in a list so it can be passed to apply functions
arg_list=list()
for (i in 1:length(indir)){
arg_list[[i]]=c(indir[i],outdir[i])
}


# Do parallel conversion
cl <- makeCluster(getOption("cl.cores", 2))
clusterExport(cl, "convert.waters")
parLapplyLB(cl,arg_list,function(x) convert.waters(x[1],x[2]))
stopCluster(cl)
}
}
