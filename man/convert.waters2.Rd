\name{convert.waters2}
\alias{convert.waters2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion of Waters' .raw data.
}
\description{
This functions converts waters' .raw files to mzXML. This function is a wrapper for masswolf that converts each "function" one at a time to a seperate mzXML file.
The function uses masswolf which needs to be in path. MassLynx also needs to be installed.
}
\usage{
convert.waters2(infiles,outdir,funcs=c(1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infiles}{.raw files (which are actually folders).}
  \item{outdir}{Folder where the mzData files should be written.}
  \item{funcs}{Vector containing the number of the functions to convert.}
}
\value{None}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
