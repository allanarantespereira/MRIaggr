\name{getFileExtention}
\title{Get file extension}
\alias{getFileExtention}
\description{
Get file extension
}
\usage{
getFileExtention(filename)
}
\arguments{
  \item{filename}{the name of the file. \emph{character}. REQUIRED.}
  }
\value{
  A \emph{list} containing :
  \itemize{
  \item\code{[[filename]]} : the filename without extension
  \item\code{[[extension]]} : the filename extension (\code{NA} if no extension)
  \item\code{[[extension.compression]]} : was the file a gzip
  }}
\keyword{functions,internal}

