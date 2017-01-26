\name{UCM}
\alias{ASAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
Function to process and analize Finger Force data and to calculate Variance across the (Good) UCM and Orthogonal to (Bad) UCM. Analysis is made across subjects
Steps:
1. Load Ramp trials and computes Enslaving Matrix
2. Plot individual finger forces
3. Load Trials from text files
4. Filter data, with Butterworth 2nd order, low pass filter, 20Hz cutoff frequency by default.
5. Computes Total Force and Force modes based on Enslaving Matrix
6. UCM analysis
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
UCM(subjects=c("S01","S02"),wd=c("/Users/Trials/"),Task="ASA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subjects} Name of subjets or folder and file name structure to search for file and analyze data,
  \item{wd} Working directory, where Folders contain each file per subject. Use Finder or File explorer to set it (use slash to separate folders), and there HAS TO BE a slash(/) at the END
  \item{Task} Name of conditions or task evaluated or file name structure to search for file,
  \item{timeR} Time (s) of Ramp task for Regression Analysis for Enslaving Matrix calculation, by default between 5.5 and 7.5 s,
  \time{timeS} Begining and End of search time (in seconds) for ASAs. It should be specified by concatenate function (i.e c(1,20)). By default starts at steady state up to the end c(6,20),
  \item{freq} Sampling frequency of data acquisition. By default is set at 1000Hz,
  \item{cutfreq} Cut off frequency. By default is set at 20Hz,
  \item{order} Filter order for data processing.
  \item{J} Jacobian Matrix, by default is [1 1 1 1]

%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Names of subjects, timeR, timeS need to be concatenated (separately)
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Angelo Bartsch, MSc%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
ASA<-UCM(subjects=c("Paula","Angelo"),wd=c("/Users/Trials/"),Task="ASA",cutfreq = 20)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ UCM }% use one of  RShowDoc("KEYWORDS")
\keyword{ ASAs }% __ONLY ONE__ keyword per line
