\name{rmvData}
\alias{rmvData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rmvData
}
\description{
Filters out all rows that contain strings present strings_to_remove in their 'Sample' or 'Target column
}
\usage{
rmvData(input_file,
        output_file,
        strings_to_remove)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_file}{
  Original .csv or .txt file that you want to filter rows out of.
}
  \item{output_file}{
  Desired name of the resulting .csv or .txt file.
}
  \item{strings_to_remove}{
  Input as single string or vector of strings. Defines the strings you want to remove rows by.
}
}
\details{
  All arguments are required.
}
\value{
  Outputs a .txt or .csv file depending on what was input. Same kind of file is output as input.
}
\author{
  Tristan Tyler Reed
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Example 1
data_in <- read.delim('myresultsinput.txt',StringsasFactors = FALSE)
head(data_in)
  Sample Target                 Cq
1   Sample1  MKI67 26.768350048071987
2   Sample1   RPS9 22.989342495730817
3   Sample1  MKI67 26.644647080591646
4   Sample1   RPS9 22.853204123658813
5   Sample2  MKI67  26.01826927010041
6   Sample2   RPS9 23.544898960731118

rmvData('myresultsinput.txt','myresultsoutput.txt','Sample1')

data_out <- read.delim('myresultsoutput.txt',StringsasFactors = FALSE)
head(data)
  Sample Target                 Cq
1   Sample2  MKI67  26.01826927010041
2   Sample2   RPS9 23.544898960731118

# Example 2
data_in <- read.delim('myresultsinput.txt',StringsasFactors = FALSE)
head(data_in)
  Sample Target                 Cq
1   Sample1  MKI67 26.768350048071987
2   Sample1   RPS9 22.989342495730817
3   Sample1  MKI67 26.644647080591646
4   Sample1   RPS9 22.853204123658813
5   Sample2  MKI67  26.01826927010041
6   Sample2   RPS9 23.544898960731118

rmvData('myresultsinput.txt','myresultsoutput.txt','MKI67')

data_out <- read.delim('myresultsoutput.txt',StringsasFactors = FALSE)
head(data)
  Sample Target                 Cq
1   Sample1   RPS9 22.989342495730817
2   Sample1   RPS9 22.853204123658813
3   Sample2   RPS9 23.544898960731118
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
