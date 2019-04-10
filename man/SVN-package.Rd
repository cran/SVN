\name{SVN-package}
\alias{SVN-package}
\alias{SVN}
\docType{package}
\title{
  Statistically validated networks
}
\description{
Statistically validated networks are built from the states of nodes and from their curious (statistically speaking) synchronization: two nodes are linked if their states are anomalously synchronized. The associated null hypothesis is that the states of nodes are Poissonian processes, i.e., independent from each other.
}
\details{
The main function SVN_links expects a matrix/data.frame/data.table of states, the line number being the index and while a given column corresponds to a given node, and  outputs the results network in an igraph format. A Multiple Hypothesis Testing correction is applied since one performs several (rather, many) tests. It is then up to the user to use some network clustering technique if needed. For the convenience of the user, the function SVN_clusters applies the infoMap or Louvain method to the output of SVN_links and returns the resulting clustering. For small enough timeseries, it may be useful to use the options(svn.memoise=TRUE) so as to avoid expensive calls to the phyper function.
}
\author{
Maintainer: Damen Challet <damien.challet@gmail.com>
}
\references{
Tumminello, M., Micciche, S., Lillo, F., Piilo, J., & Mantegna, R. N. (2011). Statistically validated networks in bipartite complex systems. PloS one, 6(3), e17994.
}
\keyword{ network, synchronization }
\seealso{
\pkg{igraph}
}
\examples{
 x=sample(c(1,0),1000,replace = TRUE)   # random vector of 0s and 1s
 xx=x; x[1]=1-x[1]   # one modifies just one element. x and xx are very similar
 x_rev=rev(x)        # x and x_rev are not
 y=sample(c(1,0),1000,replace = TRUE)   # y is another random vector of 0s and 1s
 M=cbind(x,xx,x_rev,y)                  # builds the matrix of states
 mylinks=SVN_links(M)
 print(mylinks)

 # one can compute clusters as well
 myclusters=SVN_clusters(mylinks)
}
