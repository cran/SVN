library(data.table)
library(igraph)
library(memoise)


#'@import data.table
#'@importFrom igraph graph.data.frame cluster_infomap cluster_louvain
#'@import utils
#'@import stats
#'@import memoise

link_pvalues=function(dt,adj,si,sj,alternative,link_id=""){
  ## colonnes doivent etre sorted
  m = ncol(adj)
  Ai = (adj == si)
  if(si==sj){
    Aj=Ai
  }else{
    Aj = (adj == sj)
  }
  nPQ = t(Ai) %*% (Aj)
  nP=colSums(Ai)
  nQ=colSums(Aj)
  nP = matrix(nP,m,m)
  nQ = matrix(nQ,m,m,byrow=TRUE)

  utri = upper.tri(nPQ, diag = FALSE)
  nPQ = nPQ[utri]
  nP = nP[utri]
  nQ = nQ[utri]
  i = colnames(adj)[row(utri)[utri]]
  j = colnames(adj)[col(utri)[utri]]

  overexpression= (alternative=="overexpression")

  if(getOption("svn.memoise")){
    phyper_m<-memoise(phyper)
    svn_phyper=phyper_m
  }else{
    if(exists("phyper_m")){
      forget(phyper_m)
      rm(phyper_m)
    }
    svn_phyper=phyper
  }
  pv=svn_phyper(nPQ-overexpression, nP, dt-nP, nQ, lower.tail=!overexpression)
  return(data.table(i=i,j=j,si=si,sj=sj,pv = pv))
}

all_links_pvalues = function(adj, links, dt, alternative){
  #
  si=NULL
  sj=NULL

  links = data.table(si=links[,1],sj=links[,2],link_id=links[,3]) ## transform links into a data.table

  #compute p values for all pairs of states
  links_pv = links[,link_pvalues(dt, adj,si,sj,alternative),by='link_id']
  return(links_pv)
}


.onLoad <- function(libname, pkgname){

  op<-options()
  op.svn <- list(svn.memoise=FALSE)
  toset <- !(names(op.svn) %in% names(op))
  if(any(toset)) options(op.svn[toset])
  invisible()
}

#' Statistical validated networks
#'
#' @param states_vs_t a matrix with time in lines and individual states in
#'   columns
#' @param alpha the family-wise error rate in the case of Bonferroni
#'   multiple-hypothesis correction, or the false discovery rate in the case of
#'   the FDR multiple hypothesis correction
#' @param MHT.correction the type of multiple hypothesis correction
#' @param exclude.states a vector of states to remove from states_vs_t
#' @param states.pair.types accepted value: "all", "same", "different": selects
#'   which kind of state pairs are tested: if there are two states e.g. (1,2),
#'   "all" allows (1,1), (2,2) and (1,2); "same" allows (1,1) and (2,2), while
#'   "different" allows (1,2)
#' @param alternative a string either equal to "overexpression" or "underexpression"
#'
#' @return a data.table object of the node pairs significantly correlated at the
#'   alpha level. The columns are: link_id: a unique identifier for links; i and j are the node names; si and sj are their states;  pv is the value associated to the link; threshold is the Multiple Hypothesis Testing-adjusted p-value threshold (only links with pv<=threshold are kept)
#'
##' @export
#'
#' @examples
#' x=sample(c(1,0),1000,replace = TRUE)
#' xx=x; x[1]=1-x[1]   # one modifies just one element
#' x_rev=rev(x)
#' y=sample(c(1,0),1000,replace = TRUE)
#' M=cbind(x,xx,x_rev,y)
#' mylinks=SVN_links(M)
#' print(mylinks)
#' # one can compute clusters as well
#' myclusters=SVN_clusters(mylinks)
SVN_links=function(states_vs_t, alpha=0.01, MHT.correction = "bonferroni", exclude.states=NULL, states.pair.types="all",alternative="overexpression"){
  # first, touch some variable name to please CRAN
  pv=NULL
  threshold=NULL
  si=NULL
  sj=NULL


  #if no colnames, use numbers
  if(length(colnames(states_vs_t))==0)
    colnames(states_vs_t)=1:ncol(states_vs_t)

  n=ncol(states_vs_t)
  dt=nrow(states_vs_t)

  allstates=unique(sort(states_vs_t))
  if(length(exclude.states)>0){
    allstates=allstates[!(allstates %in% exclude.states)]
  }
  #types of links
  if(states.pair.types=="all"){               # (1,1), (1,2) and (2,2)
    links=expand.grid(allstates, allstates)
  }else if(states.pair.types=="same"){        # (1,1) and (2,2)
    links=cbind(allstates,allstates)
  }else if(states.pair.types=="different"){   # (1,2) and (2,1)
    links=t(cbind(combn(allstates,2), combn(rev(allstates),2)))
  }

  #link type ids
  links=cbind(links, 'link_id'=1:nrow(links))

  numTests=nrow(links)*n*(n-1)/2

  #compute p-values for all links

  links=all_links_pvalues(states_vs_t, links, dt, alternative=alternative)

  if(MHT.correction=="bonferroni"){
    links=links[pv<alpha/numTests]
    links[,threshold:=alpha/numTests]
  }else{
    #validate links using fdr

    #order links by p-value and compare to thresholds
    ##links_fdr=links[order(links$pv),]
    setkeyv(links,cols= "pv")
    fdrt=alpha/length(links$pv)*seq_along(links$pv)

    links[,fdrt:=fdrt]

    #perform fdr
    if(tail(links$pv,1)>tail(links$fdrt,1)){
      links = links[pv<fdrt,]
    }
  }
  return(links)
}

#' Find clusters in Statistically Validated Networks
#'
#' @param links an object obtained from the SVN_links function
#' @param cluster.method a string, either "infomap" or "louvain", that selects the clustering method
#'
#' @return an igraph communities object
#' @export
#'
#' @examples  M=matrix(rbinom(200,size = 1,0.5),nrow=50)
#' mylinks=SVN_links(M)
#' # no links, then
#' print(mylinks)
#'
#' # another example
#' x=c(1,0,0,0,0,0,1,1,1,1,0,0)
#' xx=c(1,1,0,0,0,0,1,1,1,1,0,0)
#' xrev=rev(x)
#' w=sample(x,length(x))
#' M=rbind(x,xx,xrev,w)
#' mylinks=SVN_links(M)
#'
#' ## if one wants clusters as well
#' myclusters=SVN_clusters(mylinks)
#'
SVN_clusters=function(links,cluster.method="infomap"){
  i=NULL
  j=NULL
  mygraph=graph.data.frame(links[,list(i,j)],directed=FALSE)
  if(cluster.method=="infomap"){
    myclusters=cluster_infomap(mygraph,nb.trials = 10)
  }else if(cluster.method=="louvain"){
    myclusters=cluster_louvain(mygraph)
  }else{
    warning("SVN_cluster: unknown method ",cluster.method,", using infomap")
    myclusters=cluster_infomap(mygraph,nb.trials = 10)
  }
  return(myclusters)
}
