#' Initialises rsimpop
#'
#' @param seed integer - a negative number uses current time to 100th second resoultion.
#' @param b boolean - whether to force re-initiliasation to specified seed.
#'
#' @return NULL
#' @export
#' @useDynLib rsimpop initRSimPop
#' @examples
#' initSimPop(-1)
initSimPop=function(seed,bForce=FALSE){
  ##browser()
  if(exists("simpop_seed") && !bForce){
    cat("Already initialised\n")
    return(NULL)
  }
  if(seed<0){
    ##take the last 30 or so bits (use 100th second resolution)
    seed=as.integer(substr(sprintf("%-15.0f",100000*(as.numeric(Sys.time()))),start = 7,stop=15))
  }
  set.seed(seed)
  res=.C("initRSimPop",
         edges=as.integer(seed))
  simpop_seed<<-seed
  NULL
}
#' @title 
#' Simulates Evolution of Cell Population using a Continuous Time Birth-Death Model
#'
#' This is the main entry function for carrying out simulations.  The parameters (e.g. rates) are fixed for the duration of the simulation.
#' Note: the implementation uses .C - so the actual data is copied on each invocation.
#' @param tree Initial population of cells together with ancestry as an APE phylo object.  Null indicates zygote.
#' @param param List. Parameter overrides.. list(pop_size=5e4,n_sim_days=365*40,divide_rate=1/120.0,mut_per_div_rate=5.0,driver_dt=1e6,driver_mito=0.2,driver_muta=0,b_run_wright_fisher=1,b_stop_at_pop_size=0,tau_leaping_prop=0)
#' @param cfg List. Specifies the compartments together with populations size etc.
#' @param b_verbose Boolean.
#' @param b_check Boolean. Check output is valid APE phylo
#' @return simpop object
#' @export
#' @useDynLib rsimpop sim_pop2
#' @examples
#' cfg=getDefaultConfig(target_pop_size,rate=initial_division_rate,ndriver=1,basefit = fitness)
#' params=list(n_sim_days=nyears_driver_acquisition*365,            b_stop_at_pop_size=1,b_stop_if_empty=0)
#' growthphase=sim_pop(NULL,params=params,cfg)
sim_pop=function(tree,
                 params=list(),
                 cfg=list(compartment=data.frame(val=c(0,1),rate=c(-1,1/120.0),popsize=c(1,1e4)),
                          info=data.frame(val=c(0,1,1),fitness=c(-1,0,0.2))),
                 b_verbose=TRUE, 
                 b_check=FALSE)
{
  toCombine=TRUE
  if(is.null(tree)){
    toCombine=FALSE
    tree=rtree(2)
    tree$edge.length=rep(0,2)
    tree$events=data.frame(value=0:1,driverid=c(0,0),node=1:2,ts=0.0,uid=0:1)
  }
  # print(tree$events)
  tree=standardiseTree(tree)  #populate the tree object with the #divisions (ndivs), the time (tBirth) and the maximum time (maxt)
  ndrivers=tree$ndrivers
  time=tree$tBirth
  ##maxt=tree$maxt
  defaultparams=list(n_sim_days=365*40,
                     b_stop_if_empty=0,
                     b_stop_at_pop_size=0,
                     maxt=tree$maxt,
                     driver_rate_per_cell_per_day=0
  )
  fields=names(defaultparams)
  if(length(setdiff(names(params),names(defaultparams)))>0){
    badparam=setdiff(names(params),names(defaultparams))
    stop(sprintf("Unsupported params:%s",paste(badparam,collapse=",")))
  }
  for(par in names(params)){
    defaultparams[[par]]=params[[par]]
  }
  params=defaultparams
  if(b_verbose){
    for(par in names(params)){
      cat(sprintf("%s: %s\n",par,params[[par]]))
    }
  }
  
  edges=tree$edge
  nmuts=tree$edge.length
  nedge=dim(tree$edge)[1]
  if(is.null(tree$ndivs)){
    ndivs=rep(0,nedge)
  }else{
    ndivs=tree$ndivs
  }
  if(is.null(tree$ndrivers)){
    ndrivers=rep(0,nedge)
  }else{
    ndrivers=tree$ndrivers
  }
  
  # if(is.null(tree$is_combined))
  #   tree$is_combined=FALSE
  
  
  if(params[["n_sim_days"]]<=params[["maxt"]]){
    browser()
    stop("n_sim_days must be greater than max(timestamp)")
  }
  
  ntips=length(tree$tip.label)
  nevents=dim(tree$events)[1]
  
  compartment=cfg$compartment
  compartmentinfo=cfg$info
  idxd=grep("^d",colnames(compartmentinfo))
  ndriver=length(idxd)
  if(ndriver==0){
    cat("No driver columns specified\n")
    dsm=c(0)
  }else{
    dsm=as.matrix(compartmentinfo[,idxd])
  }
  # if(ndriver>1000){
  #   stop("Too many drivers: edit CellSimulation::setCompartmentInfoRecursively in C++ to support more!")
  # }
  totalpop=max(sum(compartment$popsize),length(tree$tip.label))
  MAX_SIZE=3*totalpop ##Allows for stochastic drift in population size (need to make this robust)
  
  MAX_EVENTS=ceiling(2*params[["n_sim_days"]])
  
  if(b_verbose){
    cat("MAX_EVENTS=",MAX_EVENTS,"\n")
    cat("MAX_SIZE=",MAX_SIZE,"\n")
    ## browser()
  }
  res=.C("sim_pop2",
         edges=as.integer(edges),
         ndivs=as.integer(ndivs),
         tBirth=as.double(time),
         ntips=as.integer(ntips),
         nedge=as.integer(nedge),###
         eventnode=as.integer(tree$events$node),
         eventval=as.integer(tree$events$val),
         eventdriverid=as.integer(tree$events$driverid),
         eventts=as.double(tree$events$ts),
         eventuid=as.integer(tree$events$uid),
         nevents=as.integer(nevents),###
         compartmentval=as.integer(compartment$val),
         compartmentrate=as.double(compartment$rate),
         compartmentsize=as.integer(compartment$popsize),
         ncompartment=as.integer(length(compartment$val)),
         compinfoval=as.integer(compartmentinfo$val),
         compinfofitness=as.double(compartmentinfo$fitness),
         compinfodriverstatus=as.integer(compartmentinfo$id),#as.integer(dsm),
         compinfondriver=as.integer(ndriver),
         ncomp=as.integer(length(compartmentinfo$val)),
         params=as.double(params),
         nparams=as.integer(length(params)),##
         max_size=as.integer(MAX_SIZE),
         max_events=as.integer(MAX_EVENTS),
         edgesOut=integer(2*MAX_SIZE),
         nDivsOut=integer(MAX_SIZE),
         stateOut=integer(MAX_SIZE),
         driverIDOut=integer(MAX_SIZE),
         tBirthOut=double(MAX_SIZE),
         nedgeOut=integer(1),
         nInternalNodeOut=integer(1),
         eventvalOut=integer(nevents),
         eventdriverOut=integer(nevents),
         eventtsOut=double(nevents),
         eventuidOut=integer(nevents),
         eventnodeOut=integer(nevents),
         neventOut=integer(1),
         eventTimestampOut=double(MAX_EVENTS),
         nPopSizeOut=integer(MAX_EVENTS),
         nDriverOut=integer(MAX_EVENTS),
         nEventsCount=integer(1),
         ntipsOut=integer(1),
         nCompPops=integer(length(compartmentinfo$val)),
         status=integer(1))
  if(res$status<0){
    stop("Error in C(.simpop2)")
  }
  nedge=res$nedgeOut
  cfg$info$population=res$nCompPops
  # print(neventOut)
  # print(eventuidOut)
  ##nEventsCount=res$nEventsCount
  # print(paste0('nEventsCount=',res$nEventsCount))
  res=list(edge=matrix(res$edgesOut,ncol=2)[1:nedge,],
           edge.length=res$nDivsOut[1:nedge],
           state=res$stateOut[1:nedge],
           driverid=res$driverIDOut[1:nedge],
           Nnode=res$nInternalNodeOut,
           ndivs=res$nDivsOut[1:nedge],
           tBirth=res$tBirthOut[1:nedge],
           nedge=res$nedgeOut,
           ntips=res$ntipsOut,
           timestamp=res$eventTimestampOut[1:res$nEventsCount],
           pop.size=res$nPopSizeOut[1:res$nEventsCount],
           totaldrivercount=res$nDriverOut[1:res$nEventsCount],
           ncellswithdrivers=1,
           params=params,
           events=data.frame(value=res$eventvalOut[1:res$neventOut],
                             driverid=res$eventdriverOut[1:res$neventOut],
                             node=res$eventnodeOut[1:res$neventOut],
                             ts=res$eventtsOut[1:res$neventOut],
                             uid=res$eventuidOut[1:res$neventOut]),
           cfg=cfg,
           status=res$status,
           currentEventUID=max(res$eventuidOut[1:res$neventOut])
  )
  
  class(res)="simpop"
  res=get_tree_from_simpop(res,b_check)
  if(toCombine)
  {
    res = combine_simpops(tree,res)
    # message("the 2 objects has been combined")
    # res$is_combined=FALSE
    # tree$is_combine=TRUE
  }
  res
}


#' Plot simpop population trajectory
#' @param simpop integer - a negative number uses current time to 100th second resoultion.
#'
#' @return NULL
#' @export
plot.simpop=function(simpop,...){
  T=length(simpop$timestamp)
  if(T>1e4){
    idx=seq(1,T,by = T %/% 1e4)
  }else{
    idx=1:T
  }
  if(max(simpop$timestamp)>365*2){
    plot(simpop$timestamp[idx]/365,simpop$pop.size[idx],log="y",xlab="Time(Years)",ylab="Pop Size",type="l",...)
  }else{
    plot(simpop$timestamp[idx],simpop$pop.size[idx],log="y",xlab="Time(Days)",ylab="Pop Size",type="l",...)
  }
}


standardiseTree=function(tree)
{
  if(is.null(tree$events)){
    ##Add events..
    stop("add events!")
  }
  
  nedge=dim(tree$edge)[1]
  if(is.null(tree$ndivs)){
    tree$ndivs=rep(0,nedge)
  }
  if(is.null(tree$time)){
    tree$tBirth=rep(0,nedge)
    tree$maxt=0
  }else{
    tree$maxt=max(tree$time)
  }
  tree
}

#' Converts simpop into an augmented APE phylo object of extant cells
#' @param sim input simpop simulation
#' @return ape tree
#' @export
get_tree_from_simpop=function(sim, # Gets APE ancestral tree of extant cells.
                              b_check=FALSE)
{
  if(!("simpop" %in% class(sim))){
    stop("get_tree_from_simpop: Input must be simpop object")
  }
  ttree=sim
  ttree$tip.label=sprintf("s%d",1:sim$ntips)
  class(ttree)=c("simpop","phylo")
  if(b_check)   #from here, it is supposed to check if the ttree object has the strucuture of a phylo object. When there is a problem with the structure of the object, it will return FATAL or MODERATE.
  {
    tmpf=tempfile() #creating a temporal file where the outputs of the checking process will be printed
    sink(tmpf)
    try(checkValidPhylo(ttree))  #this is the ape function that does the checking
    sink()
    chk=readLines(tmpf)
    if(any(grepl("FATAL",chk)) || any(grepl("MODERATE",chk))){
      ##browser()
      stop("invalid phylo!")
    }
  }
  ttree
}

#' Wrapper function for sub-sampling a simpop object
#' @param tree simpop object.
#' @param tips character vector - names of tips to include in sub-sample
#' @return simpop
#' @useDynLib rsimpop sub_sample
C_subsample_pop=function(tree,tips){
  if( ! tree$edge[which(tree$state==0),2] %in% tips){
    stop("tip list must contain outgroup")
  }
  ##int * edges, int * nmuts, int * ndrivers, int * nmitosis,int * ntips,int * nedge,int *edgesOut,int * nmutsOut,int * ndriversOut,int * nedgeOut
  tree=standardiseTree(tree)
  time=tree$tBirth
  edges=tree$edge
  ndivs=tree$ndivs
  nedge=dim(tree$edge)[1]
  maxt=tree$maxt
  ntips=length(tree$tip.label)### Get rid of tip label for actual sims...
  MAX_SIZE=2*nedge ## Should just be nedge
  if(is.null(tree$cfg)){
    stop("Must include cfg field to subsample tree")
  }
  cfg=tree$cfg
  nevents=length(tree$events$node)
  idxd=grep("^d",colnames(cfg$info))
  ndriver=length(idxd)
  if(ndriver<1){
    stop("Please specify driver columns in compartmentinfo")
  }
  
  res=.C("sub_sample",
         as.integer(edges),
         as.integer(ndivs),
         as.double(tree$tBirth),
         as.integer(ntips),
         as.integer(nedge),
         eventnode=as.integer(tree$events$node),
         eventval=as.integer(tree$events$val),
         eventdriverid=as.integer(tree$events$driverid),
         eventts=as.double(tree$events$ts),
         eventuid=as.integer(tree$events$uid),
         nevents=as.integer(nevents),
         compartmentval=as.integer(cfg$compartment$val),
         compartmentrate=as.double(cfg$compartment$rate),
         compartmentsize=as.integer(cfg$compartment$popsize),
         ncompartment=as.integer(length(cfg$compartment$val)),###
         compinfoval=as.integer(cfg$info$val),
         compinfofitness=as.double(cfg$info$fitness),
         compinfodriverstatus=as.integer(cfg$info$id),#as.integer(as.matrix(cfg$info[,idxd])),
         compinfondriver=as.integer(length(idxd)),
         ncomp=as.integer(length(cfg$info$val)),##
         as.integer(tips),  ##Actual tips to keep!
         as.integer(length(tips)),
         max_size=as.integer(MAX_SIZE),
         edgesOut=integer(2*MAX_SIZE),
         nDivsOut=integer(MAX_SIZE),
         stateOut=integer(MAX_SIZE),
         driverIDOut=integer(MAX_SIZE),
         tBirthOut=double(MAX_SIZE),
         nedgeOut=integer(1),
         nInternalNodeOut=integer(1),
         eventtsOut=double(nevents),
         eventvalOut=integer(nevents),
         eventdriveridOut=integer(nevents),
         eventuidOut=integer(nevents),
         eventnodeOut=integer(nevents),
         neventOut=integer(1),
         ntipsOut=integer(1),
         nCompPops=integer(length(cfg$info$val)),
         status=integer(1)
  )
  if(res$status>0){
    stop("Error call C(.sub_sample)")
  }
  
  nedge=res$nedgeOut
  ###browser()
  ##nEventsCount=res$nEventsCount
  cfg$info$population=res$nCompPops
  res=list(edge=matrix(res$edgesOut,ncol=2)[1:nedge,],
           edge.length=res$nDivsOut[1:nedge],
           state=res$stateOut[1:nedge],
           driverid=res$driverIDOut[1:nedge],
           Nnode=res$nInternalNodeOut,
           ndivs=res$nDivsOut[1:nedge],
           tBirth=res$tBirthOut[1:nedge],
           nedge=res$nedgeOut,
           ntips=res$ntipsOut,
           timestamp=tree$timestamp,
           pop.size=tree$pop.size,
           totaldrivercount=tree$totaldrivercount,
           ncellswithdrivers=1,
           params=tree$params,
           events=data.frame(value=res$eventvalOut[1:res$neventOut],
                             driverid=res$eventdriveridOut[1:res$neventOut],
                             node=res$eventnodeOut[1:res$neventOut],
                             ts=res$eventtsOut[1:res$neventOut],
                             uid=res$eventuidOut[1:res$neventOut]),
           cfg=cfg,
           currentEventUID=max(res$eventuidOut[1:res$neventOut])
           
  )
}


#' Gets subsampled tree - retaining events information appropriately
#' @param tree simpop
#' @param N integer. Number of tips to subsample (ignored if tips is specified).
#' @param tips integer. node id of tips to subsample. Must include out-group.
#' @return simpop object
#' @export
#' @examples
#' cfg=getDefaultConfig(target_pop_size =5e4,rate=0.1)
#' params=list(n_sim_days=2*365,b_stop_at_pop_size=1,b_stop_if_empty=0)
#' growthphase=sim_pop(NULL,params=params,cfg)
#' stree=get_subsampled_tree(growthphase,50)
get_subsampled_tree=function(tree,N,
                             tips=tree$edge[c(which(tree$state==0 & tree$edge[,2]<=length(tree$tip.label)),sample(which(tree$state!=0 & tree$edge[,2]<=length(tree$tip.label)),N)),2])
{
  ##tree: << A phylo object
  ##N<< Number of tips/cells to subsample
  
  N=length(tips)
  tmp=C_subsample_pop(tree,tips)
  ##browser()
  tmp$tip.label=sprintf("s%d",1:N)
  class(tmp)=c("simpop","phylo")
  checkValidPhylo(tmp)
  tmp
}


assign_celltype=function(tree, celltype, driverid, cfg=tree$cfg)
{
  #Sets the celltype for current cells (tips)
  ##<< A extended simpop phylo object
  ##<< cell type in parallel with tips NAs keep current cell type.
  ##<< Driver id (0) = no driver.
  ##<< specifies the compartment division rates etc...
  
  tree$maxt=max(tree$time)
  if(all(is.na(celltype))){
    tree$cfg=cfg
    return(tree)
  }
  
  idx=which(!is.na(celltype) & celltype!=0)
  utype=unique(celltype[idx])
  newtype=setdiff(utype,tree$events$value)
  #cat("Creating new cell types",newtype,"\n")
  if(any(newtype>0 & !(newtype %in% cfg$compartment$val))){
    stop("bad config provided")
  }
  tree$cfg=cfg
  ###browser()
  
  newevents=data.frame(value=celltype[idx],driverid=driverid[idx],node=idx,ts=rep(tree$maxt,length(idx)),uid=((1:length(idx))+tree$currentEventUID))  #(1:length(idx)+tree$currentEventUID)
  tree$events=rbind(tree$events,newevents)
  idx=which(!is.na(celltype))
  tree$state[match(idx,tree$edge[,2])]=celltype[idx]
  ##Following can be deleted
  tree$driverid[match(idx,tree$edge[,2])]=driverid[idx]
  tree$currentEventUID=max(tree$event$uid)
  tree
}

#' Gets elapsed time tree with optional mutational aquisition modelling.
#' Uses the model burden ~ sum(pois(mutatrateperdiv))+pois(backgroundrate*duration)
#' If mutrateperdivision is not provide the function returns the ultrametric elapsed time tree.
#' Note: the tree that comes out of sim_pop has branch lengths=number of symmetric divisions (i.e. edge_length=ndivs)
#' @param tree simpop phylo. output from sim_pop or get_subsampled_tree.
#' @param mutrateperdivision numeric. Mutation rate per division
#' @param backgroundrate numeric. Negative Binomial or Poisson Distributed Mutation count
#' @param odf  numeric. Overdispersion factor.
#' @return simpop/phylo - with branch lengths defined as above (sometimes stochastic)
#' @export
#'
get_elapsed_time_tree=function(tree,mutrateperdivision=NULL,backgroundrate=NULL,odf=1){
  N=length(tree$tip.label)+1
  L=length(tree$edge.length)
  TT=max(tree$timestamp)
  ##
  ##browser()
  idx.child=match(tree$edge[,2],tree$edge[,1])
  duration=ifelse(is.na(idx.child),TT-tree$tBirth,tree$tBirth[idx.child]-tree$tBirth)
  ##Set durat
  duration[which(tree$state==0)]=0
  if(!is.null(mutrateperdivision)){
    if(odf>1){
      tree$edge.length=sapply(tree$ndiv,function(n) sum(rpois(n,mutrateperdivision)))+get_nb(n=L,meanmuts=backgroundrate*duration,od=odf)
    }else{
      tree$edge.length=sapply(tree$ndiv,function(n) sum(rpois(n,mutrateperdivision)))+rpois(L,backgroundrate*duration)
    }
  }else{
    tree$edge.length=duration
  }
  tree
}



get_nb=function(meanmuts,od=2,n=1){
  ifelse(meanmuts>0,rnbinom(n,mu = meanmuts,size = meanmuts/(od-1)),0)
}


#' Concatenates simpops (maintains population size trajectory and drivercount trajectory)
#' Concatenates timestamps and population size
#' @param simpop1 simpop from initial run
#' @param simpop2 simpop from next run
#' @return simpop/phylo - with branch lengths defined as above (sometimes stochastic)
#' @export
#'
combine_simpops=function(simpop1,simpop2){
  ##cat("combine_simpops:Remember to add compartment level granularity")
  simpop2$timestamp=c(simpop1$timestamp,simpop2$timestamp)
  simpop2$pop.size=c(simpop1$pop.size,simpop2$pop.size)
  simpop2$totaldrivercount=c(simpop1$totaldrivercount,simpop2$totaldrivercount)
  simpop2
}
#' @title 
#' Gets a simple default configuration with two compartments
#' @description 
#' getDefaultConfig configures one main cell compartment using the user specified parameters and 
#' a default compartment (outgroup) which corresponds to the zygote.
#' @param targetPopSize integer - a negative number uses current time to 100th second resolution.
#' @param rate float - division rate in division per day
#' @param descr character - description of main compartment
#' @param ndriver integer - number of drivers - redundant now.
#' @param basefit float - base fitness redundant now.
#' @return list config
#' @export
getDefaultConfig=function(target_pop_size,rate,
                          descr="cellType1",ndriver=1,
                          basefit=0.1)
{
  if(!is.numeric(rate) | length(rate)>1)
    stop("per day division rate should be a numeric value.")
  if(!is.numeric(basefit) | length(basefit)>1)
    stop("base fitness should be a numeric value.")
  if(!is.numeric(ndriver) | length(ndriver)>1)
    stop("the number of driver should be an integer.")
  if(!is.numeric(target_pop_size) | length(target_pop_size)>1)
    stop("the target population size should be a scalars")
  compartment=data.frame(val=0,
                         rate=-1,
                         popsize=1,
                         desc="outgroup")
  info=data.frame(val=0,population=0,fitness=0,id=0,driver1=0)
  cfg=list(compartment=compartment,info=info)
  if(length(target_pop_size)>1){
    stop("please provide scalars")
  }
  #addCellCompartment(cfg,target_pop_size,rate,ndriver,basefit,descr)
  addCellCompartment(cfg,target_pop_size,rate,descr=descr)
}

#' Adds a new cell compartment to a simpop configuration
#' @param cfg list. simpop config
#' @param rate target population
#' @param ndriver integer. DEPRECATED
#' @param basefit numeric. DEPRECATED
#' @param descr character. Description
#' @return list simpop config
#' @export
addCellCompartment=function(cfg,population,rate,ndriver=0,basefit=0,descr=NULL){
  val=max(cfg$compartment$val)+1
  if(is.null(descr)){
    descr=sprintf("cellType%d",val)
  }
  row=cfg$info[which(cfg$info$val==0),]
  row$val=val
  row$population=0
  row$fitness=0
  row$id=0
  
  cfg$info=rbind(cfg$info,row)
  cfg$compartment=rbind(cfg$compartment,data.frame(val=val,rate=rate,popsize=population,desc=descr))
  cfg$drivers=rbind(cfg$drivers,data.frame(val=val,driver=1,fitness=0))
  cfg
}

#' Adds a differentiation events to a simpop. Can either specify relevant tips or randomly generate.
#' @param tree - simpop object.
#' @param cfg  - List. simpop config.
#' @param newCellType - character
#' @param idx - integer/vector. Index of the tip where to add the differentiation event. 
#' @param nEvents - integer. the number of differentiation events to be added
#' @param currentCompartment - integer. ID (value) of compartment where the differentiation event will be added to.
#' @return simpop object
#' @export
addDifferentiationEvents = function(tree,cfg,newCellType,idx=NULL,nEvent=-1,currentCompartment=-1)
{
  if(!newCellType %in% cfg$compartment$val){
    stop("newCellType not in cfg$compartment")
  }
  N=length(tree$tip.label)
  idxt=match(1:N,tree$edge[,2])
  if(is.null(idx))
  {
    if(nEvent<0){
      stop("Please specify nEvent")
    }
    if(currentCompartment<0){
      ##randomly pick from all non-outgroup
      idx=sample(which(tree$state[idxt]>0),size = nEvent,replace = FALSE)
    }else{
      idxx=which(tree$state[idxt]==currentCompartment)
      if(length(idxx)==0){
        stop("Choose a currentCompartment that has non-zero representation in tree!")
      }
      idx=sample(idxx,size = nEvent,replace = FALSE)
    }
  }
  else
  {
    if(nEvent<0)
      stop("Please specify nEvent")
    if(currentCompartment<0)   ##randomly pick from all non-outgroup
    {
      if(!all(idx %in% which(tree$state[idxt]>0)))
        stop("some tips indexes are out of range of non-outgroup tip indexes.")
    }
    else
    {
      if(!all(idx %in% which(tree$state[idxt]==currentCompartment)))
        stop("some tips indexes are out of range of current compartment tip indexes.")
    }
  }
  celltype=rep(NA,N)
  celltype[idx]=newCellType
  driverid=rep(0,N)
  assign_celltype(tree,celltype,driverid,cfg=cfg)  ##<< specifies the compartment division rates etc.
}

make_star_tree=function(N){
  dat=list(edge=matrix(as.integer(c(rep(N+1,N),1:N)),ncol=2),edge.length=rep(0,N),tip.label=sprintf("t%d",1:N),Nnode=as.integer(1))
  class(dat)="phylo"
  dat
}

####New driver approach
#' Adds a driver event to a current simpop
#' @param tree - simpop.
#' @param cfg  - List. simpop config.
#' @param currentCompartment - integer. ID (value) of compartment on which the driver will be added to.
#' @param fitness - numreic. the fitness associated to the driver.
#' @return 
#' simpop object
#' @export
addDriverEvent=function(tree,cfg,currentCompartment,fitness)
{
  N=length(tree$tip.label)
  idxt=match(1:N,tree$edge[,2])
  idxx=which(tree$state[idxt]==currentCompartment)
  if(length(idxx)==0){
    stop("Choose a currentCompartment that has non-zero representation in tree!")
  }
  if(abs(fitness)<1e-6){
    stop("Unable to addDriverEvent with zero fitness")
  }
  
  idx=sample(idxx,size = 1,replace = FALSE)
  celltype=rep(NA,N)
  celltype[idx]=tree$state[idxt[idx]]
  
  ##Find driver ID take the lowest with an empty population
  idxd=grep("^driver",colnames(cfg$info))
  d=-1
  for(i in 1:length(idxd)){
    if(sum(cfg$info$population[which(cfg$info$val==currentCompartment & cfg$info[,idxd[i]]==1)])==0){
      d=i
      break
    }
  }
  #browser()
  if(d<0)
  {
    ##We need to add a single new row to cfg$info
    d=length(idxd)+1
    cfg$info[[sprintf("driver%d",d)]]=0
    idxd=grep("^driver",colnames(cfg$info))
    ##Need to add a driver
    cfg$drivers=rbind(cfg$drivers,data.frame(val=currentCompartment,driver=d,fitness=fitness))
  }
  ##Drop all rows with the existing driver<d>=1
  idx.drop=which(cfg$info$val==currentCompartment & cfg$info[,idxd[d]]==1)
  if(length(idx.drop)>0){
    cfg$info=cfg$info[-idx.drop,]
  }
  olddriverid=tree$driverid[idxt[idx]]
  idxn=which(cfg$info$id==olddriverid & cfg$info$val==currentCompartment)
  if(length(idxn)!=1){
    ##We use the logic that the most recently acquired driver uniquely defines the
    ##driver "signature". i.e. if driver a is followed by b is followed by c then all
    ##cells carrying c also carry a and b.
    stop("Unexpected behaviour")
  }
  if(cfg$info$population[idxn]>1){
    ##We keep the old row and decrement the population by 1
    cfg$info$population[idxn]=cfg$info$population[idxn]-1
    ##Add a copy of the old row  with a population of one and driver<d> set to 1
    row=cfg$info[idxn,]
    row$population=1
    row[,idxd[d]]=1
    row$id=d
    cfg$info=rbind(cfg$info,row)
  }else{
    ##Old row is now redundant - so just repurpose it
    if(cfg$info$population[idxn]!=1){
      stop("Unexpected behaviour")
    }
    cfg$info$id[idxn]=d
    cfg$info[idxn,idxd[d]]=1
  }
  cfg$drivers$fitness[which(cfg$drivers$val==currentCompartment & cfg$drivers$driver==d)]=fitness
  cfg=recalculateFitnessV2(cfg,currentCompartment)
  #}
  driverid=rep(d,length(celltype))
  assign_celltype(tree,celltype,driverid,cfg=cfg)
}

recalculateFitnessV2=function(cfg,cellType){
  idx=which(cfg$info$val==cellType)
  if(length(idx)>0)
  {
    drivcol=sprintf("driver%d",cfg$drivers$driver[which(cfg$drivers$val==cellType)])
    fitness=cfg$drivers$fitness[which(cfg$drivers$val==cellType)]
    #additive
    if(length(drivcol)==1){
      cfg$info$fitness[idx]=fitness*cfg$info[idx,drivcol]
    }else{
      cfg$info$fitness[idx]=apply(cfg$info[idx,drivcol],1,function(x) sum(fitness*x))
      #multiplicative
      ##cfg$info$fitness[idx]=apply(cfg$info[idx,drivcol],1,function(x) prod((1+fitness*x))-1)
    }
  }
  cfg
}


