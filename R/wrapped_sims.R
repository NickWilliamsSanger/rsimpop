### This file contains wrapped calls to the simulator that serve as examples

#' Runs a single compartment/single driver scenario
#'
#' The driver is subject to stochastic extinction and so is repeatedly dropped in until it "takes".
#' The prevailing state when the driver is introduced is saved and reinstated with each attempted introduction.
#' @param initial_division_rate - Rate of symmetric cell division during development
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param nyears_driver_acquisition - When driver is acquired.
#' @param nyears - Total number of years to run the simulation
#' @param fitness - relative fitness advantage.  Cells carrying this divide at a rate=(1+fitness)*baserate
#' @param minprop - Minimum aberrant cell fraction at nyears for driver to be regarded as have "taken"
#' @param maxtry  - Maximum number of attempts to introduce a driver
#' @return simpop object.
#' @export
#' @examples
#' selsim=run_selection_sim(0.05,1/365,target_pop_size = 5e4,nyears = 50,fitness=0.3)
run_selection_sim=function( initial_division_rate=0.1,
                            final_division_rate=1/365,
                            target_pop_size=1e5,
                            nyears_driver_acquisition=15,
                            nyears=40,
                            fitness=0.2, ## relative fitness
                            minprop=0.001,
                            mindriver=1,
                            maxtry=40
){
  cfg=getDefaultConfig(target_pop_size,rate=initial_division_rate,ndriver=1,basefit = fitness)
  params=list(n_sim_days=nyears_driver_acquisition*365,
              b_stop_at_pop_size=1,
              b_stop_if_empty=0
  )
  growthphase=sim_pop(NULL,params=params,cfg)
  ndivkeep=0
  mdivkeep=0
  gdivkeep=0
  tree0=get_tree_from_simpop(growthphase)
  if(growthphase$status==0){
    ##We've exited before maturity and need to add driver now..
    tree1=tree0
    params[["n_sim_days"]]=nyears*365 ##years of simulation
    params[["b_stop_if_empty"]]=1

    dc=0
    tries=0
    tree1_tmp=tree1
    while(dc<max(minprop*target_pop_size,mindriver) ){
      if(tries>=maxtry){
        return(NULL)
      }
      cat("No driver found: tries=",tries,"\n")
      tries=tries+1
      tree1_tmp=addDriverEvent(tree1,tree1$cfg,1,fitness=fitness)
      print(tree1_tmp$cfg$info)
      params[["b_stop_at_pop_size"]]=1
      adult2=sim_pop(tree1_tmp,params=params,tree1_tmp$cfg)
      adult2a=combine_simpops(growthphase,adult2)
      tree2=get_tree_from_simpop(adult2a)
      params[["b_stop_at_pop_size"]]=0
      cfg=tree2$cfg
      cfg$compartment$rate[2]=final_division_rate
      cfg$compartment$popsize[2]=target_pop_size
      adult2=sim_pop(tree2,params=params,cfg)
      adult2=combine_simpops(adult2a,adult2)
      dc=adult2$cfg$info$population[3]
    }
  }else{
    gdivkeep=mean(nodeHeights(tree0)[which(tree0$edge[,2]<=length(tree0$tip.label)),2])
    cfg$compartment$rate[2]=final_division_rate
    cfg$compartment$popsize[2]=target_pop_size
    years=nyears
    params[["n_sim_days"]]=nyears_driver_acquisition*365 ##years of simulation
    params[["b_stop_at_pop_size"]]=0 ## So it doesn't stop simulating immediately
    adult1=sim_pop(tree0,params=params,cfg)
    adult1=combine_simpops(growthphase,adult1)
    tree1=get_tree_from_simpop(adult1)
    params[["n_sim_days"]]=nyears*365 ##years of simulation
    params[["b_stop_if_empty"]]=1
    dc=0
    tries=0
    tree1_tmp=tree1
    while(dc<max(minprop*target_pop_size,mindriver) ){
      if(tries>=maxtry){
        return(NULL)
      }
      cat("No driver found: tries=",tries,"\n")
      tries=tries+1
      tree1_tmp=addDriverEvent(tree1,tree1$cfg,1,fitness=fitness)
      if(TRUE){
        ndivkeep=nodeHeights(tree1_tmp)[which(tree1_tmp$edge[,2]==tree1_tmp$events$node[3]),2]
        mdivkeep=mean(nodeHeights(tree1_tmp)[which(tree1_tmp$edge[,2]<=length(tree1_tmp$tip.label)),2])
      }
      print(tree1_tmp$cfg$info)
      adult2=sim_pop(tree1_tmp,params=params,tree1_tmp$cfg)
      dc=adult2$cfg$info$population[3]
    }
    adult2=combine_simpops(adult1,adult2)
  }

  fulltree=get_tree_from_simpop(adult2)
  fulltree$tries=tries
  fulltree$ndivkeep=ndivkeep
  fulltree$mdivkeep=mdivkeep
  fulltree$gdivkeep=gdivkeep
  return(fulltree)
}

#' Runs a multiple driver scenario
#'
#' The user specifies the rate at which drivers are introduced and the distribution the selective coefficients are drawn from.
#' Note that unlike in \code{\link{run_selection_sim}} the drivers are allowed to stochastically die out.
#' The drivers at the specified rate with the gap between successive introductions exponentially distributed.
#' @param initial_division_rate - Rate of symmetric cell division during development
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param drivers_per_year - The expected number of drivers per year
#' @param nyears - Total number of years to run the simulation
#' @param bForceReseed - Whether to reseed.
#' @param offset - time offset for choosing seed.
#' @return simpop object.
#' @export
#' @examples
#' selsim=run_selection_sim(0.05,1/365,target_pop_size = 5e4,nyears = 50,fitness=0.3)
run_driver_process_sim=function( initial_division_rate=0.1,
                                  final_division_rate=1/365,
                                  target_pop_size=1e5,
                                  drivers_per_year=0.1,
                                  nyears=40,
                                  fitnessGen=function(){0},
                                  bForceReseed=FALSE,offset=0
){
  if(bForceReseed){
    delay=(offset/10 +1)
    Sys.sleep(delay)
    cat("delay=",delay,"\n")
    initSimPop(-1,bForce = TRUE)
  }
  ##driverdt=cumsum(rexp(1000,drivers_per_year/365))
  dpcpd=(drivers_per_year/365)/target_pop_size

  k=1
  nsd=nyears*365
  cfg=getDefaultConfig(1e5,rate=initial_division_rate,ndriver=1,basefit = 0)
  params=list(n_sim_days=nsd,##This is the longest that the simulation will continue
              b_stop_at_pop_size=1,
              b_stop_if_empty=0,
              driver_rate_per_cell_per_day=dpcpd
  )
  growthphase=get_tree_from_simpop(sim_pop(NULL,params=params,cfg))
  while(growthphase$status==2 || growthphase$status==0){
    k=k+1
    params=list(n_sim_days=nsd,##This is the longest that the simulation will continue
                b_stop_at_pop_size=1,
                b_stop_if_empty=0,
                driver_rate_per_cell_per_day=dpcpd
    )
    if(max(growthphase$timestamp)>params$n_sim_days){
      stop("ts>n_sim_days:unexpected bahaviour")
    }
    growthphase=addDriverEvent(growthphase,growthphase$cfg,currentCompartment = 1,fitness=fitnessGen())
    gpk=get_tree_from_simpop(sim_pop(growthphase,params=params,growthphase$cfg))
    growthphase=combine_simpops(growthphase,gpk)
  }
  tree0=get_tree_from_simpop(growthphase)
  tree0$cfg$compartment$rate[2]=final_division_rate
  tree0$cfg$compartment$popsize[2]=target_pop_size
  years=nyears
  params[["b_stop_at_pop_size"]]=0 ## So it doesn't stop simulating immediately
  if(max(tree0$timestamp)>=nyears*365){
    return(tree0)
  }
  adult=sim_pop(tree0,params=params,tree0$cfg)
  adult=combine_simpops(tree0,adult)
  tree0=get_tree_from_simpop(adult)
  if(length(which(tree0$cfg$info$fitness==0 & tree0$cfg$info$population>0))>2){
    browser()
  }
  while(TRUE){
    k=k+1
    ## Deal with end of last sim and next driver being too close together. Very rarely kicks in.
    ##params[["n_sim_days"]]=min(nyears*365,max(driverdt[k],max(tree0$timestamp)+0.1)) ##years of simulation
    if(max(tree0$timestamp)>=nyears*365){
      return(tree0)
    }
    #browser()
    tree0=addDriverEvent(tree0,tree0$cfg,1,fitness=fitnessGen())
    tmp2=tree0
    adult=sim_pop(tree0,params=params,tree0$cfg)
    adult=combine_simpops(tree0,adult)
    tree0=get_tree_from_simpop(adult)
    print(tree0$cfg$info %>% dplyr::filter(population>0))
    if(length(which(tree0$cfg$info$fitness==0 & tree0$cfg$info$population>0))>2){
      browser()
    }
    if(k>1000){
      stop("too many iterations!")
    }
  }

}

#' Runs a single compartment/single driver scenario with transient selection
#'
#' The driver is subject to stochastic extinction and so is repeatedly dropped in until it "takes".
#' The prevailing state when the driver is introduced is saved and reinstated with each attempted introduction.
#' @param initial_division_rate - Rate of symmetric cell division during development
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param nyears_driver_acquisition - When driver is acquired.
#' @param nyears_transient_end - When selective advantage is set to 0.
#' @param nyears - Total number of years to run the simulation
#' @param fitness - relative fitness advantage.  Cells carrying this divide at a rate=(1+fitness)*baserate
#' @param minprop - Minimum aberrant cell fraction at nyears for driver to be regarded as have "taken"
#' @param maxtry  - Maximum number of attempts to introduce a driver
#' @return simpop object.
#' @export
#' @examples
#' tselsim=run_transient_selection(0.05,1/365,target_pop_size = 5e4,nyears_driver_acquisition=15,
#' nyears_transient_end=30,nyears=50,fitness=0.5)
run_transient_selection=function( initial_division_rate,
                                  final_division_rate,
                                  target_pop_size=1e5,
                                  nyears_driver_acquisition=15,
                                  nyears_transient_end=30,
                                  nyears=40,
                                  fitness=0.2, ## relative fitness
                                  minprop=0.05
){
  selsim=run_selection_sim(initial_division_rate,
                                    final_division_rate,
                                    target_pop_size,
                                    nyears_driver_acquisition,
                                    nyears_transient_end,
                                    fitness,minprop)
  ##switch off
  cfg=selsim$cfg
  cfg$info$fitness[3]=0
  params=selsim$params
  params[["n_sim_days"]]=nyears*365
  params[["maxt"]]=NULL
  final=sim_pop(selsim,params=params,cfg)
  final=combine_simpops(selsim,final)
  return(get_tree_from_simpop(final))
}

#' Runs a simple neutral simulation
#'
#' @param initial_division_rate - Rate of symmetric cell division during development
#' @param final_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param target_pop_size - Size of target population
#' @param nyears - Total number of years to run the simulation
#' @return simpop object.
#' @export
#' @examples
#' neutsim=run_neutral_sim(0.05,1/365,target_pop_size = 5e4,nyears=80)
run_neutral_sim=function( initial_division_rate,
                                            final_division_rate,
                                            target_pop_size=1e5,
                                            nyears=40
){
  cfg=getDefaultConfig(target_pop_size,rate=initial_division_rate,ndriver=1,basefit = 0)
  params=list(n_sim_days=nyears*365,
              b_stop_at_pop_size=1,
              b_stop_if_empty=0
  )
  growthphase=sim_pop(NULL,params=params,cfg)
  cfg$compartment$rate[2]=final_division_rate
  cfg$compartment$popsize[2]=target_pop_size
  params[["b_stop_at_pop_size"]]=0
  adult1=sim_pop(growthphase,params=params,cfg)
  return(combine_simpops(growthphase,adult1))
}


add_driver=function(simpop,params,acceptance_threshold=0.05){
  tree1=get_tree_from_simpop(simpop)
  params[["b_stop_if_empty"]]=1
  tree1=simpop
  dc=0
  tries=0
  tree1_tmp=tree1
  while(dc/target_pop_size<0.05){
    cat("No driver found: tries=",tries,"\n")
    idx=sample(length(tree1$tip.label)-1,1)+1
    celltype=rep(NA,length(tree1$tip.label))
    celltype[idx]=-1
    tree1_tmp=assign_celltype(tree1,celltype,tree1$cfg)
    simpop2=sim_pop(tree1_tmp,params=params,tree1_tmp$cfg)

    dc=simpop2$cfg$info$population[3]
    #print(adult2$cfg$info)
    tries=tries+1
    if(tries>max_tries){
      stop("Unable to add driver.. Too many attempts")
    }
  }
  simpop=combine_simpops(simpop,simpop2)
  simpop
}

#' Runs a simple single compartment neutral simulation with a specified trajectory
#'
#' @param simpop - Rate of symmetric cell division during development
#' @param initial_division_rate - Rate of symmetric cell division once population equilibrium is reached.
#' @param trajectory - data.frame - with fields ts(timestamp in days),target_pop_size,division_rate
#' @param nyears - Total number of years to run the simulation
#' @return simpop object.
#' @export
#' @examples
#' trajectory=data.frame(ts=365*(1:80),target_pop_size=5e4+100*(1:80),division_rate=1/(2*190))
#' trajectory$target_pop_size[5:10]=2*trajectory$target_pop_size[5:10]
#' trajectory$target_pop_size[11:15]=0.2*trajectory$target_pop_size[11:15]
#' sp=run_neutral_trajectory(NULL,0.5,trajectory)
#' plot(sp)
run_neutral_trajectory=function(simpop,initial_division_rate,trajectory){
  if(any(trajectory$division_rate>1)){
    stop("Supplied division rate too high; should be less than 1.0")
  }
  if(any(trajectory$target_pop_size>20e6)){
    stop("Supplied population too high; should be less than 20e6")
  }

  if(is.null(simpop)){
    cfg=getDefaultConfig(target_pop_size = trajectory$target_pop_size[1],rate = initial_division_rate,basefit = 0,ndriver = 1)

    params=list(n_sim_days=trajectory$ts[1],##This is the longest that the simulation will continue
                b_stop_at_pop_size=1,
                b_stop_if_empty=0
    )
    sp=sim_pop(NULL,params=params,cfg,b_verbose = FALSE)
    idx=1
  }else{
    maxt=max(simpop$timestamp)
    idx=which(trajectory$ts>maxt)
    if(length(idx)>0){
      idx=idx[1]
    }else{
      idx=1
    }

    sp=simpop
  }
  params[["b_stop_at_pop_size"]]=0
  for(i in idx:(length(trajectory$ts)-1)){
    cfg=list(compartment=data.frame(val=c(0,1),rate=c(-1,trajectory$division_rate[i]),popsize=c(1,trajectory$target_pop_size[i])),
           info=sp$cfg$info)
    cfg=getDefaultConfig(target_pop_size = trajectory$target_pop_size[i],rate = trajectory$division_rate[i],basefit = 0,ndriver = 1)
    params[["n_sim_days"]]=trajectory$ts[i+1]
    spx=sim_pop(get_tree_from_simpop(sp),params=params,cfg,b_verbose = FALSE)
    sp=combine_simpops(sp,spx)
  }
  sp
}


run_2compartment_model=function(){
   ##Initialise LT-HSC
   #

}

simple_continuous_selection=function(){

}
