#require("RColorBrewer")
#require("scales")
#require("tidyverse")

#' Plots an ape phylo tree
#'
#' This is similar to APE plot.phylo but returns useful information on each branch to facilitate annotation of tree.
#' Be careful not to unintentionally supply population scale simpop
#' @param tree phylo/simpop
#' @param direction character. Only "down" is currently supported.
#' @param cex.label numeric. Size of tip labels (0 results in no tip labels)
#' @param offset numeric (0-1). Fractional offset of tip labels from branch terminus
#' @param b_do_not_plot Boolean.
#' @param lwd numeric. Line width
#' @param bars numeric. vector with named elements that should all match tip labels - adds barplot under tree.
#' @param default_edge_color colour.
#' @param ymax  numeric max y
#' @param cex.terminal.dots numeric - size of terminal dots
#' @param mar numeric vector Override default resetting of mar as in par("mar")
#' @param b.add.scale Boolean - whether to add scale axis
#' @param left.margin.prop numeric - leaves this fraction empty on left
#' @return simpop object
#' @export
#' @examples
#' test=rtree(20)
#' tree=plot_tree(test,cex.terminal=0.5,bars=sapply(sprintf("s%d",seq(2,20,2)),function(x) runif(1)))
#' node_labels(tree)
plot_tree=function(tree,direction="down",cex.label=1,offset=0,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="darkgrey",ymax=NULL,cex.terminal.dots=0,mar=NULL,b.add.scale=TRUE,left.margin.prop=0.1){
  if(is.null(mar)){
    par(mar=c(1, 1, 1, 3) + 0.1)
  }else{
    par(mar=mar)
  }
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  if(direction=="across"){
    stop("More work needed to support horizontal plotting (across) - please use down")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords

  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ymax=max(coords$a1)*1.05
    }
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  if(is.null(bars)){
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*left.margin.prop),xmax),ylim=c(0-(ymax*0.05),ymax),xlab="",ylab="")
  }else{

    plot(NULL,axes=FALSE,xlim=c(0-(xmax*left.margin.prop),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=max(coords$a1)
    ##browser()
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
      text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=tipcol,cex=cex.terminal.dots,pch=19)
    }
  }
  tree$direction=direction
  tree$top=top
  scales=c(0,1,10,50,100,200,500,1000,2000,5000)
  #scale=scales[max(which(ymax/2>scales))]
  scale=scales[max(which(ymax/5>=scales))]
  #cat("scale=",scale,"\n")
  if(b.add.scale){
    axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2)
  }
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("%s Muts",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,col = "grey")

  }
  tree
}
set_cedge=function(parent,tree){
  for(child in get_node_children(parent,tree)){
    child.idx=which(tree$edge[,2]==child)
    parent.idx=which(tree$edge[,2]==parent)
    if(length(parent.idx)==0){
      pedge=0
    }else{
      pedge=tree$cedge[parent.idx]
    }
    tree$cedge[child.idx]=pedge+tree$edge.length[child.idx];
    tree=set_cedge(child,tree)
  }
  tree
}

##Not very efficient -- recursively calculates height as average of children's height.
get_height=function(tree,node){
  N=length(tree$tip.label)
  if(node<=N){
    return(node)#tree$height[which(tree$edge[,2]==node)]=node
  }else{
    children=get_node_children(node,tree)
    return(mean(sapply(children,function(child) get_height(tree,child))))
  }
}

set_height=function(tree){
  tree$height_end=sapply(tree$edge[,2],function(i) get_height(tree,i))
  tree$height_start=tree$height[match(tree$edge[,1],tree$edge[,2])]
  N=length(tree$tip.label)
  root=N+1
  idx=which(tree$edge[,1]==root)
  tree$height_start[idx]=mean(tree$height_end[idx])
  tree
}

##horizontal elbow
elbow=function(x0,x1,y0,y1,...){
  arrows(x0=x0,y0=y0,y1=y1,length=0,...)
  arrows(x0=x0,x1=x1,y0=y1,length=0,...)
}
##vertical elbow
elbowv=function(x0,x1,y0,y1,...){
  #browser()
  arrows(x0=x0,x1=x1,y0=y0,length=0,...)
  arrows(x0=x1,y0=y0,y1=y1,length=0,...)

}

set_tree_coords=function(atree){
  ##get the cumulative distance from the root.
  tree=atree
  tree$cedge=rep(0,length(tree$edge.length))
  N=length(tree$tip.label)
  root=N+1
  tree=set_cedge(root,tree)
  tt=set_height(tree)
  atree$coords=data.frame(a0=tt$cedge-tt$edge.length,a1=tt$cedge,
                          b0=tt$height_start,b1=tt$height_end,stringsAsFactors = FALSE)
  if(!is.null(atree$color)){
    atree$coords$color=atree$color
  }
  atree
}


#' Gets all descendents of the specified node.
#'
#' TODO: Deprecate.. Does the same job as phytools::getDescendents
#' @param node id of specified node
#' @param tree phylo. Tree
#' @export
get_all_node_children=function(node,tree){
  children=tree$edge[which(tree$edge[,1]==node),2]
  offspring=children
  for(child in children){
    offspring=c(offspring,get_all_node_children(child,tree))
  }
  offspring
}

#' Gets node id of immediate descendents of the specified nodes.
#'
#' Will be empty for terminal nodes.
#' @param node id of specified node
#' @param tree phylo.
#' @export
get_node_children=function(node,tree){
  tree$edge[which(tree$edge[,1]==node),2]
}

#' Gets tip labels of tips that descend from the specified node
#' @param node id of specified node
#' @param tree phylo.
#' @export
get_samples_in_clade=function(node,tree){
  if(node<=length(tree$tip.label)){
    return(tree$tip.label[node])
  }
  tree$tip.label[intersect(get_all_node_children(node,tree),1:length(tree$tip.label))]
}

get_y_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$top-tree$coords[idx,c("a1","a0")])
}

get_x_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$coords[idx,c("b1","b0")])
}
##The following ties the PDX datastructure to the tree
get_idx_for_node=function(pdx,node){
  which(pdx$details$node==node)
}

#' Gets summary edge information specified by child node of an edge.
#'
#' @param pdx list containing a data frame "details" with an entry "node" that maps to nodes in the tree.
#' @param tree phylo. Enhanced phylo - has coords data.frame specifying the coordinates of the edges.
#' @param node integer.
#' @export
get_edge_info=function(pdx,tree,node){
  y=get_y_range(tree,node)
  x=get_x_range(tree,node)
  idx=get_idx_for_node(pdx,node)
  samples=get_samples_in_clade(node,tree)
  list(yb=y[1],yt=y[2],x=x[1],xm=x[2],idx.in.details=idx,samples=samples)
}


add_annotation=function(pdx,tree,annot_function,control=NULL){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(pdx,tree,tree$edge[i,2],control=control))
}

add_defaults=function(control,defaults,mandatory_fields=c()){
  if(is.null(control)){
    control=list()
  }
  for(field in mandatory_fields){
    if(is.null(control[[field]])){
      stop(sprintf("Required parameter %s not supplied in control list",field))
    }
  }

  for(field in names(defaults)){
    if(is.null(control[[field]])){
      control[[field]]=defaults[[field]]
    }
  }
  control
}

#' Adds node labels to an existing tree plot
#'
#' @param tree enhanced phylo returned by plot_tree
#' @param col character. Colour specification for labels
#' @param cex numeric.  Size of labels
#' @param b_include_tips boolean. Wether to include tips
#' #' @return simpop object
#' @export
#' @examples
#' test=rtree(20)
#' tree=plot_tree(test,cex.terminal=0.5,bars=sapply(sprintf("s%d",seq(2,20,2)),function(x) runif(1)))
#' node_labels(tree)
#'
node_labels=function(tree,col="blue",cex=1,b_include_tips=FALSE){
  idx=which(tree$edge[,2]>length(tree$tip.label))
  if(b_include_tips){
    idx=1:dim(tree$edge)[1]
  }else{
    idx=which(tree$edge[,2]>length(tree$tip.label))
  }
  text(tree$edge[idx,2],x=tree$coords$b1[idx],y=tree$top-tree$coords$a1[idx],col=col,cex=cex)
}

get_parents=function(node,mut_table,exclude_root=TRUE){
  idx=which(mut_table[,2]==node)
  parents=node##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=mut_table[idx,1]
    parents=c(parents,parent)
    idx=which(mut_table[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)]
  }else{
    parents
  }
}
