
### plotTree, nterminal, maxdepth by Thorsten Hothorn 

## ## utility functions for querying the number of
## ## terminal nodes and the maximal depth of (sub-)trees
nterminal <- function(node) {
    if (node$terminal) return(1)
    nl <- nterminal(node$left)
    nr <- nterminal(node$right)
    return(nl + nr)
}

maxdepth <- function(node) {
    if (node$terminal) return(1)
    nl <- maxdepth(node$left)
    nr <- maxdepth(node$right)
    return(max(c(nl, nr)) + 1)
}


plotTree <- function(node, xlim, ylim, nx, ny, 
               terminal_panel, inner_panel, edge_panel,
	       tnex = 2, drop_terminal = TRUE, debug = FALSE) {

    ### the workhorse for plotting trees

    ### set up viewport for terminal node
    if (node$terminal) {
        x <- xlim[1] + diff(xlim)/2
        y <- ylim[1] + 0.5
       
        tn_vp <- viewport(x = unit(x, "native"),
                          y = unit(y, "native") - unit(0.5, "lines"),
                          width = unit(1, "native"), 
                          height = unit(tnex, "native") - unit(1, "lines"),
			  just = c("center", "top"),
                          name = paste("Node", node$nodeID, sep = ""))
        pushViewport(tn_vp)
        if (debug)
            grid.rect(gp = gpar(lty = "dotted", col = 4))
        terminal_panel(node) 
        upViewport()
        return(NULL)
    }    

    ### number of left leafs
    nl <- nterminal(node$left)

    ### number of right leafs
    nr <- nterminal(node$right)

    ### position of inner node
    x0 <- xlim[1] + (nl / (nl + nr)) * diff(xlim)
    y0 <- max(ylim)

    ### proportion of left terminal nodes in left node
    if (node$left$terminal) {
        lf <- 1/2
    } else {
        lf <- nterminal(node$left$left) / (nterminal(node$left$left) + 
                                           nterminal(node$left$right))
    }

    ### proportion of left terminal nodes in right node
    if (node$right$terminal) {
        rf <- 1/2
    } else {
        rf <- nterminal(node$right$left) / (nterminal(node$right$left) + 
                                            nterminal(node$right$right))
    }

    ### position of left and right daugher node
    x1l <- xlim[1] + (x0 - xlim[1]) * lf
    x1r <- x0 + (xlim[2] - x0) * rf
    
    if (!drop_terminal) {
        y1l <- y1r <- y0 - 1
    } else {
        y1l <- if (node$left$terminal) tnex - 0.5 else y0 - 1
        y1r <- if (node$right$terminal) tnex - 0.5 else y0 - 1
    }

    ### draw edges
    grid.lines(x = unit(c(x0, x1l), "native"), 
               y = unit(c(y0, y1l), "native"))
    grid.lines(x = unit(c(x0, x1r), "native"), 
               y = unit(c(y0, y1r), "native"))

    ### create viewport for inner node
    in_vp <- viewport(x = unit(x0, "native"),
                      y = unit(y0, "native"),
                      width = unit(1, "native"),
                      height = unit(1, "native") - unit(1, "lines"), 
                      name = paste("Node", node$nodeID, sep = ""))
    pushViewport(in_vp)
    if (debug)
        grid.rect(gp = gpar(lty = "dotted"))
    inner_panel(node)
    upViewport()

    ps <- node$psplit
    if (ps$ordered) {
        if (!is.null(attr(ps$splitpoint, "levels"))) {
            split <- attr(ps$splitpoint, "levels")[ps$splitpoint]
        } else {
            split <- ps$splitpoint
        }
    } else {
        ### <FIXME>: always to the left? </FIXME>
        split <- attr(ps$splitpoint, "levels")[as.logical(ps$splitpoint) & (ps$table > 0)]
    }


    ### position of labels
    y1lr <- max(y1l, y1r)
    ypos <- y0 - (y0 - y1lr) * 0.5
    xlpos <- x0 - (x0 - x1l) * 0.5 * (y0 - y1lr)/(y0 - y1l)
    xrpos <- x0 - (x0 - x1r) * 0.5 * (y0 - y1lr)/(y0 - y1r)

    ### setup left label
    lsp_vp <- viewport(x = unit(xlpos, "native"),
                       y = unit(ypos, "native"),
                       width = unit(xlpos - xrpos, "native"),
                       height = unit(1, "lines"), 
                       name =  paste("lEdge", node$nodeID, sep = ""))
    pushViewport(lsp_vp)
    if (debug)
        grid.rect(gp = gpar(lty = "dotted", col = 2))
    edge_panel(split, ordered = ps$ordered, left = TRUE)
    upViewport()

    ### setup right label
    if (ps$ordered) {
        if (!is.null(attr(ps$splitpoint, "levels"))) {
            split <- attr(ps$splitpoint, "levels")[ps$splitpoint]
        } else {
            split <- ps$splitpoint
        }
    } else {
        split <- attr(ps$splitpoint, "levels")[!as.logical(ps$splitpoint) & (ps$table > 0)]
    }

    rsp_vp <- viewport(x = unit(xrpos, "native"),
                       y = unit(ypos, "native"),
                       width = unit(xlpos - xrpos, "native"),
                       height = unit(1, "lines"),
                       name =  paste("rEdge", node$nodeID, sep = ""))
    pushViewport(rsp_vp) 
    if (debug)
        grid.rect(gp = gpar(lty = "dotted", col = 2))
    edge_panel(split, ordered = ps$ordered, left = FALSE)
    upViewport()

    plotTree(node$left, c(xlim[1], x0), c(y1l, 1), nx, ny, 
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
    plotTree(node$right, c(x0, xlim[2]), c(y1r, 1), nx, ny,
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
}

##changed binary tree function for wikileaks
#author: Thomas Rusch 2010
plot.BinaryTree <- function(x, main = NULL, type = c("extended", "simple"),
                            terminal_panel = NULL, tp_args = list(),
			    inner_panel = node_inner, ip_args = list(),
                            edge_panel = edge_simple, ep_args = list(),
			    drop_terminal = (type[1] == "extended"),
			    tnex = (type[1] == "extended") + 1, 
			    newpage = TRUE,
			    pop = TRUE,
			    ...) {

    ### plot BinaryTree objects

    ### extract tree
    ptr <- x@tree
    ### total number of terminal nodes
    nx <- nterminal(ptr)
    ### maximal depth of the tree
    ny <- maxdepth(ptr)

    ### compute default settings
    type <- match.arg(type)
    if (type == "simple") {
        if (is.null(terminal_panel)) 
            terminal_panel <- node_terminal
        if (is.null(tnex)) tnex <- 1
    } else {
        if (is.null(terminal_panel))
            terminal_panel <- switch(class(response(x)[[1]])[1],
	                             "Surv" = node_surv,
                                     "factor" = node_barplot,
                                     "was_ordered" = node_barplot,
                                     "ordered" = node_barplot,
                                     node_boxplot)
        if (is.null(tnex)) tnex <- 2
    }

    ## setup newpage
    if (newpage) grid.newpage()

    ## setup root viewport
    root_vp <- viewport(layout = grid.layout(3, 3, 
    			height = unit(c(ifelse(is.null(main), 0, 3), 1, 1), 
                                      c("lines", "null", "lines")),
    			width = unit(c(1, 1, 1), 
                                     c("lines", "null", "lines"))), 
    			name = "root")       
    pushViewport(root_vp)
  
    ## viewport for main title (if any)
    if (!is.null(main)) {
        main_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1, 
                            name = "main")
        pushViewport(main_vp)
        grid.text(y = unit(1, "lines"), main, just = "center")
        upViewport()
    }

    ## setup viewport for tree
    tree_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
    			xscale = c(0, nx), yscale = c(0, ny + (tnex - 1)), 
                        name = "tree")
    pushViewport(tree_vp)

    ### setup panel functions (if necessary)
    ### the heuristic is as follows: If the first argument
    ### is `ctreeobj' than we assume a panel generating function, 
    ### otherwise the function is treated as a panel function
    if(inherits(terminal_panel, "grapcon_generator"))
      terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
    if(inherits(inner_panel, "grapcon_generator"))
      inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
    if(inherits(edge_panel, "grapcon_generator"))
      edge_panel <- do.call("edge_panel", c(list(x), as.list(ep_args)))


    if((nx <= 1 & ny <= 1)) {
      pushViewport(plotViewport(margins = rep(1.5, 4), name = paste("Node", ptr$nodeID, sep = "")))
      terminal_panel(ptr)    
    } else {
      ## call the workhorse
      plotTree(ptr,
        xlim = c(0, nx), ylim = c(0, ny - 0.5 + (tnex - 1)),
        nx = nx, ny = ny, 
        terminal_panel = terminal_panel,
        inner_panel = inner_panel,
        edge_panel = edge_panel,
        tnex = tnex,
        drop_terminal = drop_terminal,
        debug = FALSE)
    }
    upViewport()
    if (pop) popViewport() else upViewport()
}

#A number of binary tree node functions for the wikileaks paper
# we use meansddevplot3 in the end
#author: Thomas Rusch, 2010-2011
node_negbinplot_trans<- function(mobobj,
                        # col = "black",
                         col=gray(0.7),
                         pow=1/2,
                         fill = gray(0.7),
			 beside = TRUE,
		         ymax = NULL,
		         ylines = NULL,
		         widths = 1,
		         gap = NULL,
		         id = TRUE,
                         res=1)
{   

    y <- response(mobobj)[,1]
 
    if(is.null(ymax)) ymax <- max(y) * 1.1
    if(is.null(gap)) gap <- 1
    if(is.null(fill)) fill <- gray.colors(length(ylevels))
    if(is.null(ylines)) ylines <- if(beside) c(3, 2) else c(1.5, 2.5)

    ### panel function for barplots in nodes
    rval <- function(node) {
        ## parameter setup
        ysh <- rep.int(y,node$weights)
        ylevels <- 0:max(ysh)
        tmp1 <- tabulate(ysh,nbins=max(ysh))
        tmp2 <- sum(ysh==0)
        tmp <- c(tmp2,tmp1)
        pred <- tmp/sum(tmp)

	modpred <- dnbinom(ylevels,size=node$model$theta,mu=exp(node$model$coefficients))

        np <- length(pred)
	nc <- np 

	fill <- rep(fill, length.out = np)	
        widths <- rep(widths, length.out = nc)
	col <- rep(col, length.out = nc)
	ylines <- rep(ylines, length.out = 2)

	gap <- gap * sum(widths)
        yscale <- c(0, max(modpred,pred)*1.1)
        xscale <- (c(0, sum(widths) + (nc+1)*gap))^pow

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines[1], 1, ylines[2]), c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_poisplot", node$nodeID, sep = ""))

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), "n = "),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste("node_poisplot", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)

  	  xcenter <- cumsum(widths+gap) - widths/2
	  for (i in 1:np) {
            grid.rect(x = (xcenter[i])^pow, y = 0, height = pred[i], 
                      width = widths[i],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
            grid.points(x=(xcenter[i])^pow,modpred[i],pch=20,size=unit(0.3,"char"),gp=gpar(col="black",fill="black"))
          }
       grid.xspline((xcenter)^pow,modpred,default.units="native")
        woiserdenn2 <- NA
        ylev2 <- NA
          for(i in 1:floor(length(xcenter)/res))
           {
              woiserdenn2[i] <- ((xcenter)[i*res])^pow
              ylev2[i] <- ylevels[i*res]
            }
                                        #X axis 
        if(length(xcenter) > 1) grid.xaxis(at=woiserdenn2,label=FALSE)
	 grid.text(ylev2,x="", y = unit(-1, "lines"), 
                   just = c("center", "top"),
	            default.units = "native", check.overlap = TRUE)
          grid.yaxis()
       
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }
    
    return(rval)
}
class(node_negbinplot_trans) <- "grapcon_generator"

node_nbindistplot<- function(mobobj,
                         col = "black",
			 beside = TRUE,
                         ymax = NULL,
		         ylines = NULL,
		         widths = 1,
		         id = TRUE
                         )
{

    y <- response(mobobj)[,1]

     if(is.null(ymax)) ymax <- max(y)
    
    ### panel function for barplots in nodes
    rval <- function(node) {
        ## parameter setup
        ysh <- rep.int(y,node$weights)
	modpred <- ymax
        xmax <- max(ysh)
        yscale <- c(0, modpred*1.1)
        xscale <- c(0-0.1*xmax,xmax+0.1*xmax)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines[1], 1, ylines[2]), c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_meansd", node$nodeID, sep = ""))

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), "n = "),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste("node_meansd", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)
        
        distplot(ysh,type="nbinomial",pop=TRUE,newpage=FALSE,legend=FALSE,xlab="",ylab="",main=".")

        upViewport(2)
    }
    return(rval)
}
class(node_nbindistplot) <- "grapcon_generator"

node_meansddevplot3<- function(mobobj,
                         col = "black",
    			 beside = TRUE,
                         ymax = NULL,
		         ylines = NULL,
		         widths = 1,
		         id = FALSE
                         )
{

    thetis <-unlist(lapply(summary(mobobj),function(x) x$theta)) 
    meansis <- exp(coefficients(mobobj))
    devis <- unlist(lapply(summary(mobobj),function(x) x$deviance)) 
    dfis <- unlist(lapply(summary(mobobj),function(x) x$df.residual))
    fitsis <- devis/dfis
    varis <- meansis+(meansis^2/thetis)
    sdis <- sqrt(varis)

    xmax <- max(sdis)
    maxid <- which(sdis==max(sdis))
    xmin <- 0 

    if(is.null(ymax)) ymax <- max(fitsis)

    if(is.null(ylines)) ylines <- c(1.7, 2)
    
    ### panel function for barplots in nodes
    rval <- function(node) {
        ## parameter setup

        meani <- exp(node$model$coefficients)
        ti <- node$model$theta
        vari <- meani+(meani^2/ti)
        sdi <-sqrt(vari)
        modpred <- node$model$deviance/node$model$df.residual

        np <- 1
	nc <- np 

	ylines <- rep(ylines, length.out = 2)

        yscale <- c(0, ymax+0.1*1/ymax)
        xscale <- c(0-0.1*xmax,xmax+0.1*xmax)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines[1], 1, ylines[2]), c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_meansd", node$nodeID, sep = ""))

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
      	mainlab <- paste(ifelse(id, paste("R", ri, "(n = "), "n = "),	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste("node_meansd", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)
       grid.polyline(x=c(meani,meani,0,sdi),y=c(0.05,modpred+0.05,0.05,0.05),id=c(1,1,2,2),default.units="native",gp=gpar(col=col))
         grid.xaxis(at=c(round(meani,digits=2),round(sdi,digits=2)),label=FALSE)
        grid.text(round(meani,digits=2),x=unit(meani/xmax,"npc"),y=unit(-1,"lines"))
        grid.text(round(sdi,digits=2),x=unit(sdi/xmax,"npc"),y=unit(-2,"lines"))
        grid.yaxis(at=modpred+0.05,label=round(modpred,digits=2))
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }
    return(rval)
}
class(node_meansddevplot3) <- "grapcon_generator"

