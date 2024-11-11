options(stringsAsFactors = FALSE)  




AppUI <- function() {
  currentfilePath <- dirname(rstudioapi::getSourceEditorContext()$path)
  pkgname <-(pkg <- strsplit(currentfilePath, '/'))[[1]][length(pkg[[1]])-1]; pkgname 
  shiny::runApp(system.file("shiny", package=pkgname))} #
#' @name appUI
#' @title Launch interactive User Interface
#' @description  AppUI initiates in the web browser an interactive user interface of ...  This user interface enables users to easily perform nearly all standard ... functions in ... package.
#' @usage appUI()
#' @param sfPower numerical value of ... for
#' @return A user interface will be shown in users' default web browser.
#' @import shiny ggplot2 ggrepel  
#' @export 
appU <- function() { loc <- gsub('.*:', '', getAnywhere("AppUI")$where[1]) 
shiny::runApp(system.file("shiny", package=loc))  }



#' @export 
install_if_missing <- function(pkgs = NULL, bioc_pkgs = NULL) {
  if (!is.null(pkgs)) {
    missing_cran_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
    if (length(missing_cran_pkgs) > 0) {
      install.packages(missing_cran_pkgs)
      message(paste("You have successfully installed the following CRAN packages:", paste(missing_cran_pkgs, collapse = ", ")))
    } else {
      message("All specified CRAN packages are already installed.")
    }
  }
  if (!is.null(bioc_pkgs)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    missing_bioc_pkgs <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
    if (length(missing_bioc_pkgs) > 0) {
      BiocManager::install(missing_bioc_pkgs)
      message(paste("You have successfully installed the following Bioconductor packages:", paste(missing_bioc_pkgs, collapse = ", ")))
    } else {
      message("All specified Bioconductor packages are already installed.")
    }
  }
}


#' @export 
loadUrl <- function(url, downloadPath = NA, sep=c("RData"," ", "," , "\t", ";", "xls", "gsheet"), ...) {
  cat('onedrive: copy link\n googlesheet: share-> Anyone with the link\n sep: "RData", ..."xls", "gsheet"\n')
  if(!is.na(downloadPath))  { tmpFile <- downloadPath
  
  } else { tmpFile <- tempfile()  }
  url2 <- gsub("e=.*", "download=1", url)
  download.file(url2, tmpFile, mode="wb") # For windows user: mode = "wb" cf. binary
  sep <- match.arg(sep)
  if(sep == "RData") {
    print(tmpFile)
    tmpFile <-  gsub("\\\\", "/", tmpFile)
    justLoaded <- try(load(tmpFile), silent = T); 
    try(assign(justLoaded, eval(as.symbol(justLoaded)),.GlobalEnv ), silent = T);
    if(class(justLoaded)=="try-error"){ justLoaded <- try(read.delim(tmpFile, ...), silent = T); message("Need 'sep' argument, is it txt file?")  }   
  } else if(sep == "xls") {
    is.installed('readxl')
    justLoaded <- try(read_excel(tmpFile,...), silent = T)
    
  } else if(sep == "gsheet") {
    is.installed('gsheet')
    cat('gsheet should be public, click share-> Anyone with the link')
    justLoaded <- gsheet2tbl(url,...)
  } else {
    justLoaded <- try(read.delim(tmpFile, sep=sep, ...), silent = T)  
  }
  justLoaded 
}



#' @export 
browseEntrez <- function(entrezIDs) {
  for(i in entrezIDs) {
    browseURL(paste0("https://www.ncbi.nlm.nih.gov/gene/", i))
  }
}


#' @export 
peep <- function(x, boxplot = F ) { 
  if(is.null(dim(x))) { if(length(x) > 10)  { print(x[1:10]) } else { print(x) }  } else if (dim(x)[1] >=10 && dim(x)[2]>=5 ){ print(x[1:5, 1:3]); boxplot(x[1:5, 1:3]) } else {print(head(x)); boxplot(x)} }


#' @export 
normalize.q <- function(x= data.frame(matrix(sample(12, replace = T), 4)), filter.sd.quantile = 0.1, tied = c("average", "min", "max"), verbose = T ) {
  # compare to normalize.quantiles, 1. accept data.frame 2. tie control option:"average", "min", "max"  3. sd.filter 4. peep & plot & verbose...
  
  x <- x[rowSums(x)>0, ]  
  x <- x[apply(x,1,sd) >= quantile(apply(x,1,sd), filter.sd.quantile), ]  
  cat(sprintf("\nrowSums(x)==0, =<quantile(sd(row), %s) were filtered\n\n", filter.sd.quantile))
  
  tied <- match.arg(tied)  
  rank <- apply(x, 2, rank,ties.method="min"); 
  if(any(tied %in% c("average", "max"))) rank.max <- apply(x, 2, rank,ties.method="max"); 
  sorted <- apply(x, 2, sort)
  sorted.row.mean <- apply(sorted, 1, mean); 
  x2 <- apply(rank, 2, function(x) sorted.row.mean[x])
  if(any(tied %in% c("average", "max"))) x2.max <- apply(rank.max, 2, function(x) sorted.row.mean[x])
  if(tied=="average") { x2 <- (x2+x2.max)/2 } else if (tied=="max"){x2 <- x2.max } else { }
  
  if( class(x) == "data.frame") { x2 <- as.data.frame(x2); rownames(x2) <- rownames(x) }
  if(verbose) {
    op <- par(no.readonly = T); par(mfrow=c(1,2), mar=c(3,3,1,1))
    cat('Original matrix or data.frame\n'); peep(x, T)
    cat('Sort the original matrix from lowest to highest\n'); peep(rank)
    cat('Determine the ranks of original matix\n');peep(sorted)
    cat('\nCalculate the means\n\n'); peep(sorted.row.mean)
    cat('\n\nFinally substitute the means into our ranked matrix\n'); peep(x2, T)
    cat(sprintf('If the values were tied, %s is used\n\n', tied))
    par(op)
    'In the example on Wikipedia, if the values were tied, the min value is used but in the normalize.quantiles() function, the average is used'
  }
  x2
}




#' @export 
DEGs <- function(Exp, cl, adj.pval = 0.1,  logFC = 2, geomTextN=5, heatmapUpN = 25, plotDEG =T, multipleRegression=F, rowCounts=F, meanFilter=10, PDF=T, cname='temp', show_column_names=T, rect_gp = gpar(col = NA, lty = 1, lwd = 0.2)) {
  try(dev.off(), silent = T)
  
  is.installed(c('ggplot2', 'ggrepel'))
  is.installed.bioconductor(c('limma', 'ComplexHeatmap'))
  
  
  if(rowCounts) { Exp <- Exp[apply(Exp, 1, mean) > meanFilter, ]; Exp <- voom(Exp, plot = T) }
  
  if(multipleRegression) { fit <-eBayes(lmFit(Exp, model.matrix(~ .,cl))); print(topTable(fit, 2))                                     } else {
    fit <-eBayes(lmFit(Exp, model.matrix(~ .,cl[, 1, drop=F]))); print(topTable(fit, 2)) } 
  
  tT <- topTable(fit, number = dim(Exp)[1])
  tT$Gene <- rownames(tT)
  tT.up <- tT[order(tT$logFC, decreasing = T ),]; tT.down<- tT[order(tT$logFC),]
  tT.filter <- data.frame(tT[inx<-(tT$adj.P.Val<adj.pval) & (abs(tT$logFC) > logFC ), ]); print(tT.filter)
  
  if(PDF) {
    pdf(file = file.path(getwd(),sprintf("%s.pdf", cname )), width=5, height=5)  }
  if(plotDEG) {
    if(rowCounts) Exp <- Exp$E
    
    if(any(colnames(tT) == "logFC" && dim(tT.filter)[1] != 0) ) {
      require(ggplot2); require(ggrepel)
      tT$Cutoff_value <- c("Not Sig", sprintf("FDR < %s & logFC > %s", adj.pval, logFC))[as.numeric(inx)+1]
      gplot <- ggplot(tT, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(color = Cutoff_value)) + labs(title ="c") + scale_color_manual(values = c("red", "grey")) +  theme_bw(base_size = 12) + theme(legend.position = "bottom") + geom_hline(yintercept= -log10(adj.pval), linetype="dashed", color = "#FF000050") + geom_vline(xintercept= c(logFC, -logFC), linetype="dashed", color = "#FF000050") 
      
      g <- gplot + geom_text_repel(    
        data = dd <- rbind(tT.up[1:geomTextN, ], tT.down[1:geomTextN, ]),
        aes(label = Gene),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) }
    print(g)
    
    if(dim(Exp)[1] >= heatmapUpN*2 ) {
      bluered <- colorRampPalette(c("blue", "white", "red"))(256)
      stats::heatmap(Exp[rbind(tT.down[1:heatmapUpN, ],tT.up[heatmapUpN:1, ])$Gene,], col = bluered, scale = "row", main = sprintf("top%s", heatmapUpN*2), Colv = NA, Rowv=NA    )
      
      # if(HUGO) { colnames(d)[1:dim(matGS)[1]] <- .mapid(colnames(d)[1:dim(matGS)[1]]) }
      
      # colse=c("#00000020", "#000000", "#0000BF10", "#0000BF30", "#0000BF50", "#0000BF70","#0000BF90","#0000BF")
      # colTemp <- colse[as.numeric(as.factor(cl[,1]))]
      # names(colTemp ) <- cl[,1]
      # colTemp<-list(colTemp); names(colTemp) <- colnames(cl)[1] 
      
      # h <- Heatmap( t(scale(t(d<-Exp[rbind(tT.up[1:heatmapUpN, ], tT.down[heatmapUpN:1, ])$Gene,]))),  col = bluered, name="Exprs", rect_gp = rect_gp, cluster_rows = T, cluster_columns = T, show_row_names = T, show_column_names=show_column_names, row_names_gp =gpar(fontsize = 5), split = data.frame(cyl = factor(c(rep('UP', heatmapUpN), rep('DOWN', heatmapUpN)), levels=c('UP','DOWN' ))),gap = unit(1.5, "mm"), top_annotation = HeatmapAnnotation(df=cl, col=colTemp, show_annotation_name = T, annotation_name_gp = gpar(fontsize = 7), annotation_name_side ='left', annotation_height=c(1.5), use_raster = T ) ) 

      # draw(h)
    }
  }
  if(PDF) { dev.off() 
  }
  return(list(fit = fit, tT.filter= tT.filter, tT.up=tT.up, tT.down=tT.down))

  }



                                                         
#' @export                                                      
RP.custom <- function(s,FDRcutoff=.1) {
  
  is.installed.bioconductor(c('RankProd'))

  for(i in 1: length(s)) { s[[i]]$y <- as.numeric(as.factor(s[[i]]$y[,1]))-1 
  rownames(s[[i]]$x) <-  gsub("///.*", "", rownames(s[[i]]$x))
  }  
  mainTitle = sprintf("   %s and %s others", names(s)[1], (length(s)-1) )
  if(length(do.call("c", lapply(s, function(x)x$y))) >= 100) { RandomPairs = 100} else { RandomPairs = NA
  }; RandomPairs
    tt <-list(); tt.origin <-c(); tt.cl<-c()
  for( i in  1: length(s))
  { tt[[i]] <- s[[i]][[1]] 
  tt.origin <- c(tt.origin , rep(i, dim(s[[i]][[1]])[2] ))
  tt.cl <- c(tt.cl, s[[i]]$y)
  }
  ttt <- do.call(cbind, tt); dim(ttt)
  RP.adv.out <- RP.advance(ttt, tt.cl, tt.origin, logged=T, rand=123, RandomPairs = RandomPairs) 
  RP.adv.out.ind=list()
  pfp.cut.off <- FDRcutoff  
  for( i in  1: length(s)) {
    RP.adv.out.ind[[i]] <- RP.advance(s[[i]][[1]], s[[i]]$y, rep(1, length(s[[i]]$y)), RandomPairs = RandomPairs, logged=T, rand=123) 
    RP.adv.out.ind[[i]]$up <- RP.adv.out.ind[[i]]$AveFC[RP.adv.out.ind[[i]]$pfp[, 1]< pfp.cut.off, , drop=F] #  class1 < class2
    RP.adv.out.ind[[i]]$down <- RP.adv.out.ind[[i]]$AveFC[RP.adv.out.ind[[i]]$pfp[, 2]< pfp.cut.off, , drop=F] #  class1 > class2
    RP.adv.out.ind[[i]]$updown <- rbind(RP.adv.out.ind[[i]]$up, RP.adv.out.ind[[i]]$down) 
  }    # fold changes of average expressions (class1/class2). log fold-change if data has been log transformed, original fold change otherwise
  RP.adv.out$up <- RP.adv.out$AveFC[RP.adv.out$pfp[, 1]< pfp.cut.off, , drop=F] #  class1 < class2
  RP.adv.out$down <- RP.adv.out$AveFC[RP.adv.out$pfp[, 2]< pfp.cut.off, , drop=F] #  class1 > class2
  RP.adv.out.ind[[i]]$down <- RP.adv.out.ind[[i]]$AveFC[RP.adv.out.ind[[i]]$pfp[, 2]< pfp.cut.off, , drop=F]
  RP.adv.out$updown <- rbind(RP.adv.out$up, RP.adv.out$down)
  list(RP.adv.out,  RP.adv.out.ind)
}                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
                                                         
