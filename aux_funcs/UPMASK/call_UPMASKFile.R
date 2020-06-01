#! /usr/bin/Rscript


UPMASKfile <- function(filenameWithPathInput, filenameWithPathOuput, 
                      positionDataIndexes=c(1,2),
                      photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
                      photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28),
                      threshold=1, maxIter=20, starsPerClust_kmeans=50, nstarts_kmeans=50, 
            nRuns=5, runInParallel=FALSE, paralelization="multicore", 
            independent=TRUE, verbose=FALSE, autoCalibrated=FALSE, 
            considerErrors=FALSE, finalXYCut=FALSE, fileWithHeader=FALSE, 
                      nDimsToKeep=4, dimRed="PCA", scale=TRUE, sep) {
  
  # 
  # User :: Interface
  #
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Starting UPMASK...\n")
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
  # 
  # Data I/O :: Perform File Input
  # 
  if(verbose) {
    # cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Reading the input table from: ",filenameWithPathInput,"\n")
  }
  # Load the file
  ocdata_full <- read.table(filenameWithPathInput, header=fileWithHeader, sep=sep)
  # if(verbose) {
  #   cat(paste("-------------------------------------------------------------------\n"))
  # }
  
  #
  # Science :: Run the UPMASK caller function
  #
  resultsTable <- UPMASKdata(ocdata_full, 
                      positionDataIndexes=positionDataIndexes,
                        photometricDataIndexes=photometricDataIndexes,
                        photometricErrorDataIndexes=photometricErrorDataIndexes,
                        threshold=threshold, maxIter=maxIter, 
              starsPerClust_kmeans=starsPerClust_kmeans, 
              nstarts_kmeans=nstarts_kmeans, nRuns=nRuns, 
              runInParallel=runInParallel, paralelization=paralelization, 
              independent=independent, verbose=verbose, 
              considerErrors=considerErrors, finalXYCut=finalXYCut, 
                        nDimsToKeep=nDimsToKeep, dimRed=dimRed, scale=scale)

  # 
  # Data I/O :: Perform File Output
  # 
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Writing the output table at: \n\t",filenameWithPathOuput,"\n")
  }
  write.table(resultsTable, filenameWithPathOuput, sep="    ", 
              col.names=fileWithHeader, quote=FALSE, row.names=FALSE)   

  # 
  # User :: Interface
  #
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Done!\n")
    cat(paste("-------------------------------------------------------------------\n"))
  }  
}



UPMASKdata <- function(dataTable,
                      positionDataIndexes=c(1,2),
                      photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
                      photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28),
                      threshold=1, classAlgol="kmeans", maxIter=25,
                      starsPerClust_kmeans=25, nstarts_kmeans=50, nRuns=8,
                      runInParallel=FALSE, paralelization="multicore", independent=TRUE,
                      verbose=FALSE, autoCalibrated=FALSE, considerErrors=FALSE,
                      finalXYCut=FALSE, nDimsToKeep=4, dimRed="PCA", scale=TRUE) {
  
  # 
  # User :: Interface
  #
  # if(verbose) {
  #   cat(paste("-------------------------------------------------------------------\n"))
  #   cat(      " Starting UPMASK analysis!\n")
  #   cat(      " UPMASK kernels v. 1.2\n")
  #   cat(paste("-------------------------------------------------------------------\n"))
  # }
  
  # 
  # Science
  #
  # Create the smart look-up table for the random analysis
  # Unfortunatelly this must be a global variable so it can be shared among 
  # parallel processes
  stcon <- create_smartTable() 
    
  # Create a place-holder list
  pp <- list()
  
  # Run UPMASK outerloop (and innerloop, called inside the outerloop)
  if (!runInParallel) {
    if(independent) {
      for(i in 1:nRuns) {
        pp[[length(pp)+1]] <- outerLoop(dataTable, 
                                positionDataIndexes=positionDataIndexes,
                                photometricDataIndexes=photometricDataIndexes, 
                                photometricErrorDataIndexes=photometricErrorDataIndexes,
                                threshold=threshold,
                                maxIter=maxIter, plotIter=FALSE, 
                                starsPerClust_kmeans=starsPerClust_kmeans, 
                                nstarts_kmeans=nstarts_kmeans,  
                                verbose=verbose, finalXYCut=finalXYCut, 
                                autoCalibrated=autoCalibrated, 
                                considerErrors=considerErrors, run=i, 
                    smartTableDB=stcon, nDimsToKeep=nDimsToKeep, dimRed=dimRed, scale=scale)
      }
    }
  } else {
    if(paralelization=="multicore") {
      # If the user wants to run in an SMP machine, then we want to use occupy the cores!
      # library(parallel)
      # Ok, now it is just a matter of running the code using a list
      pp <- mclapply(1:nRuns, function(x) { outerLoop(dataTable, 
                                positionDataIndexes=positionDataIndexes,
                                photometricDataIndexes=photometricDataIndexes, 
                                photometricErrorDataIndexes=photometricErrorDataIndexes,
                                threshold=threshold, 
                                maxIter=maxIter, plotIter=FALSE, 
                                starsPerClust_kmeans=starsPerClust_kmeans, 
                                nstarts_kmeans=nstarts_kmeans, 
                                verbose=verbose, finalXYCut=finalXYCut, 
                                autoCalibrated=autoCalibrated, 
                                considerErrors=considerErrors, run=x, 
                                smartTableDB=stcon, nDimsToKeep=nDimsToKeep, 
                                            dimRed=dimRed, scale=scale)} )
    } else if(paralelization=="MPIcluster") {
      # if the user wants to run in a cluster
      # then we need to use MPI
      stop(" MPIcluster is not implemented yet.\n Aborting cowardly!")
    }
  }
  
  ## To do :: check if this is necessary, after all the user
  ## should perform the cuts himself
  dataTable <- performCuts(dataTable) 
  
  # Organize the results into a single data frame
  mergedResults <- data.frame(id=1:length(dataTable[,1]))
  for(i in 1:length(pp)) {
    if(length(pp[[i]])>1) {
      mergedResults <- data.frame(mergedResults, class=pp[[i]]$class)
    }
  }
  
  # Compute a frequency to output as a frequentist probability
  freq <- vector("double",length(mergedResults$id))
  nConverged <- length(mergedResults)-1
  if(nConverged > 0) {
#    for(i in 1:length(mergedResults$id)) {
    for(i in 1:nrow(mergedResults)) { 
      freq[i] <- sum(mergedResults[i,2:(nConverged+1)])/nConverged
    }
    mergedResults <- cbind(mergedResults, probability=freq)
  } else {
    #    for(i in 1:length(mergedResults$id)) {
    for(i in 1:nrow(mergedResults)) { 
      freq[i] <- 0
    }
    mergedResults <- cbind(mergedResults, probability=freq)
  }
  
  # Organize the results into a single data frame
  ocdata_out <- data.frame(dataTable, mergedResults[,2:length(mergedResults)])
  
  # Clean the global variable (argh!) storing the smart lookup table...
  #rm(smartTable, pos=globalenv())
  
  # Close the connection to the database storing the smart lookup table
  dbDisconnect(stcon)
  
  # 
  # User :: Interface
  #
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " UPMASK analysis is finished!\n")
    cat(paste("-------------------------------------------------------------------\n"))
  }
    
  # That's all folks!
  return(ocdata_out)
}

library('DBI')
library(RSQLite)
create_smartTable <- function() {
#  con <- dbConnect(RSQLite::SQLite(), ":memory:")
  con <- dbConnect(dbDriver("SQLite"), dbname = tempfile())
  dbWriteTable(con, "smartTable", data.frame(nstar=c(0), mean=c(0), sd=c(0)), row.names = F)
  return(con)
}



outerLoop <- function(ocdata_full, 
                      positionDataIndexes=c(1,2),
                      photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
                      photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28),
                      threshold=1, maxIter=25, plotIter=FALSE, 
                      verbose=FALSE, starsPerClust_kmeans=50, nstarts_kmeans=50, 
                      finalXYCut=FALSE, 
                      autoCalibrated=FALSE, considerErrors=FALSE, run=0, 
                      smartTableDB, nDimsToKeep=4, dimRed="PCA", scale=TRUE) {
  
  # The version without autoThresold is deprecated
  autoThreshold <- TRUE
  
  # Create an internal index
  idx <- which(ocdata_full$field!=0)
  ocdata_full$field[idx] <- rep(1, times=length(idx))
  
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(paste(" Points in the original data   :",length(ocdata_full[,1]),"\n"))
  }
  
  # Perform cuts in the data if necessary
  ocdata_full <- performCuts(ocdata_full)
  
  # Take errors in the data table into account if the user request
  if(considerErrors) {
    ocdata_full <- takeErrorsIntoAccount(ocdata_full, photometricDataIndexes, 
                                         photometricErrorDataIndexes)
  }
  
  if(verbose) {
    cat(paste(" Points after the data cuts    :",length(ocdata_full[,1]),"\n"))
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Starting iterations...\n")
    cat(paste(" Maximum iterations :",maxIter,"\n"))
    if(plotIter) {
      cat(paste(" You have choosen to see the iteration's plots.\n"))
    } else {
      cat(paste(" You have choosen NOT to see the iteration's plots.\n"))
    }
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
  # Start iterations
  verb <- 2
  if(!verbose) { 
    verb <- 0 
  }
  ocdata_full <- data.frame(ocdata_full, id=(1:nrow(ocdata_full))) # create an id
  ocdata_out <- ocdata_full # it start as the full list
  ocdata_res <- ocdata_full # it start as the full list
  solConv <- FALSE
  i <- 0
  for(i in 1:maxIter) {
    oldSize <- length(ocdata_res[,1])
    # Select the columns to analyse
    ocdata <- ocdata_res[,photometricDataIndexes]
    
    # Iterate the inner loop
    ocdata_res <- innerLoop(ocdata_res, ocdata, classAlgol="kmeans", 
                            autoThresholdLevel=threshold, iiter=i, plotIter=FALSE, 
                            starsPerClust_kmeans=starsPerClust_kmeans, 
                            nstarts_kmeans=nstarts_kmeans, verbosity=verb, runId=run, 
                            autoCalibrated=autoCalibrated, 
                            positionDataIndexes=positionDataIndexes, smartTableDB=smartTableDB, 
                            nDimsToKeep=nDimsToKeep, dimRed=dimRed, scale=scale)
    
    # Flag this iteration's cluster stars
    member <- data.frame(m=rep(0,length(ocdata_out[,1])))
    member$m[ocdata_res$id] <- 1
    ocdata_out <- data.frame(ocdata_out, member)
    
    # Check if the solution converged at this iteration
    if(length(ocdata_res[,1])==oldSize) {
      solConv <- TRUE
      break # not an elegant solution
    }
    
    if(length(ocdata_res[,1])==0) {
      solConv <- FALSE
      break # not an elegant solution
    }
  }
  
  # If the solution converged, then let the user know (if verbose)
  # If it did not, flag everything as non member stars
  if (solConv) {
    if(verbose) {
      cat(paste(" Convergence found at iteration,",i,"!\n"))
    }
  } else {
    if(verbose) {
        cat(" Sorry, the system never converged...\n But I am not aborting!!\n")
    }
    # Flag this iteration's cluster stars as non-members
    member <- data.frame(m=rep(0,length(ocdata_out[,1])))
    ocdata_out <- data.frame(ocdata_out, member)
  }
  
  # If the user wants to perform a final cut in order to get the spatially clustered
  # members only, then do it.
  if(finalXYCut) {
    if(verbose) {
      cat(paste("-------------------------------------------------------------------\n"))
      cat(      " Selecting the stars at the highest X-Y density levels...\n")
    }
    ocdata_out <- getStarsAtHighestDensityRegion(ocdata_out, threshold=3, verbose=FALSE)
  } else {
    ocdata_out <- data.frame(ocdata_out, finalClass=ocdata_out[,length(ocdata_out)] )   
  }
    
  # That's all folks!  
  return(data.frame(id=ocdata_out$id, class=ocdata_out$finalClass))
}

performCuts <- function(originalData) {

  # This function can be tailored to the user needs, if it is necessary to perform
  # data cuts due to missing data bands, or certain flags, etc.
  # originalData <- subset(originalData, originalData$U<30) # Cut points without U magnitude

  return(originalData)
}


getStarsAtHighestDensityRegion <- function(ocdata_out, threshold=2, posIdx=c(1,2), plotAnalysis=FALSE, verbose=FALSE) {
  # Get the cluster stars
  # oc_data_clean <- subset(ocdata_out, ocdata_out$resMclust.class==1)
  oc_data_clean <- subset(ocdata_out, ocdata_out[,ncol(ocdata_out)]==1)
  
  # Perform the Kde2d in the X-Y space
  dens_map <- kde2d(oc_data_clean[,posIdx[1]],oc_data_clean[,posIdx[2]],n=50, 
                    lims=c(range(ocdata_out[,posIdx[1]]), range(ocdata_out[,posIdx[2]])))
  
  # Compute the average density and its sd (using an iterative method)
  dens_vals <- as.vector(dens_map$z)
  stat_dens <- meanThreeSigRej(dens_vals)
  
  # Compute the average density and its sd, difference from the random fields
  # stat_dens <- analyse_randomKde2d(nfields=2000, nstars=length(oc_data_clean$x), 
  #               (max(oc_data_clean$x)-min(oc_data_clean$x)),
  #               (max(oc_data_clean$y)-min(oc_data_clean$y)), 
  #               nKde=50, showStats=FALSE, returnStats=TRUE)
  
  
  loopFlag <- TRUE
  while(loopFlag) { 
    # Create a density map with flags for the selected cluster / field regions
    if(verbose) {
      print(stat_dens)
    }
    flagged_dens_map <- dens_map
    flagged_dens_map$z[which(flagged_dens_map$z<(stat_dens$mean+threshold*stat_dens$sd))] <- 0
    flagged_dens_map$z[which(flagged_dens_map$z>=(stat_dens$mean+threshold*stat_dens$sd))] <- 1
    
    # Select the stars inside the flagged region
    oc_data_clean_tmp <- oc_data_clean
    # before anything else, reorganize the data
    dmap <- flagged_dens_map
    xvec <- rep(dmap$x,times=length(dmap$y))
    yvec <- vector("double",(length(dmap$x)*length(dmap$y)) )
    kk <- 1
    for (j in 1:length(dmap$y)) {
      yvec[kk:(kk+length(dmap$x)-1)] <- rep(dmap$y[j], length(dmap$x))
      kk <- kk+length(dmap$x)
    }
    zvec <- as.vector(dmap$z)
    flagged_region <- data.frame(x=xvec, y=yvec, z=zvec)
    # first let's get rid of the easy ones...
    flagged_region <- subset(flagged_region, flagged_region$z==1)
    dxP2 <- abs(flagged_dens_map$x[2]-flagged_dens_map$x[1])/2
    dyP2 <- abs(flagged_dens_map$y[2]-flagged_dens_map$y[1])/2
    oc_data_clean_tmp <- subset(oc_data_clean_tmp, oc_data_clean_tmp[,posIdx[1]]>=(min(flagged_region$x)-dxP2))
    oc_data_clean_tmp <- subset(oc_data_clean_tmp, oc_data_clean_tmp[,posIdx[1]]<=(max(flagged_region$x)+dxP2))
    oc_data_clean_tmp <- subset(oc_data_clean_tmp, oc_data_clean_tmp[,posIdx[2]]>=(min(flagged_region$y)-dyP2))
    oc_data_clean_tmp <- subset(oc_data_clean_tmp, oc_data_clean_tmp[,posIdx[2]]<=(max(flagged_region$y)+dyP2))
    
    loopFlag <- FALSE
    
    if(verbose) {
      print(oc_data_clean_tmp)
    } 
    
    if(length(oc_data_clean_tmp[,posIdx[1]])==0) {
      loopFlag=TRUE
      threshold<-threshold-1
      cat(paste(" WARNING: The thresholding for the final cluster selection in the X-Y density map was to high!\n"))
      cat(paste(" WARNING:     since no stars were present in the final cluster, I am lowering the threshold!\n"))
      cat(paste(" WARNING:     The NEW threshold will be",threshold,"\n"))
      if(threshold<0) {
        stop("ABORTING! Threshold for spatial clustering of member stars is less than zero!")
      }
    }
  }
  
  # ok, now check the hard ones one by one
  trueVec <- c()
  for(i in 1:length(oc_data_clean_tmp[,posIdx[1]])) { # for each star
    # check if this star is inside any of the flagged boxes
    isInside <- FALSE
    for(j in 1:length(flagged_region$z)) {
      if( (oc_data_clean_tmp[i,posIdx[1]] >= (flagged_region$x[j]-dxP2) ) &&
            (oc_data_clean_tmp[i,posIdx[1]] <= (flagged_region$x[j]+dxP2) ) &&
            (oc_data_clean_tmp[i,posIdx[2]] >= (flagged_region$y[j]-dyP2) ) &&
            (oc_data_clean_tmp[i,posIdx[2]] <= (flagged_region$y[j]+dyP2) )) {
        isInside <- TRUE
        break 
      }
    }
    if (isInside) {
      trueVec <- c(trueVec, i)
    }
  }
  oc_data_veryclean_tmp <- oc_data_clean_tmp[trueVec,]
  
  # Plot the results
  if(plotAnalysis) {
    clevels <- c((stat_dens$mean), (stat_dens$mean+1*stat_dens$sd), (stat_dens$mean+2*stat_dens$sd), 
                 (stat_dens$mean+3*stat_dens$sd), (stat_dens$mean+threshold*stat_dens$sd))    
    
    # flagged plot
    colvec <- rgb(0:1/2,0,0)
    dev.new()
    plot(ocdata_out[,posIdx[1]], ocdata_out[,posIdx[2]], pch=19, cex=0.1, col="black", type="n", main=paste("Final Open Cluster"), xlab="X", ylab="Y")
    image(flagged_dens_map, col=colvec, add=TRUE)
    contour(dens_map, levels=clevels, labels=c("mean", "1sd", "2sd","3sd", "user"), col=c("red", "red", "red", "red", "blue"), add=TRUE, labcex=1.5)
    points(oc_data_clean[,posIdx[1]], oc_data_clean[,posIdx[2]], pch=19, cex=0.2, col="white") # OC data  
    points(oc_data_clean_tmp[,posIdx[1]], oc_data_clean_tmp[,posIdx[2]], cex=1, col="blue") # OC data 
    points(oc_data_veryclean_tmp[,posIdx[1]], oc_data_veryclean_tmp[,posIdx[2]], pch=19, cex=1, col="green") # OC data  
    
    # normal plot
    collevels <- 256
    div <- 2*collevels
    colvec <- rgb(0:collevels/div,0:collevels/div,0:collevels/div)
    dev.new()
    plot(ocdata_out[,posIdx[1]], ocdata_out[,posIdx[2]], pch=19, cex=0.1, col="black", type="n", main=paste("Final Open Cluster"), xlab="X", ylab="Y")
    image(dens_map, col=colvec, add=TRUE)
    contour(dens_map, levels=clevels, labels=c("mean", "1sd", "2sd","3sd", "user"), col=c("red", "red", "red", "red", "blue"), add=TRUE, labcex=1.5)
    points(oc_data_clean[,posIdx[1]], oc_data_clean[,posIdx[2]], pch=19, cex=0.2, col="white") # OC data
    
  }
  
  retDf <- data.frame(ocdata_out, finalClass=rep(0,length(ocdata_out[,posIdx[1]])))
  
  retDf$finalClass[oc_data_veryclean_tmp$id] <- 1
  
  return(retDf) 
}


meanThreeSigRej <- function(vec, maxI=50, tolerance=0.001) {
  tvec <- vec
  newM <- mean(tvec)
  i<-0
  for(i in 1:maxI) {
    oldM <- newM
    tvec <- subset(tvec, tvec<(mean(tvec)+3*sd(tvec)))
    newM <- mean(tvec)
    newSd <- sd(tvec)
    # if the mean converged, stop the torture, please...
    if( abs(oldM-newM)/oldM < tolerance) { break }
  }
  if(i==maxI) {
    stop(" ARGH! The Iterative 3-sigma rejection reached the maximum iterations without converging...\n Aborting cowardly!")
  }
  return(data.frame(mean=newM,sd=newSd, convergenceAtIter=i))
}





library('MASS')
innerLoop <- function(ocdata_full, ocdata, classAlgol="kmeans", autoThresholdLevel=3, 
                         autoThreshold=TRUE, iiter=0, plotIter=FALSE, verbosity=1, 
                         starsPerClust_kmeans=50, nstarts_kmeans=50, runId=0, 
                         autoCalibrated=FALSE, stopIfEmpty=FALSE, 
                         positionDataIndexes=c(1,2), smartTableDB, nDimsToKeep=4, dimRed="PCA", scale=TRUE) {
  if(verbosity!=0) {
    cat(paste(" [runId:",runId,"] ITERATION:",iiter," RUNNING...\n"))
  }
  inSize <- length(ocdata_full)
  
  # Perform the dimensionality reduction step
  if(dimRed!="None") {
    if(verbosity>=1) {
      cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [1/3] Performing Dimensionality Reduction...\n"))
    }
    if(dimRed=="PCA") {
      if(verbosity>=1) {
        cat(paste(" [runId:",runId,"] ITERATION:",iiter," --       Using PCA...\n"))
      }
      ocdata_pca <- prcomp(ocdata, scale=scale, center=scale, cor=TRUE)
      ocdata_px <- predict(ocdata_pca)
    } else 
    if(dimRed=="LaplacianEigenmaps") {
      minPoints <- round(runif(1, min=5, max=starsPerClust_kmeans))
      if(starsPerClust_kmeans > 25) {
        minPoints <- round(runif(1, min=5, max=25))
      } 
      if(nrow(ocdata) <= minPoints) {
        minPoints <- round(nrow(ocdata) - 1)
      }
      if(verbosity>=1) {
        cat(paste(" [runId:",runId,"] ITERATION:",iiter," --       Using Laplacian Eigenmaps... pts =",minPoints," Data=",nrow(ocdata),"\n"))
      }
      leim <- dimRed::LaplacianEigenmaps()
      leim@stdpars$ndim <- nDimsToKeep
      leim@stdpars$knn <- minPoints
      if(scale==TRUE) {
        ocemb <- leim@fun(dimRed::dimRedData(scale(ocdata)), leim@stdpars)
      } else {
        ocemb <- leim@fun(dimRed::dimRedData(ocdata), leim@stdpars)
      }
      ocdata_px <- as.data.frame(ocemb@data@data)
    } else {
      stop("You should select PCA, LaplacianEigenmaps as the dimensionality reduction method, or implement your own. Alternatively, you can select None to ignore dimensionality reduction.")
    }
  } else {
    if(verbosity>=1) {
      cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [1/3] You have chosen not to perform dimensionality reduction...\n"))
    }
    if(scale==TRUE) {
      ocdata_px <- scale(ocdata)
    } else {
      ocdata_px <- ocdata
    }
  }
  #########
  
  # Perform the clustering step
  if (classAlgol=="kmeans") {
    if(verbosity>=1) {
      cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [2/3] Performing k-means clustering analysis in the transformed data space...\n"))
    }
    starsPerClust <- starsPerClust_kmeans
    nclust <- round(length(ocdata_full[,1])/starsPerClust)
    if(nclust > 1) {
      fit <- kmeans(ocdata_px[,1:nDimsToKeep], nclust, nstart=nstarts_kmeans, iter.max=100)
      # get cluster means
      aggregate(ocdata_px, by=list(fit$cluster), FUN=mean)
      # append cluster assignment
      ocdata_px <- data.frame(ocdata_px, resMclust.class=fit$cluster)
    } else {
      # if the number of predicted clusters is less than one, lets prevent the kmeans method from 
      # crashing the code by assigning a randomized choice as the result (this is expected to happen if
      # a very small number of stars are assigned at a certain iteration).
      ocdata_px <- data.frame(ocdata_px, resMclust.class=round(runif(length(ocdata_px[,1]), 0, 1)))
    } 
  } else {
     stop(" Error: the selected method for the clustering in the inner loop is not implemented.")
  }
  
  # Print Dimensionality reduction info and create plots -- KEPT FOR DEPURATION PURPOSES
  if(plotIter){
    dev.new()
    par(cex=0.3)
    #   plot(resMclust, data=ocdata_px[,1:4])
    pairs(ocdata_px[,1:4], pch=19, cex=0.2, col=rainbow(max(ocdata_px$resMclust.class))[ocdata_px$resMclust.class])
  }
  
  # Merge the colors, the dimensionality reduction columns and cluster classification with the original results
  ocdata_full_withColorsAndPcaCols <- data.frame(ocdata_full, ocdata, ocdata_px)
    
  # Select the classes with densities above the threshold
  # if(verbosity>=1) {
  cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [3/3] Found",max(ocdata_px$resMclust.class),"individual classes...\n"))
  #   if(verbosity==1) {
  #     cat(paste(" [runId:",runId,"] ITERATION:",iiter,"          You can get a coffee. This can take long! \n"))
  #   }
  # }
  
  vclass <- c()
  not_class <- c()
  for(i in 1:max(ocdata_px$resMclust.class)) {
    dfn <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==i)
        
    if(nrow(dfn)>2) {
      
      dif_max_mean <- kde2dForSubset(ocdata_full_withColorsAndPcaCols, setw=i, returnDistance=TRUE, showStats=FALSE, printPlots=FALSE, positionDataIndexes=positionDataIndexes)
      
      # First, get the thresholding level...
      if(autoThreshold) {
        # if(verbosity>=2) {
        #   cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- Class",i," -- Performing analysis of random fields...\n"))
        # }
        
        if(autoCalibrated) {
          # This is an experimental code for performing automatic calibration of the 
          # number of random realisations
          at <- analyse_randomKde2d_AutoCalibrated(
                 nstars=length(dfn$resMclust.class), 
                 (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])),
                 (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])), 
                 nKde=50, showStats=FALSE, returnStats=TRUE)
        } else {                    
          at <- analyse_randomKde2d_smart(
                 nfields=2000, nstars=length(dfn$resMclust.class), 
                 (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])),
                 (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])), 
                 nKde=50, showStats=FALSE, returnStats=TRUE, smartTableDB=smartTableDB)
        }
        threshold <- at$mean + autoThresholdLevel*at$sd
        if(verbosity>=2) {
          # cat(paste(" [runId:",runId,"] ITERATION:",iiter,"         Automatic threshold for class",i," spatial clustering selected at ",round(threshold,1),"above the mean density.\n"))
          # cat(paste(" [runId:",runId,"] ITERATION:",iiter,"                                 class",i," got dif_max_mean             = ",round(dif_max_mean,1),"above the mean density.\n"))
        }
      } else {
        stop(" The code without autothresolding was deprecated. \b Aborting cowardly.\n\n")
      }
      
      if(is.na(round(threshold,1))) {
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"/- PROBLEM REPORT -----------------------------------------------------------------\n"))
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"|      Class",i," has a NA value in the threshold!\n"))
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"|      probably due to its small number of stars: ", nrow(dfn),".\n"))
        kde2dForSubset(ocdata_full_withColorsAndPcaCols, setw=i, returnDistance=FALSE, showStats=TRUE, printPlots=FALSE, positionDataIndexes=positionDataIndexes)
        print(dfn)
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"\\----------------------------------------------------------------------------------\n"))
      } 
      
      if(!is.na(round(threshold,1))) {
        if(round(dif_max_mean,1) >= round(threshold,1)) {
          vclass <- c(vclass, i)
          if(verbosity>=2) {
            cat(paste(" [runId:",runId,"] ITERATION:",iiter," Class",i," ok: ",round(dif_max_mean,2),">=",round(threshold,2),"\n"))
          }
        } else {
          not_class <- c(not_class, i)
          # if(verbosity>=2) {
          #   cat(paste(" [runId:",runId,"] ITERATION:",iiter,"      --   Class",i," : WILL BE ELIMINATED!!\n"))
          # }
        }
      }
    } else {
      # if(verbosity>=2) {
      #   cat(paste(" [runId:",runId,"] ITERATION:",iiter,"   -- Class",i," -- THERE ARE TWO OR LESS STARS IN THIS CLASS! \n"))
      # }
      not_class <- c(not_class, i)
    }
  }
  
  if(length(vclass)==0) {
    cat(" No spatial clustering detected in the real space based on the clustered photometric data in the transformed data space!\n")
    #cat(paste(" Spatial thresholding at ",round(threshold,2), "\n")) ### ---- USED FOR DEPURATION PURPOSES
    if(stopIfEmpty) {
      stop(" No spatial clustering detected!\n Aborting cowardly!")
    } else {
      oc_reconst <- ocdata_full_withColorsAndPcaCols[0,]
    }
  }
  
  if(verbosity>=1) {
    if(length(vclass)!=0) {
      if(verbosity >= 2) {
      #    cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- The ",length(vclass),"selected classes from the kernel density estimation in the X-Y space are: "))
      #    print(vclass)
      # } else {
         cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- Number of classes selected from the kde analysis : ",length(vclass)))
      }
      cat("\n")
    } else {
      cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- No classes were selected from the kernel density estimation in the X-Y space.\n"))
    }
  }
  
  # Organize the data for returning to the outer loop
  if(length(vclass)>=1) {
    # First get the data of the selected objects
    for(i in 1:length(vclass)) {
      oc_tmp <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==vclass[i])
      if(i==1) {
        oc_reconst <- oc_tmp
      } else {
        oc_reconst <- rbind(oc_reconst, oc_tmp)
      }
    }
    # now the data of the not selected objects
    for(i in 1:length(not_class)) {
      not_tmp <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==not_class[i])
      if(i==1) {
        field_reconst <- not_tmp
      } else {
        field_reconst <- rbind(field_reconst, not_tmp)
      }         
    }
  }
  
  # if(verbosity!=0) {
  #   cat(paste(" [runId:",runId,"] ITERATION:",iiter," DONE!\n"))
  #   cat(paste("-------------------------------------------------------------------\n"))
  # }
  
  # What is going out of this function must be only the astronomical cluster stars only...
  return(oc_reconst[,1:inSize])    
}



kde2dForSubset <- function(df, setw=1, n=50, showStats=TRUE, printPlots=TRUE, returnDistance=FALSE, positionDataIndexes=c(1,2)) {
  # Load the 2d-kde library
  # library(MASS)
  
  # Filter the data
  dfn <- subset(df, df$resMclust.class==setw)
  
  # Create the 2d-kde
  dataX <- dfn[,positionDataIndexes[1]]
  dataY <- dfn[,positionDataIndexes[2]]
  # using the method of Sheather & Jones (1991) to select the bandwidth
  # kde2dmap <- kde2d(dataX, dataY, n=n, lims=c(range(df$x), range(df$y)), h = c(width.SJ(dataX), width.SJ(dataY)))
  # using normal bandwidth selection rule of thumb (MASS)
  kde2dmap <- kde2d(dataX, dataY, n=n, lims=c(range(df[,positionDataIndexes[1]]), range(df[,positionDataIndexes[2]])) ) 
  
  # Plot the results
  if(printPlots) {
    image(kde2dmap, main=paste("Class",setw) )
    points(dataX, dataY, pch=19, cex=0.3)
  }
  
  # Print some statistics
  if(showStats) {
    cat(paste("------ Stats - Class",setw," -----\n"))
    cat(paste("   Max dens.   : ", max(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Min dens.   : ", min(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Mean dens.  : ", mean(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Sd dens.    : ", sd(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Median dens.: ", median(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   MAD dens.   : ", mad(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Dist max from the mean (in sd): ", ((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))),"\n" ))
    cat(paste("-------------------------------------\n\n"))     
  }
  
  # Return the distance
  if(returnDistance) {
    return(((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))))
  }
}


analyse_randomKde2d_smart <- function(nfields=100, nstars, maxX, maxY, nKde=50, 
                                      showStats=FALSE, returnStats=TRUE, smartTableDB) {
  
  # Get the table from the database
  smartTable <- dbReadTable(smartTableDB, "smartTable")
  
  res <- subset(smartTable, smartTable$nstar==nstars) # only nstars matters by now...
  
  if(length(res$mean)!=0) {
    retStat <- data.frame(mean=mean(res$mean), sd=mean(res$sd))
  } else {
    # run the analysis...
    retStat <- analyse_randomKde2d(nfields, nstars, maxX, maxY, nKde, showStats, returnStats)
    # and store the results in the table, so next time you won't need to compute it again!
    dbWriteTable(smartTableDB, "smartTable", data.frame(nstar=nstars, mean=retStat$mean, sd=retStat$sd), append=TRUE, row.names = F)
  }
  
  return(retStat)
}

analyse_randomKde2d <- function(nfields=100, nstars, maxX, maxY, nKde=50, showStats=FALSE, returnStats=TRUE) {
  
  maxDistStats <- vector("double", nfields)
  
  # run the analysis
  for(i in 1:nfields) {
    maxDistStats[i] <- create_randomKde2d(nstars, maxX, maxY, nKde=nKde, returnDistance=TRUE)
  }

  # Print some statistics
  if(showStats) {
    cat(paste("------ Statistics of the sample random fields -----\n"))
    cat(paste("   Max distance    : ", max(maxDistStats),"\n"))
    cat(paste("   Min distance    : ", min(maxDistStats),"\n"))
    cat(paste("   Mean distance   : ", mean(maxDistStats),"\n"))
    cat(paste("   Sd distance     : ", sd(maxDistStats),"\n"))
    cat(paste("   Median distance : ", median(maxDistStats),"\n"))
    cat(paste("   MAD distance    : ", mad(maxDistStats),"\n"))
    cat(paste("---------------------------------------------------\n\n"))       
    hist(maxDistStats, freq=FALSE)
    lines(density(maxDistStats), col="red")
  }
  
  # Return the mean and st.dev of the distribution  
  if(returnStats) {
    return(data.frame(mean=mean(maxDistStats), sd=sd(maxDistStats)))
  }
  
}


create_randomKde2d <- function(nstars, maxX, maxY, nKde=50, printPlots=FALSE, showStats=FALSE, returnDistance=FALSE) {
  # Load the 2d-kde library
  #library(MASS)
  
  # Create a random sample in the X-Y space
  dataX <- runif(nstars, 0, maxX)
  dataY <- runif(nstars, 0, maxY)
  
  # Create the 2d-kde
  # using the method of Sheather & Jones (1991) to select the bandwidth
  # kde2dmap <- kde2d(dataX, dataY, n=nKde, lims=c(range(df$x), range(df$y)), h = c(width.SJ(dataX), width.SJ(dataY)))
  # using normal bandwidth selection rule of thumb (MASS)
  kde2dmap <- kde2d(dataX, dataY, n=nKde, lims=c(0, maxX, 0, maxY) ) 
  
  # Plot the results
  if(printPlots) {
    image(kde2dmap, main=paste("Random Field") )
    points(dataX, dataY, pch=19, cex=0.3)
  }
  
  # Print some statistics
  if(showStats) {
    cat(paste("------ Statistics of the random field -----\n"))
    cat(paste("   Max dens.   : ", max(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Min dens.   : ", min(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Mean dens.  : ", mean(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Sd dens.    : ", sd(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Median dens.: ", median(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   MAD dens.   : ", mad(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Dist max from the mean (in sd): ", ((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))),"\n" ))
    cat(paste("-------------------------------------------\n\n"))       
  }
  
  # Return the distance
  if(returnDistance) {
    return(((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))))
  }
}






















filenameWithPathInput <- "input/0.10_0.10_0.10_0.10_0.10.dat"
sep <- ","
positionDataIndexes <- c(2,3)
photometricDataIndexes <- c(8,22)
photometricErrorDataIndexes <- c(9,23)

# filenameWithPathInput <- "input/0.80_0.49_0.49_0.33_0.29.dat"
# sep <- " "
# positionDataIndexes <- c(2,3)
# photometricDataIndexes <- c(6,7)
# photometricErrorDataIndexes <- c(4,5)

filenameWithPathOuput <- file.path(getwd(), "up-RESULTS.dat")
fileWithHeader <- TRUE

nDimsToKeep <- 2
nRuns <- 50
verbose <- TRUE

# runInParallel <- FALSE, paralelization <- "multicore", independent <- TRUE, 
# threshold <- 1
# starsPerClust_kmeans <- 50
# nstarts_kmeans <- 50
# autoCalibrated <- FALSE, considerErrors <- FALSE
# dimRed="PCA", scale=TRUE
finalXYCut <- FALSE


UPMASKfile(filenameWithPathInput, filenameWithPathOuput,
positionDataIndexes, photometricDataIndexes, photometricErrorDataIndexes,
nRuns=nRuns, verbose=verbose, fileWithHeader=fileWithHeader, nDimsToKeep=nDimsToKeep,
sep=sep, finalXYCut=finalXYCut)

# Sol:
# > posIdx <- c(2,3)
# > photIdx <- c(8,20,22,24)
# > photErrIdx <- c(9,21,23,25)
# > UPMASKfile(inputFileName, outputFileName, posIdx, photIdx, photErrIdx, nRuns=50, 
# +            starsPerClust_kmeans=50, verbose=TRUE, fileWithHeader=TRUE)