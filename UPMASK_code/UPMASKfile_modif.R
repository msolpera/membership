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
  
  # if(verbose) {
  #   cat(paste("-------------------------------------------------------------------\n"))
  #   cat(      " Starting UPMASK...\n")
  #   cat(paste("-------------------------------------------------------------------\n"))
  # }
  

  if(verbose) {
    # cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Reading the input table from: ",filenameWithPathInput,"\n")
  }
  # Load the file
  ocdata_full <- read.table(filenameWithPathInput, header=fileWithHeader, sep=sep)

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


  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Writing the output table at: \n\t",filenameWithPathOuput,"\n")
  }
  write.table(resultsTable, filenameWithPathOuput, sep="    ", 
              col.names=fileWithHeader, quote=FALSE, row.names=FALSE)   

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
  
  # Create the smart look-up table for the random analysis
  # Unfortunatelly this must be a global variable so it can be shared among 
  # parallel processes
  stcon <- create_smartTable() 
    
  # Create a place-holder list
  pp <- list()
  
  # Run UPMASK outerloop (and innerloop, called inside the outerloop)
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
   
  # Close the connection to the database storing the smart lookup table
  dbDisconnect(stcon)

  # That's all folks!
  return(ocdata_out)
}


create_smartTable <- function() {
#  con <- dbConnect(RSQLite::SQLite(), ":memory:")
  con <- dbConnect(dbDriver("SQLite"), dbname = tempfile())
  dbWriteTable(con, "smartTable", data.frame(nstar=c(0), mean=c(0), sd=c(0)), row.names = F)
  return(con)
}



outerLoop <- function(ocdata_full, 
                      positionDataIndexes=c(1,2),
                      photometricDataIndexes=c(3),
                      photometricErrorDataIndexes=c(4),
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
  
  ocdata_out <- data.frame(ocdata_out, finalClass=ocdata_out[,length(ocdata_out)] )   
    
  # That's all folks!  
  return(data.frame(id=ocdata_out$id, class=ocdata_out$finalClass))
}



innerLoop <- function(ocdata_full, ocdata, classAlgol="kmeans", autoThresholdLevel=3, 
                         autoThreshold=TRUE, iiter=0, plotIter=FALSE, verbosity=1, 
                         starsPerClust_kmeans=50, nstarts_kmeans=50, runId=0, 
                         autoCalibrated=FALSE, stopIfEmpty=FALSE, 
                         positionDataIndexes=c(1,2), smartTableDB, nDimsToKeep=4,
                         dimRed="PCA", scale=TRUE) {

  cat(paste(" [runId:",runId,"] ITERATION:",iiter," RUNNING...\n"))
  inSize <- length(ocdata_full)
  # cat(paste(" ocdata_full N= ",length(ocdata_full[,1]),"\n"))
  # cat(paste(" ocdata      N= ",length(ocdata[,1]),"\n"))
  
  # Perform the dimensionality reduction step
  cat(paste(" [runId:",runId,"] ITERATION:",iiter,", Using PCA: ",nDimsToKeep," dimensions","\n"))
  ocdata_pca <- prcomp(ocdata, scale=scale, center=scale)
  ocdata_px <- predict(ocdata_pca)
  
  # Perform the clustering step
  starsPerClust <- starsPerClust_kmeans
  nclust <- round(length(ocdata_full[,1])/starsPerClust)
  # if(nclust > 1) {
  fit <- kmeans(ocdata_px[,1:nDimsToKeep], nclust, nstart=nstarts_kmeans, iter.max=100)
  # get cluster means WHAT DOES THIS LINE DO?
  aggregate(ocdata_px, by=list(fit$cluster), FUN=mean)
  # append cluster assignment
  ocdata_px <- data.frame(ocdata_px, resMclust.class=fit$cluster)
  
  # Merge the colors, the dimensionality reduction columns and cluster classification with the original results
  ocdata_full_withColorsAndPcaCols <- data.frame(ocdata_full, ocdata, ocdata_px)
    
  # Select the classes with densities above the threshold
  cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [3/3] Found",max(ocdata_px$resMclust.class),"individual classes...\n"))
  
  vclass <- c()
  not_class <- c()
  for(i in 1:max(ocdata_px$resMclust.class)) {
    dfn <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==i)
        
    if(nrow(dfn)>2) {
      
      dif_max_mean <- kde2dForSubset(ocdata_full_withColorsAndPcaCols, setw=i, returnDistance=TRUE, showStats=FALSE, printPlots=FALSE, positionDataIndexes=positionDataIndexes)
      # Analyze random fields
      at <- analyse_randomKde2d_smart(
             nfields=2000, nstars=length(dfn$resMclust.class), 
             (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])),
             (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])), 
             nKde=50, showStats=FALSE, returnStats=TRUE, smartTableDB=smartTableDB)
      threshold <- at$mean + autoThresholdLevel*at$sd
      
      if(!is.na(round(threshold,1))) {
        if(round(dif_max_mean,1) >= round(threshold,1)) {
          vclass <- c(vclass, i)
          # cat(paste(" Stats:",dif_max_mean, at$mean, at$sd))
          cat(paste(" [runId:",runId,"] ITERATION:",iiter," Class",i," ok: ",round(dif_max_mean,1),">=",round(threshold,1),"N=",nrow(dfn),"\n"))
          # cat(paste(" [runId:",runId,"] ITERATION:",iiter," Class",i,round((dif_max_mean-at$mean)/at$sd,2),"N=",nrow(dfn),"\n"))
        } else {
          not_class <- c(not_class, i)
        }
      }
    } else {
      #   cat(paste(" [runId:",runId,"] ITERATION:",iiter,"   -- Class",i," -- THERE ARE TWO OR LESS STARS IN THIS CLASS! \n"))
      not_class <- c(not_class, i)
    }
  }
  
  # Organize the data for returning to the outer loop
  # if(length(vclass)>=1) {
  # First get the data of the selected objects
  for(i in 1:length(vclass)) {
    oc_tmp <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==vclass[i])
    if(i==1) {
      oc_reconst <- oc_tmp
    } else {
      oc_reconst <- rbind(oc_reconst, oc_tmp)
    }
  }
  cat(paste("Number of classes selected from the kde analysis : ",length(vclass),"N=",nrow(oc_reconst),"\n"))

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
  # xyminmax <- list(c(range(df[,positionDataIndexes[1]])), c(range(df[,positionDataIndexes[2]])))
  # cat(paste("Cluster KDE range",xyminmax[[1]][1], xyminmax[[1]][2], xyminmax[[1]][1], xyminmax[[2]][2],"\n"))
  kde2dmap <- kde2d(dataX, dataY, n=n, lims=c(range(df[,positionDataIndexes[1]]), range(df[,positionDataIndexes[2]])) ) 
  
  # # Print some statistics
  # if(showStats) {
  #   cat(paste("------ Stats - Class",setw," -----\n"))
  #   cat(paste("   Max dens.   : ", max(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Min dens.   : ", min(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Mean dens.  : ", mean(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Sd dens.    : ", sd(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Median dens.: ", median(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   MAD dens.   : ", mad(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Dist max from the mean (in sd): ", ((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))),"\n" ))
  #   cat(paste("-------------------------------------\n\n"))     
  # }
  
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
  
  # cat(paste("Random KDE range",0, maxX, 0, maxY),"\n")
  # run the analysis
  for(i in 1:nfields) {
    maxDistStats[i] <- create_randomKde2d(nstars, maxX, maxY, nKde=nKde, returnDistance=TRUE)
  }

  # # Print some statistics
  # if(showStats) {
  #   cat(paste("------ Statistics of the sample random fields -----\n"))
  #   cat(paste("   Max distance    : ", max(maxDistStats),"\n"))
  #   cat(paste("   Min distance    : ", min(maxDistStats),"\n"))
  #   cat(paste("   Mean distance   : ", mean(maxDistStats),"\n"))
  #   cat(paste("   Sd distance     : ", sd(maxDistStats),"\n"))
  #   cat(paste("   Median distance : ", median(maxDistStats),"\n"))
  #   cat(paste("   MAD distance    : ", mad(maxDistStats),"\n"))
  #   cat(paste("---------------------------------------------------\n\n"))       
  #   hist(maxDistStats, freq=FALSE)
  #   lines(density(maxDistStats), col="red")
  # }
  
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
  
 
  # # Print some statistics
  # if(showStats) {
  #   cat(paste("------ Statistics of the random field -----\n"))
  #   cat(paste("   Max dens.   : ", max(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Min dens.   : ", min(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Mean dens.  : ", mean(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Sd dens.    : ", sd(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Median dens.: ", median(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   MAD dens.   : ", mad(as.vector(kde2dmap$z)),"\n"))
  #   cat(paste("   Dist max from the mean (in sd): ", ((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))),"\n" ))
  #   cat(paste("-------------------------------------------\n\n"))       
  # }
  
  # Return the distance
  if(returnDistance) {
    return(((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))))
  }
}


















set.seed(12345)

library('DBI')
library(RSQLite)
library('MASS')

positionDataIndexes <- c(2,3)
nRuns <- 25
starsPerClust_kmeans <- 25
verbose <- TRUE
finalXYCut <- FALSE
filenameWithPathOuput <- file.path(getwd(), "output/up-RESULTS.dat")
fileWithHeader <- TRUE

# filenameWithPathInput <- "input/0.10_0.10_0.10_0.10_0.10.dat"
# sep <- ","
# positionDataIndexes <- c(2,3)
# photometricDataIndexes <- c(8,22)
# photometricErrorDataIndexes <- c(9,23)

sep <- " "

# filenameWithPathInput <- "input/oc_12_500_3000_3.1_p019_0900_1.dat"
# photometricDataIndexes <- c(4,5,6,7,8,9)
# photometricErrorDataIndexes <- c(4,5,6,7,8,9)
# nDimsToKeep <- 4
filenameWithPathInput <- "input/78.00_1.75_5.11_3.16_49.54.dat"
photometricDataIndexes <- c(6,7)
photometricErrorDataIndexes <- c(4,5)
nDimsToKeep <- 2

# filenameWithPathInput <- "input/oc_12_500_1500_1.5_p019_0800_1.dat"
# sep <- ""
# positionDataIndexes <- c(1,2)
# photometricDataIndexes <- c(3,5,7,9,11,13,15,17)
# photometricErrorDataIndexes <- c(4,6,8,10,12,14,16,18)

# runInParallel <- FALSE, paralelization <- "multicore", independent <- TRUE, 
# threshold <- 1
# nstarts_kmeans <- 50
# autoCalibrated <- FALSE, considerErrors <- FALSE
# dimRed="PCA", scale=TRUE

UPMASKfile(filenameWithPathInput, filenameWithPathOuput,
positionDataIndexes, photometricDataIndexes, photometricErrorDataIndexes,
nRuns=nRuns, verbose=verbose, fileWithHeader=fileWithHeader, nDimsToKeep=nDimsToKeep,
sep=sep, finalXYCut=finalXYCut, starsPerClust_kmeans=starsPerClust_kmeans)

# Sol:
# > posIdx <- c(2,3)
# > photIdx <- c(8,20,22,24)
# > photErrIdx <- c(9,21,23,25)
# > UPMASKfile(inputFileName, outputFileName, posIdx, photIdx, photErrIdx, nRuns=50, 
# +            starsPerClust_kmeans=50, verbose=TRUE, fileWithHeader=TRUE)