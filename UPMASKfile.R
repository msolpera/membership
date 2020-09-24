#! /home/gabriel/miniconda3/envs/py3/lib/R/bin/Rscript

#  R package UPMASK file R/UPMASKfile.R
#  Copyright (C) 2014 Alberto Krone-Martins
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#' @title Run UPMASK in a file
#' 
#' @description \code{UPMASKfile} executes the UPMASK method using a file as an input
#' and writes another file as an output. This is a wrapper function that only reads a 
#' file into an R data frame, calls the \code{UPMASKdata} function using this data frame 
#' and the parameters passed by the user and writes the output into another file.
#' 
#' @param filenameWithPathInput a string indicating the file containing the data to run UPMASK on (with full path)
#' @param filenameWithPathOuput a string indicating the file where the output shall be written (with full path) 
#' @param positionDataIndexes an array of integers indicating the columns of the file containing the spatial position measurements
#' @param photometricDataIndexes an array of integers with the column numbers containing photometric measurements (or any other measurement to go into the PCA step)
#' @param photometricErrorDataIndexes an array of integers with the column numbers containing the errors of the photometric measurements
#' @param threshold a double indicating the thresholding level for the random field analysis
#' @param maxIter an integer the maximum amount of iterations of the outer loop before giving up convergence (usually it is not necessary to modify this)
#' @param starsPerClust_kmeans an integer with the average number of stars per k-means cluster
#' @param nstarts_kmeans an integer the amount of random re-initializations of the k-means clustering method (usually it is not necessary to modify this)
#' @param nRuns the total number of individual runs to execute the total number of outer loop runs to execute
#' @param runInParallel a boolean indicating if the code should run in parallel
#' @param paralelization a string with the type of paralilization to use. the paralelization can be: "multicore" or "MPIcluster". At this moment only "multicore" is implemented (defaults to multicore).
#' @param independent a boolean indicating if non-parallel runs should be completely independent
#' @param verbose a boolean indicating if the output to screen should be verbose
#' @param autoCalibrated a boolean indicating if the number of random field realizations for the clustering check in the position space should be autocalibrated (experimental code, defaults to FALSE).
#' @param considerErrors a boolean indicating if the errors should be taken into account
#' @param finalXYCut a boolean indicating if a final cut in the XY space should be performed (defaults to FALSE)
#' @param fileWithHeader a boolean indicating if the input file has a text header
#' @param nDimsToKeep an integer with the number of dimensions to consider (defaults to 4)
#' @param dimRed a string with the dimensionality reduction method to use (defaults to PCA. The only other options are LaplacianEigenmaps or None)
#' @param scale a boolean indicating if the data should be scaled and centered
#' 
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#'
#' @examples
#' \dontrun{
#' # Analyse a simulated open cluster using spatial and photometric data 
#' # Create strings with filenames
#' fileNameI <- "oc_12_500_1000_1.0_p019_0880_1_25km_120nR_withcolors.dat"
#' inputFileName <- system.file("extdata", fileNameI, package="UPMASK")
#' outputFileName <- file.path(tempdir(), "up-RESULTS.dat")
#' 
#' # Example of how to run UPMASK using data from a file
#' # (serious analysis require at least larger nRuns)
#' posIdx <- c(1,2)
#' photIdx <- c(3,5,7,9,11,19,21,23,25,27)
#' photErrIdx <- c(4,6,8,10,12,20,22,24,26,28)
#' UPMASKfile(inputFileName, outputFileName, posIdx, photIdx, photErrIdx, nRuns=5, 
#'            starsPerClust_kmeans=25, verbose=TRUE, fileWithHeader=TRUE)
#' 
#' # Open the resulting file to inspect the results
#' tempResults <- read.table(outputFileName, header=TRUE)
#' 
#' # Create a simple raw plot to see the results
#' pCols <- tempResults[,length(tempResults)]/max(tempResults[,length(tempResults)])
#' plot(tempResults[,1], tempResults[,2], col=rgb(0,0,0,pCols), cex=0.5, pch=19)
#' 
#' # Clean the environment
#' rm(list=c("tempResults", "inputFileName", "outputFileName", "pCols", "fileNameI"))
#' } 
#'  
#' @usage UPMASKfile(filenameWithPathInput, filenameWithPathOuput, 
#' positionDataIndexes=c(1,2), photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
#' photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28), threshold=1, 
#' maxIter=20, starsPerClust_kmeans=50, nstarts_kmeans=50, nRuns=5, 
#' runInParallel=FALSE, paralelization="multicore", independent=TRUE, verbose=FALSE, 
#' autoCalibrated=FALSE, considerErrors=FALSE, finalXYCut=FALSE, 
#' fileWithHeader=FALSE, nDimsToKeep=4, dimRed="PCA", scale=TRUE)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords misc, utilities
#' @export
#
UPMASKfile <- function(filenameWithPathInput, filenameWithPathOuput, 
					  positionDataIndexes=c(1,2),
					  photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
					  photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28),
					  threshold=1, maxIter=20, starsPerClust_kmeans=50, nstarts_kmeans=50, 
            nRuns=5, runInParallel=FALSE, paralelization="multicore", 
            independent=TRUE, verbose=FALSE, autoCalibrated=FALSE, 
            considerErrors=FALSE, finalXYCut=FALSE, fileWithHeader=FALSE, 
					  nDimsToKeep=4, dimRed="PCA", scale=TRUE) {
  
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
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Reading the input table from: \n\t",filenameWithPathInput,"\n")
  }
  # Load the file
  ocdata_full <- read.table(filenameWithPathInput, header=fileWithHeader)
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
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
  write.table(resultsTable, filenameWithPathOuput, sep="	", 
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


#  R package UPMASK file R/UPMASKdata.R
#  Copyright (C) 2014 Alberto Krone-Martins
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#' @title Run UPMASK in a data frame
#' 
#' @description \code{UPMASKdata} executes the UPMASK method on a data frame, and returns 
#' another data frame as output, including the membership analysis result as additional 
#' columns. 
#' 
#' \code{UPMASKdata} is a method for performing membership assignment in stellar 
#' clusters. The distributed code is prepared to use photometry and spatial positions, 
#' but it can take into account other types of data as well. The method is able to take 
#' into account arbitrary error models (the used must rewrite the 
#' \code{\link{takeErrorsIntoAccount}} function), and it is unsupervised, data-driven, 
#' physical-model-free and relies on as few assumptions as possible. The approach followed
#' for membership assessment is based on an iterative process, dimensionality reduction, 
#' a clustering algorithm and a kernel density estimation.
#' 
#' @param dataTable a data frame with the data to perform the analysis
#' @param positionDataIndexes an array of integers indicating the columns of the data frame containing the spatial position measurements
#' @param photometricDataIndexes an array of integers with the column numbers containing photometric measurements (or any other measurement to go into the PCA step)
#' @param photometricErrorDataIndexes an array of integers with the column numbers containing the errors of the photometric measurements
#' @param threshold a double indicating the thresholding level for the random field analysis
#' @param classAlgol a string indicating the type of clustering algorithm to consider. Only k-means is implemented at this moment (defaults to kmeans)
#' @param maxIter an integer the maximum amount of iterations of the outer loop before giving up convergence (usually it is not necessary to modify this)
#' @param starsPerClust_kmeans an integer with the average number of stars per k-means cluster
#' @param nstarts_kmeans an integer the amount of random re-initializations of the k-means clustering method (usually it is not necessary to modify this)
#' @param nRuns the total number of individual runs to execute the total number of outer loop runs to execute
#' @param runInParallel a boolean indicating if the code should run in parallel
#' @param paralelization a string with the type of paralilization to use. the paralelization can be: "multicore" or "MPIcluster". At this moment only "multicore" is implemented (defaults to multicore).
#' @param independent a boolean indicating if non-parallel runs should be completely independent
#' @param verbose a boolean indicating if the output to screen should be verbose
#' @param autoCalibrated a boolean indicating if the number of random field realizations for the clustering check in the position space should be autocalibrated (experimental code, defaults to FALSE).
#' @param considerErrors a boolean indicating if the errors should be taken into account
#' @param finalXYCut a boolean indicating if a final cut in the XY space should be performed (defaults to FALSE)
#' @param nDimsToKeep an integer with the number of dimensions to consider (defaults to 4)
#' @param dimRed a string with the dimensionality reduction method to use (defaults to PCA. The only other options are LaplacianEigenmaps or None)
#' @param scale a boolean indicating if the data should be scaled and centered
#' 
#' @return A data frame with the original data used to run the method and additional columns indicating the classification at each run, as well as a membership probability in the frequentist sense.
#' 
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#'
#' @examples
#' \dontrun{
#' # Analyse a simulated open cluster using spatial and photometric data 
#' # Load the data into a data frame
#' fileNameI <- "oc_12_500_1000_1.0_p019_0880_1_25km_120nR_withcolors.dat"
#' inputFileName <- system.file("extdata", fileNameI, package="UPMASK")
#' ocData <- read.table(inputFileName, header=TRUE)
#' 
#' # Example of how to run UPMASK using data from a data frame
#' # (serious analysis require at least larger nRuns)
#' posIdx <- c(1,2)
#' photIdx <- c(3,5,7,9,11,19,21,23,25,27)
#' photErrIdx <- c(4,6,8,10,12,20,22,24,26,28)
#' 
#' upmaskRes <- UPMASKdata(ocData, posIdx, photIdx, PhotErrIdx, nRuns=2, 
#'                         starsPerClust_kmeans=25, verbose=TRUE)
#' 
#' # Create a simple raw plot to see the results
#' pCols <- upmaskRes[,length(upmaskRes)]/max(upmaskRes[,length(upmaskRes)])
#' plot(upmaskRes[,1], upmaskRes[,2], col=rgb(0,0,0,pCols), cex=0.5, pch=19)
#' 
#' # Clean the environment
#' rm(list=c("inputFileName", "ocData", "posIdx", "photIdx", "photErrIdx", 
#'           "upmaskRes", "pCols"))
#' } 
#'  
#' @usage UPMASKdata(dataTable, positionDataIndexes=c(1,2),
#' photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
#' photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28), threshold=1, 
#' classAlgol="kmeans", maxIter=25, starsPerClust_kmeans=25, nstarts_kmeans=50, 
#' nRuns=8, runInParallel=FALSE, paralelization="multicore", independent=TRUE, 
#' verbose=FALSE, autoCalibrated=FALSE, considerErrors=FALSE, 
#' finalXYCut=FALSE, nDimsToKeep=4, dimRed="PCA", scale=TRUE)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords cluster, methods, multivariate, nonparametric
#' @import parallel
#' @export
#
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
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Starting UPMASK analysis!\n")
    cat(      " UPMASK kernels v. 1.2\n")
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
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
    # plot(resMclust, data=ocdata_px[,1:4])
    pairs(ocdata_px[,1:4], pch=19, cex=0.2, col=rainbow(max(ocdata_px$resMclust.class))[ocdata_px$resMclust.class])
  }
  
  # Merge the colors, the dimensionality reduction columns and cluster classification with the original results
  ocdata_full_withColorsAndPcaCols <- data.frame(ocdata_full, ocdata, ocdata_px)
    
  # Select the classes with densities above the threshold
  if(verbosity>=1) {
    cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [3/3] Performing comparison with random density fields for",max(ocdata_px$resMclust.class),"individual classes...\n"))
    if(verbosity==1) {
      cat(paste(" [runId:",runId,"] ITERATION:",iiter,"          You can get a coffee. This can take long! \n"))
    }
  }
  
  vclass <- c()
  not_class <- c()
  for(i in 1:max(ocdata_px$resMclust.class)) {
    dfn <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==i)
        
    if(nrow(dfn)>2) {
      
      dif_max_mean <- kde2dForSubset(ocdata_full_withColorsAndPcaCols, setw=i, returnDistance=TRUE, showStats=FALSE, printPlots=FALSE, positionDataIndexes=positionDataIndexes)
      
      # First, get the thresholding level...
      if(autoThreshold) {
        if(verbosity>=2) {
          cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- Class",i," -- Performing analysis of random fields...\n"))
        }
        
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
          cat(paste(" [runId:",runId,"] ITERATION:",iiter,"     Automatic threshold for class",i," spatial clustering selected at ",round(threshold,1),"above the mean density.\n"))
          cat(paste(" [runId:",runId,"] ITERATION:",iiter,"                             class",i," got dif_max_mean             = ",round(dif_max_mean,1),"above the mean density.\n"))
        }
      } else {
        stop(" The code without autothresolding was deprecated. \b Aborting cowardly.\n\n")
      }
      
      if(is.na(round(threshold,1))) {
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"/- PROBLEM REPORT -----------------------------------------------------------------\n"))
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"|    Class",i," has a NA value in the threshold!\n"))
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"|    probably due to its small number of stars: ", nrow(dfn),".\n"))
        kde2dForSubset(ocdata_full_withColorsAndPcaCols, setw=i, returnDistance=FALSE, showStats=TRUE, printPlots=FALSE, positionDataIndexes=positionDataIndexes)
        print(dfn)
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"\\----------------------------------------------------------------------------------\n"))
      } 
      
      if(!is.na(round(threshold,1))) {
        if(round(dif_max_mean,1) >= round(threshold,1)) {
          vclass <- c(vclass, i)
          if(verbosity>=2) {
            cat(paste(" [runId:",runId,"] ITERATION:",iiter," <<<< -- Class",i," : ok!\n"))
          }
        } else {
          not_class <- c(not_class, i)
          if(verbosity>=2) {
            cat(paste(" [runId:",runId,"] ITERATION:",iiter,"      -- Class",i," : WILL BE ELIMINATED!!\n"))
          }
        }
      }
    } else {
      if(verbosity>=2) {
        cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- Class",i," -- THERE ARE TWO OR LESS STARS IN THIS CLASS! \n"))
      }
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
         cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- The ",length(vclass),"selected classes from the kernel density estimation in the X-Y space are: "))
         print(vclass)
      } else {
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
  
  if(verbosity!=0) {
    cat(paste(" [runId:",runId,"] ITERATION:",iiter," DONE!\n"))
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
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

















































#####################################################################
cat(paste("R version: ",R.version.string))

set.seed(12345)

library('DBI')
library(RSQLite)
library('MASS')

fileNameI <- "input/189.00_0.04_1.48_17.27_7.03.dat"
# inputFileName <- system.file("extdata", fileNameI, package="UPMASK")
outputFileName <- file.path(getwd(), "output/up-RES.dat")
positionDataIndexes <- c(2,3)
photometricDataIndexes <- c(6,7)
photometricErrorDataIndexes <- c(4,5)
nDimsToKeep <- 2
nRuns <- 5
starsPerClust_kmeans <- 25
verbose <- TRUE
# verbose <- 0
finalXYCut <- FALSE
fileWithHeader <- TRUE
UPMASKfile(fileNameI, outputFileName,
positionDataIndexes, photometricDataIndexes, photometricErrorDataIndexes,
nRuns=nRuns, verbose=verbose, fileWithHeader=fileWithHeader, nDimsToKeep=nDimsToKeep,
finalXYCut=finalXYCut, starsPerClust_kmeans=starsPerClust_kmeans)