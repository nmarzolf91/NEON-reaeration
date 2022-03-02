##############################################################################################
#' @title Constant rate or slug travel time between two stations

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function determines the timestamp of the peak tracer conductance for use
#' calculating travel time between an upstream and downstream sensor set.

#' @importFrom graphics points
#' @importFrom pracma trapz
#' @importFrom stats loess.smooth

#' @param loggerDataIn User input of the R data object holding the conductivity time series
#' for a site, date, and station [dataframe]
#' @param currEventID User input of the eventID of the tracer experiment [string]
#' @param injectionType User input of the injection type either "constant" or "slug" [string]
#' @param expStartTime User input of the experiment start time, in UTC timezone [posixct]

#' @return This function returns the peak tracer timestamp [dateTime]

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, velocity, travel time, reaeration, metabolism

#' @examples
#' #Using an example file
#' #travelTime <- def.calc.travelTime(
#' #dataDir = paste(path.package("reaRate"),"inst\\extdata", sep = "\\"),
#' #currEventID = "GUIL.20150129", injectionType = "constant", bPlot = T)

#' @seealso def.calc.reaeration.R for calculating reaeration rates

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-08-03)
#     original creation
#   Kaelin M. Cawley (2021-05-25)
#     updated to handle model injection types
##############################################################################################
def.calc.peakTime <- function(
  loggerDataIn,
  currEventID,
  injectionType,
  expStartTime,
  station
){
  library(zoo)
  library(lubridate)
  station<-station
  #Trim the data for only after the experiment started
  #need a QC for duplicate measurements per time? e.g. length(loggerDataIn$dateTimeLogger) = length(unique(loggerDataIn$dateTimeLogger))
  trimTime <- ifelse(min(loggerDataIn$dateTimeLogger) < (expStartTime - 5*60), (expStartTime - 5*60), min(loggerDataIn$dateTimeLogger))
  #that trims it to the experiment start time to get rid of extra data
  loggerDataTrim <- loggerDataIn[loggerDataIn$dateTimeLogger > expStartTime,] #does not trim for the end, 
  #would only be able to do that if ran station 4 first
  #plot(loggerData$dateTimeLogger, loggerData$spCond)
  #lines(loggerDataTrim$dateTimeLogger, loggerDataTrim$spCond, col = "blue")
  
  #Create a plot where users select the range to pick the peak
  medCond <- stats::median(loggerDataTrim$spCond, na.rm =TRUE)
  stdCond <- stats::mad(loggerDataTrim$spCond, na.rm = TRUE)
  lowPlot <- ifelse(medCond-30 < 0, 0, medCond-30) #these are just for axis limits
  highPlot <- ifelse(medCond+30 > max(loggerDataTrim$spCond, na.rm = TRUE), 
                     medCond+30,
                     max(loggerDataTrim$spCond, na.rm = TRUE))
  #so the noise seems to be real but collected every 10 seconds.  It's probably not perfectly mixed by station 1
  #maybe a 5 point moving window?
  #try to add in a column with the rolling mean #load in library zoo earlier
  loggerDataTrim$rollmean<-rollmean(loggerDataTrim$spCond, 9, fill = list(rep(NA, 4), NULL, rep(NA,4)))
  
    invisible(dev.new(noRStudioGD = TRUE, width=12, height=7))
  plot(loggerDataTrim$dateTimeLogger,
       loggerDataTrim$spCond,
       cex=0.7,
       col='grey20',
       xlab = "Measurement time",
       ylab = "Specific Conductance",
       ylim = c(lowPlot, highPlot),
       main = paste(currEventID, station, sep = '\n'))  ##here is where I needed to add a title with the station number
  
  #add in a 2 minute, 10 second rolling window mean
  lines(loggerDataTrim$dateTimeLogger,
        loggerDataTrim$rollmean, col='coral2', lwd=1.5)
  
  #Have users choose if the plot has a defined peak, make the buttons closer together
  points(x = c(max(loggerDataTrim$dateTimeLogger),max(loggerDataTrim$dateTimeLogger)),
         y = c(highPlot*.9,highPlot*.75),
         col = c("darkolivegreen4", "tomato3"),
         lwd = 2,
         pch = 19,
         cex = 6)
  #title(main = paste0("Click green dot (upper lefthand) if the peak/plateau is identifiable. \nClick red dot (upper righthand) if not identifiable.\n",currEventID))
  badPlotBox <- identify(x = c(max(loggerDataTrim$dateTimeLogger),max(loggerDataTrim$dateTimeLogger)),
                         y = c(highPlot*.9,highPlot*.7),
                         n = 1,
                         tolerance = 0.4,
                         labels = c("Good", "Bad"))
  Sys.sleep(1)
  invisible(dev.off())
  #okay so this gives the user the option to cut off analysis here
  #now add an additional step where a linear regression is calculated 
  #because of the error, I do not want to fit a 3 point linear regression, even if it is smoothed

  if(length(badPlotBox) && badPlotBox==1){  #1 is good, 2 is bad
    #If things look good, move on
    invisible(dev.new(noRStudioGD = TRUE, width=12, height=7))
    plot(loggerDataTrim$dateTimeLogger,
         loggerDataTrim$spCond,
         xlab = "Measurement Number",
         ylab = "Specific Conductance",
         ylim = c(lowPlot, highPlot))
    lines(loggerDataTrim$dateTimeLogger,
          loggerDataTrim$rollmean, col='coral2', lwd=1.5)
    title(main = paste0("The plot starts at the time of the injection.\nClick right of the peak/plateau where the useful timeseries ends.\n"))
    ans <- identify(x = loggerDataTrim$dateTimeLogger,
                    y = loggerDataTrim$spCond,
                    n = 1,
                    tolerance = 0.4)
    Sys.sleep(1)
    invisible(dev.off())
    endHere <- ans #choose an end point to look for the peak
  }else{
    return(NULL)
  } #so this just skips that current event ID if no station 1 peak
  
  #Trim the loggerData to just the area specified so the user can remove the double peak
  loggerDataTrim <- loggerDataTrim[1:endHere,]
  
  #calculate the derivative using 3 minute periods
  
  loggerDataTrim$minutes<-as.numeric(difftime(loggerDataTrim$dateTimeLogger,loggerDataTrim$dateTimeLogger[1], units = "min"))
  loggerDataTrim$slope <- NA
  
  backgroundCond <- mean(loggerDataTrim$spCond[1:5], na.rm = TRUE)
  
  #Background correct the logger data
  loggerDataTrim$corrSpCond <- loggerDataTrim$spCond - backgroundCond
  
  for(i in 3:length(loggerDataTrim$dateTimeLogger)) { #also depends on the period so maybe have the user designate above
    #add that function in later
    period<-5 #number of 10 second measurements to use
    spaces<-(period-1)/2
    sub<-loggerDataTrim[(i-spaces):(i+spaces),]
    mod1<-lm(corrSpCond~minutes, data=sub)  #makes sure the slope is in minutes units
    loggerDataTrim$slope[i]<-coef(mod1)[2]
  }
  hist(loggerDataTrim$slope, breaks=100)
  #add a rolling mean?
 
  loggerDataTrim$rollmeanslope<-rollmean(loggerDataTrim$slope, 13, fill = list(rep(NA, 6), NULL, rep(NA,6)))
  
#  invisible(dev.new(noRStudioGD = TRUE, width=12, height=7))
#  plot(loggerDataTrim$dateTimeLogger,
#       loggerDataTrim$slope, cex=0.7,pch=17,col='blue',xlab = "Measurement time",ylab = "Slope",
#       #ylim = c(lowPlot, highPlot),
#       main = paste(currEventID, station, sep = '\n'))  
  
  #add in a 2 minute, 10 second rolling window mean and some boundaries
#  lines(loggerDataTrim$dateTimeLogger,
#        loggerDataTrim$rollmeanslope, col='coral2', lwd=1.5)
  #wont have the background ones - too large a window
  
#  lines(loggerDataTrim$dateTimeLogger,
#        rep(1, length(loggerDataTrim$dateTimeLogger)), col='grey50', lwd=1.5, lty='dotted')
#  lines(loggerDataTrim$dateTimeLogger,
#        rep(-1, length(loggerDataTrim$dateTimeLogger)), col='grey50', lwd=1.5, lty='dotted')
  
#   zeroset<-loggerDataTrim[loggerDataTrim$rollmeanslope > -1 & loggerDataTrim$rollmeanslope < 1,]
   #okay so the zero set tolerance for station 1 needs to be high.  Station 4 can be lower
   if(station=='Station_1') {
     zeroset<-loggerDataTrim[loggerDataTrim$rollmeanslope > -1 & loggerDataTrim$rollmeanslope < 1,] } else {
     zeroset<-loggerDataTrim[loggerDataTrim$rollmeanslope > -0.1 & loggerDataTrim$rollmeanslope < 0.1,]}
   
#   backgroundCond <- mean(loggerDataTrim$spCond[1:5], na.rm = TRUE)
   plateauset<-zeroset[zeroset$spCond>backgroundCond*1.2,]
# plateaureached<-min(zeroset$dateTimeLogger, na.rm = TRUE)
  #plateaucond<-loggerDataTrim$spCond[index]
  plateaucond<-mean(plateauset$spCond, na.rm=TRUE) 
    #back to the original dataset
 # index<-which(loggerDataTrim$dateTimeLogger==plateaureached)
  index<-min(which(loggerDataTrim$spCond >= (plateaucond - (0.03*plateaucond)) & 
                     loggerDataTrim$spCond <= (plateaucond + (0.03*plateaucond))), na.rm=TRUE)
  if(index==Inf | index == -Inf) {
    index <-min(which(loggerDataTrim$spCond >= (plateaucond - (0.05*plateaucond)) & 
                        loggerDataTrim$spCond <= (plateaucond + (0.05*plateaucond))), na.rm=TRUE) }
  if(index==Inf | index == -Inf) {
      cat('No matching plateau points within five percent of the plateau conductivity') 
    peakInfoOut <- list("station" = station,
                        "backgroundCond"=backgroundCond,
                        #       "peakArea"=areaUnderCurve,
                        "peakTime"=NA,
                        "peakCond" = NA,
                        "nominalTime" = NA,
                        "nominalCond" = NA,
                        "centroidTime"=NA,
                        #      "harmonicMeanTime"=harmonicMeanTime,
                        "endPlotTime"=loggerDataTrim$dateTimeLogger[endHere]) 
    return(peakInfoOut)}
    
    #and then trim to plateau
  if(!is.na(index)) {
    loggerDataTrim_pt<-loggerDataTrim[1:index,]
    plateaureached<-loggerDataTrim$dateTimeLogger[index]
    # plateaucond<-loggerDataTrim$spCond[index]
    nominaltime<-loggerDataTrim$dateTimeLogger[min(which(loggerDataTrim$spCond>(0.5*plateaucond)))]
    nominal<-0.5*plateaucond} else {
      cat("Warning, poor conductivity point match", currEventID)
      peakInfoOut <- list("station" = station,
                          "backgroundCond"=backgroundCond,
                          #       "peakArea"=areaUnderCurve,
                          "peakTime"=NA,
                          "peakCond" = NA,
                          "nominalTime" = NULL,
                          "nominalCond" = NA,
                          "centroidTime"=NULL,
                          #      "harmonicMeanTime"=harmonicMeanTime,
                          "endPlotTime"=loggerDataTrim$dateTimeLogger[endHere]) 
      return(peakInfoOut)}
  
  #okay now have the nominal time saved as nominal time, peak time as plateaureached, and plateau conductivity level as plateaucond
  #will need to check that this works for slug injections as well
  
  #also note - Bob's k600 Argon paper uses the nominal transport time, or time to reach 1/2 the plateau concentration
  #in the conductivity curves.  2018 argon paper
  

  
  #because I am having some issues, I think I will try trimming to where the increase starts also

  #####
  slopeindex<-min(which(loggerDataTrim_pt$slope >= 0 & 
                          loggerDataTrim_pt$corrSpCond >= (0.1*backgroundCond)), na.rm=TRUE)
  if(slopeindex == Inf | slopeindex == -Inf) {
    slopeindex <-min(which(loggerDataTrim_pt$slope >= 0 & 
                             loggerDataTrim_pt$corrSpCond >= (0.05*backgroundCond)), na.rm=TRUE) }
  if(slopeindex==Inf | slopeindex == -Inf) {
      cat('No matching points with a notable slope increase') 
    peakInfoOut <- list("station" = station,
                        "backgroundCond"=backgroundCond,
                        #       "peakArea"=areaUnderCurve,
                        "peakTime"=NA,
                        "peakCond" = NA,
                        "nominalTime" = NULL,
                        "nominalCond" = NA,
                        "centroidTime"=NULL,
                        #      "harmonicMeanTime"=harmonicMeanTime,
                        "endPlotTime"=loggerDataTrim$dateTimeLogger[endHere]) 
    return(peakInfoOut)}

  
  if(!is.na(slopeindex)) {
    loggerDataTrim_pt<-loggerDataTrim_pt[slopeindex:dim(loggerDataTrim_pt)[1],] } else {
      cat("Warning, no slope increase match", currEventID)
      peakInfoOut <- list("station" = station,
                          "backgroundCond"=backgroundCond,
                          #       "peakArea"=areaUnderCurve,
                          "peakTime"=NA,
                          "peakCond" = NA,
                          "nominalTime" = NULL,
                          "nominalCond" = NA,
                          "centroidTime"=NULL,
                          #      "harmonicMeanTime"=harmonicMeanTime,
                          "endPlotTime"=loggerDataTrim$dateTimeLogger[endHere]) 
      return(peakInfoOut)
      }
  
  slugInjTypes <- c("NaBr","model","model - slug")
  criInjTypes <- c("NaCl","model - CRI")
  #might need to convert later but I think this might work for both?
  
#  if(injectionType %in% criInjTypes){

  #for now, maybe use the first time the rolling mean of the slope is 0, but in the future program it so that it has to be within
  #some window of the 0 slope for the next few minutes?
  
 
  #Background correct the logger data
 # loggerDataTrim_pt$corrSpCond <- loggerDataTrim_pt$spCond - backgroundCond
  
  #Add the cumulative time difference column and the cumulative cond diff column
  loggerDataTrim_pt$cumTimeDiff <- as.numeric(difftime(loggerDataTrim_pt$dateTimeLogger, loggerDataTrim_pt$dateTimeLogger[1], units = "min"))
  loggerDataTrim_pt$cumCondDiff <- loggerDataTrim_pt$corrSpCond - loggerDataTrim_pt$corrSpCond[1]
  
  #Now loop through and calculate the centroid time (tc) and harmonic mean time (thm)
  areaUnderCurve <- pracma::trapz(loggerDataTrim_pt$cumTimeDiff, loggerDataTrim_pt$spCond)
  #this function calculates the trapezoidal integral....not sure where this is needed now
  
  #going to try to program the centroid and compare the times
  #first replot the spc in the chosen window just to check
  
  tc_sub <- 0
  hm <- 0
  L<-0
  Cc_sub<-0 #conductance at centroid
  
  for(i in 3:length(loggerDataTrim_pt$spCond)){
    t <- loggerDataTrim_pt$cumTimeDiff[i] #cumulative minutes passed
    dt <- loggerDataTrim_pt$cumTimeDiff[i] - loggerDataTrim_pt$cumTimeDiff[i-1] #dt in minutes
    slope2 <- (loggerDataTrim_pt$slope[i])^2 #slope was by minutes
  
    Ltoadd <- sqrt(1+slope2) * dt 
   # loggerDataTrim$tc[i] <- tcToAdd
    L <- L + Ltoadd
    tctoadd<-t * Ltoadd
    tc_sub<-tc_sub+tctoadd
    
    #trying to add a y point now
    c<-loggerDataTrim_pt$cumCondDiff[i]
    Cctoadd<-c*Ltoadd
    Cc_sub<- Cc_sub + Cctoadd
    
#    hmToAdd <- (1/t) * px * dt #1/cumulative time * time interval * cumulative time...time interval?
 #   loggerDataTrim$hm[i] <- hmToAdd
#  hm <- hm + hmToAdd
    
#    ti <- t
  }
  tc<-tc_sub/L #should be minutes
  tc #coming out to 22 minutes, looks a little too low? #but any number tweaking can change it by a few minutes...
  #will need to see how much it changes with the fourth station
  Cc<-Cc_sub/L
  Cc
  
  invisible(dev.new(noRStudioGD = TRUE, width=12, height=7))
  plot(loggerDataTrim_pt$dateTimeLogger,
       loggerDataTrim_pt$spCond,
       cex=0.7,
       col='black',
       xlab = "Measurement time",
       ylab = "SpCond",
       xlim = c(min(loggerDataTrim_pt$dateTimeLogger), max(loggerDataTrim_pt$dateTimeLogger)+120),
       ylim = c(min(loggerDataTrim_pt$spCond), max(loggerDataTrim_pt$spCond)+2),
       main = paste(currEventID, station, sep = '\n'))  
  
  #add in rolling window mean and some boundaries
  lines(loggerDataTrim_pt$dateTimeLogger,
        loggerDataTrim_pt$rollmean, col='coral2', lwd=1.5)
  
  centroidTime <- as_datetime(loggerDataTrim_pt$dateTimeLogger[1] + (tc*60))  #the time addition is in seconds, tc was in minutes
 
   points(x=as_datetime(c(nominaltime, plateaureached)), y=c(nominal, plateaucond),
         pch=17, col='blue', cex=2)
  points(centroidTime, Cc+backgroundCond, col='red', pch=17, cex=2) #better
  
  ans <- identify(loggerDataTrim_pt$dateTimeLogger,
                  loggerDataTrim_pt$spCond, n = 1, tolerance = 100, plot = F)
  invisible(dev.off())

#  Sys.sleep(1)
 # invisible(dev.off())
 
  #CURRENT STATUS: FIGURE OUT THE THM VS CENTROID TIME FORMULAS - leave this out for now
 # thm <- 1/hm

 
#  harmonicMeanTime <- loggerDataTrim$dateTimeLogger[1] + thm
  
  #Take a look at the outputs - all Kaelin's
 # plot(loggerDataTrim$dateTimeLogger, loggerDataTrim$originalSpCond)
#  lines(loggerDataTrim$dateTimeLogger, loggerDataTrim$test)
 # par(new = TRUE)
  #plot(loggerDataTrim$dateTimeLogger, loggerDataTrim$spCond, xlab = "", ylab = "", col = "blue", type = "l", axes = FALSE)
  #abline(v = peakTime, col = "cyan")
  #abline(v = tc+loggerDataTrim$dateTimeLogger[1], col = "green")
  #abline(v = thm+loggerDataTrim$dateTimeLogger[1], col = "purple")
  #graphics::legend(x = "right", legend = c("smoothed","derivative","peakTime","centroidTime","harmonicMeanTime"), lty = c(1,1,1,1), col = c("black","blue","cyan","green","purple"))
  
  peakInfoOut <- list("station" = station,
                    "backgroundCond"=backgroundCond,
               #       "peakArea"=areaUnderCurve,
                      "peakTime"=plateaureached,
                      "peakCond" = plateaucond,
                      "nominalTime" = nominaltime,
                      "nominalCond" = nominal,
                      "centroidTime"=centroidTime,
                #      "harmonicMeanTime"=harmonicMeanTime,
                      "endPlotTime"=loggerDataTrim$dateTimeLogger[endHere])
  #okay now have the nominal time saved as nominal time, peak time as plateaureached, and plateau conductivity level as plateaucond
  return(peakInfoOut)
}

