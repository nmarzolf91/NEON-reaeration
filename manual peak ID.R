##############################################################################################
#' @title Constant rate or slug travel time between two stations

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr
#' Amanda Gay DelVecchia \email{amanda.delvecchia@duke.edu} \cr

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
manualpeakTime <- function(
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
         y = c(highPlot*.95,highPlot*.9),
         col = c("darkolivegreen4", "tomato3"),
         lwd = 2,
         pch = 19,
         cex = 6)
  #title(main = paste0("Click green dot (upper lefthand) if the peak/plateau is identifiable. \nClick red dot (upper righthand) if not identifiable.\n",currEventID))
  badPlotBox <- identify(x = c(max(loggerDataTrim$dateTimeLogger),max(loggerDataTrim$dateTimeLogger)),
                         y = c(highPlot*.95,highPlot*.9),
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
    #slopeincrease<-readline(prompt='What is the minimum slope (SpC/minute) to begin integration?')
   # zerotol<-readline(prompt='What is the max tolerance range for the derivate = 0 point?')
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
  
  #mean background conductivity
  backgroundCond <- mean(loggerDataTrim$spCond[1:5], na.rm = TRUE)
  
  #Background correct the logger data
  loggerDataTrim$corrSpCond <- loggerDataTrim$spCond - backgroundCond
  
  ##how to ask the user to specify a number of 10 second measurements to use in the derivative estimate
  

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
 
  #now the slope is based on 50 seconds, and calculating a rolling mean of slope
  #13 here = 130 seconds averaged for rolling mean. 
  loggerDataTrim$rollmeanslope<-rollmean(loggerDataTrim$slope, 13, fill = list(rep(NA, 6), NULL, rep(NA,6)))
 
 
    invisible(dev.new(noRStudioGD = TRUE, width=12, height=7))
    plot(loggerDataTrim$dateTimeLogger,
         loggerDataTrim$spCond,
         xlab = "Measurement time",
         ylab = "Specific Conductance")
         #ylim = c(lowPlot, highPlot))
    lines(loggerDataTrim$dateTimeLogger,
          loggerDataTrim$rollmean, col='coral2', lwd=1.5)
    title(main = paste0("Click the peak"))
    manualpeak <- identify(x = loggerDataTrim$dateTimeLogger,
                    y = loggerDataTrim$spCond,
                    n = 1,
                    tolerance = 0.4)
    #slopeincrease<-readline(prompt='What is the minimum slope (SpC/minute) to begin integration?')
    # zerotol<-readline(prompt='What is the max tolerance range for the derivate = 0 point?')
    Sys.sleep(1)
    invisible(dev.off())
  
    loggerDataTrim_pt<-loggerDataTrim #[1:index,] #index = start of plateau
    peakreached<-loggerDataTrim$dateTimeLogger[manualpeak]
    peakcond<-loggerDataTrim$spCond[manualpeak]-backgroundCond
   
    # plateaucond<-loggerDataTrim$spCond[index]
    nominaltime<-loggerDataTrim$dateTimeLogger[min(which(loggerDataTrim$corrSpCond>(0.5*(peakcond))))]  #peakcond already corrected
    nominal<-(0.5*(peakcond)) + backgroundCond
  

  
  #Add the cumulative time difference column and the cumulative cond diff column
  loggerDataTrim_pt$cumTimeDiff <- as.numeric(difftime(loggerDataTrim_pt$dateTimeLogger, loggerDataTrim_pt$dateTimeLogger[1], units = "min"))
  loggerDataTrim_pt$cumCondDiff <- loggerDataTrim_pt$corrSpCond - loggerDataTrim_pt$corrSpCond[1]
  
#saving the area for the discharge calculation
  areaUnderCurve <- pracma::trapz(loggerDataTrim_pt$cumTimeDiff, loggerDataTrim_pt$spCond)
  #this function calculates the trapezoidal integral

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
  
  points(x=as_datetime(c(nominaltime, peakreached)), y=c(nominal, peakcond + backgroundCond),
         pch=17, col='blue', cex=2)
  
  ans <- identify(loggerDataTrim_pt$dateTimeLogger,
                  loggerDataTrim_pt$spCond, n = 1, tolerance = 100, plot = F)
  invisible(dev.off())
  
  #  Sys.sleep(1)
  # invisible(dev.off())
  
  
  peakInfoOut <- list("station" = station,
                      "backgroundCond"=backgroundCond,
                      #       "peakArea"=areaUnderCurve,
                      "peakTime"=peakreached,
                      "peakCond" = peakcond + backgroundCond,
                      "nominalTime" = nominaltime,
                      "nominalCond" = nominal,
                      'area' = areaUnderCurve,
                      "centroidTime"= NA,
                      "centroidCond" = NA,
                      #      "harmonicMeanTime"=harmonicMeanTime,
                      "endPlotTime"=loggerDataTrim$dateTimeLogger[endHere])
  #okay now have the nominal time saved as nominal time, peak time as plateaureached, and plateau conductivity level as plateaucond
  return(peakInfoOut)
} 


