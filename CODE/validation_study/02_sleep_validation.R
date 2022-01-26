

library( zoo )
library( hms )
library( data.table )
library( stringr )
library( lubridate )
library( sp )
library( stats )
library(rgeos)
library( plyr )

################# Functions #########################

## function for normalizing a vector
normalize_func <- function( x ) return( (x - mean( x, na.rm = T ) )/ sd( x, na.rm = T ) )

## function for adding age-sex class
addAgeSex<-function(df,vector){
  df$ageSex<-NA
  print("Make sure to assign this function to the dataframe you want updated. i.e. babdat<- addAgeSex(babdat,'tags')")
  for(i in 1:nrow(df)){
    if(vector[i]=="2427" ||vector[i]=="2434" || vector[i]== "2457"){
      df$ageSex[i]<-"Male Adult"
    }else{
      if(vector[i]=="2426" || vector[i]=="2432" || vector[i]=="2433" || vector[i]=="2436" || vector[i]=="2443" || vector[i]=="2452"){
        df$ageSex[i]<-"Male Subadult"
      }else{
        if(vector[i]=="2448"|| vector[i]=="2454"){
          df$ageSex[i]<-"Male Juvenile"
        }else{
          if(vector[i]=="2430"|| vector[i]=="2439"||vector[i]=="2446"||vector[i]=="2447"|| vector[i]=="2449"|| vector[i]=="2451"|| vector[i]=="2453"|| vector[i]=="2456"|| vector[i]=="2458"|| vector[i]=="2460"){
            df$ageSex[i]<-"Female Adult"
          }else{
            df$ageSex[i]<-"Female Subadult"
          }
        }
      }
    }
  }
  return(df)
}

## function for turning tag_names into colors
babcolor<-function(IDvec){
  outputVec<-c()
  for(i in 1:length(IDvec)){
    if(IDvec[i]=="2427" ||IDvec[i]=="2434" || IDvec[i]== "2457"){
      outputVec[i]<-'blue'
    }else{
      if(IDvec[i]=="2426" || IDvec[i]=="2432" || IDvec[i]=="2433" || IDvec[i]=="2436" || IDvec[i]=="2443" || IDvec[i]=="2452"){
        outputVec[i]<-'skyblue'
      }else{
        if(IDvec[i]=="2448"|| IDvec[i]=="2454"){
          outputVec[i]<-'grey'
        }else{
          if(IDvec[i]=="2430"|| IDvec[i]=="2439"||IDvec[i]=="2446"||IDvec[i]=="2447"|| IDvec[i]=="2449"|| IDvec[i]=="2451"|| IDvec[i]=="2453"|| IDvec[i]=="2456"|| IDvec[i]=="2458"|| IDvec[i]=="2460"){
            outputVec[i]<-'red'
          }else{
            outputVec[i]<-'pink'
          }
        }
      }
    }
  }
  return(outputVec)
}

## function for setting transparency of a color while plotting
transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

## function for plotting times from noon to noon. It will make correct 12:00 - 24:00 to be 00:00 - 12:00 and 00:00 - 12:00 to be 12:00 to 24:00

ts_func <- function( time_vec ){
  
  num_time <- as.numeric( as_hms( time_vec ) )
  
  corr_time <- ifelse( num_time < 12*60*60, num_time + 12*60*60, num_time - 12*60*60 )
  
  return( corr_time )
  
}

#### sleep period function

## missing_mins: this is the maximum total number of minutes of data that can be missing from a day and still have that day included in the analysis

## time_gap: this is the maximum allowable time gap between two accelerometer bursts (in seconds) that can exist in a day without removing this day from the data

## move_window: this is the size of the moving window (in minutes) used in calculating the rolling median 

## percentile: this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity 

## multiplier: this is the multiplier of the threshold value determined by the percentile above. Values below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active

## block_size: duration in minutes of the blocks of continuous inactivity that will be considered sleep

## gap_size: maximum duration between sleep blocks that will be merged

## title: if title == T, the plot will include the tag number and the night number at the top of the plot

## x_axis: if x_axis == F, the plot will be plotted without an x axis
## waso_percentile: this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity within the sleep period

## waso_multiplier: this is the multiplier of the threshold value determined by the percentile above. Values WITHIN the sleep period below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active. This threshold value may be different than the one above even when waso_percentile = percentile and waso_multiplier = multiplier, because the value of the given percentile of this variable depends on the smoothing window (see waso_window below; aka waso window might not be equal to mov_window)

## waso_window: this is the size of the moving window (in minutes) used in calculating the rolling median that will be used to find periods of wake after sleeping. A waso_window of 1 is the same as using the raw data without a rolling median

## waso_block: this is the number of consecutive minutes of inactivity needed to classify a period as sleep within the sleep period. A waso_block of 1 means that anytime the value is below the threshold, the baboon in considered sleeping and anytime the value is above the threshold the baboon is considered awake

## las: las sets las in the plotting window


sleep_per_func <- function( tag, night, missing_mins = 45, time_gap = 20*60, mov_window = 9, percentile = 0.10, multiplier = 1.125, block_size = 30, gap_size = 45, title = F, x_axis = T, plot_waso = F, waso_window = 1, waso_block = 1, waso_percentile = 0.10, waso_multiplier = 1.125, las = 1, ... ){
  
  ## save a variable denoting the total number of minutes in the day
  mins_in_day <- 60*24
  
  ## subset the data to the given tag on the given night
  night_dat <- d1[ d1$tag == tag & d1$night == night, ]
  
  ## sort the timestamps (they are probably already sorted)
  sorted_times <- sort( night_dat$local_time )
  
  ## find the time difference in seconds between each burst
  time_diffs <- as.numeric( diff( as_hms( sorted_times ) ), units = 'secs' )
  
  ## if there is more than a single burst...
  if(length(time_diffs) != 0){
    
    ## if the number of bursts exceed the minimum required number of bursts in a night (determined by missing mins) and if the gaps in the data are within the allowable gap size (determined by time_gap)...
    if( nrow( night_dat ) > ( mins_in_day - missing_mins ) & max( time_diffs ) < time_gap ){
      
      ## take the rolling median of the log VeDBA and save it as a column
      night_dat$roll_log_vedba <- rollmedian( night_dat$log_vedba, mov_window, fill = NA, align = 'center' )
      
      ## determine the threshold activity vs. inactivity threshold based on the percentile, multiplier, and the rolling median just produced
      thresh <- quantile( night_dat$roll_log_vedba, percentile, na.rm = T) * multiplier
      
      ## put the rows of the dataframe in order from noon to noon (they should already be in this order, so this should be redundant)
      night_dat <- night_dat[ order( ts_func( night_dat$local_time ) ), ]
      
      ## turn the times into numerical elements for plotting
      ts_time <- ts_func( night_dat$local_time )
      
      if( title == F ){
        ## plot the log VeDBA
        #plot( ts_time, night_dat$log_vedba, type = 'l', xlab = 'Time', ylab = '', xaxt = 'n', las = las )
        
        plot( ts_time, night_dat$log_vedba, type = 'l', xlab = 'Time', ylab = '', xaxt = 'n', las = las, ylim = c( 1.9, 7 ), ... )
        
      }else{
        ## plot the log VeDBA
        plot( ts_time, night_dat$log_vedba, type = 'l', xlab = 'Time', ylab = '', main = paste( tag, night ), xaxt = 'n', las = las, ... )
        
        
      }
      
      if( x_axis == T ){
        
        axis( 1, at = seq( 0, 60*24*60, 60*60), labels = c( as_hms( seq( 12*60*60, 60*23*60, 60*60) ), as_hms( seq( 0, 60*12*60, 60*60) ) ) ) 
        
      }
      
      title( ylab = 'log VeDBA', line = 3.9 )
      ## plot the rolling median of the log VeDBA
      #lines( ts_time, rollmedian( night_dat$log_vedba, mov_window, fill = NA, align = 'center' ), col = 'red')
      
      ## plot the threshold of the log VeDBA
      #abline( h = thresh, col = 'blue' )
      
      ### find blocks of continuous inactivity
      
      ## find the run length encoding of periods above and below the threshold
      temp <- rle(as.numeric( night_dat$roll_log_vedba < thresh ) ) 
      
      ## mark the rows that are part of runs (i.e. part of chunks that are greater than the block_size of either continuous activity or continuous inactivity )
      night_dat$runs <- as.numeric( rep( temp$lengths > block_size, times = temp$lengths ) )
      
      ## mark the rows corresponding to sleep bouts. These sleep bouts are runs of inactivity
      night_dat$sleep_bouts <- as.numeric( night_dat$roll_log_vedba < thresh & night_dat$runs == 1 )
      
      ## find when sleep bouts start and end
      diffs <- diff( c(0, night_dat$sleep_bouts ) )
      starts <- which( diffs == 1 ) [ -1 ]
      ends <- which( diffs == -1 )
      
      ## if there are any sleep bouts...
      if( length( which( diffs == 1 ) ) != 0){
        
        ## find the duration of the gaps between each sleep bout (the end of one sleep bout and the start of the next)
        gaps <- as.numeric( night_dat$local_timestamp [ starts ] - night_dat$local_timestamp [ ends[ 1: length( starts ) ] ], units = 'mins' )
        
        ## sleep bouts separated by gaps that are shorter than that specified by gap_size will be merged. Note which of these gaps are shorter than the gap_size
        inds_to_remove <- which( gaps < gap_size ) 
        
        ## if there are NO gaps between sleep bouts that are to be removed...
        if( length( inds_to_remove ) == 0 ){
          
          ## set sleep onset index to be the start of sleep bouts
          onset <- which( diffs == 1 ) 
          
          ## set waking index to be the end of sleep bouts
          wake <- ends
          
        }else{ ## if there ARE gaps between sleep bouts that are to be removed...
          
          ## set sleep onset index to be the start of sleep bouts that do not correspond to the gaps to be removed (because these will be within sleep periods, not a start of a new bout)
          onset <- which( diffs == 1 ) [ - (inds_to_remove + 1) ]
          
          ## set waking index to be the end of sleep bouts that do not correspond to the gaps to be removed
          wake <- ends [ - inds_to_remove ]
          
        }
        
        ## determine which sleep period is the longest
        per_ind <- which.max( as.numeric( night_dat$local_timestamp[ wake ] - night_dat$local_timestamp[ onset ], units = 'secs' ) )
        
        ## plot the sleep onset time and waking time on the log VeDBA plot
        abline( v = c( ts_time[ onset[ per_ind ] ], ts_time[ wake[ per_ind ] ] ), col = 'orange', lty = 3, lwd = 4 )
        
        ## if you also want to plot WASO
        if( plot_waso == T ){
          
          ## calculate the rolling median of the log VeDBA using the waso_window
          night_dat$SPT_roll_log_vedba <- rollmedian( night_dat$log_vedba, waso_window, fill = NA, align = 'center' )
          
          ## plot this rolling median
          #lines( ts_time, night_dat$SPT_roll_log_vedba, col = 'red', lty = 1 )
          
          ## calculate the threshold for sleeping and waking within the sleep period
          SPT_thresh <- quantile( night_dat$SPT_roll_log_vedba, waso_percentile, na.rm = T) * waso_multiplier
          
          ## plot the threshold
          #abline( h = SPT_thresh, col = 'blue', lty = 1, lwd = 2 )
          
          ## subset the night's data to only the sleep period
          trim_night <- night_dat[ night_dat$local_timestamp >= night_dat$local_timestamp[ onset[ per_ind ] ] & night_dat$local_timestamp <= night_dat$local_timestamp[ wake[ per_ind ] ] , ]
          
          ### find blocks of continuous inactivity
          
          ## calcuate the run length encoding
          temp <- rle(as.numeric( trim_night$SPT_roll_log_vedba < SPT_thresh ) ) 
          
          ## mark the runs of activity or inactivity
          trim_night$night_runs <- as.numeric( rep( temp$lengths >= waso_block, times = temp$lengths ) )
          
          ## mark the runs of inactivity as sleep bouts
          trim_night$night_sleep_bouts <- as.numeric( trim_night$SPT_roll_log_vedba < SPT_thresh & trim_night$night_runs == 1 )
          
          ## find the starts and ends of waking bouts
          diffs <- diff( c(1, trim_night$night_sleep_bouts ) )
          
          starts <- which( diffs == -1 )
          
          ## add back in a "- 1" at the end of this line if you wish for the start and end times to be accurate. Now I just want to make it so the polygons show up even without a border
          ends <- which( diffs == 1 )
          
          ## if there are waking bouts
          if( length( which( diffs == -1 ) ) != 0){
            
            ## if the last waking bout never ends...
            if( length( starts ) != length( ends ) ){
              
              ## make it end at the end of the sleep period
              ends <- c(ends, length( diffs) )
              
            }
            
            ## save the start times and end times of waking bouts
            starts <- ts_func( trim_night$local_time[ starts ] )
            ends <- ts_func( trim_night$local_time[ ends ] )
            
            ## plot a polygon for each distinct waking bout
            for( n in 1:length( starts ) ){
              
              polygon( x = c(starts[ n ], ends[ n ], ends[ n ], starts[ n ], starts[ n ] ), y = c( 0, 0, 10, 10, 0), col = transp('blue', .25), border = NA )
              
            }
          }
        }
        ## fill in the sleep period data frame with the sleep onset and waking time associated with the longest sleep period in the day (noon to noon)
        return( c( night_dat$local_timestamp[ onset[ per_ind ] ], night_dat$local_timestamp[ wake[ per_ind ] ] ) )
      } 
    }
  } 
}

################## Read in the d1 (accelerometer burst) data ###################

## d1 is a dataframe with a row for each minute for each baboon. Each row contains the raw (or interpolated) GPS acc burst, and several different measures calculated from bursts (like VeDBA)
d1 <- fread("DATA/validation_study/2019_full_night_and_day_data.csv")

## turn the data table into a dataframe
d1 <- as.data.frame( d1 )

## turn timestamp into POSIX element and time into character
d1$timestamp <- as.POSIXct( d1$timestamp, tz = 'UTC' )
d1$time <- as.character( d1$time )

## change times to local time here
d1$local_timestamp <- d1$timestamp + 3*60*60

## make a column for local time
d1$local_time <- str_split_fixed(d1$local_timestamp, " ", 2)[,2]

####################### Determining sleep periods with modification of Van Hees et al. 2018 method ####################################

## assign each minute of data to a given night. A night lasts from noon to noon. First, apply a time shift so that each night is a unit, and not each day
time_shift <- d1$local_timestamp - 12*60*60

## save the date of the first night of the study (the date of the night is always the date of the evening at the beginning of that night; so the first night of the study is 2012-07-31, although the data starts on 2012-08-01, because the data on that first morning is still technically part of the data for the previous night, as a night is noon to noon)
start_date <- as.Date(min(d1$local_timestamp)- 12*60*60)

## assign night as number of nights from the start of the study, with all data before the first noon representing night 1
d1$night <- as.numeric( as.Date(time_shift) - start_date + 1 )

## show how many baboon-nights there are
nrow( unique( d1[ , c( 'tag', 'night' ) ] ) )

## check where the night changes from one night to the next to see if it is at noon
d1[(diff(c( d1$night )) == 1),]


## save a variable denoting the total number of minutes in the day
mins_in_day <- 60*24

missing_mins <- 120 ## this is the maximum total number of minutes of data that can be missing from a day and still have that day included in the analysis (for sleep period time and sleep based analyses; i.e. not ave_vedba)

night_missing_mins <- 45 ## this is the maximum total number of minutes of data that can be missing from a dark period and still have that day included in the analysis

time_gap <- 20*60 ## this is the maximum allowable time gap between two accelerometer bursts (in seconds) that can exist in a noon-to-noon period without removing this noon-to-noon period from the data

mov_window <- 9 ## this is the size of the moving window (in minutes) used in calculating the rolling median of the average VeDBA

percentile <- 0.10 ## this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity (VeDBA here)

multiplier <- 1.125 ## this is the multiplier of the threshold value determined by the percentile above. Values below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active

block_size <- 30 ## duration in minutes of the blocks of continuous inactivity that will be considered sleep

gap_size <- 45 ## maximum duration between sleep blocks that will be merged

waso_percentile <- 0.10 ## this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity within the sleep period

waso_multiplier <- 1.125 ## this is the multiplier of the threshold value determined by the percentile above. Values WITHIN the sleep period below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active. This threshold value may be different than the one above even when waso_percentile = percentile and waso_multiplier = multiplier, because the value of the given percentile of this variable depends on the smoothing window (see waso_window below; aka waso window might not be equal to mov_window)

waso_window <- 1 ## this is the size of the moving window (in minutes) used in calculating the rolling median that will be used to find periods of wake after sleeping. A waso_window of 1 is the same as using the raw data without a rolling median

waso_block <- 3 ## this is the number of consecutive minutes of inactivity needed to classify a period as sleep within the sleep period. A waso_block of 1 means that anytime the value is below the threshold, the baboon in considered sleeping and anytime the value is above the threshold the baboon is considered awake

nap_block <- waso_block ## this is the number of consecutive minutes of inactivity needed to classify a period as sleep during the waking period. A nap_block of 1 means that anytime the value is below the threshold, the baboon in considered sleeping and anytime the value is above the threshold the baboon is considered awake

sep_day_night <- F ## this determines whether sleep periods and non-sleep periods are separated before finding runs of inactivity to consider as sleep

## create a vector containing the names of each baboon
tag_names <- unique( d1$tag )

## make a copy of d1. We will fill in this new dataframe with information about if the baboon was asleep in each epoch
full_dat <- d1

full_dat$sleep_per <- NA ## binary indicating whether a row belongs to the sleep period window
full_dat$sleep_bouts <- NA ## binary indicating whether the row is considered sleep, based on the waso or nap requirements
full_dat$n_bursts <- NA ## the number of bursts collected in a given noon-to-noon period (to be compared to the total number of minutes in the day). This column will indicate whether the data for a given night is insufficient to calculate the sleep period (and thus: onset, waking, SPT, sleep_eff, TST, sleep_bouts -- because this depends on a percentile of bursts' log vedba, WASO, wake_bouts, summed_VeDBA, night_VeDBA_corr, dark_TST, prev_naps, prev_day_sleep)
full_dat$max_time_diff <- NA ## the maximum difference between consecutive fixes in a given noon-to-noon period. With the previous column, this column will indicate whether the data is insufficient to calculate the sleep period (and thus: onset, waking, SPT, sleep_eff, TST, WASO, wake_bouts, summed_VeDBA, night_VeDBA_corr, prev_naps )
full_dat$dark <- NA ## a binary indicating whether a row belongs to the period of darkness (between the end of astrological evening twilight and the beginning of astrological morning twilight)
full_dat$poss_dark_bursts <- NA ## the number of potential bursts in the dark period on this night, i.e. the number of minutes between the start and end of the night
full_dat$n_dark_bursts <- NA ## the total number of bursts actually taken during the dark period on this night. This will be compared to the previous column to determine whether the data during the dark period is sufficient to calculate the dark_TST and ave_vedba
full_dat$max_dark_time_diff <- NA ## the maximum difference between consecutive fixes in a given dark period. With the previous column, this column will indicate whether the data is insufficient to calculate the dark_TST and ave_vedba
full_dat$poss_day_bursts <- NA ## the number of potential bursts in preceding day (light) period, i.e. the number of minutes between the end of morning astrological twilight (maybe this should be changed to sunrise?) and the start of evening astrological twilight (maybe this should be changed to sunset?)
full_dat$n_day_bursts <- NA ## the total number of bursts actually taken during the preceding day (light) period. This will be compared to the previous column to determine whether the data during the day period is sufficient to calculate the prev_day_sleep and prev_day_ave_vedba.
full_dat$max_day_time_diff <- NA ## the maximum difference between consecutive fixes in a given day (light) period. With the previous column, this column will indicate whether the data is insufficient to calculate the prev_day_sleep and prev_day_ave_vedba
## prev_naps depends on having sufficient data to calculate both the current SPT as well as the previous day's SPT


## create a vector containing the names of each baboon
tag_names <- unique( d1$tag )

## for each individual...
for( tag in tag_names ){
  
  ## subset the data to just this individual's data
  id_dat <- d1[ d1$tag == tag, ]
  
  ## create a vector the nights for which this individual has data
  nights <- unique( id_dat$night )
  
  ## for each night on which this individual has data
  for( night in nights ){
    
    ## subset this individual's data to just that night
    night_dat <- id_dat[ id_dat$night == night, ]
    
    ## create empty columns for the sleep period and sleep bout binary variables
    night_dat$sleep_per <- NA
    night_dat$sleep_bouts <- NA
    
    
    ## save a column of the total number of bursts for that day. This will also make it easier to remove these days from the dataframe later
    night_dat$n_bursts <- nrow( night_dat )
    
    ## sort the timestamps (they are probably already sorted)
    sorted_times <- c( '00:00:00', sort( night_dat$local_time ), '23:59:00' )
    
    ## find the time difference in seconds between each burst
    time_diffs <- as.numeric( diff( as_hms( sorted_times ) ), units = 'secs' )
    
    if( length( time_diffs ) > 0 ){ ### There is one night for one baboon with only one single burst, which is why this if statement is needed
      
      ## save a column of the maximum time difference between burst for that day (this will make it easier to pull out days with insufficient data later)
      night_dat$max_time_diff <- max( time_diffs )
      
    }else{
      
      night_dat$max_time_diff <- NA
      
    }
    
    ## take the rolling median of the log VeDBA and save it as a column
    roll_log_vedba <- rollmedian( night_dat$log_vedba, mov_window, fill = NA, align = 'center' )
    
    ## determine the threshold activity vs. inactivity threshold based on the percentile, multiplier, and the rolling median just produced
    thresh <- quantile( roll_log_vedba, percentile, na.rm = T) * multiplier
    
    ### find blocks of continuous inactivity
    
    ## find the run length encoding of periods above and below the threshold
    temp <- rle(as.numeric( roll_log_vedba < thresh ) ) 
    
    ## mark the rows that are part of runs (i.e. part of chunks that are greater than the block_size of either continuous activity or continuous inactivity )
    sleep_per_runs <- as.numeric( rep( temp$lengths > block_size, times = temp$lengths ) )
    
    ## mark the rows corresponding to sleep bouts. These sleep bouts are runs of inactivity
    sleep_per_sleep_bouts <- as.numeric( roll_log_vedba < thresh & sleep_per_runs == 1 )
    
    ## find when sleep bouts start and end
    diffs <- diff( c(0, sleep_per_sleep_bouts ) )
    starts <- which( diffs == 1 ) [ -1 ]
    ends <- which( diffs == -1 )
    
    ## if there are any sleep bouts...
    if( length( which( diffs == 1 ) ) != 0 ){
      
      ## find the duration of the gaps between each sleep bout (the end of one sleep bout and the start of the next)
      gaps <- as.numeric( night_dat$local_timestamp [ starts ] - night_dat$local_timestamp [ ends[ 1: length( starts ) ] ], units = 'mins' )
      
      ## sleep bouts separated by gaps that are shorter than that specified by gap_size will be merged. Note which of these gaps are shorter than the gap_size
      inds_to_remove <- which( gaps < gap_size ) 
      
      ## if there are NO gaps between sleep bouts that are to be removed...
      if( length( inds_to_remove ) == 0 ){
        
        ## set sleep onset index to be the start of sleep bouts
        onset <- which( diffs == 1 ) 
        
        ## set waking index to be the end of sleep bouts
        wake <- ends
        
      }else{ ## if there ARE gaps between sleep bouts that are to be removed...
        
        ## set sleep onset index to be the start of sleep bouts that do not correspond to the gaps to be removed (because these will be within sleep periods, not a start of a new bout)
        onset <- which( diffs == 1 ) [ - (inds_to_remove + 1) ]
        
        ## set waking index to be the end of sleep bouts that do not correspond to the gaps to be removed
        wake <- ends [ - inds_to_remove ]
        
      }
      
      ## determine which sleep period is the longest
      per_ind <- which.max( as.numeric( night_dat$local_timestamp[ wake ] - night_dat$local_timestamp[ onset ], units = 'secs' ) )
      
      ## fill in the sleep period data frame with the sleep onset and waking time associated with the longest sleep period in the day (noon to noon)
      
      night_dat$sleep_per <- as.numeric( night_dat$local_timestamp >= night_dat$local_timestamp[ onset[ per_ind ] ] & night_dat$local_timestamp <= night_dat$local_timestamp[ wake[ per_ind ] ] )
      
    }else{ ## if there aren't any sleep bouts, record all rows as a 0 in the sleep_per column
      
      night_dat$sleep_per <- 0
      
    }
    
    ## take the rolling median of the log VeDBA
    night_dat$roll_log_vedba <- rollmedian( night_dat$log_vedba, waso_window, fill = NA, align = 'center' )
    
    ### find blocks of continuous inactivity
    
    if( sep_day_night == T ){
      
      # first during the night. Splitting this up by day and night will allow us to use different threshold dor the duration of inactivity required to be considered sleep between day and night
      
      ## subset the night data to only the sleep period window
      trim_night <- night_dat[ night_dat$sleep_per == 1, ]
      
      ## find the run length encoding of periods above and below the threshold
      temp <- rle(as.numeric( trim_night$roll_log_vedba < thresh ) ) 
      
      ## mark which rows are part of runs of activity vs. inactivity, as determined by waso_block
      runs <- as.numeric( rep( temp$lengths >= waso_block, times = temp$lengths ) )
      
      sleep_bouts_night <- as.numeric( trim_night$roll_log_vedba < thresh & runs == 1 )
      
      ## mark which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
      trim_night$sleep_bouts <- sleep_bouts_night
      
      ## return this sleep period data back into night_dat
      night_dat[ night_dat$sleep_per == 1, ] <- trim_night
      
      if( sum( night_dat$sleep_per ) > 0 ){ ## if there are sleep periods
        
        ## subset the night data to only the times outside the sleep period window
        trim_day_before <- night_dat[ night_dat$sleep_per != 1 & night_dat$local_timestamp < min( night_dat$local_timestamp[ night_dat$sleep_per == 1 ] ), ]
        
        ## find the run length encoding of periods above and below the threshold
        temp <- rle(as.numeric( trim_day_before$roll_log_vedba < thresh ) ) 
        
        ## mark which rows are part of runs of activity vs. inactivity, as determined by nap_block
        runs <- as.numeric( rep( temp$lengths >= nap_block, times = temp$lengths ) )
        
        
        sleep_bouts_day_before <- as.numeric( trim_day_before$roll_log_vedba < thresh & runs == 1 )
        
        ## mark which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
        trim_day_before$sleep_bouts <- sleep_bouts_day_before
        
        ## return this data from outside the sleep period window back to night_data
        night_dat[ night_dat$sleep_per != 1 & night_dat$local_timestamp < min( night_dat$local_timestamp[ night_dat$sleep_per == 1 ] ), ] <- trim_day_before
        
        
        ## subset the night data to only the times outside the sleep period window
        trim_day_after <- night_dat[ night_dat$sleep_per != 1 & night_dat$local_timestamp > max( night_dat$local_timestamp[ night_dat$sleep_per == 1 ] ), ]
        
        ## find the run length encoding of periods above and below the threshold
        temp <- rle(as.numeric( trim_day_after$roll_log_vedba < thresh ) ) 
        
        ## mark which rows are part of runs of activity vs. inactivity, as determined by nap_block
        runs <- as.numeric( rep( temp$lengths >= nap_block, times = temp$lengths ) )
        
        sleep_bouts_day_after <- as.numeric( trim_day_after$roll_log_vedba < thresh & runs == 1 )
        
        ## mark which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
        trim_day_after$sleep_bouts <- sleep_bouts_day_after
        
        ## return this data from outside the sleep period window back to night_data
        night_dat[ night_dat$sleep_per != 1 & night_dat$local_timestamp > max( night_dat$local_timestamp[ night_dat$sleep_per == 1 ] ), ] <- trim_day_after
        
        
      }else{
        
        ## subset the night data to only the times outside the sleep period window
        trim_day <- night_dat[ night_dat$sleep_per != 1, ]
        
        ## find the run length encoding of periods above and below the threshold
        temp <- rle(as.numeric( trim_day$roll_log_vedba < thresh ) ) 
        
        ## mark which rows are part of runs of activity vs. inactivity, as determined by nap_block
        runs <- as.numeric( rep( temp$lengths >= nap_block, times = temp$lengths ) )
        
        sleep_bouts_day <- as.numeric( trim_day$roll_log_vedba < thresh & runs == 1 )
        
        ## mark which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
        trim_day$sleep_bouts <- sleep_bouts_day
        
        
        ## return this data from outside the sleep period window back to night_data
        night_dat[ night_dat$sleep_per != 1, ] <- trim_day
        
      }
      
    }else{
      
      ## find the run length encoding of periods above and below the threshold
      temp <- rle(as.numeric( night_dat$roll_log_vedba < thresh ) ) 
      
      ## mark which rows are part of runs of activity vs. inactivity, as determined by waso_block
      runs <- as.numeric( rep( temp$lengths >= waso_block, times = temp$lengths ) )
      
      sleep_bouts <- as.numeric( night_dat$roll_log_vedba < thresh & runs == 1 )
      
      ## make which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
      night_dat$sleep_bouts <- sleep_bouts
      
    }
    
    
    night_dat <- night_dat[ ,  names( night_dat ) != 'roll_log_vedba' ]
    
    
    ### put the night data back into full_dat
    full_dat[ full_dat$tag == tag & full_dat$night == night, ] <- night_dat
  }
  
}

pre_clean_full <- full_dat

study_nights <- min( d1$night ):max( d1$night )

day_lim_start_time <- "07:30:00"

day_lim_end_time <- "17:30:00"

sleep_per <- data.frame( tag = rep( unique( d1$tag ), each = length( study_nights ) ), night = rep( study_nights, times = length( tag ) ), onset = NA, waking = NA, SPT = NA, WASO = NA, TST = NA, sleep_eff = NA, wake_bouts = NA, summed_VeDBA = NA, night_VeDBA_corr = NA, ave_vedba = NA, dark_TST = NA, light_TST = NA, dark_sleep_eff = NA, light_sleep_eff = NA, prev_naps = NA, prev_day_sleep = NA, prev_day_sleep_lim = NA, prev_day_ave_vedba = NA, poss_dark_bursts = NA, n_dark_bursts = NA, poss_day_bursts = NA, n_day_bursts = NA, max_time_diff = NA, n_bursts= NA, max_dark_time_diff = NA, max_day_time_diff = NA )


## create empty vectors for the durations of sleep and wake bouts. We will fill these in to see if the distributions of the durations of these bouts later
sleep_durs <- c()
wake_durs <- c() 


## for each individual...
for( tag in tag_names ){
  
  ## subset the data to just this individual's data
  id_dat <- full_dat[ full_dat$tag == tag, ]
  
  ## create a vector the nights for which this individual has data
  nights <- unique( id_dat$night )
  
  ## for each night on which this individual has data
  for( night in nights ){
    
    ## subset this individual's data to just that night
    night_dat <- id_dat[ id_dat$night == night, ]
    
    ## should already be in order, but just in case
    night_dat <- night_dat[ order( night_dat$local_timestamp ), ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$n_bursts <- unique( night_dat$n_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$max_time_diff <- unique( night_dat$max_time_diff )
    
    
    SPT_dat <- night_dat[ night_dat$sleep_per == 1, ]
    
    if( nrow( SPT_dat ) > 0 ){
      
      onset <- min( SPT_dat$local_timestamp )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$onset <- onset
      
      waking <- max( SPT_dat$local_timestamp ) 
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$waking <- waking
      
      SPT <- as.numeric( waking - onset, units = 'mins' ) + 1
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$SPT <- SPT
      
      WASO <- sum( SPT_dat$sleep_bouts == 0 )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$WASO <- WASO
      
      TST <- sum( SPT_dat$sleep_bouts == 1 )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$TST <- TST
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$sleep_eff <- TST/ nrow( SPT_dat )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$summed_VeDBA <- sum( SPT_dat$log_vedba )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$night_VeDBA_corr <- sum( SPT_dat$log_vedba ) / SPT
      
      ## find the distinct sleep bouts (i.e. epochs of sleep separated by waking)
      diffs <- diff( c( 0, SPT_dat$sleep_bouts ) )
      
      ## save the number of distinct wake bouts
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$wake_bouts <- sum( diffs == -1 )
      
      ## find durations of sleep and wake bouts
      temp <- rle( SPT_dat$sleep_bouts )
      
      ## add the duration of sleep bouts to the sleep bout duration vector
      sleep_durs <- c( sleep_durs, temp$lengths[ temp$values == 1 ] )
      ## add the duration of wake bouts to the wake bout duration vector
      wake_durs <- c( wake_durs, temp$lengths[ temp$values == 0 ] )
      
      
      waking_dat <- id_dat[ id_dat$local_timestamp > sleep_per[ sleep_per$tag == tag & sleep_per$night == ( night - 1 ), ]$waking & id_dat$local_timestamp <  onset   , ]
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$prev_naps <- sum( waking_dat$sleep_bouts )
      
      
    }
    
  }
}

## for days on the later side of a noon-to-noon period with a lot of missing data, we can't have a reliable threshold for what was considered sleep in the morning, and we may be missing a lot of epochs. Therefore, remove: sleep time during the day (day defined by SPT (wake time to onset time)), sleep time during the day (day defined as prev dark period end to dark period start), sleep time during the day (defined as 0730 to 1730)

sleep_per[ paste( sleep_per$tag, sleep_per$night, sep = '_' ) %in% paste( sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ), 'tag' ], ( sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ) , 'night' ] + 1 ), sep = '_' ), c( 'prev_naps', 'prev_day_sleep', 'prev_day_sleep_lim' ) ] <- NA

## for days on the later side of a noon-to-noon period with large gaps of missing data, we can't have a reliable waking time. Therefore, remove: sleep time during the day (day defined by SPT (wake time to onset time))
sleep_per[ paste( sleep_per$tag, sleep_per$night, sep = '_' ) %in% paste( sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ), 'tag' ], ( sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ) , 'night' ] + 1 ), sep = '_' ), c( 'prev_naps' ) ] <- NA

## remove all these variable from the night, and from the days on the early side of the noon-to-noon period if the noon-to-noon period is missing a lot of data (because then we might not be able to reliably calculate the sleep VeDBA threshold, and a lot of epochs might be missing, which would skew TST and such)
sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ), c( 'onset', 'waking', 'SPT', 'sleep_eff', 'TST', 'WASO', 'wake_bouts', 'summed_VeDBA', 'night_VeDBA_corr', 'dark_TST', 'dark_sleep_eff', 'light_TST', 'light_sleep_eff', 'prev_naps', 'prev_day_sleep', 'prev_day_sleep_lim' ) ] <- NA

## remove all these variable from the night, and from the days on the early side of the noon-to-noon period (only for those depending on SPT) if the noon-to-noon period has large gaps of missing data (because then we can't reliably calculate the SPT)
sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ), c( 'onset', 'waking', 'SPT', 'sleep_eff', 'TST', 'WASO', 'wake_bouts', 'summed_VeDBA', 'light_TST', 'light_sleep_eff', 'night_VeDBA_corr', 'prev_naps' ) ] <- NA

## remove data for sleep period and sleep bouts on days when there is a lot of missing data, because we cannot reliably calculate the sleep VeDBA threshold and there may be a lot of missing epochs
full_dat[ full_dat$n_bursts < ( mins_in_day - missing_mins ), c( 'sleep_per', 'sleep_bouts' ) ] <- NA

## remove data for sleep period on days when there are large gaps of missing data, becasue we can't reliably calculate the SPT with gaps in the data
full_dat[ full_dat$max_time_diff > time_gap, 'sleep_per'  ] <- NA

## reformat sleep timestamp
sleep_per$onset <- as.POSIXct( sleep_per$onset, origin = "1970-01-01 00:00:00", tz = "UTC" )

## reformat waking timestamp
sleep_per$waking <- as.POSIXct( sleep_per$waking, origin = "1970-01-01 00:00:00", tz = "UTC" )

## make columns for just the time part of the sleep onset and waking timestamps
sleep_per$onset_time <- as_hms( sleep_per$onset )
sleep_per$waking_time <- as_hms( sleep_per$waking )


#write.csv( full_dat, 'DATA/validation_study/full_dat_2019.csv', row.names = F )

#write.csv( sleep_per, 'DATA/validation_study/sleep_per_2019.csv', row.names = F )

######## comparing focals to sleep algorithm ##########

## if the behavioral observation data was downloaded from Dryad, skip the code that immediately follows, and start running code where notified below

options( digits.secs = 6 )

full_dat <- fread( 'DATA/validation_study/full_dat_2019.csv' )

full_dat <- as.data.frame( full_dat )

full_dat$timestamp <- as.POSIXct( full_dat$timestamp, tz = 'UTC' )


## load in the original focal follow observation data #####
files <- list.files(path = "/Users/cloftus/Documents/Night_acc/sleep_validation/loopy_focal_follows_2021_09_17/", pattern = "*.csv", full.names = TRUE, recursive = TRUE)

dfs <- lapply(as.list(files), read.csv)

## make a list of new dfs which are going to have the timestamps of behaviors based on the frame timestamps when they actually occur
new_dfs <- vector( mode = 'list', length = length( dfs ) )

for( i in 1:length( dfs ) ){
  
  temp_dat <- dfs[[ i ]] 
  
  vid_name <- as.character( unique( temp_dat$Scoring ) )
  
  vid_date <- str_split_fixed( vid_name, '_', 2 )[ ,1 ]
  
  ## read in frame information
  frame_info <- read.table( paste0( "Z:\\baboon\\archive\\processed\\video\\thermal\\2019_summer\\cliff_data\\mp4_files\\viewpoint_1\\T1020\\", vid_date, '\\', vid_name, '.txt' ) )
  
  
  ## use the frame information to associate a timestamp with each frame
  names( frame_info ) <- c( "file", "frame_info" )
  
  frame_timestamps_temp <- str_split_fixed( frame_info$frame_info, "_", 2 )[, 2 ]
  
  frame_timestamps <- str_split_fixed( frame_timestamps_temp, ".tiff", 2 )[, 1 ]
  
  ## reformat to the typical timestamp format by parsing the string (using the format field of as.POSIX can get me part way there to the reformatting, but causes problems with the milliseconds)
  
  final_timestamps <- paste0( substring( frame_timestamps, 1, 4 ), "-", substring( frame_timestamps, 5, 6 ), "-", substring( frame_timestamps, 7, 8 ), " ", substring( frame_timestamps, 10, 11 ), ":", substring( frame_timestamps, 12, 13 ), ":", substring( frame_timestamps, 14, 15 ), ".", substring( frame_timestamps, 16, 21 ) )
  
  ## associate each frame number in the tracklets dataframe with the timestamp of that frame
  behav_start_timestamps <- final_timestamps[ ifelse( temp_dat$Start_Frame == length( final_timestamps ), temp_dat$Start_Frame, (temp_dat$Start_Frame + 1 ) ) ]
  
  behav_stop_timestamps <- final_timestamps[ ifelse( temp_dat$Stop_Frame == length( final_timestamps ), temp_dat$Stop_Frame, (temp_dat$Stop_Frame + 1 ) ) ]
  
  temp_dat$start_timestamp <- as.POSIXct( behav_start_timestamps, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS" ) 
  
  temp_dat$stop_timestamp <- as.POSIXct( behav_stop_timestamps, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS" ) 
  
  new_dfs[[ i ]] <- temp_dat
  
  
}

full_focal <- ldply( new_dfs, rbind )


collar_dat <- read.csv( "DATA/validation_study/tag_metadata.csv" )

names( collar_dat ) <- c( 'collar_id', 'collar_color', 'battery' )

collar_dat <- collar_dat[ !is.na( collar_dat$collar_id ), ]

collar_dat$name <- paste( collar_dat$collar_color, collar_dat$battery, sep = '_' )

fuller_focal <- merge( x = full_focal, y = collar_dat[ , c( 'name', 'collar_id' ) ], by.x = 'Subject', by.y = 'name', all.x = T, all.y = F, sort = F )

# write.csv( fuller_focal, "DATA/validation_study/fuller_focal.csv", row.names = F )

fuller_focal <- read.csv( "DATA/validation_study/fuller_focal.csv" )

fuller_focal$start_timestamp <- as.POSIXct( fuller_focal$start_timestamp, tz = 'UTC' )

fuller_focal$stop_timestamp <- as.POSIXct( fuller_focal$stop_timestamp, tz = 'UTC' )

fuller_focal$duration <- fuller_focal$stop_timestamp - fuller_focal$start_timestamp

full_trim_focal <- fuller_focal[ fuller_focal$Value %in% c( 'Active', 'Alert', 'Unalert' ), ]


pub_dat <- full_trim_focal[ full_trim_focal$collar_id %in% c( '2428', '2433', '2434', '2436', '2441', '2447', '2450' ), ] ## these tags were determined by running this whole code without trimming to just these tags, and running the validation below, and seeing which tags had sleep data and behavioral observation data that were collected synchronously

pub_dat <- pub_dat[ , c( 'collar_id', 'Subject', 'Scoring', 'start_timestamp', 'stop_timestamp', 'Value', 'duration' ) ]

write.csv( pub_dat, 'DATA/validation_study/2019_Papio_anubis_behavioral_scoring_Loftus_et_al_Dryad.csv', row.names = F )

### If downloading the behavioral data from Dryad, run start running code here

full_dat <- fread( 'DATA/validation_study/full_dat_2019.csv' ) ## must run the sleep classification algorithm at the top of this script to produce full_dat

full_dat <- as.data.frame( full_dat )

full_dat$timestamp <- as.POSIXct( full_dat$timestamp, tz = "UTC" )

pub_dat <- read.csv( 'DATA/validation_study/2019_Papio_anubis_behavioral_scoring_Loftus_et_al_Dryad.csv' )

pub_dat$start_timestamp <- as.POSIXct( pub_dat$start_timestamp, tz = 'UTC' )
pub_dat$stop_timestamp <- as.POSIXct( pub_dat$stop_timestamp, tz = 'UTC' )

accuracies <- c()
inds <- c()

new_full_dat <- full_dat

new_full_dat$behavior <- NA

timestamp_correction <- 16 ## the amount of time (in seconds) that has to be added to the video time to get to the GPS time

tags <- as.character( unique( pub_dat$collar_id ) )

for( focal_tag in tags ){
  
  trim_focal <- pub_dat[ pub_dat$collar_id == focal_tag & pub_dat$Value %in% c( 'Active', 'Alert', 'Unalert' ), ]
  
  trim_focal <- trim_focal[ order( trim_focal$start_timestamp ),]

  trim_focal$start_timestamp_corr <- trim_focal$start_timestamp + timestamp_correction
  
  trim_focal$stop_timestamp_corr <- trim_focal$stop_timestamp + timestamp_correction
  
  trim_focal <- trim_focal[ order( trim_focal$start_timestamp_corr), ]
  
  sec_focal_dat <- data.frame( tag = focal_tag, timestamp = seq( floor_date( min( trim_focal$start_timestamp_corr ), unit = "min" ), ceiling_date( max( trim_focal$stop_timestamp_corr ), unit = "min" ), by = '1 secs' ), behavior = NA )
  
  for( i in 1:nrow( trim_focal ) ){
    
    start_ind <- min( which( sec_focal_dat$timestamp > trim_focal$start_timestamp_corr[ i ] ) )
    
    stop_ind <- max( which( sec_focal_dat$timestamp < ( trim_focal$stop_timestamp_corr[ i ] - 1 ) ) )
    
    if( start_ind < stop_ind ){
      
      sec_focal_dat$behavior[ start_ind:stop_ind ] <- trim_focal$Value[ i ]
      
    }
    
  }
  
  sec_merge <- merge( x = full_dat[ full_dat$tag == focal_tag, ], y = sec_focal_dat, by.x = c( 'tag', 'timestamp' ), by.y = c( 'tag', 'timestamp' ), all.x = T, all.y = F, sort = F )
  
  sec_merge <- sec_merge[ order( sec_merge$timestamp ), ]
  
  new_full_dat[ new_full_dat$tag == focal_tag, ] <- sec_merge

}

valid_dat <- new_full_dat[ !is.na( new_full_dat$sleep_bouts ) & !is.na( new_full_dat$behavior ), ]

unique( valid_dat$tag )

confusion_matrix <- table( new_full_dat$sleep_bouts, new_full_dat$behavior )

accur <- ( confusion_matrix[ '0', 'Active' ] + confusion_matrix[ '0', 'Alert' ] + confusion_matrix[ '1', 'Unalert' ] ) / sum( confusion_matrix )

print( accur )

