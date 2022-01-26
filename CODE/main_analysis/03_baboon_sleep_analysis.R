
library( sjPlot )
library( zoo )
library( hms )
library( data.table )
library( stringr )
library( lubridate )
library( suncalc )
library( LaplacesDemon )
library( dplyr )
library( purrr )
library( multcomp )
library( nlme )
library(tidyr) 
library(lmerTest)
library( sp )
library( rgdal )
library( stats )
library(rgeos)
library( entropy )
library( reshape2 )
library( plyr )
library( rstan )
library( brms )
library( glmmTMB )
library( plotrix )
library( psych )
library(ggplot2)
library( mgcv )

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

lags <- function( sleep_per_df = sleep_per,  vars = c( "sleep_eff", "TST" ), n = 1 ){
  
  homeo <- sleep_per_df
  
  homeo$ind <- 1:nrow( homeo )
  
  tag_names <- unique( homeo$tag )
  
  old_homeo <- homeo
  
  
  for( var in vars ){
    
    for( i in 1: n ){
      
      homeo[ , ncol( homeo ) + 1 ] <- NA
      
    }
    
  }
  
  for( tag in tag_names ){
    
    id_dat <- old_homeo[ old_homeo$tag == tag, ]
    
    
    for( var in vars ){
      
      for( i in 1: n ){
        
        col <- id_dat[ , names( id_dat ) == var ]
        
        id_dat[ , ncol( id_dat ) + 1 ] <-  data.table::shift( col, i )
        
        names( id_dat )[ ncol( id_dat ) ] <- paste( 'prev', var, i, sep = '_' )
        
      }
      
    }
    
    homeo[ homeo$tag == tag, ] <- id_dat
    
  }
  
  naming <- paste( rep( paste( 'prev', vars, sep = '_' ), each = n ), 1:n, sep = '_' )
  
  names( homeo )[ ( ncol( homeo ) - length( naming ) + 1 ):ncol( homeo ) ] <- naming
  
  
  for( var in vars ){
    
    cols <- which( grepl( paste( 'prev', var, sep = '_' ), names( homeo ) ) )
    
    if( var %in% c( 'sleep_eff', 'ave_vedba', 'prev_day_ave_ved' ) ){
      
      homeo[ , ncol( homeo ) + 1 ] <- apply( data.frame( homeo[ , cols ] ) , 1, FUN = mean )
      
    }else{
      
      homeo[ , ncol( homeo ) + 1 ] <- apply( data.frame( homeo[ , cols ] ), 1, FUN = sum )
      
    }
  }
  
  naming_2 <- paste( 'prev', vars, 'cum', sep = "_" )
  
  names( homeo )[ ( ncol( homeo ) - length( naming_2 ) + 1 ):ncol( homeo ) ] <- naming_2
  
  
  return( homeo )
}



cum_lags <- function( sleep_per_df = sleep_per, vars = c( 'TST', 'prev_naps'), n = 1, inc_earliest_day = F ){
  
  print( "the previous day variable must be the second variable of the two in the vector" )
  
  homeo <- sleep_per_df
  
  homeo$ind <- 1:nrow( homeo )
  
  tag_names <- unique( homeo$tag )
  
  old_homeo <- homeo
  
  for( var in vars ){
    
    for( i in 1: n ){
      
      homeo[ , ncol( homeo ) + 1 ] <- NA
      
    }
    
  }
  
  if( inc_earliest_day == F ){
    
    homeo <- homeo[ ,  - ncol( homeo ) ]
    
  }
  
  for( tag in tag_names ){
    
    id_dat <- old_homeo[ old_homeo$tag == tag, ]
    
    
    for( var in vars ){
      
      for( i in 1: n ){
        
        col <- id_dat[ , names( id_dat ) == var ]
        
        id_dat[ , ncol( id_dat ) + 1 ] <-  shift( col, i )
        
        names( id_dat )[ ncol( id_dat ) ] <- paste( 'prev', var, i, sep = '_' )
        
      }
      
    }
    
    if( inc_earliest_day == F ){
      
      id_dat <- id_dat[ ,  - ncol( id_dat ) ]
      
    }
    
    homeo[ homeo$tag == tag, ] <- id_dat
    
  }
  
  naming <- paste( rep( paste( 'prev', vars, sep = '_' ), each = n ), 1:n, sep = '_' )
  
  
  if( inc_earliest_day == F ){
    
    naming <- naming[ - length( naming ) ]
    
  }
  
  names( homeo )[ ( ncol( homeo ) - length( naming ) + 1 ):ncol( homeo ) ] <- naming
  
  
  cols_of_int <- c( vars[ 2 ], naming )
  
  to_cum <- homeo[ , names( homeo ) %in% cols_of_int ]
  
  if( vars[ 1 ] == 'TST' | vars[ 1 ] == 'dark_TST'){
    
    homeo[ , ncol( homeo ) + 1 ] <- apply( to_cum, 1, FUN = sum )
    
  }else{
    
    homeo[ , ncol( homeo ) + 1 ] <- apply( to_cum, 1, FUN = mean )
    
  }
  
  
  names( homeo )[ ncol( homeo ) ] <- paste( 'prev', vars[ 1 ], 'cum', sep = "_" )
  
  return( homeo )
  
}



spat_disc_func <- function( n ){ ## n is the meters at which to spatially discretize
  
  ## instantiate a column to fill in with step lengths
  bdat$spat.disc.dist_1 <- NA
  
  for(id in tag_names){ ## for each individual
    
    ## print the tag name for progress indicator
    print(id)
    
    ## subset the GPS data to just this individual's data
    idDat <- bdat[bdat$id == id,]
    
    ## for each day of the individual's data...
    for(d in unique(idDat$day)){    ## for each day of the individual's data...
      ## subset the data to just that day
      dayDat <- idDat[idDat$day == d,]
      
      ## instantiate the column for step lengths (redundant)
      dayDat$spat.disc.dist_1 <- NA
      
      ## save their first GPS point
      temp.y <- dayDat$y[1] 
      temp.x <- dayDat$x[1]
      
      ## set the temporary counter to 1
      temp.count <- 1
      
      ## set the distance traveled since previous passage time to 0
      dist <- 0
      
      ## set the row counter to 0
      i <- 0
      
      ## set the total step length distance to 0
      totalDist <- 0
      
      while( i <= nrow( dayDat ) ){ ## while the row counter has not surpassed the last row of the data frame
        
        
        while( dist < n ){ ## while the individual has not yet reached the first passage radius since its previous passage
          
          
          if( i == nrow( dayDat ) ){ ## if the row counter is at the last row 
            
            break ## break the loop to continue below, where we will fill the step length for this last row with a 0
          }
          
          if( i == 0 || i == temp.count ){ ## if we are at the first GPS fix of the day or if the baboon has just made a passage
            
            i <- i + 1 ## advance the row counter
            
            dist <- sqrt( ( dayDat$y[ i ] - temp.y )**2 + ( dayDat$x[ i ] - temp.x )**2 ) ## calculate the travel distance between the current GPS and the previous location of first passage
            
          }else{ ## if we are neither at the first fix nor at the position where the baboon has just made a passage
            
            ## insert a 0 for the baboons step length (because it did not make a passage in the previous attempt)
            dayDat$spat.disc.dist_1[ i ] <- 0
            
            ## fill the latitude and xgitude at this row with the latitude where the baboon previously made a passage, because this it's "functional" location
            dayDat$y[ i ] <- dayDat$y[ temp.count ]
            
            dayDat$x[ i ] <- dayDat$x[ temp.count ]
            
            ## advance the row counter
            i <- i + 1
            
            ## calculate the distance from this new row to the baboon's functional location (the location at which it last made a passage)
            dist<-sqrt((dayDat$y[i]-temp.y)**2 + (dayDat$x[i]-temp.x)**2) 
            
          }
          
        }
        
        ## at this point we have made it out of the loop above. Either because we are at the end of the dataframe or because the distance from the previous passage has reached the necessary threshold (aka the baboon just made a new passage)
        if( dist < n ){ ## if the distance is less than the threshold (only would occur here if we are at the end of the dataframe)
          
          ## set the step length at this row to 0 and break out of the lop
          dayDat$spat.disc.dist_1[ i ] <- 0
          
          break
          
        }else{ ## if the baboon did indeed just make a passage
          
          ## save its step length as the distance between this passage and previous passage
          dayDat$spat.disc.dist_1[ i ] <- dist
          
          ## save the latitude and longitude of this location of passage (to compare subsequent fixes to, in order to find the next passage)
          temp.y <- dayDat$y[ i ]
          
          temp.x <- dayDat$x[ i ]
          
          ## save the row of this passage as temp.count
          temp.count <- i
          
          ## reset the distance from previous passage to 0
          dist <- 0
        }
        
      }
      ## put this day of data back into the individual's full dataframe
      idDat[idDat$day == d,]$spat.disc.dist_1 <- dayDat$spat.disc.dist_1
    }
    ## put this individual's data back into the complete dataframe
    bdat[bdat$id == id,]$spat.disc.dist_1 <- idDat$spat.disc.dist_1
  }
  ## return the complete dataframe, with the new discretized step length column added
  return(bdat$spat.disc.dist_1)
}


lonlat_to_utm <- function( df, lon_name = 'lon', lat_name = 'lat', crs_utm = CRS("+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # this is for Panama: "+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
)){
  library( sp )
  library( rgdal )
  crs_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  df_sp <- SpatialPointsDataFrame(coords = df[, c(lon_name, lat_name)], proj4string = crs_longlat, data = df)
  
  df_sp <- spTransform(df_sp, crs_utm)
  
  df <- as.data.frame(df_sp)
  
  names(df) [ names(df) == paste( lon_name, 1, sep = '.') ] <- 'x'
  names(df) [ names(df) == paste( lat_name, 1, sep = '.') ] <- 'y'
  
  return( df )
}




################## Read in the d1 (accelerometer burst) data ###################

## d1 is a dataframe with a row for each minute for each baboon. Each row contains the raw (or interpolated) acc burst, and several different measures calculated from bursts (like VeDBA)
d1 <- fread("DATA/main_data/full_night_and_day_data.csv")

## turn the data table into a dataframe
d1 <- as.data.frame( d1 )

## turn timestamp into POSIX element and time into character
d1$timestamp <- as.POSIXct( d1$timestamp, tz = 'UTC' )
d1$time <- as.character( d1$time )

## change times to local time here
d1$local_timestamp <- d1$timestamp + 3*60*60

## make a column for local time
d1$local_time <- str_split_fixed(d1$local_timestamp, " ", 2)[,2]

#################### Read in the GPS data #######################

## from Movebank, I downloaded the GPS data that is annotated with some environmental data. This line is reading that GPS data in
bdat <- fread( "DATA/main_data/Collective movement in wild baboons.csv", fill = T )

## turning bdat into a dataframe from a data table
bdat <- as.data.frame( bdat )

## keeping only the relevant columns
bdat <- bdat[, c('tag-local-identifier','timestamp','location-lat','location-long' ) ]

## rename the GPS columns
names(bdat) <- c('id','timestamp','lat','lon')

## make the timestamp into a POSIXct element
bdat$timestamp <- as.POSIXct(bdat$timestamp, tz = 'UTC' )

## make a column named day that reprsents the day from the start of the study period, with the first day being day 1
bdat$day <- as.numeric(as.Date(bdat$timestamp) - as.Date(min(bdat$timestamp)) + 1)

## make a column named time that will have just the time component, and not the date, of the timestamp
bdat$time <- str_split_fixed( bdat$timestamp, " ",2)[,2]

## remove rows of the dataframe that don't have a fix
bdat <- bdat[!is.na(bdat$lat),]

## find the average longitude and latitude of the study GPS points. These will be used to find the times of sunset and sunrise
ave_lon <- mean( bdat$lon )
ave_lat <- mean( bdat$lat )

## look at how many 'baboon days' there are in the data
nrow( unique( bdat[ , c( 'id', 'day' ) ] ) )


############## Make a dataframe of sunrise and sunset times ###############

## make an empty dataframe with each date of the study. For each date, we will fill in when sunset occurred, when the dark period began (the end of evening astronomical twilight), when the dark period ended the following morning (the beginning of morning astronomical twilight), and when sunrise the following morning occured
sun_dat <- data.frame( date = c( ( min( as.Date( d1$local_timestamp ) ) - 2:1 ), unique( as.Date( d1$local_timestamp ) ) ), sunset = NA, night_start = NA, night_end = NA, sunrise = NA )

## fill in sunset time and night start time (dark period start time) on each date with the getsunlighttimes function 
sun_dat[, c( 'sunset', 'night_start' ) ] <- getSunlightTimes( date = sun_dat$date, lon = ave_lon, lat = ave_lat )[, c( 'sunset', 'night' ) ]

## fill in rise time and night end time (dark period end time) on each date with those from the following date with the getsunlighttimes function. The reason we are using the following date is because we want to know when a night ended, which happens on the date following the start of that night 
sun_dat[, c( 'sunrise', 'night_end' ) ] <- getSunlightTimes( date = ( sun_dat$date + 1 ), lon = ave_lon, lat = ave_lat )[, c( 'sunrise', 'nightEnd' ) ]

## put sun data in local time
sun_dat[ , 2:5 ] <- sun_dat[ , 2:5 ] + 3*60*60

## save the date of the first night of the study (the date of the night is always the date of the evening at the beginning of that night; so the first night of the study is 2012-07-31, although the data starts on 2012-08-01, because the data on that first morning is still technically part of the data for the previous night, as a night is noon to noon)
start_date <- as.Date(min(d1$local_timestamp)- 12*60*60)

## make a column for night that matches with the night column in bdat and d1. This is the nights from the beginning of the study period, with the first night of the study period being night 1
sun_dat$night <- as.numeric( sun_dat$date - start_date ) + 1


############## Make a dataframe of moon phases ###############

## save the unique dates associated with the local timestamps in d1
dates <- c( ( min( as.Date( d1$local_timestamp ) ) - 2:1 ), unique( as.Date( d1$local_timestamp ) ) )

## using the getmoonillumination function, get information about the moon on the study dates
moon_dat <- getMoonIllumination( date = c( min( dates ) - 1, dates, max( dates ) + 1 ) )

## make a column for night that matches the night column in bdat and d1. Because the values in this dataframe give the moon information at 00:00:00 for each particular date, the date in the date column actually corresponds to one night earlier than that date
moon_dat$night <- as.numeric( moon_dat$date - start_date )


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

frag_block <- 2 ## this is the number of minutes of waking that need to be consecutive to be considered a wake bout during the night (other epochs of wake that do not meet this criterion will still be considered wake for WASO and wake_bouts, but not frag_wake_bouts)

sep_day_night <- F ## this determines whether sleep periods and non-sleep periods are separated before finding runs of inactivity to consider as sleep

## shows the time (as well as one previous time and one later time) where a minute is skipped. This shows that throughout the data, there is never a burst that occurs at 14:59:00
sort( unique( d1$time ) ) [ which( diff( as_hms( sort( unique( d1$time ) ) ) ) != as_hms( '00:01:00' ) ) + -1:1 ] 

## again confirms that every minute is represented in the data except for one (can tell this by comparing this number to the minutes_in_day variable above)
length( unique(d1$time) )

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
    
    
    
    dark_start <- sun_dat$night_start[ sun_dat$night == night ]
    
    dark_end <- sun_dat$night_end[ sun_dat$night == night ] 
    
    night_dat$dark <- as.numeric( night_dat$local_timestamp >= dark_start & night_dat$local_timestamp <= dark_end )
    
    night_dat$poss_dark_bursts <- floor( as.numeric( dark_end - dark_start, units = 'mins' ) )
    
    night_dat$n_dark_bursts <- sum( night_dat$local_timestamp >= dark_start & night_dat$local_timestamp <= dark_end )
    
    ## sort the timestamps (they are probably already sorted)
    sorted_times <- c( dark_start, sort( night_dat[ night_dat$dark == 1, 'local_timestamp'] ), dark_end )
    
    ## find the time difference in seconds between each burst
    time_diffs <- as.numeric( diff( sorted_times ), units = 'secs' )
    
    if( length( time_diffs ) != 0 ){
      ## save the maximum time difference during the night
      night_dat$max_dark_time_diff <- max( time_diffs )  ## so if this variable is NA, that means that there are no bursts taken from the nighttime period on this day for this individual
    }else{
      
      night_dat$max_dark_time_diff <- NA
    }
    
    
    
    first_light <- sun_dat$night_end[ sun_dat$night == ( night - 1 ) ]
    
    night_dat$poss_day_bursts <- floor( as.numeric( dark_start - first_light, units = 'mins' ) )
    
    day_dat <- id_dat[ id_dat$local_timestamp > first_light & id_dat$local_timestamp < dark_start, ]
    
    night_dat$n_day_bursts <- nrow( day_dat )
    
    
    ## sort the timestamps (they are probably already sorted)
    sorted_times <- c( first_light, sort( day_dat$local_timestamp ), dark_start )
    
    ## find the time difference in seconds between each burst
    time_diffs <- as.numeric( diff( sorted_times ), units = 'secs' )
    
    if( length( time_diffs ) != 0 ){
      ## save the maximum time difference during the night
      night_dat$max_day_time_diff <- max( time_diffs )  ## so if this variable is NA, that means that there are no bursts taken from the nighttime period on this day for this individual
    }else{
      
      night_dat$max_day_time_diff <- NA
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

sleep_per <- data.frame( tag = rep( unique( d1$tag ), each = length( study_nights ) ), night = rep( study_nights, times = length( tag ) ), onset = NA, waking = NA, SPT = NA, WASO = NA, TST = NA, sleep_eff = NA, wake_bouts = NA, frag_wake_bouts = NA, summed_VeDBA = NA, night_VeDBA_corr = NA, ave_vedba = NA, dark_TST = NA, light_TST = NA, dark_sleep_eff = NA, light_sleep_eff = NA, prev_naps = NA, prev_day_sleep = NA, prev_day_sleep_lim = NA, prev_day_ave_vedba = NA, poss_dark_bursts = NA, n_dark_bursts = NA, poss_day_bursts = NA, n_day_bursts = NA, max_time_diff = NA, n_bursts= NA, max_dark_time_diff = NA, max_day_time_diff = NA )


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
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$poss_dark_bursts <- unique( night_dat$poss_dark_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$n_dark_bursts <- unique( night_dat$n_dark_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$max_dark_time_diff <- unique( night_dat$max_dark_time_diff )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$poss_day_bursts <- unique( night_dat$poss_day_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$n_day_bursts <- unique( night_dat$n_day_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$max_day_time_diff <- unique( night_dat$max_day_time_diff )
    
    
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
      
      temp <- rle( SPT_dat$sleep_bouts )
      
      runs <- as.numeric( rep( temp$lengths >= frag_block, times = temp$lengths ) )
      
      frag_wake_bouts <- as.numeric( SPT_dat$sleep_bouts == 0 & runs == 1 )
      
      diffs <- diff( c( 1, frag_wake_bouts ) )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$frag_wake_bouts <- sum( diffs == 1 )
      
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
    
    
    
    
    dark_dat <- night_dat[ night_dat$dark == 1, ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$ave_vedba <- mean( dark_dat$log_vedba )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$dark_TST <- sum( dark_dat$sleep_bouts == 1 )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$dark_sleep_eff <- sum( dark_dat$sleep_bouts == 1 ) / nrow( dark_dat )
    
    light_dat <- night_dat[ night_dat$dark == 0 & night_dat$sleep_per == 1, ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$light_TST <- sum( light_dat$sleep_bouts == 1 )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$light_sleep_eff <- sum( light_dat$sleep_bouts == 1 ) / nrow( dark_dat )
    
    
    first_light <- sun_dat$night_end[ sun_dat$night == ( night - 1 ) ]
    
    last_light <- sun_dat$night_start[ sun_dat$night == night ]
    
    day_dat <- id_dat[ id_dat$local_timestamp > first_light & id_dat$local_timestamp < last_light, ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$prev_day_sleep <- sum( day_dat$sleep_bouts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$prev_day_ave_vedba <- mean( day_dat$log_vedba )
    
    day_lim_start <- as.POSIXct( paste( str_split_fixed( first_light, " ", 2 )[ , 1 ], day_lim_start_time ), tz = 'UTC' )
    
    day_lim_end <- as.POSIXct( paste( str_split_fixed( last_light, " ", 2 )[ , 1 ], day_lim_end_time ), tz = 'UTC' )
    
    day_lim <- id_dat[ id_dat$local_timestamp > day_lim_start & id_dat$local_timestamp < day_lim_end, ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$prev_day_sleep_lim <- sum( day_lim$sleep_bouts )
    
  }
}


########### Add in moon phase data ###########

sleep_per <- merge( x = sleep_per, y = moon_dat, by = 'night', all.x = T, all.y = F, sort = F )

sum( !is.na( sleep_per$SPT ) )

### check number of nights for which sleep period was calculated and inspect those for which no sleep period was calculated ###

sleep_per_nona <- sleep_per[ !is.na( sleep_per$SPT ), ]

nrow( sleep_per_nona )

left_out <- unique( d1[ , c( 'tag', 'night' ) ] )[ !paste( unique( d1[ , c( 'tag', 'night' ) ] )$tag, unique( d1[ , c( 'tag', 'night' ) ] )$night ) %in% paste( sleep_per_nona$tag, sleep_per_nona$night ), ]


for( i in 1:nrow( left_out ) ){
  
  tag_night_dat <- d1[ d1$tag == left_out$tag[ i ] & d1$night == left_out$night[ i ], ]
  
  plot( tag_night_dat$local_timestamp, tag_night_dat$log_vedba )
  
} 

nrow( left_out )

# just what an example of what the plot of log_vedba from a full night of data should look like
tag_night_dat <- d1[ d1$tag == 2426 & d1$night == 5, ]

## just as an example of what the previous plots are supposed to look like
plot( tag_night_dat$local_timestamp, tag_night_dat$log_vedba )


##### Add in temperature  ######

## first make a data frame for every hour of the study to import into Movebank

hours <- seq( floor_date( min( full_dat$timestamp ), 'hour' ), ceiling_date( max( full_dat$timestamp ), 'hour' ), by = '1 hour' )

bdat$night <- as.numeric( as.Date( bdat$timestamp + 3*60*60 - 12*60*60 ) - as.Date(min(bdat$timestamp) + 3*60*60 - 12*60*60) + 1 )

morn_lats <- bdat$lat[ bdat$time <= '03:15:00' & bdat$night <= 21 ]

morn_lons <- bdat$lon[ bdat$time <= '03:15:00' & bdat$night <= 21 ]

env_data <- data.frame( `location-lat` = median( morn_lats ), `location-long` = median( morn_lons ), timestamp = hours )

env_data$timestamp <- paste( env_data$timestamp, '.000', sep = '' )

write.csv( env_data, "DATA/main_data/env_data.csv", row.names = F )

### import the above data frame into Movebank and use their Env-DATA functionality to add the temperature data, and then load the resulting data frame back in. This resulting dataframe is also available on Dryad
env_data <- read.csv( 'DATA/main_data/env_data.csv-6150899038464587825.csv' )

names( env_data ) <- c( 'lat', 'lon', 'timestamp', 'temperature' )

env_data$timestamp <- as.POSIXct( env_data$timestamp, tz = 'UTC' )

env_data$local_timestamp <- env_data$timestamp + 3*60*60

time_shift <- env_data$local_timestamp - 12*60*60

start_date <- as.Date( min( env_data$local_timestamp ) - 12*60*60 )

env_data$night <- as.numeric( as.Date( time_shift ) - start_date + 1 )

env_data$local_time <- str_split_fixed( env_data$local_timestamp, ' ', 2 )[ , 2 ]

env_sleep_dat <- data.frame( night = unique( env_data$night ) )


for( night in unique( env_data$night ) ){
  
  night_dat <- env_data[ env_data$night == night, ]
  
  env_sleep_dat[ env_sleep_dat$night == night, 'temperature' ] <- min( night_dat$temperature )
  
}

sleep_per <- merge( x = sleep_per, y = env_sleep_dat, by = 'night', all.x = T, all.y = F, sort = F )

############# Cleaning the dataframes of data on nights with insufficient data ################

## for days on the later side of a noon-to-noon period with a lot of missing data, we can't have a reliable threshold for what was considered sleep in the morning, and we may be missing a lot of epochs. Therefore, remove: sleep time during the day (day defined by SPT (wake time to onset time)), sleep time during the day (day defined as prev dark period end to dark period start), sleep time during the day (defined as 0730 to 1730)

sleep_per[ paste( sleep_per$tag, sleep_per$night, sep = '_' ) %in% paste( sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ), 'tag' ], ( sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ) , 'night' ] + 1 ), sep = '_' ), c( 'prev_naps', 'prev_day_sleep', 'prev_day_sleep_lim' ) ] <- NA

## for days on the later side of a noon-to-noon period with large gaps of missing data, we can't have a reliable waking time. Therefore, remove: sleep time during the day (day defined by SPT (wake time to onset time))
sleep_per[ paste( sleep_per$tag, sleep_per$night, sep = '_' ) %in% paste( sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ), 'tag' ], ( sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ) , 'night' ] + 1 ), sep = '_' ), c( 'prev_naps' ) ] <- NA

## remove all these variable from the night, and from the days on the early side of the noon-to-noon period if the noon-to-noon period is missing a lot of data (because then we might not be able to reliably calculate the sleep VeDBA threshold, and a lot of epochs might be missing, which would skew TST and such)
sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ), c( 'onset', 'waking', 'SPT', 'sleep_eff', 'TST', 'WASO', 'wake_bouts', 'summed_VeDBA', 'night_VeDBA_corr', 'dark_TST', 'dark_sleep_eff', 'light_TST', 'light_sleep_eff', 'prev_naps', 'prev_day_sleep', 'prev_day_sleep_lim' ) ] <- NA

sleep_per_nona <- sleep_per[ !is.na( sleep_per$SPT ), ]

nrow( sleep_per_nona )


## remove all these variable from the night, and from the days on the early side of the noon-to-noon period (only for those depending on SPT) if the noon-to-noon period has large gaps of missing data (because then we can't reliably calculate the SPT)
sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ), c( 'onset', 'waking', 'SPT', 'sleep_eff', 'TST', 'WASO', 'wake_bouts', 'summed_VeDBA', 'light_TST', 'light_sleep_eff', 'night_VeDBA_corr', 'prev_naps' ) ] <- NA

sleep_per_nona <- sleep_per[ !is.na( sleep_per$SPT ), ]

nrow( sleep_per_nona )

## remove data for sleep period and sleep bouts on days when there is a lot of missing data, because we cannot reliably calculate the sleep VeDBA threshold and there may be a lot of missing epochs
full_dat[ full_dat$n_bursts < ( mins_in_day - missing_mins ), c( 'sleep_per', 'sleep_bouts' ) ] <- NA

## remove data for sleep period on days when there are large gaps of missing data, because we can't reliably calculate the SPT with gaps in the data
full_dat[ full_dat$max_time_diff > time_gap, 'sleep_per'  ] <- NA

write.csv( full_dat, 'DATA/main_data/full_dat.csv', row.names = F )


## reformat sleep timestamp
sleep_per$onset <- as.POSIXct( sleep_per$onset, origin = "1970-01-01 00:00:00", tz = "UTC" )

## reformat waking timestamp
sleep_per$waking <- as.POSIXct( sleep_per$waking, origin = "1970-01-01 00:00:00", tz = "UTC" )

## make columns for just the time part of the sleep onset and waking timestamps
sleep_per$onset_time <- as_hms( sleep_per$onset )
sleep_per$waking_time <- as_hms( sleep_per$waking )


############ Basic plots ##########

sleep_trim <- sleep_per[ !is.na( sleep_per$SPT ), ]

### Fig 1 ###

sun_trim <- sun_dat[ sun_dat$night %in% unique( full_dat$night ), ]

ave_sunset <- as_hms( mean( as.numeric( as_hms( sun_trim$sunset ) ) ) )
ave_sunset <- ts_func( ave_sunset )

ave_sunrise <- as_hms( mean( as.numeric( as_hms( sun_trim$sunrise ) ) ) )
ave_sunrise <- ts_func( ave_sunrise )

ave_night_start <- as_hms( mean( as.numeric( as_hms( sun_trim$night_start ) ) ) )
ave_night_start <- ts_func( ave_night_start )

ave_night_end <- as_hms( mean( as.numeric( as_hms( sun_trim$night_end ) ) ) )
ave_night_end <- ts_func( ave_night_end )

par(mar=c(6,6,4,4))
layout(matrix(1:3, ncol = 1), widths = 1, heights = c(2,2,2.5), respect = FALSE)


## Fig 1A

par(mar = c(0, 5, 0.5, 1 ) )

sleep_per_func( tag = 2430, night = 3, x_axis = F, plot_waso = T, title = F, waso_block = waso_block, lwd = 0.5 )

# legend( x = 67000, y = 7.6, col = c( 'orange', transp( 'blue', 0.25 ) ), lty = c( 2, 0 ), pch = c( NA, 15 ), pt.cex = 2, legend = c( 'Sleep period time window', 'Wake after sleep onset'), bty = 'n', cex = 0.88 )

#legend( x = 67000, y = 6.9, col = c( 'orange', 'blue', transp( 'blue', 0.25 ) ), lty = c( 2, 0, 0 ), pch = c( NA, NA, 15 ), pt.cex = 2, legend = c( 'Sleep period time window', '', 'Wake after sleep onset'), bty = 'n', cex = 0.88 ) ##y was 7.6


#legend( 'topright', col = c( 'orange', transp( 'blue', 0.25 ) ), lty = c( 02, 0 ), pch = c( NA, 15 ), pt.cex = 2, legend = c( '', ''), bty = 'n' )


polygon( x = c( ave_sunset, ave_sunrise, ave_sunrise, ave_sunset ), y = c( -10, -10, 10, 10), col = transp('grey', .25), border = NA )

polygon( x = c( ave_night_start, ave_night_end, ave_night_end, ave_night_start ), y = c( -10, -10, 10, 10), col = transp('grey', .25), border = NA )


## Fig 1A


par(mar = c(0, 5, 0, 1))

ave_ved <- aggregate( d1$log_vedba, by = list( d1$local_time ), FUN = mean, na.rm = T ) 
names( ave_ved ) <- c( 'local_time', 'ave_VeDBA')

ave_ved <- ave_ved[ order( ts_func( ave_ved$local_time ) ), ]


se_ved <- aggregate( d1$log_vedba, by = list( d1$local_time ), FUN = std.error, na.rm = T ) 
names( se_ved ) <- c( 'local_time', 'se_VeDBA')

se_ved <- se_ved[ order( ts_func( se_ved$local_time ) ), ]

se_ved$se_VeDBA[ is.na( se_ved$se_VeDBA) ] <- 0

x_vals <- ts_func( ave_ved$local_time )

plot( x_vals , ave_ved$ave_VeDBA, type = 'l', xaxt = 'n', xlab = 'Time', ylab = '', las = 1, lwd = 0.5  )

title( ylab = "log VeDBA", line = 3.9 )


low_se <- ave_ved$ave_VeDBA - se_ved$se_VeDBA
high_se <- ave_ved$ave_VeDBA + se_ved$se_VeDBA


polygon( x = c( ( x_vals ) ,rev( x_vals ), x_vals[ 1 ] ), y = c( low_se, rev( high_se ), low_se[ 1 ] ), col=transp( "red" , 0.5 ), border=NA ) # fills area between the curves

polygon( x = c( ave_sunset, ave_sunrise, ave_sunrise, ave_sunset ), y = c( -10, -10, 10, 10), col = transp('grey', .25), border = NA )

polygon( x = c( ave_night_start, ave_night_end, ave_night_end, ave_night_start ), y = c( -10, -10, 10, 10), col = transp('grey', .25), border = NA )


## Fig 1C

par(mar = c(4.1, 5, 0, 1))

## make a density plot of the onset time and waking time
onset_dens <- density( ts_func( sleep_trim$onset_time ) )

wake_dens <- density(  ts_func( sleep_trim$waking_time ) )

waso_times <- full_dat[ !is.na( full_dat$sleep_per ) & full_dat$sleep_per == 1 & full_dat$sleep_bouts == 0, 'local_time' ]

waso_dens <- density( ts_func( waso_times ) )

x_range <- c( 0, 24*60*60 )

y_range <- c( min( onset_dens$y, wake_dens$y ), max( onset_dens$y, wake_dens$y ) )

plot( onset_dens$x, onset_dens$y, type = 'l', ylab = '', xlab = 'Time', xaxt = 'n', yaxt = "n", xlim = x_range, ylim = y_range, lty = 2, las = 1 )

#axis( 2, at = seq( 0, 0.0004, by = 0.0002 ), las = 1 )

title( ylab = "Probability Density", line = 3.9 )

lines( wake_dens$x, wake_dens$y, type = 'l', ylab = 'Density', xlab = 'Waking Time', xaxt = 'n', lty = 3 )


#legend( x = 'topright', cex = 0.9, col = c( 'red', 'blue' ), lty = c( 1, 1 ), legend = c( 'Sleep onset time', 'Waking time' ), bty = 'n' )

polygon( x = c( ave_sunset, ave_sunrise, ave_sunrise, ave_sunset ), y = c( -10, -10, 10, 10), col = transp('grey', .25), border = NA )

polygon( x = c( ave_night_start, ave_night_end, ave_night_end, ave_night_start ), y = c( -10, -10, 10, 10), col = transp('grey', .25), border = NA )

axis( 1, at = seq( 0, 60*24*60, 3*60*60), labels = c( as_hms( seq( 12*60*60, 60*23*60, 3*60*60) ), as_hms( seq( 0, 60*12*60, 3*60*60) ) ) ) 

axis( 1, at = seq( 0, 60*24*60, 1*60*60), labels = F ) 

axis( 2, at = c( 0.0, 0.0002, 0.0004 ), labels = c( '0.0', '0.0002', '0.0004' ), las = 2 )

#### End of Fig


dev.off()


############## Sleep efficiency and TST during sleep period as a function of sleep during previous waking period  ( ALSO: sleep during the next waking period as a function of sleep quality at night)  ############

nap_dat <- sleep_per

nap_dat$next_naps <- NA
nap_dat$next_day_sleep <- NA
nap_dat$next_day_sleep_lim <- NA
nap_dat$next_day_ave_vedba <- NA

tag_names <- unique( nap_dat$tag )

for( tag in tag_names ){
  
  id_dat <- nap_dat[ nap_dat$tag == tag, ]
  
  id_dat$next_naps <- shift( id_dat$prev_naps, -1 )
  
  id_dat$next_day_sleep <- shift( id_dat$prev_day_sleep, -1 )
  
  id_dat$next_day_sleep_lim <- shift( id_dat$prev_day_sleep_lim, -1 )
  
  id_dat$next_day_ave_vedba <- shift( id_dat$prev_day_ave_vedba, -1 )
  
  nap_dat[ nap_dat$tag == tag, ] <- id_dat
  
}

write.csv( nap_dat, "DATA/main_data/nap_dat.csv", row.names = F )

######## Setting up full_dat for later epoch by epoch analysis(sleep efficiency as a function of time from the start of the sleep period) ##########


full_dat_sleep <- full_dat[ full_dat$sleep_per == 1 & !is.na( full_dat$sleep_per ) , ]

full_dat_sleep$time_from_sp_start <- NA

full_dat_sleep$prop_from_sp_start <- NA

tag_names <- unique( full_dat_sleep$tag )

for( tag in tag_names ){
  
  id_dat <- full_dat_sleep[ full_dat_sleep$tag == tag, ]
  
  nights <- unique( id_dat$night )
  
  for( night in nights ){
    
    night_dat <- id_dat[ id_dat$night == night, ]
    
    sleep_per_start <- min( night_dat$local_timestamp )
    
    sleep_per_end <- max( night_dat$local_timestamp )
    
    sleep_per_dur <- as.numeric( sleep_per_end - sleep_per_start, units = 'mins' )
    
    night_dat$time_from_sp_start <- as.numeric( night_dat$local_timestamp - sleep_per_start, units = 'mins' )
    
    night_dat$prop_from_sp_start <- night_dat$time_from_sp_start / sleep_per_dur
    
    id_dat[ id_dat$night == night, ] <- night_dat
    
  }
  
  full_dat_sleep[ full_dat_sleep$tag == tag, ] <- id_dat
  
}


write.csv( full_dat_sleep, "DATA/main_data/full_dat_sleep.csv", row.names = F )


######################### Effect of physical exertion on sleep quality #########################

## add UTM coordinates to the GPS data
bdat <- lonlat_to_utm( bdat )

## save the names of the individuals in the dataframe
tag_names <- unique( bdat$id )

####### Calculate a spatially discretized travel distance ##########

bdat$spat.disc.dist_5 <- spat_disc_func( 5 )

## create a duplicate of the full dataframe, from which we will remove insufficient data
t_bdat <- bdat

## remove any data on tag 2457 before night 18, because his GPS error was high before this night
t_bdat <- t_bdat[-which(t_bdat$id == '2457' & t_bdat$day < 18),]

## for each individual 
for(i in unique( t_bdat$id ) ){
  
  ## for each day
  for(j in unique( t_bdat$day ) ){
    
    ## subset the data to this individual on this day
    newDF <- t_bdat[ t_bdat$id == i & t_bdat$day == j, ]
    
    ## if this individual has data on this day
    if(nrow(newDF) != 0){
      
      ## if the individual's data cuts off before 14:00:00 UTC (i.e. 5 PM local time) or if it starts after 04:30:00 UTC (i.e. 7:30 AM local time)
      if(!is.na(max(newDF$time)) & max(newDF$time) < "14:00:00" || (!is.na(max(newDF$time)) & min(newDF$time) > "04:30:00" )){  ### Change the time here to be more or less restrictive on when a collar has to continue working for to be included in that day
        
        ## remove this individual's data on this day from the full dataframe
        t_bdat<-t_bdat[ -which( t_bdat$id ==i & t_bdat$day == j ), ]
        
      }
    }
  }
}

nrow( unique( t_bdat[ , c( 'id', 'day' ) ] ) )

## sum the step lengths for each day to calculate a daily travel distance
df1 <- aggregate( t_bdat[ , 'spat.disc.dist_5' ], by = list( t_bdat$id, t_bdat$day ), FUN = sum, na.rm = T )

## record the number of GPS fixes for each individual for each day
df2 <- aggregate(t_bdat$timestamp,by = list(t_bdat$id, t_bdat$day), FUN = length)

## record the time of each individual's first and last GPS fix on each day
df3 <- aggregate( t_bdat$timestamp, by = list( t_bdat$id, t_bdat$day ), FUN = min )
df4 <- aggregate( t_bdat$timestamp, by = list( t_bdat$id, t_bdat$day ), FUN = max )

## combine this data into one dataframe
dat <- cbind( df1, df2[ , 3 ], df3[ , 3 ], df4[ , 3 ] )

## rename the columns of the dataframe
names(dat) <- c( 'id', 'day', 'spat.disc_5', 'n', 'start_time', 'end_time' )

## put the data in order by day
dat <- dat[ order( dat$day ), ]

## write the travel distance data to a csv
write.csv( dat, "DATA/main_data/trav_dists_trim.csv", row.names = F )

### read travel distance data back into R
dat <- read.csv( "DATA/main_data/trav_dists_trim.csv" )

## order the data by day
dat <- dat[ order( dat$day ), ]

## make a column for night. This will correspond to the night following this day, which we will use to merge this dataframe with the sleep period dataframe
dat$night <- dat$day + 1

## merge the sleep period data frame with the travel distance dataframe
full_sleep <- merge( x = sleep_per, y = dat[ , c('id', 'night', 'spat.disc_5') ], by.x = c( 'tag', 'night' ), by.y = c( 'id', 'night' ), all.x = T, all.y = T, drop = F, sort = F )

full_sleep <- full_sleep[ order( full_sleep$night ), ]
full_sleep <- full_sleep[ order( full_sleep$tag ), ]

write.csv( full_sleep, "DATA/main_data/full_sleep.csv", row.names = F )


full_sleep <- read.csv( "DATA/main_data/full_sleep.csv" )


############ Effect of sleep location on sleep efficiency and TST ####################

## save a shapefile of the morning GPS fixes so we can visualize them in ArcGIS. This well help us trace the sleep trees in ArcGIS by showing us which trees the baboons actually slept in

mornings <- bdat[ bdat$time < '03:15:00', ]

mornings$ind <- 1:nrow( mornings )

mornings <- mornings[ mornings$ind %% 60 == 0, ]

crs_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

sp_mornings <- SpatialPointsDataFrame(coords = mornings[ , c('lon','lat')], proj4string = crs_longlat, data = mornings ) 

writeOGR( sp_mornings, 'DATA/main_data/morning_GPS', driver="ESRI Shapefile", layer='morning_GPS')

## load the sleep trees that were traced in ArcGIS
sleep_trees <- readOGR(dsn = "DATA/main_data", layer = "sleep_trees")

crs_utm <- CRS("+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" )

sleep_trees_utm <- spTransform( sleep_trees, crs_utm )

tree_dat <- as.data.frame( unique( bdat[ , c( 'id', 'day' ) ] ) )

names( tree_dat ) <- c( 'tag', 'day' )

tree_dat$tree <- NA
tree_dat$x <- NA
tree_dat$y <- NA

tag_names <- as.character( unique( tree_dat$tag ) )

morning_time_thresh <- '03:15:00'

cl_bdat <-  bdat

baboon_mornings <- 0

baboon_mornings_1 <- 0

near_not_in <- 0

total_assigned <- 0

not_in <- 0

total_night_with_dat <- 0

for( tag in tag_names ){
  
  id_dat <- cl_bdat[ cl_bdat$id == tag, ]
  
  days <- as.character( unique( id_dat$day ) )
  
  for( day in days ){
    
    day_dat <- id_dat[ id_dat$day == day, ]
    
    morning_fixes <- day_dat[ day_dat$time < morning_time_thresh, ]
    
    baboon_mornings_1 <- baboon_mornings_1 + 1
    
    if( !nrow( morning_fixes ) < 10 ){
      
      baboon_mornings <- baboon_mornings + 1
      
      total_night_with_dat <- total_night_with_dat + 1
      
      first_fixes <- morning_fixes[ 1:10, ]
      
      x <- median( first_fixes$x )
      y <- median( first_fixes$y )
      
      
      tree_dat[ tree_dat$tag == tag & tree_dat$day == day, c( 'x', 'y' ) ] <- c( x, y ) 
      
      
      first_fix <- data.frame( tag = unique( day_dat$id ),  x = x, y = y )
      
      sp_first_fix <- SpatialPointsDataFrame( coords = first_fix[, c( 'x', 'y' ) ], proj4string = crs_utm, data = first_fix ) 
      
      if( min( gDistance( sp_first_fix, sleep_trees_utm, byid = T ) ) < 10 ){
        
        total_assigned <- total_assigned + 1
        
        if( min( gDistance( sp_first_fix, sleep_trees_utm, byid = T ) ) != 0 ){
          
          near_not_in <- near_not_in + 1
          
        }
        
        min_ind <- which.min( gDistance( sp_first_fix, sleep_trees_utm, byid = T ) )
        
        tree_dat[ tree_dat$tag == tag & tree_dat$day == day, 'tree' ] <- min_ind
        
      }else{
        
        print( paste( tag, day, 'not within 10m of tree' ) )
        
        not_in <- not_in + 1 
      }
      
      
    }else{
      
      print( paste( tag, day, 'does not have 10 morning fixes' ) )
      
    }
    
    
  }
  
}

baboon_mornings

baboon_mornings_1

total_night_with_dat

not_in

not_in / total_night_with_dat

total_assigned

near_not_in

near_not_in / total_assigned

# the baboons slept at a different sleep site on day 22, 23, and 24
tree_dat$site <- ifelse( tree_dat$day %in% c( 22, 23, 24 ), 2, 1 )


## This seems like it is a day off but it is actually right, because really what we want is the morning locations to match up to the data from the night (and day) before
tree_dat$night <- tree_dat$day

tree_dat <- tree_dat[ , names( tree_dat ) != "day" ]


fuller_sleep <- merge( x = full_sleep, y = tree_dat, by.x = c( 'tag', 'night' ), by.y = c( 'tag', 'night' ), sort = F )

write.csv( fuller_sleep, "DATA/main_data/fuller_sleep.csv", row.names = F )


fuller_sleep <- read.csv( "DATA/main_data/fuller_sleep.csv" )


######################## Site fidelity ############################


# remove individuals with 3 days of data or less, because we can't say anything about consistency in sleep tree choice when we only have data from a few nights
tag_names <- as.character( unique( fuller_sleep$tag ) )

rem_vec <- c()

for( tag in tag_names ){
  
  id_dat <- fuller_sleep[ fuller_sleep$tag == tag, ]
  
  nights <- unique( id_dat$night[ !is.na( id_dat$tree ) ] )
  
  if( length( nights ) < 4 ){
    
    rem_vec <- c( rem_vec, tag )
    
  }
}

## see when many of the tags drop out
aggregate( bdat$id, by = list( bdat$day ), FUN = function( x ) length( unique( x ) ) )

## trim to the first 15 days, otherwise there aren't enough other baboons with which to swap identities (and this is going to inflate the number of randomizations that end up exactly the same as the real data, which is going to reduce the power of the randomization test)
short_dat <- fuller_sleep[ fuller_sleep$night < 15 & !is.na( fuller_sleep$tree ), ]


set.seed( 160 )

## set the number of randomizations
n_rands <- 1000

## save the names of the individuals
tag_names <- as.character( unique( short_dat$tag ) )

## save the nights on which there is data
nights <- as.character( unique( short_dat$night ) )

## create an empty dataframe. The first column will be the name of each individual. The second column will be the empirical value of entropy of their tree choice. Then there will be a column for the entropy of their tree choice in each randomization. 
ent_dat <- as.data.frame( matrix( nrow = length( tag_names ), ncol = 2 + n_rands ) )

## assign the names of the columns
names( ent_dat ) <- c( "tag", "emp", paste0( "rand_", 1:n_rands ) )

## fill the first column with the tag names
ent_dat$tag <- tag_names

## initialize an empty matrix that we will fill out with when we run the randomizations

sleep_tree_arr <- array( 0, dim = c( length( tag_names ), length( unique( short_dat$tree ) ), n_rands + 1 ), dimnames = list( tag_names, sort( unique( short_dat$tree ) ), c( 'emp', paste0( "rand_", 1:n_rands ) ) ) )


## for each randomization
for( n in 1:n_rands ){
  
  ## print randomization number for progress monitoring
  if( n %% 100 == 0 ){
    print( n )
  } 
  
  ## if on the first run
  if( n == 1 ){
  
    # calculate the empirical entropy for each individual
    
    ## for each individual
    for( tag in tag_names ){
      
      ## subset the datat to just this individual's data
      id_dat <- short_dat[ short_dat$tag == tag, ]
      
      ## create a table of counts of nights that the baboon slept in each tree
      counts <- table( id_dat$tree )
      
      ## find the row of the sleep tree array associated with this individual
      row_ind <- which( rownames( sleep_tree_arr ) == tag )
      
      ## find the columns of the sleep tree array associated with the trees this individual slept in
      col_inds <- match( names( counts ), colnames( sleep_tree_arr ) )
      
      ## record the number of times the individual slept in each of the sleep trees (empirical) in the sleep tree array
      sleep_tree_arr[ row_ind, col_inds, 1 ] <- counts
      
      ## record the entropy of this sleep tree choice in the entropy dataframe
      ent_dat[ ent_dat$tag == tag, 'emp' ] <- entropy( counts )
      
    }
  }
  
  
  ## create a duplicate of short_dat which we will manipulate for the randomization
  rand_dat <- short_dat
  
  ## for each night
  for( night in nights ){
    
    ## subset the data to just the data on this night
    night_dat <- rand_dat[ rand_dat$night == night, ]
    
    ## shuffle the identities associated with each row of the dataframe
    night_dat$tag <- sample( night_dat$tag )
    
    ## reinsert this night of data back into the dataframe
    rand_dat[ rand_dat$night == night, ] <- night_dat
    
  }
  
  ## for each individual
  for( tag in tag_names ){
    
    ## subset the data to just this individual's data
    id_dat <- rand_dat[ rand_dat$tag == tag, ]
    
    ## create a table of counts of nights that the baboon slept in each tree
    counts <- table( id_dat$tree )
    
    ## find the row of the sleep tree array associated with this individual
    row_ind <- which( rownames( sleep_tree_arr ) == tag )
    
    ## find the columns of the sleep tree array associated with the trees this individual slept in
    col_inds <- match( names( counts ), colnames( sleep_tree_arr ) )
    
    ## record the number of times the individual slept in each of the sleep trees (randomized) in the sleep tree array
    sleep_tree_arr[ row_ind, col_inds, n + 1 ] <- counts
    
    ## record the entropy of this randomized "sleep tree choice" in the entropy dataframe
    ent_dat[ ent_dat$tag == tag, n + 2 ] <- entropy( counts )
    
  }
  
}

sleep_tree_arr <- sleep_tree_arr[ !rownames( sleep_tree_arr ) %in% rem_vec, , ]

ent_dat <- ent_dat[ !ent_dat$tag %in% rem_vec, ]

# show the empirical and random distributions for the popultaion
## assign the empirical values of entropies to a vector
emp_dist <- ent_dat$emp

## assign the values of entropies from randomizations to a vector
rand_dist <- unlist( ent_dat[ , 3:ncol( ent_dat )] )

## create a density object of empirical entropies
emp_dens <- density( emp_dist )

## create a density object of randomized entropies
rand_dens <- density( rand_dist )

## find the min and max of the x and y values of empirical and randomized values of entropies. These will be the limits of our plot
xs <- c( min( emp_dens$x, rand_dens$x ), max( emp_dens$x, rand_dens$x ) )
ys <- c( min( emp_dens$y, rand_dens$y ), max( emp_dens$y, rand_dens$y ) )

## plot the distribution of empirical entropies
plot( emp_dens$x, emp_dens$y, xlim = xs, ylim = ys, type = 'l', col = 'red', ylab = 'Probability Density', xlab = 'Entropy' )

plot( 1, type = 'n', xlim = xs, ylim = ys, ylab = 'Probability Density', xlab = 'Entropy',  col = 'red' )
axis(1, at = seq(-0.5, 2.5, by = 0.5), labels = rep( "", 7 ) )
axis(2, at = seq(0, 2, by = 0.5), labels = rep("", 5))

lines( emp_dens$x, emp_dens$y, col = 'red' )
## plot the distribution of randomized entropies
lines( rand_dens$x, rand_dens$y )

## add a legend
legend( 'topleft', cex = 1, legend = c( 'empirical data', "", 'permuted data' ), col = c( 'red', "white", 'black' ), lty = 1, bty = 'n' , text.col = 'black' )

## ks-test shows that the empirical distribution and random distributions are significantly different
km_mod <- ks.test( emp_dist, rand_dist, alternative = 'greater' )

km_mod

####################### Testing the effect of sleeping in preferred tree on sleep quality #######################

## extract the random sleep tree choices into their own array
rand_arr <- sleep_tree_arr[ , , -1 ]

## take the average of each randomization, so we understand, on average, how much each individual is expected to use each tree
rand_ave <- apply( rand_arr, MARGIN = c( 1, 2 ), FUN = mean )

## extract the empirical sleep tree choices into their own array
emp_arr <- sleep_tree_arr[ , , 1 ]

## subtract the "expected use" from the "observed use" of trees and divide it by the number of observations (i.e. the number of nights they have data on) for each individual. This gives an idea of how much an individual "prefers" each tree, corrected for the amount of data so that individuals with many nights and individuals with fewer nights have the same capacity to "prefer" a tree 
pref_arr <- ( emp_arr - rand_ave ) / rowSums( emp_arr )


## make the preference array into a dataframe with a column for individual, a column for tree, and a column for the individual's preference for that tree
tree_pref <- adply( pref_arr, c( 1, 2 ) )

## name the columns of the dataframe
names( tree_pref ) <- c( 'tag', 'tree', 'pref_score' )

## merge preference data into fuller_sleep
fullest_sleep <- merge( x = fuller_sleep, y = tree_pref, by.x = c( 'tag', 'tree' ), by.y = c( 'tag', 'tree' ), all.x = T, all.y = F, sort = F )

fullest_sleep <- fullest_sleep[ order( fullest_sleep$night ), ]
fullest_sleep <- fullest_sleep[ order( fullest_sleep$tag ), ]

write.csv( fullest_sleep, "DATA/main_data/fullest_sleep.csv", row.names = F )

fullest_sleep <- read.csv( "DATA/main_data/fullest_sleep.csv" )

##################### Adding age and sex into the dataframe ######################

final_sleep <- addAgeSex( fullest_sleep, fullest_sleep$tag )

final_sleep$sex <- str_split_fixed( final_sleep$ageSex, ' ', 2 )[ , 1 ]

final_sleep$age <- str_split_fixed( final_sleep$ageSex, ' ', 2 )[ , 2 ]


write.csv( final_sleep, "DATA/main_data/final_sleep.csv", row.names = F )

final_sleep <- read.csv( "DATA/main_data/final_sleep.csv" )

################## Does the number of other individuals in the tree influence sleep? #####################

nights <- unique( final_sleep$night )

trees <- sort( unique( final_sleep$tree ) ) ## sort also drops the NAs

social_tree_dat <- data.frame( night = rep( nights, each = length( trees ) ), tree = rep( trees, times = length( nights ) ), total_num = NA )

for( night in nights ){
  
  nigh_dat <- final_sleep[ final_sleep$night == night, ]
  
  for( tree in trees ){
    
    tre_dat <- nigh_dat[ nigh_dat$tree == tree & !is.na( nigh_dat$tree ), ]
    
    if( nrow( tre_dat ) != 0 ){
      
      social_tree_dat[ social_tree_dat$night == night & social_tree_dat$tree == tree, 'total_num' ] <- nrow( tre_dat )
      
    }
  }
}

finaler_sleep <- merge( x = final_sleep, y = social_tree_dat, by = c( 'tree', 'night' ), all.x = T, all.y = F, sort = F )

num_tracked <- aggregate( finaler_sleep$tree,  by = list( finaler_sleep$night ), FUN = function( x ) sum( !is.na( x ) ) )

names( num_tracked ) <- c( 'night', 'n_tracked' )

finalest_sleep <- merge( x = finaler_sleep, y = num_tracked, by = 'night', all.x = T, all.y = F, sort = F )

finalest_sleep$total_num_corr <- finalest_sleep$total_num / finalest_sleep$n_tracked


### adding in daytime VeDBA data ###

ved_dat <- read.csv( "DATA/main_data/vedba_mean_2012.csv" )

ved_dat$timestamp <- as.POSIXct( ved_dat$timestamp, format = "%d/%m/%Y %H:%M:%S", tz = "UTC" )

ved_dat$local_timestamp <- ved_dat$timestamp + 3*60*60

ved_dat$local_time <- str_split_fixed( ved_dat$local_timestamp, " ", 2 )[ , 2 ]

ved_dat <- ved_dat[ ved_dat$local_time > "06:00:00" & ved_dat$local_time < "18:00:00", ]

## save the date of the first night of the study (the date of the night is always the date of the evening at the beginning of that night; so the first night of the study is 2012-07-31, although the data starts on 2012-08-01, because the data on that first morning is still technically part of the data for the previous night, as a night is noon to noon)
start_date <- as.Date( min( d1$local_timestamp )- 12*60*60 )

## make a column for night that matches with the night column in bdat and d1. This is the nights from the beginning of the study period, with the first night of the study period being night 1
ved_dat$night <- as.numeric( as.Date( ved_dat$timestamp ) - start_date ) + 1

head( ved_dat  )

agg_ved <- aggregate( ved_dat$vedba, by = list( ved_dat$tag, ved_dat$night ), FUN = sum, na.rm = T )

names( agg_ved ) <- c( "tag", "night", "vedba" )

finalest_sleep <- merge( x = finalest_sleep, y = agg_ved, by = c( "tag", "night" ), all.x = T, all.y = F, sort = T )

write.csv( finalest_sleep, "DATA/main_data/finalest_sleep.csv", row.names = F )

finalest_sleep <- read.csv( "DATA/main_data/finalest_sleep.csv" )

################# CHECKPOINT 1 ###########################

#save.image( "DATA/main_data/sleep_analysis_pub_code_pre_mods.RData" )

load( "DATA/main_data/sleep_analysis_pub_code_pre_mods.RData" )

finalest_sleep <- finalest_sleep[ order( finalest_sleep$night ), ]
finalest_sleep <- finalest_sleep[ order( finalest_sleep$tag ), ]

## sleep descriptors

mean( finalest_sleep$SPT, na.rm = T ) / 60

sd( finalest_sleep$SPT/60, na.rm = T ) / sqrt( sum( !is.na( finalest_sleep$SPT/60 ) ) )

mean( finalest_sleep$TST, na.rm = T ) / 60

sd( finalest_sleep$TST/60, na.rm = T ) / sqrt( sum( !is.na( finalest_sleep$TST/60 ) ) )

mean( finalest_sleep$sleep_eff, na.rm = T )

sd( finalest_sleep$sleep_eff, na.rm = T ) / sqrt( sum( !is.na( finalest_sleep$sleep_eff ) ) )

finalest_sleep$fragG2 <- finalest_sleep$frag_wake_bouts / ( finalest_sleep$TST / 60 ) 

mean( finalest_sleep$fragG2, na.rm = T )

sd( finalest_sleep$fragG2, na.rm = T ) / sqrt( sum( !is.na( finalest_sleep$fragG2 ) ) )

des_dat <- merge( x = finalest_sleep, y = sun_dat, by = 'night', all.x = T, all.y = F, sort = F )

des_dat$before_dusk <- as.numeric( des_dat$night_start - as.POSIXct( des_dat$onset, tz = 'UTC' ), units = 'mins' )

mean( des_dat$before_dusk, na.rm = T )

sd( des_dat$before_dusk, na.rm = T ) / sqrt( sum( !is.na( des_dat$before_dusk ) ) )


des_dat$after_dawn <- as.numeric(  as.POSIXct( des_dat$waking, tz = 'UTC' ) - des_dat$night_end, units = 'mins' )

mean( des_dat$after_dawn, na.rm = T )

sd( des_dat$after_dawn, na.rm = T ) / sqrt( sum( !is.na( des_dat$after_dawn ) ) )


## reformat sleep timestamp
finalest_sleep$onset <- as.POSIXct( finalest_sleep$onset, origin = "1970-01-01 00:00:00", tz = "UTC" )

## reformat waking timestamp
finalest_sleep$waking <- as.POSIXct( finalest_sleep$waking, origin = "1970-01-01 00:00:00", tz = "UTC" )

finalest_sleep$onset_time <- as_hms( finalest_sleep$onset )
finalest_sleep$waking_time <- as_hms( finalest_sleep$waking )

finalest_sleep$onset_time_num <- as.numeric( finalest_sleep$onset_time )
finalest_sleep$waking_time_num <- as.numeric( finalest_sleep$waking_time )

## correlations between sleep metrics
corr.test( finalest_sleep[ , c( 'TST', 'onset_time_num', 'waking_time_num', 'SPT', 'sleep_eff', 'fragG2' ) ] )


######### preparing variables for model  ########

finalest_sleep$SPT <- finalest_sleep$SPT / 60

total_mod_dat_temp <- lags( sleep_per_df = finalest_sleep, vars = c( 'TST', 'fragG2', 'sleep_eff' ), n = 1 )

ave_TSTs <- aggregate( total_mod_dat_temp[ , c( 'TST', 'fragG2' ) ], by = list( total_mod_dat_temp$tag ), FUN = mean, na.rm = T )

names( ave_TSTs ) <- c( 'tag', 'ave_TST', 'ave_fragG2' )

total_mod_dat <- merge( x = total_mod_dat_temp, y = ave_TSTs, by = 'tag', all.x = T, all.y = F, sort = F )

total_mod_dat$TST <- total_mod_dat$TST / 60

total_mod_dat$spat.disc_5 <- total_mod_dat$spat.disc_5 / 1000

total_mod_dat$prev_TST_diff <- total_mod_dat$prev_TST_cum - total_mod_dat$ave_TST

total_mod_dat$prev_TST_diff_std <- normalize_func( total_mod_dat$prev_TST_diff )

total_mod_dat$prev_fragG2_diff <- total_mod_dat$prev_fragG2_cum - total_mod_dat$ave_fragG2 

total_mod_dat$prev_fragG2_diff_std <- normalize_func( total_mod_dat$prev_fragG2_diff )

total_mod_dat$fragG2_std <- normalize_func( total_mod_dat$fragG2 )

total_mod_dat$TST_std <- normalize_func( total_mod_dat$TST )

total_mod_dat$spat.disc_5_std <- normalize_func( total_mod_dat$spat.disc_5 )

total_mod_dat$pref_score_std <- normalize_func( total_mod_dat$pref_score )

total_mod_dat$prev_day_sleep_lim_std <- normalize_func( total_mod_dat$prev_day_sleep_lim )

total_mod_dat$fraction_std <- normalize_func( total_mod_dat$fraction )

total_mod_dat$total_num_corr_std <- normalize_func( total_mod_dat$total_num_corr )

total_mod_dat$temperature_std <- normalize_func( total_mod_dat$temperature )

total_mod_dat$vedba_std <- normalize_func( total_mod_dat$vedba )

total_mod_dat$night_num <- as.numeric( total_mod_dat$night )
total_mod_dat$night <- as.factor( total_mod_dat$night )
total_mod_dat$tree <- as.factor( total_mod_dat$tree )
total_mod_dat$age <- as.factor( total_mod_dat$age )
total_mod_dat$sex <- as.factor( total_mod_dat$sex )
total_mod_dat$tag <- as.factor( total_mod_dat$tag)
total_mod_dat$site <- as.factor( total_mod_dat$site )
total_mod_dat$onset_time_num <- as.numeric( as_hms( as.character( total_mod_dat$onset_time ) ) )
total_mod_dat$onset_time_num_std <- normalize_func( total_mod_dat$onset_time_num )
total_mod_dat$waking_time_num <- as.numeric( as_hms( as.character( total_mod_dat$waking_time ) ) )
total_mod_dat$waking_time_num_std <- normalize_func( total_mod_dat$waking_time_num )

TST_model_dat <- total_mod_dat[ !is.na( total_mod_dat$TST ) & !is.na( total_mod_dat$total_num_corr ) & !is.na( total_mod_dat$tree ) & !is.na( total_mod_dat$spat.disc_5 ) & !is.na( total_mod_dat$pref_score ) & !is.na( total_mod_dat$prev_day_sleep_lim ) & !is.na( total_mod_dat$prev_TST_diff ) & !is.na( total_mod_dat$prev_fragG2_diff ) & total_mod_dat$night_num < 21, ]

options(mc.cores = parallel::detectCores() )


priors <- c(prior(normal(0,2), class = "b", coef = pref_score_std),
            
            prior(normal(0,2), class = "b", coef = fraction_std),
            
            prior(normal(0,2), class = "b", coef = total_num_corr_std),
            
            prior(normal(0,2), class = "b", coef = spat.disc_5_std),
            
            prior(normal(0,2), class = "b", coef = prev_day_sleep_lim_std),
            
            prior(normal(0,2), class = "b", coef = prev_TST_diff_std),
            
            prior(normal(0,2), class = "b", coef = prev_fragG2_diff_std),
            
            prior(normal(0,2), class = "b", coef = temperature_std),
            
            prior(normal(0,2), class = "b", coef = sexMale),
            
            prior(normal(0,2), class = "b", coef = ageJuvenile),
            
            prior(normal(0,2), class = "b", coef = ageSubadult),
            
            prior(normal(0,2), class = "sd", group = night),
            
            prior(normal(0,2), class = "sd", group = tag ),
            
            prior(normal(0,2), class = "b", coef = tree2 ),
            prior(normal(0,2), class = "b", coef = tree3 ),
            prior(normal(0,2), class = "b", coef = tree4 ),
            prior(normal(0,2), class = "b", coef = tree5 ),
            prior(normal(0,2), class = "b", coef = tree6 ),
            prior(normal(0,2), class = "b", coef = tree7 ),
            prior(normal(0,2), class = "b", coef = tree8 ),
            prior(normal(0,2), class = "b", coef = tree10 ),
            prior(normal(0,2), class = "b", coef = tree11 )
            
)



####### FINAL MODELS ##########

## Full TST model

TST_model_dat <- total_mod_dat[ !is.na( total_mod_dat$TST ) & !is.na( total_mod_dat$total_num_corr ) & !is.na( total_mod_dat$tree ) & !is.na( total_mod_dat$spat.disc_5 ) & !is.na( total_mod_dat$pref_score ) & !is.na( total_mod_dat$prev_day_sleep_lim ) & !is.na( total_mod_dat$prev_TST_diff ) & !is.na( total_mod_dat$prev_fragG2_diff ) & total_mod_dat$night_num < 21, ]

hist( TST_model_dat$TST, main = "Histogram of total sleep time", xlab = "Total sleep time (hours)" ) 

Bayes_mod_TST <- brm( TST_std ~ spat.disc_5_std + prev_day_sleep_lim_std + prev_TST_diff_std + prev_fragG2_diff_std + pref_score_std + total_num_corr_std + temperature_std + fraction_std + age + sex + tree + ( 1| tag ) + ( 1| night ), data = TST_model_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_mod_TST, file = "RESULTS/models/TST_std.rds" )

Bayes_mod_TST <- readRDS( file = "RESULTS/models/TST_std.rds" )

summary( Bayes_mod_TST )

tab_model( Bayes_mod_TST )

effects <- fixef( Bayes_mod_TST )[ - c( 1, 13:nrow( fixef( Bayes_mod_TST ) ) ), ]

par( bg = 'white' )

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n' )

axis( 1, seq( -1.5, 1, by = .5 ) )
abline( v = 0, lty = 2 )

for( i in 1:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16 )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i  )
  
  text( x = -1.8, y = i, labels = rownames( effects )[ i ]  )
  
}

### for presentation purposes, I will run the model without standardizing all variables

TST_model_dat <- total_mod_dat[ !is.na( total_mod_dat$TST ) & !is.na( total_mod_dat$total_num_corr ) & !is.na( total_mod_dat$tree ) & !is.na( total_mod_dat$spat.disc_5 ) & !is.na( total_mod_dat$pref_score ) & !is.na( total_mod_dat$prev_day_sleep_lim ) & !is.na( total_mod_dat$prev_TST_diff ) & !is.na( total_mod_dat$prev_fragG2_diff ) & total_mod_dat$night_num < 21, ]

priors_unstd <- c(prior(normal(0,2), class = "b", coef = pref_score ),
                  
                  prior(normal(0,2), class = "b", coef = fraction ),
                  
                  prior(normal(0,2), class = "b", coef = total_num_corr),
                  
                  prior(normal(0,2), class = "b", coef = spat.disc_5 ),
                  
                  prior(normal(0,2), class = "b", coef = prev_day_sleep_lim ),
                  
                  prior(normal(0,2), class = "b", coef = prev_TST_diff ),
                  
                  prior(normal(0,2), class = "b", coef = temperature ),
                  
                  prior(normal(0,2), class = "b", coef = sexMale ),
                  
                  prior(normal(0,2), class = "b", coef = ageJuvenile),
                  
                  prior(normal(0,2), class = "b", coef = ageSubadult),
                  
                  prior(normal(0,2), class = "sd", group = night),
                  
                  prior(normal(0,2), class = "sd", group = tag ),
                  
                  prior(normal(0,2), class = "b", coef = tree2 ),
                  prior(normal(0,2), class = "b", coef = tree3 ),
                  prior(normal(0,2), class = "b", coef = tree4 ),
                  prior(normal(0,2), class = "b", coef = tree5 ),
                  prior(normal(0,2), class = "b", coef = tree6 ),
                  prior(normal(0,2), class = "b", coef = tree7 ),
                  prior(normal(0,2), class = "b", coef = tree8 ),
                  prior(normal(0,2), class = "b", coef = tree10 ),
                  prior(normal(0,2), class = "b", coef = tree11 )
                  
)

Bayes_mod_TST_unstd <- brm( TST ~ spat.disc_5 + prev_day_sleep_lim + prev_TST_diff + prev_fragG2_diff + pref_score + total_num_corr + temperature + fraction + age + sex + tree + ( 1| tag ) + ( 1| night ), data = TST_model_dat, family= "gaussian", iter = 5000, prior = priors_unstd, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_mod_TST_unstd, file = "RESULTS/models/TST_unstd.rds" )

Bayes_mod_TST_unstd <- readRDS( file = "RESULTS/models/TST_unstd.rds" )

summary( Bayes_mod_TST_unstd )

tab_model( Bayes_mod_TST_unstd  )

plot( conditional_effects( Bayes_mod_TST_unstd, "tree" ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )

plot( conditional_effects( Bayes_mod_TST_unstd, "tree" ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] 

plot( conditional_effects( Bayes_mod_TST_unstd, "prev_day_sleep_lim", spaghetti = F ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )[[ 1 ]] + ylim( 8.2, 10 )

plot( conditional_effects( Bayes_mod_TST_unstd, "prev_day_sleep_lim", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( TST_model_dat$TST )  ) + geom_point( data = TST_model_dat, aes( x = prev_day_sleep_lim, y = TST ), alpha = 0.3, inherit.aes = F )


plot( conditional_effects( Bayes_mod_TST_unstd, "prev_TST_diff", spaghetti = F ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )[[ 1 ]] + ylim( 8.2, 10 )

plot( conditional_effects( Bayes_mod_TST_unstd, "prev_TST_diff", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( TST_model_dat$TST )  ) + geom_point( data = TST_model_dat, aes( x = prev_TST_diff, y = TST ), alpha = 0.3, inherit.aes = F )

plot( conditional_effects( Bayes_mod_TST_unstd, "spat.disc_5", spaghetti = F ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )[[ 1 ]] + ylim( 8.2, 10 )

plot( conditional_effects( Bayes_mod_TST_unstd, "spat.disc_5", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( TST_model_dat$TST )  ) + geom_point( data = TST_model_dat, aes( x = spat.disc_5, y = TST ), alpha = 0.3, inherit.aes = F )

plot( conditional_effects( Bayes_mod_TST_unstd, "prev_fragG2_diff", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( TST_model_dat$TST )  ) + geom_point( data = TST_model_dat, aes( x = prev_fragG2_diff, y = TST ), alpha = 0.3, inherit.aes = F )


plot( conditional_effects( Bayes_mod_TST_unstd, "pref_score", spaghetti = F ),  line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )

plot( conditional_effects( Bayes_mod_TST_unstd, "pref_score", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( TST_model_dat$TST )  ) + geom_point( data = TST_model_dat, aes( x = pref_score, y = TST ), alpha = 0.3, inherit.aes = F )


plot( conditional_effects( Bayes_mod_TST_unstd, "total_num_corr", spaghetti = F ),  line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme(  text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )

plot( conditional_effects( Bayes_mod_TST_unstd, "total_num_corr", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( TST_model_dat$TST )  ) + geom_point( data = TST_model_dat, aes( x = total_num_corr, y = TST ), alpha = 0.3, inherit.aes = F )


## this is what the theme should be if you want a white background:
# plot( conditional_effects( Bayes_mod_TST_unstd, "total_num_corr", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )


#### TST model with vedba instead of travel distance 


TST_model_dat <- total_mod_dat[ !is.na( total_mod_dat$TST ) & !is.na( total_mod_dat$total_num_corr ) & !is.na( total_mod_dat$tree ) & !is.na( total_mod_dat$vedba ) & !is.na( total_mod_dat$pref_score ) & !is.na( total_mod_dat$prev_day_sleep_lim ) & !is.na( total_mod_dat$prev_TST_diff ) & !is.na( total_mod_dat$prev_fragG2_diff ) & total_mod_dat$night_num < 21, ]

options(mc.cores = parallel::detectCores() )


priors <- c(prior(normal(0,2), class = "b", coef = pref_score_std),
            
            prior(normal(0,2), class = "b", coef = fraction_std),
            
            prior(normal(0,2), class = "b", coef = total_num_corr_std),
            
            prior(normal(0,2), class = "b", coef = vedba_std),
            
            prior(normal(0,2), class = "b", coef = prev_day_sleep_lim_std),
            
            prior(normal(0,2), class = "b", coef = prev_TST_diff_std),
            
            prior(normal(0,2), class = "b", coef = prev_fragG2_diff_std),
            
            prior(normal(0,2), class = "b", coef = temperature_std),
            
            prior(normal(0,2), class = "b", coef = sexMale),
            
            prior(normal(0,2), class = "b", coef = ageJuvenile),
            
            prior(normal(0,2), class = "b", coef = ageSubadult),
            
            prior(normal(0,2), class = "sd", group = night),
            
            prior(normal(0,2), class = "sd", group = tag ),
            
            prior(normal(0,2), class = "b", coef = tree2 ),
            prior(normal(0,2), class = "b", coef = tree3 ),
            prior(normal(0,2), class = "b", coef = tree4 ),
            prior(normal(0,2), class = "b", coef = tree5 ),
            prior(normal(0,2), class = "b", coef = tree6 ),
            prior(normal(0,2), class = "b", coef = tree7 ),
            prior(normal(0,2), class = "b", coef = tree8 ),
            prior(normal(0,2), class = "b", coef = tree10 ),
            prior(normal(0,2), class = "b", coef = tree11 )
            
)

hist( TST_model_dat$TST, main = "Histogram of total sleep time", xlab = "Total sleep time (hours)" ) 

Bayes_mod_ved_TST <- brm( TST_std ~ vedba_std + prev_day_sleep_lim_std + prev_TST_diff_std + prev_fragG2_diff_std + pref_score_std + total_num_corr_std + temperature_std + fraction_std  + age + sex + tree + ( 1| tag ) + ( 1| night ), data = TST_model_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_mod_ved_TST, file = "RESULTS/models/ved_TST.rds" )

Bayes_mod_ved_TST <- readRDS( file = "RESULTS/models/ved_TST.rds" )

tab_model( Bayes_mod_ved_TST )

effects <- fixef( Bayes_mod_ved_TST )[ -1, ]

par( bg = "black" )

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white'  )

axis( 1, seq( -1.5, 1, by = .5 ), col = 'white', col.ticks = 'white', col.axis = 'white')
abline( v = 0, lty = 2, col = 'white' )

for( i in 1:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16, col = 'white' )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i, col = 'white'  )
  
  text( x = -1.8, y = i, labels = rownames( effects )[ i ], col = 'white'  )
  
}



### for presentation purposes, I will run the model without standardizing all variables

priors_unstd <- c(prior(normal(0,2), class = "b", coef = pref_score ),
                  
                  prior(normal(0,2), class = "b", coef = fraction ),
                  
                  prior(normal(0,2), class = "b", coef = total_num_corr),
                  
                  prior(normal(0,2), class = "b", coef = vedba ),
                  
                  prior(normal(0,2), class = "b", coef = prev_day_sleep_lim ),
                  
                  prior(normal(0,2), class = "b", coef = prev_TST_diff ),
                  
                  prior(normal(0,2), class = "b", coef = temperature ),
                  
                  prior(normal(0,2), class = "b", coef = sexMale ),
                  
                  prior(normal(0,2), class = "b", coef = ageJuvenile),
                  
                  prior(normal(0,2), class = "b", coef = ageSubadult),
                  
                  prior(normal(0,2), class = "sd", group = night),
                  
                  prior(normal(0,2), class = "sd", group = tag ),
                  
                  prior(normal(0,2), class = "b", coef = tree2 ),
                  prior(normal(0,2), class = "b", coef = tree3 ),
                  prior(normal(0,2), class = "b", coef = tree4 ),
                  prior(normal(0,2), class = "b", coef = tree5 ),
                  prior(normal(0,2), class = "b", coef = tree6 ),
                  prior(normal(0,2), class = "b", coef = tree7 ),
                  prior(normal(0,2), class = "b", coef = tree8 ),
                  prior(normal(0,2), class = "b", coef = tree10 ),
                  prior(normal(0,2), class = "b", coef = tree11 )
                  
)

Bayes_mod_ved_TST_unstd <- brm( TST ~ vedba + prev_day_sleep_lim + prev_TST_diff + prev_fragG2_diff + pref_score + total_num_corr + temperature + fraction + age + sex + tree + ( 1| tag ) + ( 1| night ), data = TST_model_dat, family= "gaussian", iter = 5000, prior = priors_unstd, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_mod_ved_TST_unstd, file = "RESULTS/models/ved_TST_unstd.rds" )

Bayes_mod_ved_TST_unstd <- readRDS( file = "RESULTS/models/ved_TST_unstd.rds" )

summary( Bayes_mod_ved_TST_unstd )

tab_model( Bayes_mod_ved_TST_unstd )

##### napping next day mod

nap_dat <- read.csv( "DATA/main_data/nap_dat.csv" )

nap_dat$TST <- nap_dat$TST/60
nap_dat$TST_std <- normalize_func( nap_dat$TST )

nap_dat$SPT <- nap_dat$SPT / 60

nap_dat$fragG2 <- nap_dat$frag_wake_bouts / nap_dat$TST

nap_dat$fragG2_std <- normalize_func( nap_dat$fragG2 )

nap_dat$next_day_sleep_lim_std <- normalize_func( nap_dat$next_day_sleep_lim )

nap_dat$night_num <- as.numeric( nap_dat$night )

nap_dat$night <- as.factor( nap_dat$night )

nap_mod_dat <- nap_dat[ !is.na( nap_dat$TST ) & !is.na( nap_dat$fragG2 ) & !is.na( nap_dat$next_day_sleep_lim ) & nap_dat$night_num < 21, ]

hist( nap_mod_dat$next_day_sleep_lim, main = "Histogram of time spent napping", xlab = "Time spent napping (mins)"  )

priors <- c(      prior(normal(0,2), class = "b", coef = TST_std),
                  
                  prior(normal(0,2), class = "b", coef = fragG2_std),
                  
                  prior(normal(0,2), class = "sd", group = night),
                  
                  prior(normal(0,2), class = "sd", group = tag ))


nap_mod <- brm( next_day_sleep_lim_std ~  ( 1 | tag ) + ( 1 | night ) + TST_std + fragG2_std, data = nap_mod_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = nap_mod, file = "RESULTS/models/nap_std.rds" )

nap_mod <- readRDS( file = "RESULTS/models/nap_std.rds" )

effects <- fixef( nap_mod )[ -1, ]

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white'  )

axis( 1, seq( -1.5, 1, by = .5 ), col = 'white', col.ticks = 'white', col.axis = 'white')
abline( v = 0, lty = 2, col = 'white' )


for( i in 1:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16, col = 'white' )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i, col = 'white' )
  
  text( x = -1.5, y = i, labels = rownames( effects )[ i ], col = 'white' )
  
}

tab_model( nap_mod )

## and again, without standardizing, for presentation purposes


priors_unstd <- c(      prior(normal(0,2), class = "b", coef = TST),
                        
                        prior(normal(0,2), class = "b", coef = fragG2),
                        
                        prior(normal(0,2), class = "sd", group = night),
                        
                        prior(normal(0,2), class = "sd", group = tag ))


nap_mod_unstd <- brm( next_day_sleep_lim ~  ( 1 | tag ) + ( 1 | night ) + TST + fragG2, data = nap_mod_dat, family= "gaussian", iter = 5000, prior = priors_unstd, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = nap_mod_unstd, file = "RESULTS/models/nap_unstd.rds" )

nap_mod_unstd <- readRDS( file = "RESULTS/models/nap_unstd.rds" )

tab_model( nap_mod_unstd )

plot( conditional_effects( nap_mod_unstd, "TST", spaghetti = F ),  line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )





#### breaking down site to the different night conditions 


priors <- c(prior(normal(0,2), class = "b", coef = fraction_std),
            
            prior(normal(0,2), class = "b", coef = spat.disc_5_std),
            
            prior(normal(0,2), class = "b", coef = prev_day_sleep_lim_std),
            
            prior(normal(0,2), class = "b", coef = prev_TST_diff_std),
            
            prior(normal(0,2), class = "b", coef = prev_fragG2_diff_std),
            
            prior(normal(0,2), class = "b", coef = temperature_std),
            
            prior(normal(0,2), class = "b", coef = sexMale),
            
            prior(normal(0,2), class = "b", coef = ageJuvenile),
            
            prior(normal(0,2), class = "b", coef = ageSubadult),
            
            prior(normal(0,2), class = "sd", group = tag ),
            
            prior(normal(0,2), class = "b", coef = cond_night1),
            prior(normal(0,2), class = "b", coef = cond_night2),
            prior(normal(0,2), class = "b", coef = cond_night3),
            prior(normal(0,2), class = "b", coef = cond_night4),
            prior(normal(0,2), class = "b", coef = cond_night5)
            
            
)

total_mod_dat$cond_night <- as.factor( ifelse( total_mod_dat$night_num == 21, 1, ifelse( total_mod_dat$night_num == 22, 2, ifelse( total_mod_dat$night_num == 23, 3, ifelse( total_mod_dat$night_num == 24, 4, ifelse( total_mod_dat$night_num > 24, 5, 0 ) ) ) ) ) )

night_cond_mod_TST_dat <- total_mod_dat[ !is.na( total_mod_dat$TST ) & !is.na( total_mod_dat$spat.disc_5 ) & !is.na( total_mod_dat$site ) & !is.na( total_mod_dat$prev_day_sleep_lim ) & !is.na( total_mod_dat$prev_TST_diff ) & !is.na( total_mod_dat$prev_fragG2_diff_std ) , ]

hist( night_cond_mod_TST_dat$TST, main = "Histogram of total sleep time", xlab = "Total sleep time (hours)" ) 

Bayes_night_cond_mod_TST <- brm( TST_std ~ cond_night + age + sex + spat.disc_5_std + prev_day_sleep_lim_std + prev_TST_diff_std + prev_fragG2_diff_std + temperature_std + fraction_std + ( 1| tag ), data = night_cond_mod_TST_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_night_cond_mod_TST, file = "RESULTS/models/night_cond_std.rds" )

Bayes_night_cond_mod_TST <- readRDS( file = "RESULTS/models/night_cond_std.rds" )

summary( Bayes_night_cond_mod_TST )

effects <- fixef( Bayes_night_cond_mod_TST )

tab_model( Bayes_night_cond_mod_TST )

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white'  )

axis( 1, seq( -1.5, 1, by = .5 ), col = 'white', col.ticks = 'white', col.axis = 'white')
abline( v = 0, lty = 2, col = 'white' )

for( i in 1:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16, col = 'white' )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i, col = 'white' )
  
  text( x = -1.5, y = i, labels = rownames( effects )[ i ], col = 'white' )
  
}


## rerunning without standardization. For presentation purposes

priors <- c(prior(normal(0,2), class = "b", coef = fraction),
            
            prior(normal(0,2), class = "b", coef = spat.disc_5),
            
            prior(normal(0,2), class = "b", coef = prev_day_sleep_lim),
            
            prior(normal(0,2), class = "b", coef = prev_TST_diff),
            
            prior(normal(0,2), class = "b", coef = prev_fragG2_diff),
            
            prior(normal(0,2), class = "b", coef = temperature),
            
            prior(normal(0,2), class = "b", coef = sexMale),
            
            prior(normal(0,2), class = "b", coef = ageJuvenile),
            
            prior(normal(0,2), class = "b", coef = ageSubadult),
            
            prior(normal(0,2), class = "sd", group = tag ),
            
            prior(normal(0,2), class = "b", coef = cond_night2),
            prior(normal(0,2), class = "b", coef = cond_night3),
            prior(normal(0,2), class = "b", coef = cond_night4),
            prior(normal(0,2), class = "b", coef = cond_night5)
            
            
)

Bayes_night_cond_mod_TST_unstd <- brm( TST ~ cond_night + age + sex + spat.disc_5 + prev_day_sleep_lim + prev_TST_diff + prev_fragG2_diff + temperature + fraction + ( 1| tag ) + ( 1| night ), data = night_cond_mod_TST_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_night_cond_mod_TST_unstd, file = "RESULTS/models/night_cond_unstd.rds" )

Bayes_night_cond_mod_TST_unstd <- readRDS( file = "RESULTS/models/night_cond_unstd.rds" )

plot( conditional_effects( Bayes_night_cond_mod_TST_unstd, effects = c( 'cond_night' ), spaghetti = F ),  line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )

plot( conditional_effects( Bayes_night_cond_mod_TST_unstd, effects = c( 'cond_night' ), spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )


tab_model( Bayes_night_cond_mod_TST_unstd )

## plot the raw data of the nights ### Fig 3

# Basic violin plot
p <- ggplot( total_mod_dat, aes( x = night, y = TST, color = site, fill = site ) ) + 
  geom_violin( trim = FALSE )

p1 <- p+ scale_color_manual(values=c("green", "blue" ) ) + scale_fill_manual(values=c("green", "blue" ) ) + theme_classic() + theme( legend.position = 'none' ) 

p1 + scale_x_discrete(labels= c('', 'Aug-01', rep( '', 7 ), 'Aug-09', rep( '', 7 ), 'Aug-17', rep( '', 7 ), 'Aug-25', rep( '', 7 ), 'Sep-02', rep( '', 7 ), 'blah' ) )



######### Models of fragG2  #############

fragG2_model_dat <- total_mod_dat[ !is.na( total_mod_dat$fragG2 ) & !is.na( total_mod_dat$total_num_corr ) & !is.na( total_mod_dat$tree ) & !is.na( total_mod_dat$spat.disc_5 ) & !is.na( total_mod_dat$pref_score ) & !is.na( total_mod_dat$prev_day_sleep_lim ) & !is.na( total_mod_dat$prev_TST_diff ) & !is.na( total_mod_dat$prev_fragG2_diff ) & total_mod_dat$night_num < 21, ]

options(mc.cores = parallel::detectCores() )

priors <- c(prior(normal(0,2), class = "b", coef = pref_score_std),
            
            prior(normal(0,2), class = "b", coef = fraction_std),
            
            prior(normal(0,2), class = "b", coef = total_num_corr_std),
            
            prior(normal(0,2), class = "b", coef = spat.disc_5_std),
            
            prior(normal(0,2), class = "b", coef = prev_day_sleep_lim_std),
            
            prior(normal(0,2), class = "b", coef = prev_TST_diff_std),

            prior(normal(0,2), class = "b", coef = prev_fragG2_diff_std),
            
            prior(normal(0,2), class = "b", coef = temperature_std),
            
            prior(normal(0,2), class = "b", coef = sexMale),
            
            prior(normal(0,2), class = "b", coef = ageJuvenile),
            
            prior(normal(0,2), class = "b", coef = ageSubadult),
            
            prior(normal(0,2), class = "sd", group = night),
            
            prior(normal(0,2), class = "sd", group = tag ),
            
            prior(normal(0,2), class = "b", coef = tree2 ),
            prior(normal(0,2), class = "b", coef = tree3 ),
            prior(normal(0,2), class = "b", coef = tree4 ),
            prior(normal(0,2), class = "b", coef = tree5 ),
            prior(normal(0,2), class = "b", coef = tree6 ),
            prior(normal(0,2), class = "b", coef = tree7 ),
            prior(normal(0,2), class = "b", coef = tree8 ),
            prior(normal(0,2), class = "b", coef = tree10 ),
            prior(normal(0,2), class = "b", coef = tree11 )
            
)



####### FINAL MODELS ##########

## fragG2 model

hist( fragG2_model_dat$fragG2, main = "Histogram of sleep fragmentation", xlab = "Sleep fragmentation (wake bouts / hour of sleep)", col = 'white' ) 

Bayes_mod_fragG2 <- brm( fragG2_std ~ spat.disc_5_std + prev_day_sleep_lim_std + prev_TST_diff_std + prev_fragG2_diff_std + pref_score_std + total_num_corr_std + temperature_std + fraction_std + age + sex + tree + ( 1| tag ) + ( 1| night ), data = fragG2_model_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )
 
saveRDS( object = Bayes_mod_fragG2, file = "RESULTS/models/fragG2_std.rds" )

Bayes_mod_fragG2 <- readRDS( file = "RESULTS/models/fragG2_std.rds" )

summary( Bayes_mod_fragG2 )

tab_model( Bayes_mod_fragG2 )

effects <- fixef( Bayes_mod_fragG2 )[ - c( 1, 13:nrow( fixef( Bayes_mod_TST ) ) ), ]

par( bg = 'white' )

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n' )

axis( 1, seq( -1.5, 1, by = .5 ) )
abline( v = 0, lty = 2 )

for( i in 1:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16 )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i  )
  
  text( x = -1.8, y = i, labels = rownames( effects )[ i ]  )
  
}

### for presentation purposes, I will run the model without standardizing all variables

priors_unstd <- c(prior(normal(0,2), class = "b", coef = pref_score ),
                  
                  prior(normal(0,2), class = "b", coef = fraction ),
                  
                  prior(normal(0,2), class = "b", coef = total_num_corr),
                  
                  prior(normal(0,2), class = "b", coef = spat.disc_5 ),
                  
                  prior(normal(0,2), class = "b", coef = prev_day_sleep_lim ),
                  
                  prior(normal(0,2), class = "b", coef = prev_TST_diff ),
                  
                  prior(normal(0,2), class = "b", coef = temperature ),
                  
                  prior(normal(0,2), class = "b", coef = sexMale ),
                  
                  prior(normal(0,2), class = "b", coef = ageJuvenile),
                  
                  prior(normal(0,2), class = "b", coef = ageSubadult),
                  
                  prior(normal(0,2), class = "sd", group = night),
                  
                  prior(normal(0,2), class = "sd", group = tag ),
                  
                  prior(normal(0,2), class = "b", coef = tree2 ),
                  prior(normal(0,2), class = "b", coef = tree3 ),
                  prior(normal(0,2), class = "b", coef = tree4 ),
                  prior(normal(0,2), class = "b", coef = tree5 ),
                  prior(normal(0,2), class = "b", coef = tree6 ),
                  prior(normal(0,2), class = "b", coef = tree7 ),
                  prior(normal(0,2), class = "b", coef = tree8 ),
                  prior(normal(0,2), class = "b", coef = tree10 ),
                  prior(normal(0,2), class = "b", coef = tree11 )
                  
)

Bayes_mod_fragG2_unstd <- brm( fragG2 ~ spat.disc_5 + prev_day_sleep_lim + prev_TST_diff + prev_fragG2_diff + pref_score + total_num_corr + temperature + fraction + age + sex + tree + ( 1| tag ) + ( 1| night ), data = fragG2_model_dat, family= "gaussian", iter = 5000, prior = priors_unstd, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_mod_fragG2_unstd, file = "RESULTS/models/fragG2_unstd.rds" )

Bayes_mod_fragG2_unstd <- readRDS( file = "RESULTS/models/fragG2_unstd.rds" )

summary( Bayes_mod_fragG2_unstd )

tab_model( Bayes_mod_fragG2_unstd  )

plot( conditional_effects( Bayes_mod_fragG2_unstd, "tree" ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] 

plot( conditional_effects( Bayes_mod_fragG2_unstd, "tree" ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )


plot( conditional_effects( Bayes_mod_fragG2_unstd, "prev_day_sleep_lim", spaghetti = F ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )[[ 1 ]] 

plot( conditional_effects( Bayes_mod_fragG2_unstd, "prev_day_sleep_lim", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( fragG2_model_dat$fragG2 )  ) + geom_point( data = fragG2_model_dat, aes( x = prev_day_sleep_lim, y = fragG2 ), alpha = 0.3, inherit.aes = F )



plot( conditional_effects( Bayes_mod_fragG2_unstd, "prev_TST_diff", spaghetti = F ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )[[ 1 ]]

plot( conditional_effects( Bayes_mod_fragG2_unstd, "prev_TST_diff", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]]

plot( conditional_effects( Bayes_mod_fragG2_unstd, "prev_TST_diff", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( fragG2_model_dat$fragG2 )  ) + geom_point( data = fragG2_model_dat, aes( x = prev_TST_diff, y = fragG2 ), alpha = 0.3, inherit.aes = F )

plot( conditional_effects( Bayes_mod_fragG2_unstd, "prev_fragG2_diff", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( fragG2_model_dat$fragG2 )  ) + geom_point( data = fragG2_model_dat, aes( x = prev_fragG2_diff, y = fragG2 ), alpha = 0.3, inherit.aes = F )



plot( conditional_effects( Bayes_mod_fragG2_unstd, "spat.disc_5", spaghetti = F ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )[[ 1 ]]

plot( conditional_effects( Bayes_mod_fragG2_unstd, "spat.disc_5", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]]

plot( conditional_effects( Bayes_mod_fragG2_unstd, "spat.disc_5", spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + ylim( range( fragG2_model_dat$fragG2 )  ) + geom_point( data = fragG2_model_dat, aes( x = spat.disc_5, y = fragG2 ), alpha = 0.3, inherit.aes = F )



plot( conditional_effects( Bayes_mod_fragG2_unstd, "pref_score", spaghetti = F ),  line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme( text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )



plot( conditional_effects( Bayes_mod_fragG2_unstd, "total_num_corr", spaghetti = F ),  line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme(  text = element_text( size = 25 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )

#### with vedba instead of travel distance 
fragG2_model_dat <- total_mod_dat[ !is.na( total_mod_dat$fragG2 ) & !is.na( total_mod_dat$total_num_corr ) & !is.na( total_mod_dat$tree ) & !is.na( total_mod_dat$vedba ) & !is.na( total_mod_dat$pref_score ) & !is.na( total_mod_dat$prev_day_sleep_lim ) & !is.na( total_mod_dat$prev_TST_diff ) & !is.na( total_mod_dat$prev_fragG2_diff ) & total_mod_dat$night_num < 21, ]

options(mc.cores = parallel::detectCores() )

priors <- c(prior(normal(0,2), class = "b", coef = pref_score_std),
            
            prior(normal(0,2), class = "b", coef = fraction_std),
            
            prior(normal(0,2), class = "b", coef = total_num_corr_std),
            
            prior(normal(0,2), class = "b", coef = vedba_std),
            
            prior(normal(0,2), class = "b", coef = prev_day_sleep_lim_std),
            
            prior(normal(0,2), class = "b", coef = prev_TST_diff_std),
            
            prior(normal(0,2), class = "b", coef = prev_fragG2_diff_std),
            
            prior(normal(0,2), class = "b", coef = temperature_std),
            
            prior(normal(0,2), class = "b", coef = sexMale),
            
            prior(normal(0,2), class = "b", coef = ageJuvenile),
            
            prior(normal(0,2), class = "b", coef = ageSubadult),
            
            prior(normal(0,2), class = "sd", group = night),
            
            prior(normal(0,2), class = "sd", group = tag ),
            
            prior(normal(0,2), class = "b", coef = tree2 ),
            prior(normal(0,2), class = "b", coef = tree3 ),
            prior(normal(0,2), class = "b", coef = tree4 ),
            prior(normal(0,2), class = "b", coef = tree5 ),
            prior(normal(0,2), class = "b", coef = tree6 ),
            prior(normal(0,2), class = "b", coef = tree7 ),
            prior(normal(0,2), class = "b", coef = tree8 ),
            prior(normal(0,2), class = "b", coef = tree10 ),
            prior(normal(0,2), class = "b", coef = tree11 )
            
)

hist( fragG2_model_dat$fragG2, main = "Histogram of total sleep time", xlab = "Total sleep time (hours)" ) 

Bayes_mod_ved_fragG2 <- brm( fragG2_std ~ vedba_std + prev_day_sleep_lim_std + prev_TST_diff_std + prev_fragG2_diff_std + pref_score_std + total_num_corr_std + temperature_std + fraction_std  + age + sex + tree + ( 1| tag ) + ( 1| night ), data = fragG2_model_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_mod_ved_fragG2, file = "RESULTS/models/ved_fragG2.rds" )

Bayes_mod_ved_fragG2 <- readRDS( file = "RESULTS/models/ved_fragG2.rds" )

tab_model( Bayes_mod_ved_fragG2 )

effects <- fixef( Bayes_mod_ved_fragG2 )[ -1, ]

par( bg = 'black' )

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white'  )

axis( 1, seq( -1.5, 1, by = .5 ), col = 'white', col.ticks = 'white', col.axis = 'white')
abline( v = 0, lty = 2, col = 'white' )

for( i in 1:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16, col = 'white' )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i, col = 'white'  )
  
  text( x = -1.8, y = i, labels = rownames( effects )[ i ], col = 'white'  )
  
}

### for presentation purposes, I will run the model without standardizing all variables

priors_unstd <- c(prior(normal(0,2), class = "b", coef = pref_score ),
                  
                  prior(normal(0,2), class = "b", coef = fraction ),
                  
                  prior(normal(0,2), class = "b", coef = total_num_corr),
                  
                  prior(normal(0,2), class = "b", coef = vedba ),
                  
                  prior(normal(0,2), class = "b", coef = prev_day_sleep_lim ),
                  
                  prior(normal(0,2), class = "b", coef = prev_TST_diff ),
                  
                  prior(normal(0,2), class = "b", coef = temperature ),
                  
                  prior(normal(0,2), class = "b", coef = sexMale ),
                  
                  prior(normal(0,2), class = "b", coef = ageJuvenile),
                  
                  prior(normal(0,2), class = "b", coef = ageSubadult),
                  
                  prior(normal(0,2), class = "sd", group = night),
                  
                  prior(normal(0,2), class = "sd", group = tag ),
                  
                  prior(normal(0,2), class = "b", coef = tree2 ),
                  prior(normal(0,2), class = "b", coef = tree3 ),
                  prior(normal(0,2), class = "b", coef = tree4 ),
                  prior(normal(0,2), class = "b", coef = tree5 ),
                  prior(normal(0,2), class = "b", coef = tree6 ),
                  prior(normal(0,2), class = "b", coef = tree7 ),
                  prior(normal(0,2), class = "b", coef = tree8 ),
                  prior(normal(0,2), class = "b", coef = tree10 ),
                  prior(normal(0,2), class = "b", coef = tree11 )
                  
)

Bayes_mod_ved_fragG2_unstd <- brm( fragG2 ~ vedba + prev_day_sleep_lim + prev_TST_diff + pref_score + total_num_corr + temperature + fraction + age + sex + tree + ( 1| tag ) + ( 1| night ), data = fragG2_model_dat, family= "gaussian", iter = 5000, prior = priors_unstd, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_mod_ved_fragG2_unstd, file = "RESULTS/models/ved_fragG2_unstd.rds" )

Bayes_mod_ved_fragG2_unstd <- readRDS( file = "RESULTS/models/ved_fragG2_unstd.rds" )

summary( Bayes_mod_ved_fragG2_unstd )

tab_model( Bayes_mod_ved_fragG2_unstd )

#### breaking down site to the different night conditions

priors <- c(prior(normal(0,2), class = "b", coef = fraction_std),
            
            prior(normal(0,2), class = "b", coef = spat.disc_5_std),
            
            prior(normal(0,2), class = "b", coef = prev_day_sleep_lim_std),
            
            prior(normal(0,2), class = "b", coef = prev_TST_diff_std),
            
            prior(normal(0,2), class = "b", coef = prev_fragG2_diff_std),
            
            prior(normal(0,2), class = "b", coef = temperature_std),
            
            prior(normal(0,2), class = "b", coef = sexMale),
            
            prior(normal(0,2), class = "b", coef = ageJuvenile),
            
            prior(normal(0,2), class = "b", coef = ageSubadult),
            
            prior(normal(0,2), class = "sd", group = tag ),
            
            prior(normal(0,2), class = "b", coef = cond_night1),
            prior(normal(0,2), class = "b", coef = cond_night2),
            prior(normal(0,2), class = "b", coef = cond_night3),
            prior(normal(0,2), class = "b", coef = cond_night4),
            prior(normal(0,2), class = "b", coef = cond_night5)
            
            
)

total_mod_dat$cond_night <- as.factor( ifelse( total_mod_dat$night_num == 21, 1, ifelse( total_mod_dat$night_num == 22, 2, ifelse( total_mod_dat$night_num == 23, 3, ifelse( total_mod_dat$night_num == 24, 4, ifelse( total_mod_dat$night_num > 24, 5, 0 ) ) ) ) ) )

night_cond_mod_fragG2_dat <- total_mod_dat[ !is.na( total_mod_dat$fragG2 ) & !is.na( total_mod_dat$temperature ) & !is.na( total_mod_dat$spat.disc_5 ) & !is.na( total_mod_dat$cond_night ) & !is.na( total_mod_dat$prev_day_sleep_lim ) & !is.na( total_mod_dat$prev_TST_diff ) & !is.na( total_mod_dat$prev_fragG2_diff ) , ]


hist( night_cond_mod_fragG2_dat$fragG2, main = "Histogram of sleep fragmentation", xlab = "Sleep fragmentation (wake bouts / hour of sleep)", col = 'white' ) 

Bayes_night_cond_mod_fragG2 <- brm( fragG2_std ~ cond_night + age + sex + spat.disc_5_std + prev_day_sleep_lim_std + prev_TST_diff_std + prev_fragG2_diff_std + temperature_std + fraction_std + ( 1| tag ), data = night_cond_mod_fragG2_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_night_cond_mod_fragG2, file = "RESULTS/models/night_cond_fragG2_std.rds" )

Bayes_night_cond_mod_fragG2 <- readRDS( file = "RESULTS/models/night_cond_fragG2_std.rds" )

summary( Bayes_night_cond_mod_fragG2 )

effects <- fixef( Bayes_night_cond_mod_fragG2 )

tab_model( Bayes_night_cond_mod_fragG2 )

par( bg = 'black' )

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white'  )

axis( 1, seq( -1.5, 1, by = .5 ), col = 'white', col.ticks = 'white', col.axis = 'white')
abline( v = 0, lty = 2, col = 'white' )

for( i in 2:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16, col = 'white' )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i, col = 'white' )
  
  text( x = -1.5, y = i, labels = rownames( effects )[ i ], col = 'white' )
  
}


## rerunning without standardization. For presentation purposes

priors <- c(prior(normal(0,2), class = "b", coef = fraction),
            
            prior(normal(0,2), class = "b", coef = spat.disc_5),
            
            prior(normal(0,2), class = "b", coef = prev_day_sleep_lim),
            
            prior(normal(0,2), class = "b", coef = prev_TST_diff),
            
            prior(normal(0,2), class = "b", coef = temperature),
            
            prior(normal(0,2), class = "b", coef = sexMale),
            
            prior(normal(0,2), class = "b", coef = ageJuvenile),
            
            prior(normal(0,2), class = "b", coef = ageSubadult),
            
            prior(normal(0,2), class = "sd", group = tag ),
            
            prior(normal(0,2), class = "b", coef = cond_night2),
            prior(normal(0,2), class = "b", coef = cond_night3),
            prior(normal(0,2), class = "b", coef = cond_night4),
            prior(normal(0,2), class = "b", coef = cond_night5)
            
            
)

Bayes_night_cond_mod_fragG2_unstd <- brm( fragG2 ~ cond_night + age + sex + spat.disc_5 + prev_day_sleep_lim + prev_TST_diff + prev_fragG2_diff + temperature + fraction + ( 1| tag ) + ( 1| night ), data = night_cond_mod_fragG2_dat, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = Bayes_night_cond_mod_fragG2_unstd, file = "RESULTS/models/night_cond_fragG2_unstd.rds" )

Bayes_night_cond_mod_fragG2_unstd <- readRDS( file = "RESULTS/models/night_cond_fragG2_unstd.rds" )

plot( conditional_effects( Bayes_night_cond_mod_fragG2_unstd, effects = c( 'cond_night' ), spaghetti = F ),  line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )

plot( conditional_effects( Bayes_night_cond_mod_fragG2_unstd, effects = c( 'cond_night' ), spaghetti = F ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )


tab_model( Bayes_night_cond_mod_fragG2_unstd )

######################### CHECKPOINT 1.5 ############################


# save.image( "DATA/main_data/sleep_analysis_pub_code_part_mods.RData" )

load( "DATA/main_data/sleep_analysis_pub_code_part_mods.RData" )

# rm( list = c( "d1", "cl_bdat", "bdat", "full_sleep", "fuller_sleep", "fuller_sleep_temp", "morning", "opp_dat", "pre_clean_full", "t_bdat", "tag_night_dat" ) )


### gamm with splines instead of quadratic

epoch_model_1 <- gamm( sleep_bouts~s( prop_from_sp_start ), data = full_dat_sleep, family = 'binomial', random = list( tag = ~1, night = ~1 ),
                       correlation = corAR1() )

summary( epoch_model_1$gam )

par(bg = 'black')

plot( epoch_model_1$gam, all.terms = TRUE, bty = "l", las = 1, col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white', ylab = '', xlab = '' )
axis(1, col = 'white', at = seq(0, 1, by = 0.2), labels = rep( "", 6 ) )
axis(2, col = 'white', at = seq(-1, 0.4, by = 0.2), labels = rep("", 8 ) )

######################### CHECKPOINT 2 ############################

#save.image( "DATA/main_data/sleep_analysis_pub_code_post_mods.RData" )

load( "DATA/main_data/sleep_analysis_pub_code_post_mods.RData" )

#### Tables for supplementary information #####

tag_meta <- data.frame( tag = sort( as.character( unique( bdat$id ) ) ), start_GPS = NA, start_ACC = NA, end_GPS = NA, end_ACC = NA )

for( i in 1:nrow( tag_meta ) ){
  
  tag_GPS <- bdat[ bdat$id == tag_meta$tag[ i ], ]
  
  tag_ACC <- d1[ d1$tag == tag_meta$tag[ i ], ]
  
  tag_meta$start_GPS[ i ] <- min( tag_GPS$timestamp )
  
  tag_meta$end_GPS[ i ] <- max( tag_GPS$timestamp )
  
  tag_meta$start_ACC[ i ] <- min( tag_ACC$timestamp )
  
  tag_meta$end_ACC[ i ] <- max( tag_ACC$timestamp )
}


for( i in 2:5 ){
  
  tag_meta[ , i ] <- as.POSIXct( tag_meta[ , i], origin = '1970-01-01', tz = 'UTC' )
  
  tag_meta[ , i ] <- as.Date( tag_meta[ , i] )
  
}


IDs <- read.csv( "DATA/main_data/IDs.csv" )


merged_dat <- merge( x = IDs, y = tag_meta, by.x = 'Collar', by.y = 'tag', all.x = F, all.y = T, sort = F )

keep_dat <- merged_dat[ , c( 'Collar', 'Sex', 'Age', 'Weight', 'Day.Captured', names( tag_meta )[ -1 ] ) ]

keep_dat$Day.Captured <- as.Date( keep_dat$Day.Captured, '%d/%m/%Y' )

keep_dat$Day.Captured <- paste0( "20", str_split_fixed( keep_dat$Day.Captured, "00", 2 )[ , 2 ] )

write.csv( keep_dat, "DATA/main_data/tag_meta.csv", row.names = F )


## TST difference between most preferred tree and least preferred tree

multiplier <- fixef( Bayes_mod_TST_unstd )[ 'pref_score', 'Estimate' ]

tag_names <- unique( total_mod_dat$tag )

diff_vec <- c() 

for( tag in tag_names ){
  
  id_dat <- total_mod_dat[ total_mod_dat$tag == tag, ]
  
  if( sum( !is.na( id_dat$pref_score ) ) != 0 ){
    
    min_score <- min( id_dat$pref_score, na.rm = T )
    
    max_score <- max( id_dat$pref_score, na.rm = T )
    
    sleep_dff <- ( max_score * multiplier - min_score * multiplier)
    
    diff_vec <- c( diff_vec, sleep_dff )
  }
  
}

mean( diff_vec )

max( diff_vec )*60
