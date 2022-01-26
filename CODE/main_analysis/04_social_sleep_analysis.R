
library( zoo )
library( data.table )
library( stringr )
library( brms )
library( sjPlot )
library( ggplot2 )

## function for setting transparency of a color while plotting
transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

## function for normalizing a vector
normalize_func <- function( x ) return( (x - mean( x, na.rm = T ) )/ sd( x, na.rm = T ) )

na.max <- function( x ){
  
  if( sum( !is.na( x ) ) == 0 ){
    
    return( NA )
  }else{
    
    return( max( x, na.rm = T ) )
  }
}


doing_same <- function( x ){
  
  if( sum( !is.na( x ) ) == 0 ){
    
    return( NA )
    
  }else{
    
    s1 <- sum( x == 0, na.rm = T ) / sum( !is.na( x ) )
    
    s2 <- max( s1, 1 - s1 )
    
    
    return( s2 )
  }
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

## apply sleep classification to the night-time period

waso_percentile <- 0.10 ## this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity within the sleep period

waso_multiplier <- 1.125 ## this is the multiplier of the threshold value determined by the percentile above. Values WITHIN the sleep period below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active. This threshold value may be different than the one above even when waso_percentile = percentile and waso_multiplier = multiplier, because the value of the given percentile of this variable depends on the smoothing window (see waso_window below; aka waso window might not be equal to mov_window)

waso_window <- 1 ## this is the size of the moving window (in minutes) used in calculating the rolling median that will be used to find periods of wake after sleeping. A waso_window of 1 is the same as using the raw data without a rolling median

waso_block <- 3 ## this is the number of consecutive minutes of inactivity needed to classify a period as sleep within the sleep period. A waso_block of 1 means that anytime the value is below the threshold, the baboon in considered sleeping and anytime the value is above the threshold the baboon is considered awake


night_start <- '21:00:00'

night_end <- '05:00:00'

tag_names <- as.character( unique( d1$tag ) )

timestamps <- seq( min( d1$local_timestamp) , max( d1$local_timestamp), by = '1 min' )

total_dat <- as.data.frame( matrix( NA, nrow = length( timestamps ), ncol = 1 + length( tag_names ) ) )

names( total_dat ) <- c( 'local_timestamp', tag_names )

total_dat$local_timestamp <- timestamps

total_dat$local_time <- str_split_fixed( total_dat$local_timestamp, ' ', 2 )[ , 2 ]

total_dat$night <- as.numeric( as.Date( total_dat$local_timestamp - 12*60*60 ) - min( as.Date( total_dat$local_timestamp - 12*60*60 ) ) + 1 )

tag_names <- as.character( unique( d1$tag ) )

for( tag in tag_names ){
  
  id_dat <- d1[ d1$tag == tag, ]
  
  nights <- as.character( unique( id_dat$night ) )
  
  for( night in nights ){
    
    night_dat <- id_dat[ id_dat$night == night, ]
    
    night_dat$roll_log_vedba <- rollmedian( night_dat$log_vedba, waso_window, fill = NA, align = 'center' )
    
    thresh <- quantile( night_dat$roll_log_vedba, waso_percentile, na.rm = T) * waso_multiplier
    
    temp <- rle(as.numeric( night_dat$roll_log_vedba < thresh ) ) 
    
    night_dat$runs <- as.numeric( rep( temp$lengths >= waso_block, times = temp$lengths ) )
    
    night_dat$awake <- as.numeric( !( night_dat$roll_log_vedba < thresh & night_dat$runs == 1 ) )
    total_dat[ match( night_dat$local_timestamp, total_dat$local_timestamp ), tag ] <- night_dat$awake
    
  }
}

### this loop is going to remove partial nights of data at the end of an individual's data so that we are only randomizing full nights of data
for( tag in tag_names ){
  
  max_ind <- max( which( !is.na( total_dat[ , tag ] ) ) )
  
  max_time <- total_dat[ max_ind, 'local_time' ]
  
  max_night <- total_dat[ max_ind, 'night' ]
  
  if( max_time < night_end | max_time > night_start ){
    
    
    total_dat[ total_dat$night == max_night, tag ] <- NA
    
  }
  
}

sent_dat <- total_dat[ total_dat$local_time > night_start | total_dat$local_time < night_end, ]

## remove the first night because there is only one individual on this night that actually has data
sent_dat <- sent_dat[ sent_dat$night != 1, ]

## calculating individual and collective vigilance ##

coll_vig_mean <- c()

coll_vig_raw <- c()

ind_vig <- c()

for( nigh in 1:length( unique( sent_dat$night ) ) ){
  
  nigh_dat <- sent_dat[ sent_dat$night == unique( sent_dat$night )[ nigh ], ]
  
  coll_vig_raw <- c( coll_vig_raw, sum( apply( nigh_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max ), na.rm = T ) )
  
  coll_vig_mean <- c( coll_vig_mean, mean( apply( nigh_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max ), na.rm = T ) )
  
  for( tag in tag_names ){
    
    t_vec <- nigh_dat[ , tag ]
  
    if( sum( !is.na( t_vec ) ) != 0 ){
      
      ind_vig <- c( ind_vig, sum( t_vec, na.rm =  T ) )
      
    }
  }
}

mean( coll_vig_raw )

sd( coll_vig_raw ) / sqrt( length( coll_vig_raw ) )

mean( coll_vig_mean )

mean( ind_vig )

sd( ind_vig ) / sqrt( length( ind_vig ) )


## Test whether at least one individual in the group is awake more often than expected by chance (sentinel hypothesis), by shifting each individual's sleep-wake data within each given night by a random amount and comparing that to the empirical data

num_rands <- 1000

shift_results_sent <- data.frame( matrix( NA, nrow = nrow( sent_dat ), ncol = 3 + num_rands ) ) 

names( shift_results_sent ) <- c( 'local_timestamp', 'night', 'emp', paste( 'rand', 1:num_rands, sep = '_' ) )

shift_results_sent$local_timestamp <- sent_dat$local_timestamp

shift_results_sent$night <- sent_dat$night

shift_results_sent$emp <- apply( sent_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max )

emp_sent_prop <- mean( apply( sent_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max ), na.rm = T )


## now randomize

rand_sent_vec <- c()

nights <- unique( sent_dat$night )

for( n in 1:num_rands ){
  
  print( n )
  
  rand_dat <- sent_dat
  
  for( night in nights ){
    
    night_dat <- rand_dat[ rand_dat$night == night, ]
    
    for( tag in tag_names ){
      
      shift <- sample( 1: nrow( night_dat ), 1 )
      
      new_inds <- 1:nrow( night_dat ) + shift
      new_inds <- new_inds %% nrow( night_dat ) + 1
      
      
      night_dat[ , tag ] <- night_dat[ , tag ][ new_inds ]
      
    }
    
    rand_dat[ rand_dat$night == night, ] <- night_dat
    
  }
  
  shift_results_sent[ , paste( 'rand', n, sep = '_' ) ] <- apply( rand_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max )
  
  sent_prop <- mean( apply( rand_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max ), na.rm = T )
  
  rand_sent_vec <- c( rand_sent_vec, sent_prop )
  
}

## find the mean of each column (each column represents one randomization and the first column is the empirical)
means <- apply( shift_results_sent[ , -( 1:2 ) ], 2, mean, na.rm = T )

## pull out the mean of the empirical
emp <- means[ 1 ]

## save a vector with the mean of each randomization
rand_vec <- means[ -1 ]

## the proportion of random that are as extreme or more extreme than the empirical (aka the p-value)
mean( emp >= rand_vec )

dens <- density( rand_vec )

dens_x <- dens$x
dens_y <- dens$y

par(bg = 'black')

plot( dens_x, dens_y, ylim = range( dens_y ), xlim = range( c( dens_x ) ), type = 'l', lwd = 2, xlab = 'Proportion of time at least one individual is awake', ylab = 'Probability density', bty = 'l', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )
axis(1, col = 'white', at = seq(0.865, 0.88, by = 0.005), labels = rep( "", 4 ) )
axis(2, col = 'white', at = seq(0, 150, by = 50), labels = rep("", 4))

polygon( x = dens_x, y = dens_y, col = transp( 'white', 0.3) )


plot( dens_x, dens_y, ylim = range( dens_y ), xlim = range( c( dens_x, emp_sent_prop ) ), type = 'l', lwd = 2, xlab = 'Proportion of time at least one individual is awake', ylab = 'Probability density', bty = 'l', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )
axis(1, col = 'white', at = seq(0.83, 0.88, by = 0.01), labels = rep( "", 6 ) )
axis(2, col = 'white', at = seq(0, 150, by = 50), labels = rep("", 4))


polygon( x = dens_x, y = dens_y, col = transp( 'white', 0.3) )

segments( x0 = emp_sent_prop, x1 = emp_sent_prop, y0 = 0, y1 = max( dens_y), lwd = 2, col = 'red', lty = 2)


##### Test the sentinel hypothesis by permuting sleep-wake data between different nights, rather than by applying random shifts to the data within the night #####

num_rands <- 1000

swap_results_sent <- data.frame( matrix( NA, nrow = nrow( sent_dat ), ncol = 3 + num_rands ) ) 

names( swap_results_sent ) <- c( 'local_timestamp', 'night', 'emp', paste( 'rand', 1:num_rands, sep = '_' ) )

swap_results_sent$local_timestamp <- sent_dat$local_timestamp

swap_results_sent$night <- sent_dat$night

swap_results_sent$emp <- apply( sent_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max )

emp_sent_prop <- mean( apply( sent_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max ), na.rm = T )


start_date <- as.Date( min( d1$local_timestamp - 12*60*60 ) )


## now randomize

rand_sent_vec <- c()

for( n in 1:num_rands ){
  
  print( n )
  
  rand_dat <- sent_dat
  
  for( tag in tag_names ){
    
    if( sum( !is.na( sent_dat[ , tag ] ) ) != 0 ){
      
      min_ind <- min( which( !is.na( sent_dat[ , tag ] ) ) )
      max_ind <- max( which( !is.na( sent_dat[ , tag ] ) ) )
      
      sub <- sent_dat[ min_ind:max_ind, ]
      
      nights_avail <- sort( unique( sub$night ) )
      
      if( length( nights_avail ) != 1 ){
        
        new_nights <- sample( nights_avail )
        
        look_up <- data.frame( original = nights_avail, new = new_nights )
        
        temp <- merge( x = sub, y = look_up, by.x = 'night', by.y = 'original' )
        
        ## first add the new date based on the new randomized timestamp
        temp$new_timestamp <- start_date + temp$new - 1
        
        ## add a day to the rows that happen after midnight (when the date should increase by one)=
        temp$new_timestamp[ temp$local_time < '12:00:00' ] <- temp$new_timestamp[ temp$local_time < '12:00:00' ] + 1
        
        ## add the time to our new randomized date
        temp$new_timestamp <- as.POSIXct( paste( temp$new_timestamp, temp$local_time ), tz = 'UTC' )
        
        ## reorder the vector of sleep to put the new timestamps in order
        temp$new_dat <- temp[ match( temp$local_timestamp, temp$new_timestamp ), tag ]
        
        rand_dat[ min_ind:max_ind, tag ] <- temp$new_dat
        
      }
      
    }
    
  }
  
  
  swap_results_sent[ , paste( 'rand', n, sep = '_' ) ] <- apply( rand_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max )
  
  sent_prop <- mean( apply( rand_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = na.max ), na.rm = T )
  
  rand_sent_vec <- c( rand_sent_vec, sent_prop )
  
}

## find the mean of each column (each column represents one randomization and the first column is the empirical)
means <- apply( swap_results_sent[ , -( 1:2 ) ], 2, mean, na.rm = T )

## pull out the mean of the empirical
emp <- means[ 1 ]

## save a vector with the mean of each randomization
rand_vec <- means[ -1 ]

## the proportion of random that are as extreme or more extreme than the empirical (aka the p-value)
mean( emp >= rand_vec )

dens <- density( rand_vec )

dens_x <- dens$x
dens_y <- dens$y

plot( dens_x, dens_y, ylim = range( dens_y ), xlim = range( c( dens_x, emp_sent_prop ) ), type = 'l', lwd = 2, xlab = 'Prop of time at leats one individual is awake', ylab = 'Probability density', bty = 'l' )

polygon( x = dens_x, y = dens_y, col = transp( 'black', 0.3) )

segments( x0 = emp_sent_prop, x1 = emp_sent_prop, y0 = 0, y1 = max( dens_y), lwd = 2, col = 'red', lty = 2)




plot( dens_x, dens_y, ylim = range( dens_y ), xlim = range( c( dens_x ) ), type = 'l', lwd = 2, xlab = 'Proportion of time at least one individual is awake', ylab = 'Probability density', bty = 'l', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )
axis(1, col = 'white', at = seq(0.865, 0.88, by = 0.005), labels = rep( "", 4 ) )
axis(2, col = 'white', at = seq(0, 150, by = 50), labels = rep("", 4))

polygon( x = dens_x, y = dens_y, col = transp( 'white', 0.3) )


plot( dens_x, dens_y, ylim = range( dens_y ), xlim = range( c( dens_x, emp_sent_prop ) ), type = 'l', lwd = 2, xlab = 'Proportion of time at least one individual is awake', ylab = 'Probability density', bty = 'l', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )
axis(1, col = 'white', at = seq(0.83, 0.88, by = 0.01), labels = rep( "", 6 ) )
axis(2, col = 'white', at = seq(0, 150, by = 50), labels = rep("", 4))


polygon( x = dens_x, y = dens_y, col = transp( 'white', 0.3) )

segments( x0 = emp_sent_prop, x1 = emp_sent_prop, y0 = 0, y1 = max( dens_y), lwd = 2, col = 'red', lty = 2)



##### Testing for group-level synchronization in sleep-wake state ###

## first comparing empirical data to randomizations produced by shifting each individual's sleep-wake data by a random amount within each given night

sync_dat <- total_dat[ total_dat$local_time > night_start | total_dat$local_time < night_end, ]

## remove the first night because there is only one individual on this night that actually has data
sync_dat <- sync_dat[ sync_dat$night != 1, ]

num_rands <- 1000

shift_results_sync <- data.frame( matrix( NA, nrow = nrow( sync_dat ), ncol = 3 + num_rands ) ) 

names( shift_results_sync ) <- c( 'local_timestamp', 'night', 'emp', paste( 'rand', 1:num_rands, sep = '_' ) )

shift_results_sync$local_timestamp <- sync_dat$local_timestamp

shift_results_sync$night <- sync_dat$night

shift_results_sync$emp <- apply( sync_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = doing_same )

emp_sync_prop <- mean( apply( sync_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = doing_same ), na.rm = T )


## now randomize

rand_sync_vec <- c()

nights <- unique( sync_dat$night )

for( n in 1:num_rands ){
  
  print( n )
  
  rand_dat <- sync_dat
  
  for( night in nights ){
    
    night_dat <- rand_dat[ rand_dat$night == night, ]
    
    for( tag in tag_names ){
      
      shift <- sample( 1: nrow( night_dat ), 1 )
      
      new_inds <- 1:nrow( night_dat ) + shift
      new_inds <- new_inds %% nrow( night_dat ) + 1
      
      
      night_dat[ , tag ] <- night_dat[ , tag ][ new_inds ]
      
    }
    
    rand_dat[ rand_dat$night == night, ] <- night_dat
    
  }
  
  shift_results_sync[ , paste( 'rand', n, sep = '_' ) ] <- apply( rand_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = doing_same )
  
  sync_prop <- mean( apply( rand_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = doing_same ), na.rm = T )
  
  rand_sync_vec <- c( rand_sync_vec, sync_prop )
  
}

## find the mean of each column (each column represents one randomization and the first column is the empirical)
means <- apply( shift_results_sync[ , -( 1:2 ) ], 2, mean, na.rm = T )

## pull out the mean of the empirical
emp <- means[ 1 ]

## save a vector with the mean of each randomization
rand_vec <- means[ -1 ]

## the proportion of random that are as extreme or more extreme than the empirical (aka the p-value)
mean( emp >= rand_vec )

dens <- density( rand_vec )

dens_x <- dens$x
dens_y <- dens$y

par( bg = "white" )

plot( dens_x, dens_y, ylim = range( dens_y ), xlim = range( c( dens_x, emp ) ), type = 'l', lwd = 2, xlab = 'Mean proportion of group synchronized', ylab = 'Probability density', bty = 'l' )

polygon( x = dens_x, y = dens_y, col = transp( 'black', 0.3) )

segments( x0 = emp, x1 = emp, y0 = 0, y1 = max( dens_y), lwd = 2, col = 'red', lty = 2)

## now compare the empirical group-level synchronization in sleep-wake state to a randomized dataset produced by permuting data across nights, rather than by shifting data wtihin a given night

sync_dat <- total_dat[ total_dat$local_time > night_start | total_dat$local_time < night_end, ]

## remove the first night because there is only one individual on this night that actually has data
sync_dat <- sync_dat[ sync_dat$night != 1, ]

num_rands <- 1000

swap_results_sync <- data.frame( matrix( NA, nrow = nrow( sync_dat ), ncol = 3 + num_rands ) ) 

names( swap_results_sync ) <- c( 'local_timestamp', 'night', 'emp', paste( 'rand', 1:num_rands, sep = '_' ) )

swap_results_sync$local_timestamp <- sync_dat$local_timestamp

swap_results_sync$night <- sync_dat$night

swap_results_sync$emp <- apply( sync_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = doing_same )

emp_sync_prop <- mean( apply( sync_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = doing_same ), na.rm = T )

start_date <- as.Date( min( d1$local_timestamp - 12*60*60 ) )

## now randomize

rand_sync_vec <- c()

for( n in 1:num_rands ){
  
  print( n )
  
  rand_dat <- sync_dat
  
  for( tag in tag_names ){
    
    if( sum( !is.na( sync_dat[ , tag ] ) ) != 0 ){
      
      min_ind <- min( which( !is.na( sync_dat[ , tag ] ) ) )
      max_ind <- max( which( !is.na( sync_dat[ , tag ] ) ) )
      
      sub <- sync_dat[ min_ind:max_ind, ]
      
      nights_avail <- sort( unique( sub$night ) )
      
      if( length( nights_avail ) != 1 ){
        
        new_nights <- sample( nights_avail )
        
        look_up <- data.frame( original = nights_avail, new = new_nights )
        
        temp <- merge( x = sub, y = look_up, by.x = 'night', by.y = 'original' )
        
        ## first add the new date based on the new randomized
        temp$new_timestamp <- start_date + temp$new - 1
        
        ## add a day to the rows that happen after midnight (when the date should increase by one)
        temp$new_timestamp[ temp$local_time < '12:00:00' ] <- temp$new_timestamp[ temp$local_time < '12:00:00' ] + 1
        
        ## add the time to our new randomized date
        temp$new_timestamp <- as.POSIXct( paste( temp$new_timestamp, temp$local_time ), tz = 'UTC' )
        
        ## reorder the vector of sleep to put the new timestamps in order
        temp$new_dat <- temp[ match( temp$local_timestamp, temp$new_timestamp ), tag ]
        
        rand_dat[ min_ind:max_ind, tag ] <- temp$new_dat
        
      }
      
    }
    
  }
  
  
  swap_results_sync[ , paste( 'rand', n, sep = '_' ) ] <- apply( rand_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = doing_same )
  
  sync_prop <- mean( apply( rand_dat[ , 2 : ( length( tag_names ) + 1 ) ], 1, FUN = doing_same ), na.rm = T )
  
  rand_sync_vec <- c( rand_sync_vec, sync_prop )
  
}


## find the mean of each column (each column represents one randomization and the first column is the empirical)
means <- apply( swap_results_sync[ , -( 1:2 ) ], 2, mean, na.rm = T )

## pull out the mean of the empirical
emp <- means[ 1 ]

## save a vector with the mean of each randomization
rand_vec <- means[ -1 ]

## the proportion of random that are as extreme or more extreme than the empirical (aka the p-value)
mean( emp >= rand_vec )

dens <- density( rand_vec )

dens_x <- dens$x
dens_y <- dens$y

par(bg = 'black')

plot( dens_x, dens_y, ylim = range( dens_y ), xlim = range( c( dens_x ) ), type = 'l', lwd = 2, xlab = 'Mean proportion of group synchronized', ylab = 'Probability density', bty = 'l', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )
axis(1, col = 'white', at = seq(0.8366, 0.8380, by = 0.0002), labels = rep( "", 8 ) )
axis(2, col = 'white', at = seq( 0, 2000, by = 500 ), labels = rep("", 5))

polygon( x = dens_x, y = dens_y, col = transp( 'white', 0.3) )


plot( dens_x, dens_y, ylim = range( dens_y ), xlim = range( c( dens_x, emp_sync_prop ) ), type = 'l', lwd = 2, xlab = 'Mean proportion of group synchronized', ylab = 'Probability density', bty = 'l', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )
axis(1, col = 'white', at = seq(0.837, 0.840, by = 0.001), labels = rep( "", 4 ) )
axis(2, col = 'white', at = seq( 0, 2000, by = 500 ), labels = rep("", 5))


polygon( x = dens_x, y = dens_y, col = transp( 'white', 0.3) )

segments( x0 = emp_sync_prop, x1 = emp_sync_prop, y0 = 0, y1 = max( dens_y), lwd = 2, col = 'red', lty = 2)


##### Do individuals synchronize more when in the same tree? #############

tag1 <- c()
tag2 <- c()
nights <- c()

tag_names <- as.character( unique( d1$tag ) )

for( a in 1: ( length( tag_names ) - 1 ) ){
  
  for( b in ( a + 1 ) : length( tag_names ) ){
    
    for( ni in unique( sent_dat$night ) ){
      
      id1 <- tag_names[ a ]
      
      id2 <- tag_names[ b ]
      
      night_dat <- sent_dat[ sent_dat$night == ni, ]
      
      if( sum( !is.na( night_dat[ , id1 ] ) ) != 0  & sum( !is.na( night_dat[ , id2 ] ) ) != 0  ){
        
        tag1 <- c( tag1, id1 )
        
        tag2 <- c( tag2, id2 )
        
        nights <- c( nights, ni )
      }
    }
  }
}

dy_sync <- data.frame( tag1 = tag1, tag2 = tag2, night = as.character( nights ) )

dy_sync$tag1 <- as.character( dy_sync$tag1 )

dy_sync$tag2 <- as.character( dy_sync$tag2 )

dy_sync$sync <- NA

for( i in 1:nrow( dy_sync ) ){
  
  id1 <- as.character( dy_sync$tag1[ i ] )
  
  id2 <- dy_sync$tag2[ i ]
  
  ni <- dy_sync$night[ i ]
  
  dy_sync$sync[ i ] <- mean( sent_dat[ sent_dat$night == ni , id1 ] == sent_dat[ sent_dat$night == ni , id2 ], na.rm = T )
  
}


final_sleep <- read.csv( "DATA/main_data/final_sleep.csv" )

trimmed <- final_sleep[ , c( 'tag', 'night', 'tree' ) ]

trimmed$tree <- as.character( trimmed$tree )

temp_dy_sync <- merge( x = dy_sync, y = trimmed, by.x = c( 'tag1', 'night' ), by.y = c( 'tag', 'night' ), all.x = T, all.y = F, sort = F ) 

names( temp_dy_sync )[ ncol( temp_dy_sync ) ] <- 'tree1'

final_dy_sync <- merge( x = temp_dy_sync, y = trimmed, by.x = c( 'tag2', 'night' ), by.y = c( 'tag', 'night' ), all.x = T, all.y = F, sort = F ) 

names( final_dy_sync )[ ncol( final_dy_sync ) ] <- 'tree2'

final_dy_sync$same_tree <- as.factor( as.character( as.numeric( final_dy_sync$tree1 == final_dy_sync$tree2 ) ) )

final_dy_sync$dy_name <- NA

for( i in 1: nrow( final_dy_sync ) ){
  
  final_dy_sync$dy_name[ i ] <- paste( sort( c( final_dy_sync[ i, "tag1" ], final_dy_sync[ i, "tag2" ] ) )[ 1 ], sort( c( final_dy_sync[ i, "tag1" ], final_dy_sync[ i, "tag2" ] ) )[ 2 ], sep = "_" )
  
}


hist( final_dy_sync$sync, main = 'Histogram of dyadic synchronization', xlab = 'Proportion of minutes dyad members exhibit same behavior')

final_dy_sync$sync_std <- normalize_func( final_dy_sync$sync )

options(mc.cores = parallel::detectCores() )

mod_dat <- final_dy_sync[ !is.na( final_dy_sync$same_tree ) & !is.na( final_dy_sync$sync ), ]

model <- brm( sync_std ~ same_tree + ( 1 | night ) + ( 1 | tag1 ) + ( 1 | tag2 ) + ( 1 | dy_name ), data = mod_dat )

saveRDS( model, "RESULTS/models/sync_mod.rds" )

model <- readRDS( "RESULTS/models/sync_mod.rds" )

summary( model )

tab_model( model )


model_unstd <- brm( sync ~ same_tree + ( 1 | night ) + ( 1 | tag1 ) + ( 1 | tag2 ) + ( 1 | dy_name ), data = mod_dat )

saveRDS( model_unstd, "RESULTS/models/sync_mod_unstd.rds" )

model_unstd <- readRDS( "RESULTS/models/sync_mod_unstd.rds" )

summary( model_unstd )

tab_model( model_unstd )

plot( conditional_effects( model_unstd, "same_tree", spaghetti = F ), line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )

final_dy_sync$same_tree_jitter <- jitter( as.numeric( as.character( final_dy_sync$same_tree ) ), factor = 2 ) + 1

plot( conditional_effects( model_unstd, "same_tree", spaghetti = F ), line_args = list( colour = "blue", size = 2 ), errorbar_args = list( colour = "blue", size = 1 ), cat_args = list( colour = "blue", size = 2 ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white", colour = "white"), panel.background = element_rect(fill = "white" ), axis.line = element_line(colour = "black"), axis.text = element_text(size = rel(0.8) ) ) )[[ 1 ]] + geom_point( data = final_dy_sync, aes( x = same_tree_jitter, y = sync ), alpha = 0.1, inherit.aes = F )

      