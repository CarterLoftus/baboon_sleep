
library(readr)
library(lubridate)
library(plyr)
library(stringr)
library(rethinking)
library(igraph)
library(data.table)
library(cluster)

#### Functions
babcolor<-function(IDvec){
  outputVec<-c()
  for(i in 1:length(IDvec)){
    if( IDvec[i] %in% c("2427", "2434", "2457") ){
      outputVec[i] <- 'blue'
    }else{
      if( IDvec[i] %in% c("2426", "2432", "2433", "2436", "2443", "2452") ){
        outputVec[i] <- 'skyblue'
      }else{
        if( IDvec[i] %in% c("2448", "2454") ){
          outputVec[i] <- 'grey'
        }else{
          if(IDvec[i] %in% c("2430", "2439", "2446", "2447", "2449", "2451", "2453", "2456", "2458", "2460") ){
            outputVec[i] <- 'red'
          }else{
            outputVec[i] <- 'pink'
          }
        }
      }
    }
  }
  return(outputVec)
}


dy_acc <- function(vect, win_size = 7){
  
  pad_size <- win_size/2 - 0.5
  
  padded <- unlist( c(rep(NA, pad_size), vect, rep(NA, pad_size)) )
  acc_vec <- rep(NA, length = length( vect ) )
  
  ## sliding window
  for(i in 1:length(vect)){
    win <- padded[i:(i+(2*pad_size))] ## subset the window
    m_ave <- mean( win, na.rm = T ) ## take the average over the window
    acc_comp <- vect[ i ] - m_ave ## finds the difference between the static component (mean) and the actual value. This is the dynamic component of the acceleration at this time point
    acc_vec[i] <- acc_comp 
  }
  
  return( unlist( acc_vec) )
}

#### End of functions

d1 <- fread("DATA/main_data/all_burst_acc.csv")

d1 <- as.data.frame(d1)

names(d1) <- c('tag','timestamp','eobs_accelerations_raw')

d1$timestamp <- as.POSIXct( d1$timestamp, tz = 'UTC' )

## What data do we have?
plot(as.numeric(as.factor(d1$tag))~d1$timestamp,cex=0.3,pch=16,main="Overnight acceleromter bursts",xlab="",xaxt='n',yaxt='n',ylab="ID",col=babcolor(d1$tag))
axis(2,at=1:length(unique(d1$tag)),labels=sort(unique(d1$tag)),las=1,cex=0.3)
axis.POSIXct(1,at=seq(min(d1$timestamp),max(d1$timestamp),by="1 day"),las=2)

d1$night <- as.Date( d1$timestamp - 12*60*60) 

tag_names <- as.character( unique( d1$tag ) )

d1$time <- str_split_fixed(d1$timestamp, " ", 2)[,2]

d2 <- as.data.frame(str_split(d1$eobs_accelerations_raw, " ", simplify = T))

for(i in 1:ncol(d2)){
  d2[,i] <- as.numeric(as.character((d2[,i])))
}

names(d2) <- paste(rep(c("x","y","z"),ncol(d2)/3),rep(1:(ncol(d2)/3), each = 3), sep = '')

d2$timestamp <- d1$timestamp
d2$tag <- d1$tag

d2[!complete.cases(d2),]

inds <- complete.cases(d2)
d1 <- d1[ inds ,]
d2 <- d2[ inds ,]

names(d2) <- c( paste( rep(c("x","y","z"), ncol(d2)/3), rep(1:(ncol(d2)/3), each = 3), sep = ''), 'timestamp')

x_d <- d2[,grepl('x',names(d2))]
head(x_d)

y_d <- d2[,grepl('y',names(d2))]
head(y_d)

z_d <- d2[,grepl('z',names(d2))]
head(z_d)

d1$ave_vedba <- apply( sqrt( apply( x_d, 1, FUN = function(x) abs( dy_acc( x ) ) )**2 + apply( y_d, 1, FUN = function(x) abs( dy_acc( x ) ) )**2 + apply( z_d, 1, FUN = function(x) abs( dy_acc( x ) ) )**2) , 2, FUN = mean )
d1$log_vedba <- log( d1$ave_vedba )

write.csv(d1, "DATA/main_data/full_night_and_day_data.csv", row.names = F)

