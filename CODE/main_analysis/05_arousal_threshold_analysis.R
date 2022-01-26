
library(brms)
library( stats )
library( sjPlot )

na_sum <- function( x ){
  
  if( sum( is.na( x ) ) == length( x ) ){
    
    return( NA )
    
  }else{
    
    return( sum( x, na.rm = T ) )
  }
}

load( "DATA/main_data/sleep_analysis_pub_code_pre_mods.RData" )

full_dat_trim <- full_dat[ !is.na( full_dat$sleep_bouts ) & full_dat$night > 1 & full_dat$night < 15 & ( full_dat$local_time > "21:00:00" | full_dat$local_time < "05:00:00" ), ]

full_tree_dat <- merge( x = full_dat_trim, y = finalest_sleep[ , c( 'tag', 'night', 'tree' ) ], by = c( 'tag', 'night' ), all.x = T, all.y = F, sort = F )

full_tree_dat <- full_tree_dat[ order( full_tree_dat$local_timestamp ), ]
full_tree_dat <- full_tree_dat[ order( full_tree_dat$tag ), ]

full_tree_dat$others_in_tree_awake <- NA
full_tree_dat$awake_prev_minute <- NA

for( i in 1:nrow( full_tree_dat ) ){
  
  if( i %% 1000 == 0) print( i / nrow( full_tree_dat ) )
  
  others_awake <- ( 1 - full_tree_dat[ full_tree_dat$tag != full_tree_dat$tag[ i ] & full_tree_dat$local_timestamp == ( full_tree_dat$local_timestamp[ i ] - 60 ) & full_tree_dat$tree == full_tree_dat$tree[ i ], 'sleep_bouts' ] )
  
  if( length( others_awake ) > 0 ){
    
    full_tree_dat$others_in_tree_awake[ i ] <- ifelse( na_sum( others_awake ) > 0, 1, 0 )
    
  }
  
  if( length( 1 - ( full_tree_dat[ full_tree_dat$tag == full_tree_dat$tag[ i ] & full_tree_dat$local_timestamp == ( full_tree_dat$local_timestamp[ i ] - 60 ), 'sleep_bouts' ] ) ) > 0 ){
    
    full_tree_dat$awake_prev_minute[ i ] <- ( 1 - ( full_tree_dat[ full_tree_dat$tag == full_tree_dat$tag[ i ] & full_tree_dat$local_timestamp == ( full_tree_dat$local_timestamp[ i ] - 60 ), 'sleep_bouts' ] ) )
    
  }
  
}


finalest_sleep$frag_wake_bouts[ is.na( finalest_sleep$TST ) ] <- NA

finalest_sleep$SPT <- finalest_sleep$SPT / 60

finalest_sleep$fragG2 <- finalest_sleep$frag_wake_bouts / ( finalest_sleep$TST / 60 ) 

total_mod_dat_temp <- lags( sleep_per_df = finalest_sleep, vars = c( 'TST', 'fragG2' ), n = 1 )

ave_TSTs <- aggregate( total_mod_dat_temp[ , c( 'TST', 'fragG2' ) ], by = list( total_mod_dat_temp$tag ), FUN = mean, na.rm = T )

names( ave_TSTs ) <- c( 'tag', 'ave_TST', 'ave_fragG2' )

total_mod_dat <- merge( x = total_mod_dat_temp, y = ave_TSTs, by = 'tag', all.x = T, all.y = F, sort = F )

total_mod_dat$prev_TST_diff <- total_mod_dat$prev_TST_cum - total_mod_dat$ave_TST

total_mod_dat$prev_TST_diff_std <- normalize_func( total_mod_dat$prev_TST_diff )

total_mod_dat$prev_fragG2_diff <- total_mod_dat$prev_fragG2_cum - total_mod_dat$ave_fragG2 

total_mod_dat$prev_fragG2_diff_std <- normalize_func( total_mod_dat$prev_fragG2_diff )


full_tree_TST <- merge( x = full_tree_dat, y = total_mod_dat[ , c( 'tag', 'night', 'prev_TST_diff', 'prev_TST_diff_std', 'prev_fragG2_diff', 'prev_fragG2_diff_std' ) ], by = c( 'tag', 'night' ), all.x = T, all.y = F, sort = F )

full_tree_TST$curr_awake <- 1 - full_tree_TST$sleep_bouts

full_tree_TST <- full_tree_TST[ order( full_tree_TST$local_timestamp ), ]
full_tree_TST <- full_tree_TST[ order( full_tree_TST$tag ), ]

#write.csv( full_tree_TST, "DATA/main_data/full_tree_TST.csv", row.names = F )

full_tree_TST <- read.csv( "DATA/main_data/full_tree_TST.csv" )

full_tree_TST$curr_awake <- as.factor( full_tree_TST$curr_awake )
full_tree_TST$others_in_tree_awake <- as.factor( full_tree_TST$others_in_tree_awake )


was_asleep_dat <- full_tree_TST[ which( full_tree_TST$awake_prev_minute == 0 ), ] 

options( mc.cores = parallel::detectCores() )

priors <- c( prior(normal(0,2), class = "b", coef = others_in_tree_awake1),
             
             prior(normal(0,2), class = "b", coef = prev_TST_diff_std ),
             
             prior(normal(0,2), class = "b", coef = prev_fragG2_diff_std ),
             
             prior(normal(0,2), class = "b", coef = others_in_tree_awake1:prev_TST_diff_std),
             
             prior(normal(0,2), class = "b", coef = others_in_tree_awake1:prev_fragG2_diff_std)
             
)


arousal_thresh_mod <- brm( curr_awake ~ others_in_tree_awake + prev_TST_diff_std + prev_fragG2_diff_std +  others_in_tree_awake:prev_TST_diff_std + others_in_tree_awake:prev_fragG2_diff_std + ( 1 | night ) + ( 1 | tag ), data = was_asleep_dat, family= "bernoulli", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( object = full_mod, file = "RESULTS/models/arousal_thresh_mod.rds" )

arousal_thresh_mod <- readRDS( "RESULTS/models/arousal_thresh_mod.rds" )

summary( arousal_thresh_mod )

tab_model( arousal_thresh_mod )

