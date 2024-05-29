# Laura Tibbs Cortes
# Aug 21, 2020

# Sub-functions for running CERES, based on Xianran's code (sub_funcs_20200707.r) but also based on my own previous
# work with running my vitamins through CERES

# Updated July 14 2020 to allow for more environmental indices

# Make ptt_ptr file using function from source:
 # this function calls sub-functions to calculate the day length for each day the plants were growing, up until our searching_daps limit,
 # then, find the weather from the 6 nearest weather stations for those dates and extract average high and low daily temperature
 # and use to calculate GDDs and other environmental parameters
Compile_PTT_PTR_local_GHCN <-  function(exp_dir, env_meta_info, exp_s_year, exp_e_year, searching_daps,ptt_ptr_file, t_base, t_max1, t_max2) {

  # make directories:
 sp_env_dir <- paste(exp_dir, 'envs/', sep = '');     if (!dir.exists(sp_env_dir))  { dir.create(sp_env_dir)};
 # GHCN is global historical climatology network
 sp_ghcn_dir <- paste(sp_env_dir, 'ghcn/', sep = ''); if (!dir.exists(sp_ghcn_dir)) { dir.create(sp_ghcn_dir)};
 sp_navy_dir <- paste(sp_env_dir, 'Geoshpere/', sep = ''); if (!dir.exists(sp_navy_dir)) { dir.create(sp_navy_dir)};
 
 # set latitude and longitude range of data
 lat_range <- range(env_meta_info$lat, na.rm = T); lon_range <- range(env_meta_info$lon, na.rm = T);
 
 # pull the weather stations with info relevant to your environments
 local_ghcn_st_file <- paste('all_ghcn_stations', sep = '');  
 local_ghcn_target_st_file <- paste(sp_ghcn_dir, '0target_ghcn_stations', sep = '');
 if (!file.exists(local_ghcn_target_st_file)) {
   
   # if you don't have a pre-made list of all ghcn stations, then get that info
   if (!file.exists(local_ghcn_st_file)) { ghcn_all_sts <- ghcnd_stations() # this function returns data for all GHCN weather stations available
   } else { ghcn_all_sts <- read.csv(local_ghcn_st_file);}
   
   # filter the ghcn info to only include those within your begin and end year, latitude, and longitude
   # only keep the max and min temperature
   STS_ghcn <- dplyr::filter(ghcn_all_sts, first_year <= exp_s_year & last_year >= exp_e_year 
                                             & between(latitude, lat_range[1] - 2, lat_range[2] + 2) 
                                             & between(longitude,lon_range[1] - 2, lon_range[2] + 2)
                                             & element %in% c("TMAX", "TMIN")
                                             );
                                           
   write.csv(STS_ghcn, local_ghcn_target_st_file  ) # write out file
 }
 STS_ghcn <- read.csv(local_ghcn_target_st_file); 
 
 st_ids <- c();
 for (e_i in 1:nrow(env_meta_info)) { 
   
   # for each environment, pull the metadata:
  env_code <- env_meta_info$env_code[e_i];
  field_lat <- env_meta_info$lat[e_i]; field_lon <- env_meta_info$lon[e_i]; 
  planting_date <- env_meta_info$PlantingDate[e_i];
  planting_year <- year(planting_date);

  # set the dates to inspect for this environment (from planting until your set end date)
  DAPs <- seq(as.Date(planting_date), length.out = searching_daps, by = "day"); # make a list of all the dates this field was growing, up until the specified searching_daps limit
  ending_date <- DAPs[searching_daps]; 
  ending_year <- year(ending_date);
  
  # Make file containing daylength if it does not already exist
  # Using DayLength_from_Equation function
  local_dl_file <- paste(sp_navy_dir, env_code, '_DL_', searching_daps, 'DAPs', sep = '')
  if (!file.exists(local_dl_file)) {
    DayLength_from_Equation(planting_year, ending_year, 
                                                            field_lat, field_lon, 
                                                            local_dl_file, DAPs ) 
    }

  # find closest 6 weather stations within 50 km of the field location
  field_sts <- meteo_distance(STS_ghcn, lat = field_lat, lon = field_lon, radius = 50, limit = 6);
  
  # save the IDs of those weather stations
  st_ids <- append(st_ids, unique(as.vector(field_sts$id)));
 }
 
 # keep only unique weather stations (in case environments near one another)
 st_ids <- unique(st_ids);
 
 # work through all years of the overall experiment to read in weather data:
 for (e_y in 1:(exp_e_year - exp_s_year + 1)) {
   gz_ghcn_file <- paste(ghcn_year_dir, exp_s_year + e_y - 1, '.csv.gz', sep = ''); # set weather file name

   # read in weather data, and filter to keep only the weather stations we want AND only keep the max and min temp data from them
   gz_ghcn <- fread(gz_ghcn_file, select = c(1:4),showProgress = FALSE);
   gz_ghcn <- gz_ghcn[gz_ghcn$V1 %in% st_ids & gz_ghcn$V3 %in% c('TMAX', 'TMIN', 'PRCP')][,1:4] # keep precip plus temp 
   
#    # definitions of weather units:
#               PRCP = Precipitation (tenths of mm)
#    	   SNOW = Snowfall (mm)
# 	   SNWD = Snow depth (mm)
#            TMAX = Maximum temperature (tenths of degrees C)
#            TMIN = Minimum temperature (tenths of degrees C)
# 	   
# 	   The other elements are:
# 	   
# 	   ACMC = Average cloudiness midnight to midnight from 30-second 
# 	          ceilometer data (percent)
# 	   ACMH = Average cloudiness midnight to midnight from 
# 	          manual observations (percent)
#            ACSC = Average cloudiness sunrise to sunset from 30-second 
# 	          ceilometer data (percent)
# 	   ACSH = Average cloudiness sunrise to sunset from manual 
# 	          observations (percent)
#            AWDR = Average daily wind direction (degrees)
# 	   AWND = Average daily wind speed (tenths of meters per second)
# 	   DAEV = Number of days included in the multiday evaporation
# 	          total (MDEV)
# 	   DAPR = Number of days included in the multiday precipiation 
# 	          total (MDPR)
#            DASF = Number of days included in the multiday snowfall 
# 	          total (MDSF)		  
# 	   DATN = Number of days included in the multiday minimum temperature 
# 	         (MDTN)
# 	   DATX = Number of days included in the multiday maximum temperature 
# 	          (MDTX)
#            DAWM = Number of days included in the multiday wind movement
# 	          (MDWM)
# 	   DWPR = Number of days with non-zero precipitation included in 
# 	          multiday precipitation total (MDPR)
# 	   EVAP = Evaporation of water from evaporation pan (tenths of mm)
# 	   FMTM = Time of fastest mile or fastest 1-minute wind 
# 	          (hours and minutes, i.e., HHMM)
# 	   FRGB = Base of frozen ground layer (cm)
# 	   FRGT = Top of frozen ground layer (cm)
# 	   FRTH = Thickness of frozen ground layer (cm)
# 	   GAHT = Difference between river and gauge height (cm)
# 	   MDEV = Multiday evaporation total (tenths of mm; use with DAEV)
# 	   MDPR = Multiday precipitation total (tenths of mm; use with DAPR and 
# 	          DWPR, if available)
# 	   MDSF = Multiday snowfall total 
# 	   MDTN = Multiday minimum temperature (tenths of degrees C; use with 
# 	          DATN)
# 	   MDTX = Multiday maximum temperature (tenths of degress C; use with 
# 	          DATX)
# 	   MDWM = Multiday wind movement (km)
#            MNPN = Daily minimum temperature of water in an evaporation pan 
# 	         (tenths of degrees C)
#            MXPN = Daily maximum temperature of water in an evaporation pan 
# 	         (tenths of degrees C)
# 	   PGTM = Peak gust time (hours and minutes, i.e., HHMM)
# 	   PSUN = Daily percent of possible sunshine (percent)
# 	   SN*# = Minimum soil temperature (tenths of degrees C)
# 	          where * corresponds to a code
# 	          for ground cover and # corresponds to a code for soil 
# 		  depth.  
# 		  
# 		  Ground cover codes include the following:
# 		  0 = unknown
# 		  1 = grass
# 		  2 = fallow
# 		  3 = bare ground
# 		  4 = brome grass
# 		  5 = sod
# 		  6 = straw multch
# 		  7 = grass muck
# 		  8 = bare muck
# 		  
# 		  Depth codes include the following:
# 		  1 = 5 cm
# 		  2 = 10 cm
# 		  3 = 20 cm
# 		  4 = 50 cm
# 		  5 = 100 cm
# 		  6 = 150 cm
# 		  7 = 180 cm
# 		  
# 	   SX*# = Maximum soil temperature (tenths of degrees C) 
# 	          where * corresponds to a code for ground cover 
# 		  and # corresponds to a code for soil depth. 
# 		  See SN*# for ground cover and depth codes. 
#            TAVG = Average temperature (tenths of degrees C)
# 	          [Note that TAVG from source 'S' corresponds
# 		   to an average for the period ending at
# 		   2400 UTC rather than local midnight]
#            THIC = Thickness of ice on water (tenths of mm)	
#  	   TOBS = Temperature at the time of observation (tenths of degrees C)
# 	   TSUN = Daily total sunshine (minutes)
# 	   WDF1 = Direction of fastest 1-minute wind (degrees)
# 	   WDF2 = Direction of fastest 2-minute wind (degrees)
# 	   WDF5 = Direction of fastest 5-second wind (degrees)
# 	   WDFG = Direction of peak wind gust (degrees)
# 	   WDFI = Direction of highest instantaneous wind (degrees)
# 	   WDFM = Fastest mile wind direction (degrees)
#            WDMV = 24-hour wind movement (km)	   
#            WESD = Water equivalent of snow on the ground (tenths of mm)
# 	   WESF = Water equivalent of snowfall (tenths of mm)
# 	   WSF1 = Fastest 1-minute wind speed (tenths of meters per second)
# 	   WSF2 = Fastest 2-minute wind speed (tenths of meters per second)
# 	   WSF5 = Fastest 5-second wind speed (tenths of meters per second)
# 	   WSFG = Peak gust wind speed (tenths of meters per second)
# 	   WSFI = Highest instantaneous wind speed (tenths of meters per second)
# 	   WSFM = Fastest mile wind speed (tenths of meters per second)
# 	   WT** = Weather Type where ** has one of the following values:
# 	   
#                   01 = Fog, ice fog, or freezing fog (may include heavy fog)
#                   02 = Heavy fog or heaving freezing fog (not always 
# 		       distinquished from fog)
#                   03 = Thunder
#                   04 = Ice pellets, sleet, snow pellets, or small hail 
#                   05 = Hail (may include small hail)
#                   06 = Glaze or rime 
#                   07 = Dust, volcanic ash, blowing dust, blowing sand, or 
# 		       blowing obstruction
#                   08 = Smoke or haze 
#                   09 = Blowing or drifting snow
#                   10 = Tornado, waterspout, or funnel cloud 
#                   11 = High or damaging winds
#                   12 = Blowing spray
#                   13 = Mist
#                   14 = Drizzle
#                   15 = Freezing drizzle 
#                   16 = Rain (may include freezing rain, drizzle, and
# 		       freezing drizzle) 
#                   17 = Freezing rain 
#                   18 = Snow, snow pellets, snow grains, or ice crystals
#                   19 = Unknown source of precipitation 
#                   21 = Ground fog 
#                   22 = Ice fog or freezing fog
# 		  
#             WV** = Weather in the Vicinity where ** has one of the following 
# 	           values:
# 		   
# 		   01 = Fog, ice fog, or freezing fog (may include heavy fog)
# 		   03 = Thunder
# 		   07 = Ash, dust, sand, or other blowing obstruction
# 		   18 = Snow or ice crystals
# 		   20 = Rain or snow shower
   
   # save this data
   if (e_y == 1) {ghcn_TM <- gz_ghcn} else {ghcn_TM <- rbind(ghcn_TM, gz_ghcn)}
 }
  # ghcn_TM <- ghcn_TM[,V4:=ifelse(V4==9999|V4==-9999, NA, round(32 + 9 * V4/50, 3))]; ## convert to degrees F; noaa record is tenths of degrees C
 ghcn_TM <- ghcn_TM %>%
   mutate(V4=ifelse(V3=="PRCP", V4, ifelse(V4==9999|V4==-9999, NA, round(32 + 9 * V4/50, 3))))
 ghcn_TM[, V2 := ymd(V2)]; # convert charcter date to ymd format
  
 PTT_PTR <- matrix(ncol = 8, nrow = 0); 
 for (e_i in 1:nrow(env_meta_info)) { 
   
  # for each environment, pull the metadata:
  env_code <- env_meta_info$env_code[e_i];
  field_lat <- env_meta_info$lat[e_i]; field_lon <- env_meta_info$lon[e_i]; 
  planting_date <- env_meta_info$PlantingDate[e_i];
  planting_year <- year(planting_date);
  
  # set the dates to inspect for this environment (from planting until your set end date); see above for details
  DAPs <- seq(as.Date(planting_date), length.out = searching_daps, by = "day");
  ending_date <- DAPs[searching_daps]; 
  ending_year <- year(ending_date);
  
  # read in day length file made earlier
  local_dl_file <- paste(sp_navy_dir, env_code, '_DL_', searching_daps, 'DAPs', sep = '')
  DL <- read.table(local_dl_file, header = T, sep = "\t", stringsAsFactors = F);

  # set name of temperature file and make it
  local_tm_file <- paste(sp_ghcn_dir, env_code, '_TM_', searching_daps, 'DAPs', sep = '' );

  # this subfunction calculates the daily max and min temperatures, averaged over the nearest weather stations, for each environment
  # if (!file.exists(local_tm_file)) {

    TM_from_local_NOAA(ghcn_TM %>% filter(V3 %in% c("TMAX", "TMIN")), 
                                                       field_lat, field_lon, 
                                                       local_tm_file, STS_ghcn, DAPs)
    # } ;

  TM_0 <- read.table(local_tm_file, header = T, sep = "\t", stringsAsFactors = F); # read in temperature file created by above function

  PRCP_0 <- PRCP_from_local_NOAA(ghcn_TM %>% filter(V3 == "PRCP"), field_lat, field_lon, 
                                                       STS_ghcn, DAPs)
  PRCP_0$date <- as.character(PRCP_0$date)
  
  # merge day length with temperature data AND precip--
  DL_TM <- merge(DL, TM_0, all.x = T) ;
  DL_TM <- left_join(DL_TM, PRCP_0, by="date")
   
  # fill in missing temperatures
  DL_TM$TMAX <- Fill_Missning_TM(DL_TM$TMAX, env);
  DL_TM$TMIN <- Fill_Missning_TM(DL_TM$TMIN, env);
  DL_TM$precip <- Fill_Missning_TM(DL_TM$precip, env);
  
  # save max and min temperatures as vectors
  Tmax <- DL_TM$TMAX; Tmin <- DL_TM$TMIN; precip.vector <- DL_TM$precip
  
  # calculate GDD with sub-function:
  GDDs <- Adjusting_TM_4_GDD(Tmax, Tmin, Haun_threshold, t_base, t_max1, t_max2);
  
  # Calculate other environmental parameters:
  PTTs <- round(GDDs * DL_TM$DL, 4); # PTT=GDDxDL
  PTRs <- round(GDDs / DL_TM$DL, 4); # PTR=GDD/DL
  PTD1s <- round((DL_TM$TMAX - DL_TM$TMIN) * DL_TM$DL, 4); #PTD1=(TMAX-TMIN)xDL
  PTD2s <- round((DL_TM$TMAX - DL_TM$TMIN) / DL_TM$DL, 4); #PTD2=(TMAX-TMIN)/DL
  PTSs <- round(((DL_TM$TMAX^2) - (DL_TM$TMIN^2)) * (DL_TM$DL^2), 4); #PTS=(TMAX^2 - TMIN^2)*(DL^2)
  avg.temp <- (Tmax+Tmin)/2
  DTR <- Tmax-Tmin
  
  
  env_codes <- rep(env_code, searching_daps); # make vector of the environmental code of the length of number of days of interest
  
  # make data frame to hold this data, and output
  t_df <- data.frame(env_code = env_codes, date = DAPs, TMAX = DL_TM$TMAX, TMIN = DL_TM$TMIN, DL = DL_TM$DL, GDD = GDDs, PTT = PTTs, 
                     PTR = PTRs, PTD1 = PTD1s, PTD2 = PTD2s, PTS = PTSs, PRECIP=precip.vector, avg.temp=avg.temp)
  PTT_PTR <- rbind(PTT_PTR, t_df);
 }
 
 write.table(PTT_PTR, file = ptt_ptr_file, sep = "\t", quote = F, row.name = F)
 invisible(gc());
}

# Expand to include more environment parameters! Don't want to break the old one in the process.
Compile_PTT_PTR_local_GHCN_extend <-  function(exp_dir, env_meta_info, exp_s_year, exp_e_year, searching_daps,ptt_ptr_file, t_base, t_max1, t_max2, petfile) {

  # make directories:
 sp_env_dir <- paste(exp_dir, 'envs/', sep = '');     if (!dir.exists(sp_env_dir))  { dir.create(sp_env_dir)};
 # GHCN is global historical climatology network
 sp_ghcn_dir <- paste(sp_env_dir, 'ghcn/', sep = ''); if (!dir.exists(sp_ghcn_dir)) { dir.create(sp_ghcn_dir)};
 sp_navy_dir <- paste(sp_env_dir, 'Geoshpere/', sep = ''); if (!dir.exists(sp_navy_dir)) { dir.create(sp_navy_dir)};
 
 # set latitude and longitude range of data
 lat_range <- range(env_meta_info$lat, na.rm = T); lon_range <- range(env_meta_info$lon, na.rm = T);
 
 # pull the weather stations with info relevant to your environments
 local_ghcn_st_file <- paste('all_ghcn_stations', sep = '');  
 local_ghcn_target_st_file <- paste(sp_ghcn_dir, '0target_ghcn_stations', sep = '');
 if (!file.exists(local_ghcn_target_st_file)) {
   
   # if you don't have a pre-made list of all ghcn stations, then get that info
   if (!file.exists(local_ghcn_st_file)) { 
     ghcn_all_sts <- ghcnd_stations() # this function returns data for all GHCN weather stations available
   } else { 
     ghcn_all_sts <- read.csv(local_ghcn_st_file);
     }
   
   # filter the ghcn info to only include those within your begin and end year, latitude, and longitude
   # only keep the max and min temperature
   STS_ghcn <- dplyr::filter(ghcn_all_sts, first_year <= exp_s_year & last_year >= exp_e_year 
                                             & between(latitude, lat_range[1] - 2, lat_range[2] + 2) 
                                             & between(longitude,lon_range[1] - 2, lon_range[2] + 2)
                                             & element %in% c("TMAX", "TMIN", "PRCP")
                                             );
                                           
   write.csv(STS_ghcn, local_ghcn_target_st_file  ) # write out file
 }
 STS_ghcn <- read.csv(local_ghcn_target_st_file); 
 
 st_ids <- c();
 for (e_i in 1:nrow(env_meta_info)) { 
   
   # for each environment, pull the metadata:
  env_code <- env_meta_info$env_code[e_i];
  field_lat <- env_meta_info$lat[e_i]; field_lon <- env_meta_info$lon[e_i]; 
  planting_date <- env_meta_info$PlantingDate[e_i];
  planting_year <- year(planting_date);

  # set the dates to inspect for this environment (from planting until your set end date)
  DAPs <- seq(as.Date(planting_date), length.out = searching_daps, by = "day"); # make a list of all the dates this field was growing, up until the specified searching_daps limit
  ending_date <- DAPs[searching_daps]; 
  ending_year <- year(ending_date);
  
  # Make file containing daylength if it does not already exist
  # Using DayLength_from_Equation function
  local_dl_file <- paste(sp_navy_dir, env_code, '_DL_', searching_daps, 'DAPs', sep = '')
  # if (!file.exists(local_dl_file)) {
    DayLength_from_Equation(planting_year, ending_year, 
                                                            field_lat, field_lon, 
                                                            local_dl_file, DAPs ) 
    # }

  # find closest 6 weather stations within 50 km of the field location
  field_sts <- meteo_distance(STS_ghcn, lat = field_lat, lon = field_lon, radius = 50, limit = 6);
  
  # save the IDs of those weather stations
  st_ids <- append(st_ids, unique(as.vector(field_sts$id)));
 }
 
 # keep only unique weather stations (in case environments near one another)
 st_ids <- unique(st_ids);
 
 # work through all years of the overall experiment to read in weather data:
 for (e_y in 1:(exp_e_year - exp_s_year + 1)) {
   gz_ghcn_file <- paste(ghcn_year_dir, exp_s_year + e_y - 1, '.csv.gz', sep = ''); # set weather file name

   # read in weather data, and filter to keep only the weather stations we want AND only keep the max and min temp data from them
   gz_ghcn <- fread(gz_ghcn_file, select = c(1:4),showProgress = FALSE);
   gz_ghcn <- gz_ghcn[gz_ghcn$V1 %in% st_ids & gz_ghcn$V3 %in% c('TMAX', 'TMIN', 'PRCP')][,1:4] # keep precip plus temp 
   
#    # definitions of weather units:
#               PRCP = Precipitation (tenths of mm)
#    	   SNOW = Snowfall (mm)
# 	   SNWD = Snow depth (mm)
#            TMAX = Maximum temperature (tenths of degrees C)
#            TMIN = Minimum temperature (tenths of degrees C)
# 	   
# 	   The other elements are:
# 	   
# 	   ACMC = Average cloudiness midnight to midnight from 30-second 
# 	          ceilometer data (percent)
# 	   ACMH = Average cloudiness midnight to midnight from 
# 	          manual observations (percent)
#            ACSC = Average cloudiness sunrise to sunset from 30-second 
# 	          ceilometer data (percent)
# 	   ACSH = Average cloudiness sunrise to sunset from manual 
# 	          observations (percent)
#            AWDR = Average daily wind direction (degrees)
# 	   AWND = Average daily wind speed (tenths of meters per second)
# 	   DAEV = Number of days included in the multiday evaporation
# 	          total (MDEV)
# 	   DAPR = Number of days included in the multiday precipiation 
# 	          total (MDPR)
#            DASF = Number of days included in the multiday snowfall 
# 	          total (MDSF)		  
# 	   DATN = Number of days included in the multiday minimum temperature 
# 	         (MDTN)
# 	   DATX = Number of days included in the multiday maximum temperature 
# 	          (MDTX)
#            DAWM = Number of days included in the multiday wind movement
# 	          (MDWM)
# 	   DWPR = Number of days with non-zero precipitation included in 
# 	          multiday precipitation total (MDPR)
# 	   EVAP = Evaporation of water from evaporation pan (tenths of mm)
# 	   FMTM = Time of fastest mile or fastest 1-minute wind 
# 	          (hours and minutes, i.e., HHMM)
# 	   FRGB = Base of frozen ground layer (cm)
# 	   FRGT = Top of frozen ground layer (cm)
# 	   FRTH = Thickness of frozen ground layer (cm)
# 	   GAHT = Difference between river and gauge height (cm)
# 	   MDEV = Multiday evaporation total (tenths of mm; use with DAEV)
# 	   MDPR = Multiday precipitation total (tenths of mm; use with DAPR and 
# 	          DWPR, if available)
# 	   MDSF = Multiday snowfall total 
# 	   MDTN = Multiday minimum temperature (tenths of degrees C; use with 
# 	          DATN)
# 	   MDTX = Multiday maximum temperature (tenths of degress C; use with 
# 	          DATX)
# 	   MDWM = Multiday wind movement (km)
#            MNPN = Daily minimum temperature of water in an evaporation pan 
# 	         (tenths of degrees C)
#            MXPN = Daily maximum temperature of water in an evaporation pan 
# 	         (tenths of degrees C)
# 	   PGTM = Peak gust time (hours and minutes, i.e., HHMM)
# 	   PSUN = Daily percent of possible sunshine (percent)
# 	   SN*# = Minimum soil temperature (tenths of degrees C)
# 	          where * corresponds to a code
# 	          for ground cover and # corresponds to a code for soil 
# 		  depth.  
# 		  
# 		  Ground cover codes include the following:
# 		  0 = unknown
# 		  1 = grass
# 		  2 = fallow
# 		  3 = bare ground
# 		  4 = brome grass
# 		  5 = sod
# 		  6 = straw multch
# 		  7 = grass muck
# 		  8 = bare muck
# 		  
# 		  Depth codes include the following:
# 		  1 = 5 cm
# 		  2 = 10 cm
# 		  3 = 20 cm
# 		  4 = 50 cm
# 		  5 = 100 cm
# 		  6 = 150 cm
# 		  7 = 180 cm
# 		  
# 	   SX*# = Maximum soil temperature (tenths of degrees C) 
# 	          where * corresponds to a code for ground cover 
# 		  and # corresponds to a code for soil depth. 
# 		  See SN*# for ground cover and depth codes. 
#            TAVG = Average temperature (tenths of degrees C)
# 	          [Note that TAVG from source 'S' corresponds
# 		   to an average for the period ending at
# 		   2400 UTC rather than local midnight]
#            THIC = Thickness of ice on water (tenths of mm)	
#  	   TOBS = Temperature at the time of observation (tenths of degrees C)
# 	   TSUN = Daily total sunshine (minutes)
# 	   WDF1 = Direction of fastest 1-minute wind (degrees)
# 	   WDF2 = Direction of fastest 2-minute wind (degrees)
# 	   WDF5 = Direction of fastest 5-second wind (degrees)
# 	   WDFG = Direction of peak wind gust (degrees)
# 	   WDFI = Direction of highest instantaneous wind (degrees)
# 	   WDFM = Fastest mile wind direction (degrees)
#            WDMV = 24-hour wind movement (km)	   
#            WESD = Water equivalent of snow on the ground (tenths of mm)
# 	   WESF = Water equivalent of snowfall (tenths of mm)
# 	   WSF1 = Fastest 1-minute wind speed (tenths of meters per second)
# 	   WSF2 = Fastest 2-minute wind speed (tenths of meters per second)
# 	   WSF5 = Fastest 5-second wind speed (tenths of meters per second)
# 	   WSFG = Peak gust wind speed (tenths of meters per second)
# 	   WSFI = Highest instantaneous wind speed (tenths of meters per second)
# 	   WSFM = Fastest mile wind speed (tenths of meters per second)
# 	   WT** = Weather Type where ** has one of the following values:
# 	   
#                   01 = Fog, ice fog, or freezing fog (may include heavy fog)
#                   02 = Heavy fog or heaving freezing fog (not always 
# 		       distinquished from fog)
#                   03 = Thunder
#                   04 = Ice pellets, sleet, snow pellets, or small hail 
#                   05 = Hail (may include small hail)
#                   06 = Glaze or rime 
#                   07 = Dust, volcanic ash, blowing dust, blowing sand, or 
# 		       blowing obstruction
#                   08 = Smoke or haze 
#                   09 = Blowing or drifting snow
#                   10 = Tornado, waterspout, or funnel cloud 
#                   11 = High or damaging winds
#                   12 = Blowing spray
#                   13 = Mist
#                   14 = Drizzle
#                   15 = Freezing drizzle 
#                   16 = Rain (may include freezing rain, drizzle, and
# 		       freezing drizzle) 
#                   17 = Freezing rain 
#                   18 = Snow, snow pellets, snow grains, or ice crystals
#                   19 = Unknown source of precipitation 
#                   21 = Ground fog 
#                   22 = Ice fog or freezing fog
# 		  
#             WV** = Weather in the Vicinity where ** has one of the following 
# 	           values:
# 		   
# 		   01 = Fog, ice fog, or freezing fog (may include heavy fog)
# 		   03 = Thunder
# 		   07 = Ash, dust, sand, or other blowing obstruction
# 		   18 = Snow or ice crystals
# 		   20 = Rain or snow shower
   
   # save this data
   if (e_y == 1) {ghcn_TM <- gz_ghcn} else {ghcn_TM <- rbind(ghcn_TM, gz_ghcn)}
 }
  # ghcn_TM <- ghcn_TM[,V4:=ifelse(V4==9999|V4==-9999, NA, round(32 + 9 * V4/50, 3))]; ## convert to degrees F; noaa record is tenths of degrees C
 ghcn_TM <- ghcn_TM %>%
   mutate(V4=ifelse(V3=="PRCP", V4, ifelse(V4==9999|V4==-9999, NA, round(32 + 9 * V4/50, 3))))
 ghcn_TM[, V2 := ymd(V2)]; # convert charcter date to ymd format
 
 # read in PET (potential evapotranspiration) data:
 pet <- read_csv(petfile)
  
 PTT_PTR <- matrix(ncol = 8, nrow = 0); 
 for (e_i in 1:nrow(env_meta_info)) { 
   
  # for each environment, pull the metadata:
  env_code <- env_meta_info$env_code[e_i];
  field_lat <- env_meta_info$lat[e_i]; field_lon <- env_meta_info$lon[e_i]; 
  planting_date <- env_meta_info$PlantingDate[e_i];
  planting_year <- year(planting_date);
  
  # set the dates to inspect for this environment (from planting until your set end date); see above for details
  DAPs <- seq(as.Date(planting_date), length.out = searching_daps, by = "day");
  ending_date <- DAPs[searching_daps]; 
  ending_year <- year(ending_date);
  
  # read in day length file made earlier
  local_dl_file <- paste(sp_navy_dir, env_code, '_DL_', searching_daps, 'DAPs', sep = '')
  DL <- read.table(local_dl_file, header = T, sep = "\t", stringsAsFactors = F);

  # set name of temperature file and make it
  local_tm_file <- paste(sp_ghcn_dir, env_code, '_TM_', searching_daps, 'DAPs', sep = '' );
 #  local_tm_file <- paste(sp_isd_dir, env_code, '_TM_', searching_daps, 'DAPs', sep = '' );
  
  # this subfunction calculates the daily max and min temperatures, averaged over the nearest weather stations, for each environment
  # if (!file.exists(local_tm_file)) {

    TM_from_local_NOAA(ghcn_TM %>% filter(V3 %in% c("TMAX", "TMIN")), 
                                                       field_lat, field_lon, 
                                                       local_tm_file, STS_ghcn, DAPs)
    # } ;

  TM_0 <- read.table(local_tm_file, header = T, sep = "\t", stringsAsFactors = F); # read in temperature file created by above function

  PRCP_0 <- PRCP_from_local_NOAA(ghcn_TM %>% filter(V3 == "PRCP"), field_lat, field_lon, 
                                                       STS_ghcn, DAPs)
  PRCP_0$date <- as.character(PRCP_0$date)
  
  # merge day length with temperature data AND precip--
  DL_TM <- merge(DL, TM_0, all.x = T) ;
  DL_TM <- left_join(DL_TM, PRCP_0, by="date")
   
  # fill in missing temperatures
  DL_TM$TMAX <- Fill_Missning_TM(DL_TM$TMAX, env);
  DL_TM$TMIN <- Fill_Missning_TM(DL_TM$TMIN, env);
  DL_TM$precip <- Fill_Missning_TM(DL_TM$precip, env);
  
  # pull PET and join:
  current.env <- env_code
  current.pet <- pet%>%
    filter(env_code==current.env) %>%
    dplyr::select(-env_code)
  current.pet$date <- as.character(current.pet$date)
  DL_TM <- left_join(DL_TM, current.pet, by="date")
  
  # save max and min temperatures as vectors
  Tmax <- DL_TM$TMAX; Tmin <- DL_TM$TMIN; precip.vector <- DL_TM$precip
  
  # calculate GDD with sub-function:
  GDDs <- Adjusting_TM_4_GDD(Tmax, Tmin, Haun_threshold, t_base, t_max1, t_max2);
  
  # Calculate other environmental parameters:
  PTTs <- round(GDDs * DL_TM$DL, 4); # PTT=GDDxDL
  PTRs <- round(GDDs / DL_TM$DL, 4); # PTR=GDD/DL
  PTD1s <- round((DL_TM$TMAX - DL_TM$TMIN) * DL_TM$DL, 4); #PTD1=(TMAX-TMIN)xDL
  PTD2s <- round((DL_TM$TMAX - DL_TM$TMIN) / DL_TM$DL, 4); #PTD2=(TMAX-TMIN)/DL
  PTSs <- round(((DL_TM$TMAX^2) - (DL_TM$TMIN^2)) * (DL_TM$DL^2), 4); #PTS=(TMAX^2 - TMIN^2)*(DL^2)
  avg.temp <- (Tmax+Tmin)/2
  DTR <- Tmax-Tmin # diurnal temperature range
  H20.balance <- (DL_TM$precip)/10 - (DL_TM$PET)/100 # precip in 10th of mm, PET in 100th of mm so correct here
  env_codes <- rep(env_code, searching_daps); # make vector of the environmental code of the length of number of days of interest
  
  # make data frame to hold this data, and output
  t_df <- data.frame(env_code = env_codes, date = DAPs, TMAX = DL_TM$TMAX, TMIN = DL_TM$TMIN, DL = DL_TM$DL, GDD = GDDs, PTT = PTTs, 
                     PTR = PTRs, PTD1 = PTD1s, PTD2 = PTD2s, PTS = PTSs, PRECIP=precip.vector, DTR=DTR, PET=DL_TM$PET, H20.balance=H20.balance)
  PTT_PTR <- rbind(PTT_PTR, t_df);
 }
 
 write.table(PTT_PTR, file = ptt_ptr_file, sep = "\t", quote = F, row.name = F)
 invisible(gc());
}

# Function used in Compile_PTT_PTR_local_GHCN:
# Function to calculate day length for each day of the growing season
DayLength_from_Equation <- function(Y1, Y2, lat, lon, local_file, days) {
  DLs <- wholeyear_daylength_from_function(Y1, lat, lon); # actually go calculate the day lengths
  
  # if ending year was later than planting year (that is, if planted in fall/winter and harvested next spring),
  # you need to get the day length for that next year too
  if (Y2 > Y1) {
    DL_2 <- wholeyear_daylength_from_function(Y2, lat, lon); 
    DLs <- rbind(DLs, DL_2);
  }
  
  # keep only the day lengths in the dates you're interested in (planting to DAP limit), and output
  DL_window <- DLs[DLs$date %in% days, ];
  write.table(DL_window, file = local_file, sep = "\t", row.name = F, quote = F)
}

# sub-function used in Compile_PTT_PTR_local_GHCN: 
# calculate day lengths for a given location in a given year
wholeyear_daylength_from_function <- function(Y, lat_dec, lon_dec) {
  # set beginning and ending dates
  d1 <- paste(Y, '-1-1', sep = ''); 
  dL <- paste(Y, '-12-31', sep = '');
  
  # make a vector of the names of all the dates of that year:
  Ds <- seq(as.Date(d1), as.Date(dL), by = "days");
  
  # calculate the daylength for a given latitude and DAY of year (not date)
  DL <- round(daylength(lat_dec, 1:length(Ds)), 3);   
  
  # make dataframe with day length by date, and output it
  DL_df <- data.frame(date = Ds, DL = DL)
 return (DL_df);
}

# calculate GDD, including adjusting for the limit max and min temp and for the Haun threshold
# used in both Compile_PTT_PTR_local_GHCN and Compile_PTT_PTR
Adjusting_TM_4_GDD <- function(Tmax, Tmin, threshold, t_base, t_max1, t_max2) {
 Tmax[Tmax < t_base] <- t_base; Tmin[Tmin < t_base] <- t_base; # if the temperatures are BELOW the base temp, set it as the base temp
 if (Tmax[1] > t_max2) {Tmax[1] <- t_max2}; # if starting temperatures are above the second (Haun) threshold,
 if (Tmin[1] > t_max2) {Tmin[1] <- t_max2}; # set them as that threshold
 if (threshold > 0) { # IF the Haun threshold exists for this species
   gdd_cum <- (Tmax[1] + Tmin[1]) / 2 - t_base; # GDD formula, for first day of season
   for (i in 2:length(Tmax)) { # for each day of the season:
    if (gdd_cum < threshold) { t_max0 <- t_max2} else {t_max0 <- t_max1}; # set the max temp as appropriate, 
    if (Tmax[i] > t_max0) { Tmax[i] <- t_max0 }; #depending on whether cumulative GDD have surpassed Haun threshold
    if (Tmin[i] > t_max0) { Tmin[i] <- t_max0 };
    gdd_cum <- gdd_cum + (Tmax[i] + Tmin[i]) / 2 - t_base; # add day's results to cumulative GDD
    }    
  } else { # or without the Haun threshold, limit the daily temps 
     Tmax[Tmax > t_max1] <- t_max1; Tmin[Tmin > t_max1] <- t_max1;
     }
 gdds <- round((Tmax + Tmin) / 2 - t_base, 4) # calculate GDD!
 return (gdds)
} 

# used in Compile_PTT_PTR_local_GHCN
TM_from_local_NOAA <- function(TM_0,  lat, lon,  local_file, sts_ghcn, daps) {

  # pull nearest 6 weather stations within 50km of field location
  field_sts <- meteo_distance(sts_ghcn, lat = lat, lon = lon, radius = 50, limit = 6);

  # pull the weather station IDs:
  st_ids <- unique(as.vector(field_sts$id));
  T_mean <-  data.frame(date = daps); # make data frame currently holding all dates in growing season through DAPs limit
 
  TM <- TM_0[TM_0$V2 %in% daps][,1:4]; # pull temperature data for the dates of interest
 
  # extract daily max temperature data for each weather station
  TMAX <- TM[TM$V1 %in% st_ids & TM$V3 == 'TMAX'];
  # and find the average daily max temp data for each day, averaged over the different weather stations
  TMAX_mean <- TMAX[, .(TMAX = mean(V4)), by = V2];
  # format and add info to T_mean dataframe:
  colnames(TMAX_mean)[1] <- 'date'; 
  T_mean  <- merge(T_mean , TMAX_mean, all.x = T)
  
  # repeat for the minimum temperatures:
  TMIN <- TM[TM$V1 %in% st_ids & TM$V3 == 'TMIN'];
  TMIN_mean <- TMIN[, .(TMIN = mean(V4)), by = V2];
  colnames(TMIN_mean)[1] <- 'date';
  T_mean <- merge(T_mean, TMIN_mean, all.x = T);
   
  # write out result
  write.table(T_mean, file = local_file, sep = "\t", row.name = F, quote = F) ;
}

# used in Compile_PTT_PTR_local_GHCN
# make this work for Precip also
PRCP_from_local_NOAA <- function(PRCP_0, lat, lon, sts_ghcn, daps) {
  # pull nearest 6 weather stations within 50km of field location
  field_sts <- meteo_distance(sts_ghcn, lat = lat, lon = lon, radius = 50, limit = 6);

  # pull the weather station IDs:
  st_ids <- unique(as.vector(field_sts$id));
  T_mean <-  data.frame(date = daps); # make data frame currently holding all dates in growing season through DAPs limit
 
  PRCP <- PRCP_0[PRCP_0$V2 %in% daps][,1:4]; # pull temperature data for the dates of interest
 
  # extract daily precip data for each weather station
  precip <- PRCP[PRCP$V1 %in% st_ids & PRCP$V3 %in% c('PRCP')];
  # and find the average daily max temp data for each day, averaged over the different weather stations
  precip_mean <- precip[, .(precip = mean(V4)), by = V2];
  # format and add info to T_mean dataframe:
  colnames(precip_mean)[1] <- 'date'; 
  T_mean  <- merge(T_mean , precip_mean, all.x = T)
  return(T_mean)
}

# used in Compile_PTT_PTR_local_GHCN
Fill_Missning_TM <- function(Tx, env) {
 NA_inds <- which(is.na(Tx)); # identify dates with missing max temp data
 A_inds <- which(!is.na(Tx)) # and all other dates
 for (NA_ind in NA_inds ) {
   if (NA_ind == 1) {  # if you are missing the first day's max temp data,
     Tx[NA_ind] <- Tx[A_inds[1]] # just assign the next day's data to it
     } else if (NA_ind >= max(A_inds)) { # or, if you are missing the last day's temp data, 
       Tx[NA_ind] <- Tx[max(A_inds)] # just assign the previous day's data to it
       } else { # otherwise, 
         pre_ind <- A_inds[max(which(A_inds < NA_ind))]; # find the nearest previous non-missing day
         suff_ind <- A_inds[min(which(A_inds > NA_ind))]; # and the nearest following non-missing day
         Tx[NA_ind] <- mean(Tx[c(pre_ind, suff_ind)] ); # and average the two
        }
 }
 return(Tx);
}

# Function to search for the best environmental parameter and timeframe to predict phenotype
Exhaustive_search_full <- function(env_mean_trait, env_paras, searching_daps, exp_pop_dir,
                              FTdaps, trait, dap_x, dap_y, LOO, min.window.size) {

  # make file names
 pop_cor_file <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_', LOO, 'LOO_cor_whole_pop', sep = '');
 pop_pval_file <- paste(exp_pop_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_', LOO, 'LOO_P_whole_pop', sep = '');
 

 # if (!file.exists(pop_cor_file)) {
  dap_win <- (searching_daps) * (searching_daps)  / 2; # how many search windows do you need to iterate over? It's half of a square, so days^2 /2

  # make matrix to hold correlations and their p values
    pop_cors_matrix <- matrix(ncol = ((ncol(env_paras)-2)*2)+3, nrow = dap_win * 1);
    pop_pval_matrix <- matrix(ncol = ((ncol(env_paras)-2))+3, nrow = dap_win * 1);
    colnames(pop_cors_matrix) <-c("Day_x", "Day_y", "window",
        paste0("R_", colnames(env_paras)[-c(1:2)]),
        paste0("nR_", colnames(env_paras)[-c(1:2)]))
    colnames(pop_pval_matrix)<-c("Day_x", "Day_y", "window",
        paste0("pval_", colnames(env_paras)[-c(1:2)]))

  # # make matrix to hold correlations and their p values, with no negative versions of the param!
  #   pop_cors_matrix <- matrix(ncol = ((ncol(env_paras)-2))+3, nrow = dap_win * 1);
  #   pop_pval_matrix <- matrix(ncol = ((ncol(env_paras)-2))+3, nrow = dap_win * 1);
  #   colnames(pop_cors_matrix) <-c("Day_x", "Day_y", "window",
  #       paste0("R_", colnames(env_paras)[-c(1:2)]))
  #       # ,
  #       # paste0("nR_", colnames(env_paras)[-c(1:2)]))
  #   colnames(pop_pval_matrix)<-c("Day_x", "Day_y", "window",
  #       paste0("pval_", colnames(env_paras)[-c(1:2)]))
  #     # c("pop_code", 'Day_x', 'Day_y', 'window', 'R_DL', 'R_GDD', 'R_PTT', 'R_PTR', 'R_PTD', 'R_PTD2', 'R_PTS', 'nR_DL', 'nR_GDD', 'nR_PTT', 'nR_PTR', 'nR_PTD', 'nR_PTD2', 'nR_PTS');
    
    # loop through the season in different window sizes:
    n <- 0; # initialize iterator variable
      for (d1 in 1:(dap_y - min.window.size)) { # the +/- X here prevents it from choosing the very last day as the start date,
        for (d2 in (d1 + min.window.size):dap_y) { # or from choosing a day too close to the start date as the window end date
         days <- c(d1:d2); # make list of days in current window
         env_facts_matrix <- matrix(nrow = nrow(env_mean_trait), ncol = ncol(env_paras)-2);
         for (e_i in 1:nrow(env_mean_trait)) { # loop through the environments
           e <- env_mean_trait$env_code[e_i]; # set current environment code
           env_para <- subset(env_paras, env_paras$env_code == e); # pull all environmental info for this environ
           # env_mean <- colMeans(env_para[days, (c(8, 9, 10, 11, 12, 13, 14) - 3)]); # find mean for these values in this environ: ### DL, GDD, PTT, PTR, PTD, PTD2, PTS
           env_mean <- colMeans(env_para[days, c(3:ncol(env_para))])
           env_facts_matrix[e_i,] <- env_mean; # put the mean values of the parameters in this environment into a matrix
         }
         n <- n + 1; # increment iterator
         ### leave one environment out and get the median correlation
         Ymean_envPara <- cbind(env_facts_matrix, env_mean_trait$meanY); # add mean observed trait value column to your environmental parameters matrix
         rs <- c(); # save r (correlation) values
         ps <- c(); # need to keep p values  now too
         if (LOO == 0) { # if you are NOT leaving-one-out when calculating correlation:
          for (k in c(1:(ncol(env_paras)-2))) {
           rs[k] <- round(cor(Ymean_envPara[,ncol(Ymean_envPara)], Ymean_envPara[,k]), digits = 4) # find correlation between each environmental parameter and the environment mean
           ps[k] <- round(-log(cor.test(Ymean_envPara[,ncol(Ymean_envPara)], Ymean_envPara[,k])$p.value, 10), digits = 4) # save the -log10 p values for the correlation
           }
         
           } else { # if leaving one out when calculating correlation:
          loo_rs_matrix <- matrix(nrow = nrow(Ymean_envPara)+ 0, ncol = ncol(env_paras)-2); # make matrix to hold correlations
          loo_ps_matrix <- matrix(nrow = nrow(Ymean_envPara)+ 0, ncol = ncol(env_paras)-2); # make matrix to hold p val
          for (k in c(1:(ncol(env_paras)-2))) { ## loop through environment parameters
           for (e_x in c(1:nrow(Ymean_envPara))) { # loop through the environments
             t_matrix <- Ymean_envPara[-e_x,]; # drop one environment from your matrix
             loo_rs_matrix[e_x, k] <- round(cor(t_matrix[,ncol(Ymean_envPara)], t_matrix[,k]), digits = 4) # find correlation between trait value and environmental param, leaving out one environ
             loo_ps_matrix[e_x, k] <- round(-log(cor.test(Ymean_envPara[,ncol(Ymean_envPara)], t_matrix[,k])$p.value, 10), digits = 4) # find p value of correlation between trait value and environmental param, leaving out one environ
            }
           }
          rs <- apply(loo_rs_matrix, 2, median); # get median of the correlations for leaving one out to get final correlation between environmental param and trait, by environ param
          ps <- apply(loo_ps_matrix, 2, median); # get median of the pval for leaving one out to get final pval between environmental param and trait, by environ param
         }
         pop_cors_matrix[n, ] <- c(d1, d2, d2 - d1, rs, 0 - rs); # add window's info to the matrix AND add the NEGATIVE correlations too (using 0-rs)
         # pop_cors_matrix[n, ] <- c(d1, d2, d2 - d1, rs); # add window's info to the matrix WITHOUT the NEGATIVE correlations 
         pop_pval_matrix[n,] <- c(d1,d2,d2-d1, ps) # keep the p values

         }
      }
    pop_cors_matrix <- pop_cors_matrix[1:n,]
    pop_pval_matrix <- pop_pval_matrix[1:n,]
    write.table(pop_cors_matrix, file = pop_cor_file, sep = "\t", row.names = F, quote = F);
    write.table(pop_pval_matrix, file = pop_pval_file, sep = "\t", row.names = F, quote = F);
    
    # }
}

# Plot results of exhaustive search
exhaustive_plot <- function(pop_cor_file, pop_pval_file, env_paras, dap_x, dap_y,
                            FTdaps, type=c("p.only", "r.only", "both") # use type to choose which type of plot to print
                            ) {
  # read data back in
 pop_cors <- read.table(pop_cor_file, header = T, sep = "\t");
 pop_pvals <- read.table(pop_pval_file, header=T, sep="\t")

 if (type=="r.only") {
   # make a pdf to show all of the environmental parameters and their optimal windows:
 pdf(paste0(exp_pop_dir, "/MaxR_", trait, ".pdf"),
     width= (ncol(env_paras)-2),height= 2,
     pointsize=6)
 layout(matrix(c(1:((ncol(env_paras)-2)*2)), 2, (ncol(env_paras)-2), byrow = T))
 
 for (k in c(1:((ncol(env_paras)-2)*2))) { # loop through the environmental parameters, plus all of them again but with color scale switched
   # pop_cor_0 <- subset(pop_cors, pop_cors$pop_code == p); # keep only the window correlations for those in the desired pop
   pop_cor_0 <- pop_cors
   pop_cor <- pop_cor_0[,c(1:3, k + 3)]; # pull correlations for the environ parameter of interest
   colnames(pop_cor)[4] <- 'R';
   pop_cor <- pop_cor[order(pop_cor$R),]; # sort by correlation
   
   xs <- pop_cor$Day_x;  ys <-  pop_cor$Day_y; # save window start and end dates, in order of correlation
   mid_R <- median(pop_cor$R); # find median correlation across windows
   
   # set colors: 
   cell_col <- floor(pop_cor$R * 12) + 13; ### the same color scale
   pop_cor$cell_col <- cell_col; 
 
   max_R <- pop_cor[which.max(abs(pop_cor$R))[1], ]; # pull the maximum absolute value of correlation found
   
   # set parameters--
   par(mar = c(0.5, 1.0, 1, 0.5) , # set margins on all sides
       mgp = c(0.05, 0.1, 0), # set margin lines for axis title, axis labels, axis line
       tck = -0.01, # length of tick marks as fraction of plotting region
       bty = "n"); # character string, suppressing box drawing
   plot(-50, -50, xlim = c(0, dap_x), ylim = c(0, dap_y), col = "white",  xlab = '', 
        xaxt = "n", yaxt = "n", ylab = 'Days after planting', bty = "n", cex.lab = 1);
   arrows(-1, 10, -1, dap_y - 10, length = 0.05, angle = 15, lwd = .5,  col = "grey59");
   mtext(c(1, 50, 100, dap_y), side = 2, at = c(1,50, 100, dap_y), line = -1, cex = .6)
  
   rect(xs - 0.5, ys - 0.5, xs + 0.5, ys + 0.5, col = col_palette[pop_cor$cell_col], border = "NA")
   rect(max(pop_cor$Day_x) - 0.5, max(pop_cor$Day_y) - 0.5, max(pop_cor$Day_x) + 0.5, max(pop_cor$Day_y) + 0.5, border = "NA", col = "white", lwd = 0.001)

   arrows(10, dap_y + 4, dap_x - 10, dap_y + 4, angle = 15, length = 0.05, lwd = .5, col = "grey59")
   mtext("Days after planting", side = 3, at = dap_x / 2, line = -0.1, cex = .6)
   mtext(c(1, 50, 100, dap_y), side = 3, at = c(1, 50, 100, dap_y), line = -1.1, cex = .6)
   arrows(max_R$Day_x + 4,  max_R$Day_y - 4,  max_R$Day_x,  max_R$Day_y, length = 0.05, angle = 15, lwd = .5, col = "grey59")
   
   # Make legend:
   box_ys <- seq(1, 50, by = 2); box_xs <- rep(dap_x - 15, 25); 
   rect(box_xs - .5 * 2, box_ys - 0.5 * 2, box_xs + 0.5 * 2, box_ys + 0.5 * 2, border = "NA", col = col_palette)
   text(dap_x - 10 - 5, 52, 'r', cex = .5);
   r_lab_top <- 1; r_lab_mid <- 0; r_lab_bottom <- -1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", max_R$R), sep = '');
   # if (k >(ncol(env_paras)-2)) { # because we want the first SEVEN charts to have 1 at the top, and next 7 at the BOTTOM
   #   r_lab_top <- -1; r_lab_bottom <- 1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", 0 - max_R$R), sep = '');
   #   }
   legend(max_R$Day_x - 4 , max_R$Day_y - 4 , c(paste( max_R$Day_x, ' to ', max_R$Day_y, ' DAP ', colnames(pop_cor_0)[k+3], sep = ''), max_r_lab),  cex = .6, bty = "n");
   text(dap_x - 10 + 3, 50, r_lab_top, cex = .5)
   text(dap_x - 10 + 3, 27, r_lab_mid, cex = .5);
   text(dap_x - 10 + 3, 1,  r_lab_bottom, cex = .5)
   
   # Make box-and-whisker plot of the phenotypes
   boxplot(FTdaps,   at = 145,  add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
   boxplot(FTdaps,   at = 1, horizontal = T, add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
   # text(mean(FTdaps), 5, 'Days to anthesis', cex = .5)
   # text(mean(FTdaps, na.rm = T), 10, paste('Trait: ', trait, sep = ''), cex = 2)
   # title(paste('Trait: ', trait, sep = ''), adj=0)
 }
 dev.off()
 }
 else if (type=="p.only") {
   # make a pdf to show all of the environmental parameters and their optimal windows:
 pdf(paste0(exp_pop_dir, "/MaxP_", trait, ".pdf"),
     width= (ncol(env_paras)-2),height= 1,
     pointsize=6)
 layout(matrix(c(1:((ncol(env_paras)-2))), 1, (ncol(env_paras)-2), byrow = T))
 
 for (k in c(1:((ncol(env_paras)-2)))) { # loop through the environmental parameters
   # pop_cor_0 <- subset(pop_cors, pop_cors$pop_code == p); # keep only the window correlations for those in the desired pop
   pop_pval_0 <- pop_pvals
   pop_pval <- pop_pval_0[,c(1:3, k + 3)]; # pull correlations for the environ parameter of interest
   colnames(pop_pval)[4] <- 'P';
   pop_pval <- pop_pval[order(pop_pval$P),]; # sort by correlation
   
   xs <- pop_pval$Day_x;  ys <-  pop_pval$Day_y; # save window start and end dates, in order of correlation
   mid_pval <- median(pop_pval$P); # find median correlation across windows
   max_P <- pop_pval[which.max(abs(pop_pval$P))[1], ]; # pull the maximum absolute value of correlation found

   # set colors: 
   cell_col <-(floor(pop_pval$P/(max_P$P/12)))+13
     # unique((floor(pop_pval$P/(max_P$P/12)))+13)
   # cell_col <- floor(pop_pval$P * 12) +13; ### the same color scale
   pop_pval$cell_col <- cell_col; 
 
   
   # set parameters--
   par(mar = c(0.5, 1.0, 1, 0.5) , # set margins on all sides
       mgp = c(0.05, 0.1, 0), # set margin lines for axis title, axis labels, axis line
       tck = -0.01, # length of tick marks as fraction of plotting region
       bty = "n"); # character string, suppressing box drawing
   plot(-50, -50, xlim = c(0, dap_x), ylim = c(0, dap_y), col = "white",  xlab = '', 
        xaxt = "n", yaxt = "n", ylab = 'Days after planting', bty = "n", cex.lab = 1);
   arrows(-1, 10, -1, dap_y - 10, length = 0.05, angle = 15, lwd = .5,  col = "grey59");
   mtext(c(1, 50, 100, dap_y), side = 2, at = c(1,50, 100, dap_y), line = -1, cex = .6)
  
   rect(xs - 0.5, ys - 0.5, xs + 0.5, ys + 0.5, col = col_palette[pop_pval$cell_col], border = "NA")
   rect(max(pop_pval$Day_x) - 0.5, max(pop_pval$Day_y) - 0.5,
        max(pop_pval$Day_x) + 0.5, max(pop_pval$Day_y) + 0.5, border = "NA", col = "white", lwd = 0.001)

   arrows(10, dap_y + 4, dap_x - 10, dap_y + 4, angle = 15, length = 0.05, lwd = .5, col = "grey59")
   mtext("Days after planting", side = 3, at = dap_x / 2, line = -0.1, cex = .6)
   mtext(c(1, 50, 100, dap_y), side = 3, at = c(1, 50, 100, dap_y), line = -1.1, cex = .6)
   arrows(max_P$Day_x + 4,  max_P$Day_y - 4,  max_P$Day_x,  max_P$Day_y, length = 0.05, angle = 15, lwd = .5, col = "grey59")
   
   # Make legend:
   box_ys <- seq(1, 50, by = 2); box_xs <- rep(dap_x - 15, 25); 
   rect(box_xs - .5 * 2, box_ys - 0.5 * 2, box_xs + 0.5 * 2, box_ys + 0.5 * 2, border = "NA", col = col_palette)
   text(dap_x - 10 - 5, 52, '-logP', cex = .5);
   p_lab_top <- max_P$P;  p_lab_bottom <- min(pop_pval$P); p_lab_mid <- (p_lab_top + p_lab_bottom)/2; max_p_lab <- paste( '-logP = ', sprintf( "%.3f", max_P$P), sep = '');
   # if (k >(ncol(env_paras)-2)) { # because we want the first SEVEN charts to have 1 at the top, and next 7 at the BOTTOM
   #   r_lab_top <- -1; r_lab_bottom <- 1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", 0 - max_R$R), sep = '');
   #   }
   legend(max_P$Day_x - 4 , max_P$Day_y - 4 , c(paste( max_P$Day_x, ' to ', max_P$Day_y, ' DAP ',  
                                                       colnames(pop_pval_0)[k+3],sep = ''), max_p_lab),  cex = .6, bty = "n");
   text(dap_x - 10 + 3, 50, p_lab_top, cex = .5)
   text(dap_x - 10 + 3, 27, p_lab_mid, cex = .5);
   text(dap_x - 10 + 3, 1,  p_lab_bottom, cex = .5)
   
   # Make box-and-whisker plot of the phenotypes
   boxplot(FTdaps,   at = 145,  add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
   boxplot(FTdaps,   at = 1, horizontal = T, add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
   # text(mean(FTdaps), 5, 'Days to anthesis', cex = .5)
   # text(mean(FTdaps, na.rm = T), 10, paste('Trait: ', trait, sep = ''), cex = 2)
   # title(paste('Trait: ', trait, sep = ''), adj=0)
 }
 dev.off()
 }
 else if (type == "both") {
   pdf(paste0(exp_pop_dir, "/MaxRandP_", trait, ".pdf"),
     width= (ncol(env_paras)-2),height= 2,
     pointsize=6)
 layout(matrix(c(1:((ncol(env_paras)-2)*2)), 2, (ncol(env_paras)-2), byrow = T))
 
 for (k in c(1:((ncol(env_paras)-2)*2))) { # loop through the environmental parameters, plus all of them again but with color scale switched
   
   # First, pull r plot info:
   
   # pop_cor_0 <- subset(pop_cors, pop_cors$pop_code == p); # keep only the window correlations for those in the desired pop
   pop_cor_0 <- pop_cors
   pop_cor <- pop_cor_0[,c(1:3, k + 3)]; # pull correlations for the environ parameter of interest
   colnames(pop_cor)[4] <- 'R';
   pop_cor <- pop_cor[order(pop_cor$R),]; # sort by correlation
   
   xs <- pop_cor$Day_x;  ys <-  pop_cor$Day_y; # save window start and end dates, in order of correlation
   mid_R <- median(pop_cor$R); # find median correlation across windows
   
   # set colors: 
   cell_col <- floor(pop_cor$R * 12) + 13; ### the same color scale
   pop_cor$cell_col <- cell_col; 
 
   max_R <- pop_cor[which.max(abs(pop_cor$R))[1], ]; # pull the maximum absolute value of correlation found
   
   # Now, pull P plot info:
   
   # make new iterator for pval info
   if (k %in% c(1:((ncol(env_paras)-2)))) {
     l <- k 
   } else {
     l <- k-(ncol(env_paras)-2)
   }
   
   pop_pval_0 <- pop_pvals
   pop_pval <- pop_pval_0[,c(1:3, l + 3)]; # pull correlations for the environ parameter of interest
   colnames(pop_pval)[4] <- 'P';
   pop_pval <- pop_pval[order(pop_pval$P),]; # sort by correlation
   
   pxs <- pop_pval$Day_x;  pys <-  pop_pval$Day_y; # save window start and end dates, in order of correlation
   mid_pval <- median(pop_pval$P); # find median correlation across windows
   max_P <- pop_pval[which.max(abs(pop_pval$P))[1], ]; # pull the maximum absolute value of correlation found

   # set colors: 
   pcell_col <-(floor(pop_pval$P/(max_P$P/12)))+13
     # unique((floor(pop_pval$P/(max_P$P/12)))+13)
   # cell_col <- floor(pop_pval$P * 12) +13; ### the same color scale
   pop_pval$pcell_col <- pcell_col; 
   
# Make the figure:
      
   # set parameters--
   par(mar = c(0.5, 1.0, 1, 0.5) , # set margins on all sides
       mgp = c(0.05, 0.1, 0), # set margin lines for axis title, axis labels, axis line
       tck = -0.01, # length of tick marks as fraction of plotting region
       bty = "n"); # character string, suppressing box drawing
   plot(-50, -50, xlim = c(0, dap_x), ylim = c(0, dap_y), col = "white",  xlab = '', 
        xaxt = "n", yaxt = "n", ylab = 'Days after planting', bty = "n", cex.lab = 1);
   arrows(-1, 10, -1, dap_y - 10, length = 0.05, angle = 15, lwd = .5,  col = "grey59");
   mtext(c(1, 50, 100, dap_y), side = 2, at = c(1,50, 100, dap_y), line = -1, cex = .6)
  
   rect(xs - 0.5, ys - 0.5, xs + 0.5, ys + 0.5, col = col_palette[pop_cor$cell_col], border = "NA") # r values
   rect(pys-0.5, pxs-0.5, pys+0.5, pxs+0.5, col=col_palette[pop_pval$pcell_col], border="NA") # p values
   
   rect(max(pop_cor$Day_x) - 0.5, max(pop_cor$Day_y) - 0.5, max(pop_cor$Day_x) + 0.5, max(pop_cor$Day_y) + 0.5, border = "NA", col = "white", lwd = 0.001)

   arrows(10, dap_y + 4, dap_x - 10, dap_y + 4, angle = 15, length = 0.05, lwd = .5, col = "grey59")
   mtext("Days after planting", side = 3, at = dap_x / 2, line = -0.1, cex = .6)
   mtext(c(1, 50, 100, dap_y), side = 3, at = c(1, 50, 100, dap_y), line = -1.1, cex = .6)
   
   # make R arrow, then P arrow
   arrows(max_R$Day_x + 4,  max_R$Day_y - 4,  max_R$Day_x,  max_R$Day_y, length = 0.05, angle = 15, lwd = .5, col = "grey59")
   arrows(max_P$Day_y + 4,  max_P$Day_x - 4,  max_P$Day_y,  max_P$Day_x, length = 0.05, angle = 15, lwd = .5, col = "grey59")

   # Make legend (r):
   box_ys <- seq(1, 50, by = 2); box_xs <- rep(dap_x - 15, 25); 
   rect(box_xs - .5 * 2, box_ys - 0.5 * 2, box_xs + 0.5 * 2, box_ys + 0.5 * 2, border = "NA", col = col_palette)
   text(dap_x - 20, 55, 'r', cex = .5);
   r_lab_top <- 1; r_lab_mid <- 0; r_lab_bottom <- -1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", max_R$R), sep = '');
   # if (k >(ncol(env_paras)-2)) { # because we want the first SEVEN charts to have 1 at the top, and next 7 at the BOTTOM
   #   r_lab_top <- -1; r_lab_bottom <- 1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", 0 - max_R$R), sep = '');
   #   }
   legend(max_R$Day_x - 4 , max_R$Day_y - 4 , c(paste( max_R$Day_x, ' to ', max_R$Day_y, ' DAP ', colnames(pop_cor_0)[k+3], sep = ''), max_r_lab),  cex = .6, bty = "n");
   text(dap_x - 20, 50, r_lab_top, cex = .5)
   text(dap_x - 20, 27, r_lab_mid, cex = .5);
   text(dap_x - 20, 1,  r_lab_bottom, cex = .5)
   
   # Make legend (p)
   # box_pys <- seq(1, 50, by = 2); box_pxs <- rep(dap_x, 25); 
   # rect(box_pxs - .5 * 2, box_pys - 0.5 * 2, box_pxs + 0.5 * 2, box_pys + 0.5 * 2, border = "NA", col = col_palette)
   text(dap_x - 10, 55, '-logP', cex = .5);
   p_lab_top <- max_P$P;  p_lab_bottom <- min(pop_pval$P); p_lab_mid <- (p_lab_top + p_lab_bottom)/2; max_p_lab <- paste( '-logP = ', sprintf( "%.3f", max_P$P), sep = '');
   legend(max_P$Day_y - 4 , max_P$Day_x - 4 , c(paste( max_P$Day_x, ' to ', max_P$Day_y, ' DAP ',  
                                                       colnames(pop_pval_0)[l+3],sep = ''), max_p_lab),  cex = .6, bty = "n");
   text(dap_x - 7 , 50, p_lab_top, cex = .5)
   text(dap_x - 7 , 27, p_lab_mid, cex = .5);
   text(dap_x - 7, 1,  p_lab_bottom, cex = .5)
   # Make box-and-whisker plot of the phenotypes
   # boxplot(FTdaps,   at = 145,  add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
   # boxplot(FTdaps,   at = 1, horizontal = T, add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
   # text(mean(FTdaps), 5, 'Days to anthesis', cex = .5)
   # text(mean(FTdaps, na.rm = T), 10, paste('Trait: ', trait, sep = ''), cex = 2)
   # title(paste('Trait: ', trait, sep = ''), adj=0)
 }
 dev.off()
 }
}

# Plot environmental parameter means vs trait mean 
Plot_Trait_mean_envParas <- function( env_mean_trait, env_paras, d1, d2, trait, exp_trait_dir, env_cols){
  days <- c(d1:d2); # make window
  env_facts_matrix <- matrix(nrow = nrow(env_mean_trait), ncol = ncol(env_paras)-1); # make matrix to hold info for each environment
  for (e_i in 1:nrow(env_mean_trait)) { # loop through the environments
    e <- env_mean_trait$env_code[e_i]; # pull current environment code
    env_para <- subset(env_paras, env_paras$env_code == e); # pull the environmental parameters for that environment
    # env_mean <- colMeans(env_para[days, c(8, 9, 10, 11, 12, 13, 14) - 3]); # find mean parameter values for across environments ### DL, GDD, PTT, PTR, PTS
    env_mean <- colMeans(env_para[days, c(3:ncol(env_para))])

    env_facts_matrix[e_i,] <- c(env_mean_trait$meanY[e_i], round(env_mean, 4) ); # put the mean Y and the environmental parameter means into a matrix
  }
  colnames(env_facts_matrix) <- c( 'meanY', colnames(env_paras)[-c(1:2)])
                                   # 'DL', 'GDD', 'PTT', 'PTR', 'PTD', 'PTD2', 'PTS'); # add column names
  envMeanPara_file <- paste(exp_trait_dir, trait, '_envMean', kPara_Name, d1, '_', d2, sep = ''); # make file name
  envMeanPara <- merge(env_mean_trait, env_facts_matrix); # combine mean Y with environmental param means into one table
  write.table(envMeanPara, file = envMeanPara_file, sep = "\t", row.names = F, quote = F); # and output the table

  # Make a pdf figure:
  pdf_file <- paste(exp_trait_dir, trait, 'Mean_', nrow(env_mean_trait), 'Env', kPara_Name, '.pdf', sep = ''); 
  pdf(pdf_file,
      width= ncol(env_paras)-2,height= 1,
      pointsize=6)

  layout(matrix(c(1:(ncol(env_paras)-2)), ncol = (ncol(env_paras)-2)));
  
  for (i in 1:(ncol(env_paras)-2)) { # loop through environmental parameters:
   par(mar = c(2.5, 2.0, 1, 0.5) , mgp = c(0.7, 0.01, 0), tck = -0.01, family = "mono");
    
    #plot environmental parameter (x) vs phenotype (Y), using mean of each for the environment
   plot(env_facts_matrix[, i + 1], env_facts_matrix[,1],
        xlab = colnames(env_paras)[i + 2], ylab = paste(trait, ' mean', sep = ''),  pch = 19, col = env_cols);
    abline(lm(env_facts_matrix[,1] ~ env_facts_matrix[, i + 1]), lty = 2); # add regression line
    
    # add correlation label
   r1 <- round(cor(env_facts_matrix[,1] , env_facts_matrix[, i + 1]), 3);
   legend("bottom", paste('r = ', r1, sep = ''), bty = "n")
   legend_p <- "topleft";  if (r1 < 0) { legend_p <- "topright"};
   
   # add legend with environmental codes
   if (i == 1) { legend(legend_p, env_mean_trait$env_code, pch = 19, col = env_cols, bty = "n" )};
  }
  dev.off()
}

# regress each line/genotype to the environmental mean and parameter
Slope_Intercept <- function(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, exp_trait, line_codes, exp_pop_dir, PTT_PTR_ind, filter.less.than) {
 maxR_win <- c(maxR_dap1:maxR_dap2); # pull the window with the max R value
 prdM <- env_mean_trait;
 kPara <- c();
 for (e_i in 1:nrow(env_mean_trait)) { # loop through the environments
   envParas <- subset(PTT_PTR, PTT_PTR$env_code == env_mean_trait$env_code[e_i]); # pull the PTT etc info for this environment
      envPara <- mean(envParas[maxR_win, PTT_PTR_ind]); # find mean PTT or other info for the environment over the course of the growing season
   kPara[e_i] <- envPara; # save mean environmental parameter for this environment
 }
 prdM$kPara <- kPara;
 lm_ab_matrix <- matrix(ncol = 6,
                        nrow = length(line_codes));
 for (l in 1: length(line_codes)) { # for each line/taxa:
   l_trait <- subset(exp_trait, exp_trait$ril_code == line_codes[l]); # pull observed values for this line in all environments
   if(nrow(l_trait) >= filter.less.than) { # if this line was grown in at least 5 environments (number based on Xianran's code):
     l_trait <- merge(l_trait, prdM); # add optimal environmental parameter and meanY info
     lm_Mean <- lm(Yobs ~ meanY, data = l_trait); # linear model with Y=observedY and X = meanY
     lm_Para <- lm(Yobs ~ kPara, data = l_trait); # linear model with Y=observedY and X = environmental mean of chosen parameter
     
     # Get the predicted Y (what we're calling "intercept") from the models with meanY (lm_Mean) and environmental param (lm_Para) as the X values:
     a_Mean <- as.vector(round(predict(lm_Mean, # use `predict` to predict the expected Y=observedY given this X=meanY
                                       data.frame(meanY = mean(prdM$meanY))), # make a dataframe holding the overall mean phenotype (Y) value
                               4)); 
     a_Para <- as.vector(round(predict(lm_Para, data.frame(kPara = mean(prdM$kPara))), 4)); ## adjusted by the population mean

     # Pull slopes (b)
     b_Mean <- as.vector(round(lm_Mean$coefficient[2], 4)); # pull slope from the lm_Mean model
     # b_Para <- as.vector(round(lm_Para$coefficient[2],4)); # pull slope from the lm_Para model
     # EDITED Nov 1 2022 -- for trait KPR, slopes are < 0.0001 so getting all shrunk to zero!!
     # So, DON'T ROUND
     b_Para <- as.vector(lm_Para$coefficient[2]); # pull slope from the lm_Para model
     
     a_Para_ori <- as.vector(round(lm_Para$coefficient[1],4)); # pull intercept from lm_Para model
     
     # store data:
     lm_ab_matrix[l,] <- c(line_codes[l], a_Mean, b_Mean, a_Para, a_Para_ori, b_Para);
   }
 }
 lm_ab_matrix <- lm_ab_matrix[!is.na(lm_ab_matrix[,2]),];
 colnames(lm_ab_matrix) <- c('line_codes',
                             'Intcp_mean', # = a_Mean = predicted Y ("intercept") when using pheno mean as predictor
                             'Slope_mean', # = b_Mean = slope when using pheno mean as predictor
                             'Intcp_para_adj', # = a_Para = predicted Y ("intercept") when using env param as predictor
                             'Intcp_para', # = a_Para_ori = mathematical intercept when using env param as predictor
                             'Slope_para'); # = b_Para = slope when using env param as predictor
 out_file <- paste(exp_pop_dir, 'Intcp_Slope', sep = '');
 write.table(lm_ab_matrix, file = out_file, sep = "\t", quote = F, row.names = F)

} 

# leave-one-out cross-validation
LOOCV <- function(maxR1, maxR2, env_mean_trait, PTT_PTR, PTT_PTR_ind, exp_trait, obs_prd_file, filter.less.than, line_codes, output=TRUE) {
 maxR_win <- c(maxR1:maxR2); # create optimal window
 prdM <- env_mean_trait; # assign mean phenotype by environment to new object
 maxR_envPara <- matrix(ncol = 2, nrow = nrow(env_mean_trait));
 kPara <- c();
 for (e_i in 1:nrow(env_mean_trait)) { # loop through the environments
   envParas <- subset(PTT_PTR, PTT_PTR$env_code == env_mean_trait$env_code[e_i]); # pull PTT etc info for relevant environ
   envPara <- mean(envParas[maxR_win, PTT_PTR_ind]); # find mean environmental parameter for optimal window + optimal environ param
   kPara[e_i] <- envPara; # save this -- the mean value of the chosen parameter in this environment
 }
 prdM$kPara <- kPara; 
 obs_prd_m <- matrix(0, ncol = 6 ,nrow = nrow(exp_trait));
 
 n <- 0; 
 
 for (l in line_codes) { # loop through the taxa/lines:
   l_trait <- subset(exp_trait, exp_trait$ril_code == l); # pull phenotype in all environs for given line
   ril_data <- merge(prdM, l_trait,  all.x = T); # add this pheno info to the table with environmental param and mean Y per environ
     if (length(which(!is.na(ril_data$Yobs))) >= filter.less.than) { # make sure at least 5 observations of this line
     for (e in 1:nrow(ril_data)) { # loop through environments
       obs_trait <- ril_data$Yobs[e]; # pull observed phenotype for current environment, for current line
       if (!is.na(obs_trait)) {
         trn <- ril_data[-e,]; # drop current environment
         l_mean <- mean(trn$Yobs, na.rm = T); # find mean of remaining observations

         # Predict the current environment using model based on other environments
         prd_trait_mean  <- round(predict( lm(Yobs ~ meanY, data = trn), ril_data[e,]), digit = 3); # model with Y=observedY and X=meanY
         prd_trait_kpara <- round(predict( lm(Yobs ~ kPara, data = trn), ril_data[e,]), digit = 3); # model with Y=observedY and X=environmental parameter
         n <- n + 1; # increment iterator
         obs_prd_m[n,] <- c(ril_data$env_code[e], 
                            # ril_data$pop_code[e], 
                            l, prd_trait_mean, prd_trait_kpara, obs_trait, l_mean);

       }
     }
   }
 }


 
 obs_prd_m <- obs_prd_m[1:n,]
 colnames(obs_prd_m) <- c('env_code', 
                          # 'pop_code', 
                          'ril_code', # aka line name
                          'Prd_trait_mean', # predicted value for trait in this environment (using model based on other environments' data) using environmental mean as predictor
                          'Prd_trait_kPara', # same as above, but with environmental parameter as predictor
                          'Obs_trait', # true observed value
                          'Line_mean'); # mean of all the other environments

 if(output) {write.table(obs_prd_m, file = obs_prd_file, sep = "\t", quote = F);
   return(prdM)
   } else {return(list(prdM=prdM, obs_prd=obs_prd_m))}
}

# Plot the observed vs predicted results using several prediction methods
Plot_prediction_result <- function(obs_prd_file, all_env_code, prdM, kPara_Name, forecast_png_file, trait=trait,
                                   print.legend=T, plot.A=T, plot.B=T, plot.C=T, plot.D=T, path=T, save.output=T) {
  if (path==T) { # if you were provided a path, read it in
     Obs_Prd_m <- read.table(obs_prd_file, sep = "\t", header = T); # pull in LOO prediction results
  } else {Obs_Prd_m <- obs_prd_file} #otherwise, just use the file as the obs_prd_m
 
 Obs_Prd_m <- Obs_Prd_m[!is.na(Obs_Prd_m$Obs_trait),]; # drop missing observations
 prd_env <- as.vector(unique(Obs_Prd_m$env_code)); # pull environment names
 env_rs <- matrix(ncol = 3, nrow = length(prd_env));
 
 # # redudant--
 # for (e_i in 1:length(prd_env)) { # loop through environments
 #   env_obs_prd <- subset(Obs_Prd_m, Obs_Prd_m$env_code == prd_env[e_i]); # pull data for current environment
 #   if (nrow(env_obs_prd) > 0) { # make sure data is available
 #    env_rs[e_i,] <- c( sprintf( "%.2f", cor(env_obs_prd$Prd_trait_mean, env_obs_prd$Obs_trait, 
 #                                            use = "complete.obs")), # corr btw value predicted with pheno mean, and observed
 #                       sprintf( "%.2f", cor(env_obs_prd$Prd_trait_kPara, env_obs_prd$Obs_trait,
 #                                            use = "complete.obs")), # corr btw value predicted with PTT etc, and observed
 #                       sprintf( "%.2f", cor(env_obs_prd$Line_mean, env_obs_prd$Obs_trait, use = "complete.obs"))); # corr btw mean value for line in all other environs, and observed
 #   }
 #    
 # }
 
 # xy_lim <- range(Obs_Prd_m[,3:5],na.rm = T) # set limits for figure
 
 # Make figure:
 if(save.output==T) {
    png(forecast_png_file, width = 4/1.5, height= 4/1.5,pointsize=6, units = "in", res = 600)
 }
 if(sum(plot.A, plot.B, plot.C, plot.D)>1){
    layout(matrix(c(1:4), 2, 2, byrow = T)); # only need layout if plotting multiple panels
 }
   obs_prd_m <- Obs_Prd_m

 if(plot.A==T){
 # Plot A: observed vs predicted by environmental mean
  # obs_prd_m <- subset(Obs_Prd_m, Obs_Prd_m$pop_code == p); # keep only desired pop
 xy_limA <- range(obs_prd_m[,c(3,5)], na.rm=T)
  par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
  # plot(obs_prd_m[,5], obs_prd_m[,3], # plot value predicted with pheno mean, vs observed
  plot(obs_prd_m[,3], obs_prd_m[,5], # plot value predicted with pheno mean, vs observed, with obs on y axis
       col = env_cols[match(obs_prd_m[,1], all_env_codes)], pch = 19, cex = .4, 
       xlab = paste('Predicted ', trait, ' by envMean', sep = ''), 
       ylab = paste('Observed ', trait, '', sep = ''), 
       xlim = xy_limA, ylim = xy_limA);
  abline(a = 0, b = 1, lty = 2, col = "gray59"); # add y=x line
  r1 <- sprintf("%.2f", cor(obs_prd_m[,3], obs_prd_m[,5], use = "complete.obs")); # add correlation to plot
  legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r1)), bty = "n")
  mtext('A', side = 3, at = xy_limA[1], line = 0.1, cex = .8); # label panel A
 }
 
 if(plot.B==T) {
  # Plot B: value predicted by environmental parameter (PTT etc) vs observed
   xy_limB <- range(obs_prd_m[,c(4,5)], na.rm=T)
  
  all_env_r <- obs_prd_m %>% # pull the correlations for each environment individually
      group_by(env_code) %>%
      mutate(env.cor=sprintf("%.2f", cor(Prd_trait_kPara, Obs_trait, use = "complete.obs"))) %>%
      dplyr::select(env_code, env.cor) %>%
      distinct%>%
      mutate(legend=paste0(env_code, " r=", env.cor))
  
  par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
  plot(obs_prd_m$Prd_trait_kPara, obs_prd_m$Obs_trait, # plot value predicted with environmental parameter (PTT etc), vs observed
       col = env_cols[match(obs_prd_m$env_code, all_env_codes)], pch = 19, cex = .4, 
       xlab = paste('Predicted ', trait, ' by ', kPara_Name, sep = ''), 
       ylab = paste('Observed ', trait, '', sep = ''), 
       xlim = xy_limB, ylim = xy_limB);
  abline(a = 0, b = 1, lty = 2, col = "gray59");
  mtext('B', side = 3, at = xy_limB[1], line = 0.1, cex = .8);
  r2 <- sprintf("%.2f", cor(obs_prd_m[,4], obs_prd_m[,5], use = "complete.obs")); # add correlation
  legend("top", legend= substitute(paste("overall ", italic('r'), " = ", R1), list(R1 = r2)), bty = "n")
  # legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r2)), bty = "n")
  # add legend with environmental codes
   if (print.legend) {
     legend(x=min(xy_limB), y=max(xy_limB), all_env_r$legend, pch = 19, col = env_cols, bty = "n" )
   }
 }
  if(plot.C==T) {
  # Plot C: environmental parameter vs mean observation in each environment
     xy_limC <- range(c(prdM$meanY), na.rm=T)
  par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
  plot(
    # prdM$PTR,
       prdM$kPara,
       prdM$meanY, # environmental parameter vs mean observation in each environment
       col = env_cols[match(prdM$env_code, all_env_codes)],  
       ylim = xy_limC, pch = 19, cex = .65, xlab = kPara_Name, ylab = 'Observed population mean');
  mtext(prdM$env_code, side = 1,
        at = prdM$kPara,
        # at=prdM$PTR,
        las = 2, line = -2, cex = .6 )
  abline(lm(prdM$meanY ~ prdM$kPara)) # add regression line
  # abline(lm(prdM$meanY ~ prdM$PTR)) # add regression line
  r2 <- sprintf("%.2f", cor(prdM$meanY, 
                            prdM$kPara)); # calculate correlation
                            # prdM$PTR)); # calculate correlation
legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r2)), bty = "n")
  mtext('C', side = 3, at = min(prdM$kPara), line = 0.1, cex = .8);
  }
 
  if(plot.D==T) {
  # Plot D: mean value for line in all other environs, vs observed
   xy_limD <- range(obs_prd_m[,c(5,6)], na.rm=T)
  par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
  plot(obs_prd_m[,6], obs_prd_m[,5],
       col = env_cols[match(obs_prd_m[,1], all_env_codes)], pch = 19, cex = .4, 
       ylab = paste('Observed ', trait, sep = ''), xlab = paste('Predicted ', trait, ' by BLUE', sep = ''), 
       xlim = xy_limD, ylim = xy_limD);
  abline(a = 0, b = 1, lty = 2, col = "gray59");
  r1 <- sprintf("%.2f", cor(obs_prd_m[,6], obs_prd_m[,5], use = "complete.obs"));
  legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r1)), bty = "n")
  mtext('D', side = 3, at = xy_limD[1], line = 0.1, cex = .8);
  }
 if(save.output==T) {dev.off()}

}



# function for UETG (Untested Environment, Tested Genotype) 1 to 2 prediction:
UETG.1to2.function <- function(maxR_dap1, maxR_dap2, PTT_PTR_ind, env_mean_trait_trn, 
                               env_mean_trait_test, PTT_PTR, exp_trait, test_env, filter.less.than) {
  maxR_win <- c(maxR_dap1:maxR_dap2); # set window
  prdM <- rbind(env_mean_trait_trn, env_mean_trait_test); 
  kPara <- c();

  # pull environmental parameter values for each environment
  for (e_i in 1:nrow(prdM)) {
    envParas <- subset(PTT_PTR, PTT_PTR$env_code == prdM$env_code[e_i]);
    envPara <- mean(envParas[maxR_win, PTT_PTR_ind]);
    kPara[e_i] <- envPara;
  }
  prdM$kPara <- kPara;
  
  obs_prd_m <- matrix(0, ncol = 6,nrow = nrow(exp_trait));
  n1 <- 1; n2 <- 0;
  for (l in line_codes) { # loop through lines
    l_trait <- subset(exp_trait, exp_trait$ril_code == l);# pull phenotype in all environs for given line
    ril_data <- merge(prdM, l_trait,  all.x = T); # add this pheno info to the table with environmental param and mean Y per environ
    trn <- ril_data[!(ril_data$env_code %in% test_env), ]; # pull environments for training (not your testing env)
    prd <- ril_data[(ril_data$env_code %in% test_env),]; # pull the  testing environment
     if (sum(!is.na(trn$Yobs)) >= filter.less.than-1 & sum(!is.na(prd$Yobs)) >= 1 ) { # make sure at least filter.less.than-1 training observations (filter.less.than-1), and one to predict
      obs_trait <- prd$Yobs; # pull observed value in testing set

      prd_trait_mean  <- round(predict( lm(Yobs ~ meanY, data = trn), prd), digit = 3); # predict the value in the testing set, using the testing set data. AND--use mean Y value as predictor
      prd_trait_kpara <- round(predict( lm(Yobs ~ kPara, data = trn), prd), digit = 3); # same as above, but use the environmental parameter as the predictor

      # save output
      obs_prd_t <- cbind(prd$env_code, rep(l, nrow(prd)), 
                         prd_trait_mean, prd_trait_kpara, prd$Yobs,
                         rep(mean(trn$Yobs, na.rm = T), nrow(prd)));
      n2 <- n1 + nrow(prd) - 1;
      obs_prd_m[n1:n2,] <- obs_prd_t;
      n1 <- n2 + 1;
    }
  }
  # output data: see LOOCV for more details
  obs_prd_m <- obs_prd_m[1:n2,]
  colnames(obs_prd_m) <- c('env_code', 'ril_code', 'Prd_trait_mean', 'Prd_trait_kPara', 'Obs_trait', 'Line_mean');
  # write.table(obs_prd_m, file = obs_prd_file, sep = "\t", quote = F, row.name = F);
  return(list(obs_prd_m, prdM));
  
}

# Plot result of 1 to 3 (TEUG = Tested Environment, Untested Genotype) prediction; 
# can also apply to other prediction scenarios
Plot_TEUG_result <- function(obs_prd_file, all_env_codes, kPara_Name, forecast_png_file, trait=trait, print.legend=T, path=T, save.output=T) {
  if (path==T) { # if you were provided a path, read it in
     Obs_Prd_m <- read.table(obs_prd_file, sep = "\t", header = T); # pull in LOO prediction results
  } else {Obs_Prd_m <- obs_prd_file} #otherwise, just use the file as the obs_prd_m
 
 # Obs_Prd_m <- Obs_Prd_m[!is.na(Obs_Prd_m$Obs_trait),]; # drop missing observations
 Obs_Prd_m <- Obs_Prd_m[!is.na(Obs_Prd_m[,3]),]; # drop missing observations
 prd_env <- as.vector(sort(unique(Obs_Prd_m$env_code))); # pull environment names
 env_rs <- matrix(ncol = 1, nrow = length(prd_env));
 
 for (e_i in 1:length(prd_env)) { # loop through environments
 # for (e_i in 1:2) { # loop through environments
   env_obs_prd <- subset(Obs_Prd_m, Obs_Prd_m$env_code == prd_env[e_i]); # pull data for current environment
   if (nrow(env_obs_prd) > 0) { # make sure data is available
    env_rs[e_i,] <- c( 
      # sprintf( "%.2f", cor(env_obs_prd[,3], env_obs_prd[,5], 
      #                                       use = "complete.obs")), # corr btw value predicted with pheno mean, and observed
                       sprintf( "%.2f", cor(env_obs_prd[,3], env_obs_prd[,7],
                                            use = "complete.obs")))
                           # sprintf( "%.2f", cor(env_obs_prd[,5], env_obs_prd[,4],
                           #                  use = "complete.obs")))
                       # , # corr btw value predicted with PTT etc, and observed
                       # sprintf( "%.2f", cor(env_obs_prd[,6], env_obs_prd[,5], use = "complete.obs"))); # corr btw mean value for line in all other environs, and observed
   }
    
 }
 
 # xy_lim <- range(Obs_Prd_m[,3:5],na.rm = T) # set limits for figure
  obs_prd_m <- Obs_Prd_m

 # Make figure:
  if (save.output==T) {
 png(forecast_png_file, width = 4/1.5, height= 4/1.5,pointsize=6, units = "in", res = 600)
 # layout(matrix(c(1:4), 2, 2, byrow = T));
  }
 
 
  # Plot B: value predicted by environmental parameter (PTT etc) vs observed
   xy_limB <- range(obs_prd_m[,c(3,7)], na.rm=T)
  
  all_env_r <- obs_prd_m %>% # pull the correlations for each environment individually
      group_by(env_code) %>%
      mutate(env.cor=sprintf("%.2f", cor(y.hat, get(colnames(obs_prd_m)[3]), use = "complete.obs"))) %>%
      select(env_code, env.cor) %>%
      distinct%>%
      mutate(legend=paste0(env_code, " r=", env.cor)) %>%
    arrange(env_code)
  
  par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
  plot(obs_prd_m[[7]], obs_prd_m[[3]], # plot value predicted with environmental parameter (PTT etc), vs observed
  # plot(obs_prd_m[,7], obs_prd_m[,3], # NOTE: this is original version for CERIS-JGRA, I changed to above for FR-gBLUP
      # col = env_cols[match(obs_prd_m[,1], all_env_codes)], pch = 19, cex = .4,  # also for CERIS-JGRA
      col = env_cols[match(obs_prd_m[[1]], all_env_codes)], pch = 19, cex = .4, # for FR-gBLUP
      xlab = paste('Predicted ', trait, ' by ', kPara_Name, sep = ''), 
       ylab = paste('Observed ', trait, '', sep = ''), 
       xlim = xy_limB, ylim = xy_limB);
  abline(a = 0, b = 1, lty = 2, col = "gray59");
  # mtext('B', side = 3, at = xy_limB[1], line = 0.1, cex = .8);
  r2 <- sprintf("%.2f", cor(obs_prd_m[,3], obs_prd_m[,7], use = "complete.obs")); # add correlation
  # r2 <- sprintf("%.2f", cor(obs_prd_m[[3]], obs_prd_m[[7]], use = "complete.obs")); # add correlation
  legend("top", legend= substitute(paste("overall ", italic('r'), " = ", R1), list(R1 = r2)), bty = "n")
  # legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r2)), bty = "n")
  # add legend with environmental codes
   if (print.legend) {
     legend(x=min(xy_limB), y=max(xy_limB), all_env_r$legend, pch = 19, col = env_cols, bty = "n" )
     }
 
   if (save.output==T) {
 dev.off()
}
}