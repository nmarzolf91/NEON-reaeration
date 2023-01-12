# plot conductivity time-series from S1 and S4


plot_reaeration_cond <- function() {
  dir <- 'data/raw/Reaeration field and lab collection/'
  
  sites <- list.files(dir)
  
  for(s in 1:length(sites)){
    site <- sites[s]
    
    # subset to only S1 and S4 conductivity probes
    locs <- c(glue(site, '.AOS.reaeration.station.01'),
              glue(site, '.AOS.reaeration.station.04'))
    
    cond_dat <- read_feather(glue(dir, site,'/rea_conductivityFieldData.feather')) 
    
    cond_dat <- cond_dat %>% 
      arrange(dateTimeLogger) %>% 
      filter(namedLocation %in% locs) %>% 
      mutate(eventID = paste(siteID, lubridate::date(startDate), sep = '_'),
             sens_pos = str_split_fixed(hoboSampleID, '_',n=3)[,2])
    
    for(d in 1:length(unique(cond_dat$eventID))){
      
      event <- unique(cond_dat$eventID)[d]
      
      g <- cond_dat %>% 
        dplyr::filter(eventID %in% event) %>% 
        ggplot(.,
               aes(x = dateTimeLogger,
                   y = lowRangeSpCondNonlinear))+
        geom_point()+
        facet_wrap(sens_pos ~ .,
                   scales = 'free')+
        labs(x = 'Timestamp',
             y = 'Sp Cond (us/cm)')+
        ggtitle(event)
      
      ggsave(plot = g, 
             filename = glue('figures/ts_cond/{site}/{event}.png'),
             width = 6, height = 3)
  
    } # end for loop: event
  } # end for loop: site
} # end function

