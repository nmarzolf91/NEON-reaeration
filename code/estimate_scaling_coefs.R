estimate_scaling_coefs <- function(dir = 'data/derived/hydraulics') {
  
  files <- list.files(dir, full.names = TRUE)
  
  out <- matrix(ncol = 18,
                nrow = length(files))
  colnames(out) <- c('site', 'a', 'a_se', 'b', 'b_se', 'c', 'c_se', 'd', 'd_se', 'e', 'e_se',
                     'f', 'f_se','r2_width', 'r2_vel', 'r2_depth', 'prod_coefs', 'sum_coefs')

  for(i in 1:length(files)) {
    
    dat <- readr::read_csv(files[i]) %>% 
      filter(mean_velocity_ms > 0)
    
    site <- str_split(files[i], '/')[[1]][4] %>% 
      substr(.,1,4)
    
    lm_depth <- lm(data = dat %>% 
                     filter(is.finite(log(mean_depth_m)),
                            is.finite(log(calcQ_m3s))
                            ),
                   log(mean_depth_m) ~ log(calcQ_m3s))
    c <- as.numeric(summary(lm_depth)$coefficients[1])
    c_se <- as.numeric(summary(lm_depth)$coefficients[3])
    d <- as.numeric(summary(lm_depth)$coefficients[2])
    d_se <- as.numeric(summary(lm_depth)$coefficients[4])
    r2_depth <- summary(lm_depth)$adj.r.squared
    
    lm_width <- lm(data = dat %>% 
                     filter(is.finite(log(mean_width_m)),
                            is.finite(log(calcQ_m3s))
                     ),
                   log(mean_width_m) ~ log(calcQ_m3s))
    a <- as.numeric(summary(lm_width)$coefficients[1])
    a_se <- as.numeric(summary(lm_width)$coefficients[3])
    b <- as.numeric(summary(lm_width)$coefficients[2])
    b_se <- as.numeric(summary(lm_width)$coefficients[4])
    r2_width <- summary(lm_width)$adj.r.squared
    
    lm_vel <- lm(data = dat %>% 
                   filter(is.finite(log(mean_velocity_ms)),
                          is.finite(log(calcQ_m3s))
                   ),
                 log(mean_velocity_ms) ~ log(calcQ_m3s))
    e <- as.numeric(summary(lm_vel)$coefficients[1])
    e_se <- as.numeric(summary(lm_vel)$coefficients[3])
    f <- as.numeric(summary(lm_vel)$coefficients[2])
    f_se <- as.numeric(summary(lm_vel)$coefficients[4])
    r2_vel <- summary(lm_vel)$adj.r.squared
    
    sum_coefs <- b+d+f
    prod_coefs <- a*c*e
    
    out[i,1] <- site
    out[i,2] <- a
    out[i,3] <- a_se
    out[i,4] <- b
    out[i,5] <- b_se
    out[i,6] <- c
    out[i,7] <- c_se
    out[i,8] <- d
    out[i,9] <- d_se
    out[i,10] <- e
    out[i,11] <- e_se
    out[i,12] <- f
    out[i,13] <- f_se
    out[i,14] <- r2_width
    out[i,15] <- r2_vel
    out[i,16] <- r2_depth
    out[i,17] <- prod_coefs
    out[i,18] <- sum_coefs
    
    plot <- dat %>% 
      select(-est_mean_width_m) %>% 
      pivot_longer(cols = 2:4,
                   names_to = 'meas',
                   values_to = 'val') %>% 
      ggplot(.,
             aes(x = log(calcQ_m3s),
                 y = log(val)))+
      geom_point()+
      geom_smooth(method = 'lm')+
      # scale_x_log10()+
      # scale_y_log10()+
      facet_grid(meas ~ ., 
                 switch = 'both')+
      theme_bw()+
      labs(x = 'Discharge (m3/s)',
           y = 'Width (m), Velocity (m/s), Depth (m)')+
      ggtitle(site)
    
    fig_dir <- 'figures'
    if(!dir.exists(fig_dir))
      dir.create(fig_dir)
    
    fig_save_dir <- 'figures/hydraulic_scaling'
    if(!dir.exists(fig_save_dir))
      dir.create(fig_save_dir)
    
    ggsave(plot = plot,
           glue('figures/hydraulic_scaling/{site}_hydraulicScaling.png'),
           height = 6, width = 4)
    
  } # end for loop
  
  readr::write_csv(format(data.frame(out), digits = 3),
                   'data/derived/scaling_coefs/NEON_site_scaling_coefs.csv')
  
  return(out)
  
} # end function
