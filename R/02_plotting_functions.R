
#' plot_fitsim
#'
#' Function that produces Figure to visualize model fit.
#'
#' @param fit modelfit, result of "sampling_model"
#' @param data data to which model was fit, vector with value per week
#' @param params list of parameters used to fit model
#'
#' @return
#' @export
#'
#' @examples
plot_fitsim = function(fit, data, params) {

  # posterior predictive check
  dat_ = data.frame(variable = "I_t_simulated", time = 1:length(data), null = "", median = data )
  pred_ = fit$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  post_ = fit$samples_prior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  g1 = ggplot() +
    geom_ribbon(data=post_,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Prior predictive check"),alpha=.3) +
    geom_ribbon(data=pred_,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior predictive check"),alpha=.3) +
    geom_line(data=pred_,aes(x=as.numeric(time),y=median),colour=col_post) +
    geom_point(data=dat_,aes(x=as.numeric(time),y=median,colour="Simulated data")) +
    scale_fill_manual(values=c("Prior predictive check"=col_prior,
                               "Posterior predictive check"=col_post)) +
    scale_colour_manual(values=c("Simulated data"=col_simul)) +
    #scale_y_continuous(trans='log10') +
    labs(x="Time (weeks)",y="Incidence cases",colour=NULL,fill=NULL) +
    theme(legend.position = "bottom", text=element_text(size=16))

  if(grepl(",",dat_$time[1])) {
    table_age_sex = tibble(group=as.character(1:10),
                           age=rep(paste0("age_group_",1:5),2),
                           sex=rep(paste0("sex_",1:2),each=5))
    dat_ = tidyr::separate(dat_,time,",",into=c("time","group")) %>% dplyr::left_join(table_age_sex,by="group")
    pred_ = tidyr::separate(pred_,time,",",into=c("time","group")) %>% dplyr::left_join(table_age_sex,by="group")
    post_ = tidyr::separate(post_,time,",",into=c("time","group")) %>% dplyr::left_join(table_age_sex,by="group")
    g1 = ggplot() +
      geom_ribbon(data=post_,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Prior predictive check"),alpha=.3) +
      geom_ribbon(data=pred_,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior predictive check"),alpha=.3) +
      geom_line(data=pred_,aes(x=as.numeric(time),y=median),colour=col_post) +
      geom_point(data=dat_,aes(x=as.numeric(time),y=median,colour="Simulated data")) +
      facet_grid(sex~age) +
      scale_fill_manual(values=c("Prior predictive check"=col_prior,
                                 "Posterior predictive check"=col_post)) +
      scale_colour_manual(values=c("Simulated data"=col_simul)) +
      labs(x="Time (weeks)",y="Prevalent cases",colour=NULL,fill=NULL) +
      theme(legend.position = "bottom", text=element_text(size=16))
  }


  # compare prior, data and posterior
  params$gamma = 1./(params$generation_time*(1-params$fraction_pre))

  R0 = (inv.logit(params$a[1])*params$contact)/params$gamma

  sim_ = as_tibble(rbind( "R0" = R0, "I0"=params$I0, "theta"=params$theta, "p_detect2" = params$p_detect2))
  colnames(sim_) <- "median"

  sim_ = sim_ %>%
    dplyr::mutate(variable=median,
                  type="Chosen value",
                  median=median,q5=median,q95=median,
                  variable = c("R0", "I0", "theta", "p_detect2" ))

  sim_ = sim_ %>% mutate_at( "variable", as.character )

  post_ = fit$samples_posterior$summary(c("R0","I0_raw", "theta",  "p_detect2"
  )) %>%
    dplyr::mutate(type="Posterior distribution", median = median, q5 = q5, q95=q95)

  prior_ = fit$samples_prior$summary(c("R0","I0_raw", "theta",  "p_detect2"))  %>%
    dplyr::mutate(type="Prior distribution")

  dat_g2 <- dplyr::bind_rows(sim_,post_,prior_) %>%
    dplyr::mutate(type=factor(type,levels=c("Posterior distribution","Chosen value","Prior distribution")),
                  variable=gsub("_raw","",variable))

  g2 =  ggplot( dat_g2 ) +
    geom_pointrange(aes(x=type,y=median,ymin=q5,ymax=q95,colour=type,shape=type)) +
    facet_wrap(~variable,scales="free",ncol=3) +
    scale_shape_manual(values=c(16,15,95),guide="none") +
    scale_colour_manual(values=c(col_post,col_simul,col_prior),guide="none") +
    coord_flip() +
    labs(x=NULL,y="Parameter values (median and 90% CrI)")

  full_plot <- cowplot::plot_grid(g1,g2,nrow=2,labels=LETTERS)

  full_plot
  return(full_plot)
}


#' plot_comparesim
#'
#' Function that produces Figure to visualize and compare three different model fits.
#'
#' @param fit1 modelfit 1 of 3, result of "sampling_model"
#' @param fit2 modelfit 2 of 3, result of "sampling_model"
#' @param fit3 modelfit 3 of 3, result of "sampling_model"
#' @param data data to which model was fit, vector with value per week
#'
#' @return plot
#' @export
#'
#' @examples
plot_comparesim = function(fit1, fit2, fit3 , data) {
  # posterior predictive check
  dat_ = fit1$samples_posterior$summary(c("I_t_simulated")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  dat_$date = data$date

  pred_1 = fit1$samples_prior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  post_1 = fit1$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  pred_1$date = data$date
  post_1$date = data$date

  pred_2 = fit2$samples_prior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  post_2 = fit2$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  pred_2$date = data$date
  post_2$date = data$date

  pred_3 = fit3$samples_prior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  post_3 = fit3$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  pred_3$date = data$date
  post_3$date = data$date

  g1 = ggplot() +
    #geom_ribbon(data=pred_1,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Prior predictive check"),alpha=.3) +
    geom_ribbon(data=post_1,aes(x=date,ymin=q5,ymax=q95,fill="Posterior predictive check"),alpha=.3) +
    geom_line(data=post_1,aes(x=date,y=median),colour=col_post) +
    geom_point(data=dat_,aes(x=date,y=median,colour="Case data Bern")) +
    scale_fill_manual(values=c(#"Prior predictive check"=col_prior,
      "Posterior predictive check"=col_post)) +
    scale_colour_manual(values=c("Case data Bern"=col_simul)) +
    #scale_y_continuous(trans='log10') +
    labs(x="Time (weeks)",y="# pos. tests",colour=NULL,fill=NULL, title="Poisson model") +
    theme(legend.position = "bottom",text=element_text(size=16))

  g2 = ggplot() +
    #geom_ribbon(data=pred_2,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Prior predictive check"),alpha=.3) +
    geom_ribbon(data=post_2,aes(x=date,ymin=q5,ymax=q95,fill="Posterior predictive check"),alpha=.3) +
    geom_line(data=post_2,aes(x=date,y=median),colour=col_post) +
    geom_point(data=dat_,aes(x=date,y=median,colour="Case data Bern")) +
    scale_fill_manual(values=c(#"Prior predictive check"=col_prior,
      "Posterior predictive check"=col_post)) +
    scale_colour_manual(values=c("Case data Bern"=col_simul)) +
    #scale_y_continuous(trans='log10') +
    labs(x="Time (weeks)",y="# pos. tests",colour=NULL,fill=NULL, title = " Quasi Poisson model") +
    theme(legend.position = "bottom", text=element_text(size=16))

  g3 = ggplot() +
    #geom_ribbon(data=pred_3,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Prior predictive check"),alpha=.3) +
    geom_ribbon(data=post_3,aes(x=date,ymin=q5,ymax=q95,fill="Posterior predictive check"),alpha=.3) +
    geom_line(data=post_3,aes(x=date,y=median),colour=col_post) +
    geom_point(data=dat_,aes(x=date,y=median,colour="Case data Bern")) +
    scale_fill_manual(values=c(#"Prior predictive check"=col_prior,
      "Posterior predictive check"=col_post)) +
    scale_colour_manual(values=c("Case data Bern"=col_simul)) +
    #scale_y_continuous(trans='log10') +
    labs(x="Time (weeks)",y="# pos. tests",colour=NULL,fill=NULL, title="Negative binomial model") +
    theme(legend.position = "bottom",text=element_text(size=16))

  # plot the fit of the probability of transmission over time
  prob_p = fit1$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  prob_qp = fit2$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  prob_nb = fit3$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))

  prob_plot <- ggplot() +
    geom_ribbon(data=prob_p,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Poisson model"),alpha=.3) +
    geom_ribbon(data=prob_qp,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Quasi Poisson model"),alpha=.3) +
    geom_ribbon(data=prob_nb,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Negative Binomial model"),alpha=.3) +
    geom_line(data=prob_p,aes(x=as.numeric(time),y=median,colour="Poisson model")) +
    geom_line(data=prob_qp,aes(x=as.numeric(time),y=median,colour="Quasi Poisson model")) +
    geom_line(data=prob_nb,aes(x=as.numeric(time),y=median,colour="Negative Binomial model")) +
    scale_fill_manual(values=c("Poisson model"="Orange",
                               "Quasi Poisson model"="Red",
                               "Negative Binomial model"="blue")) +
    scale_colour_manual(values=c("Poisson model"="Orange",
                                 "Quasi Poisson model"="Red",
                                 "Negative Binomial model"="blue")) +
    labs(x="Time (weeks)",y="Probability of transmission upon contact",colour=NULL,fill=NULL) +
    theme(legend.position = "bottom",text=element_text(size=16))

  legend <- get_legend(g1 + theme(legend.box.margin = margin(0, 0, 0, 12), legend.position = "bottom") )

  prow <- cowplot::plot_grid(g1 + theme(legend.position="none"),
                             g2 + theme(legend.position="none"),
                             g3 + theme(legend.position="none"), nrow=3,labels=LETTERS)

  full_plot = plot_grid(prow, legend, nrow =2, rel_widths = c(3, .4), rel_heights = c(3,.2))

  return(list(full_plot, prob_plot))
}

#' plot_fitsim_strat_realdata
#'
#' Function that produces Figure to visualize stratified model fit to data from Geneva.
#'
#' @param fit modelfit of stratified model, result of "sampling_model"
#' @param data data to which model was fit, array with value per week and age-class
#' @param cum_inc vector of cumulative incidence over time
#' @param seroprev seroprevalence data used for model fit
#'
#' @return Figure
#' @export
#'
#' @examples
plot_fitsim_strat_realdata = function(fit, data, cum_inc, seroprev) {

  # posterior predictive check
  dat_ = data.frame(variable = "I_t_simulated", time = 1:dim(data)[1], null = "", median = data )
  post_ = fit$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) %>%
    tidyr::separate(step,",",into=c("age_group","time"))
  prior_ = fit$samples_prior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) %>%
    tidyr::separate(step,",",into=c("age_group","time"))

  post_$date = rep(dat_$median.date, each = 3)
  prior_$date = rep(dat_$median.date, each = 3)
  dat_$date = as.Date(dat_$median.date)

  g1 = ggplot() +
    geom_rect(aes(xmin = as.Date("2020-04-06"), xmax = as.Date("2020-05-10"), ymin = -Inf, ymax = Inf),
              alpha = I(.5), fill = I("lightblue")) +
    geom_rect(aes(xmin = as.Date("2020-06-02"), xmax = as.Date("2020-06-30"), ymin = -Inf, ymax = Inf),
              alpha = I(.5), fill = I("lightblue")) +
    geom_rect(aes(xmin = as.Date("2020-11-20"), xmax = as.Date("2020-12-31"), ymin = -Inf, ymax = Inf),
              alpha = I(.5), fill = I("lightblue")) +
    geom_point(data=dat_ ,aes(x=date,y=median.cases_age1, colour="0-19")) +
    geom_point(data=dat_ ,aes(x=date,y=median.cases_age2, colour="20-59")) +
    geom_point(data=dat_ ,aes(x=date,y=median.cases_age3, colour="60+")) +
    #geom_ribbon(data=post_ %>% filter(age_group==2), aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Prior predictive check"),alpha=.3) +
    geom_ribbon(data=post_ %>% filter(age_group==1), aes(x=date,ymin=q5,ymax=q95,fill="0-19"),alpha=.3) +
    geom_line(data=post_ %>% filter(age_group==1),aes(x=date,y=median,colour="0-19")) +
    geom_ribbon(data=post_ %>% filter(age_group==2), aes(x=date,ymin=q5,ymax=q95,fill="20-59"),alpha=.3) +
    geom_line(data=post_ %>% filter(age_group==2),aes(x=date,y=median,colour="20-59")) +
    geom_ribbon(data=post_ %>% filter(age_group==3), aes(x=date,ymin=q5,ymax=q95,fill="60+"),alpha=.3) +
    geom_line(data=post_ %>% filter(age_group==3),aes(x=date,y=median,colour="60+")) +
    scale_fill_manual(values=c("0-19"=col_simul,
                               "20-59"="orange",
                               "60+"="red")) +
    scale_colour_manual(values=c("0-19"=col_simul,
                                 "20-59"="orange",
                                 "60+"="red" )) +
    #scale_y_continuous(trans='log10') +
    labs(x="",y="PCR confirmed cases per week",colour=NULL,fill=NULL) +
    theme(legend.position = "bottom", text=element_text(size=16))+
    scale_x_date(date_labels = "%b %Y")

  post_prob = fit$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) %>%
    tidyr::separate(step,",",into=c("age_group","time"))
  prior_prob = fit$samples_prior$summary(c("prob_infection")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) %>%
    tidyr::separate(step,",",into=c("age_group","time"))

  post_prob$date = as.Date(rep(dat_$median.date, each = 3))
  prior_prob$date = as.Date(rep(dat_$median.date, each = 3))

  g1prob = ggplot() +
    geom_rect(aes(xmin = as.Date("2020-03-13"), xmax = as.Date("2020-04-04"), ymin = -Inf, ymax = Inf),
              alpha = I(.5), fill = I("lightblue"))  +
    geom_rect(aes(xmin = as.Date("2020-10-18"), xmax = as.Date("2020-12-31"), ymin = -Inf, ymax = Inf),
              alpha = I(.5), fill = I("lightblue")) +
    geom_ribbon(data=post_prob %>% filter(age_group==1),aes(x=date, ymin=q5,ymax=q95,fill="0-19"),alpha=.3) +
    geom_ribbon(data=post_prob %>% filter(age_group==2),aes(x=date, ymin=q5,ymax=q95,fill="20-59"),alpha=.3) +
    geom_ribbon(data=post_prob %>% filter(age_group==3),aes(x=date, ymin=q5,ymax=q95,fill="60+"),alpha=.3) +
    #geom_ribbon(data=pred_,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior predictive check"),alpha=.3) +
    geom_line(data=post_prob%>% filter(age_group==1),aes(x=date,y=median,colour="0-19")) +
    geom_line(data=post_prob%>% filter(age_group==2),aes(x=date,y=median,colour="20-59")) +
    geom_line(data=post_prob%>% filter(age_group==3),aes(x=date,y=median,colour="60+")) +
    scale_fill_manual(values=c("0-19"=col_simul,
                               "20-59"="orange",
                               "60+"="red")) +
    scale_colour_manual(values=c("0-19"=col_simul,
                                 "20-59"="orange",
                                 "60+"="red" )) +
    labs(x="",y="rho",colour=NULL,fill=NULL) +
    theme(legend.position = "bottom",text=element_text(size=16))+
    scale_x_date(date_labels = "%b %Y")

  post_var = fit$samples_posterior$summary(c("I0_raw", "theta","p_detect1", "p_detect2", "p_detect3" )) %>%
    dplyr::mutate(type="Posterior distribution", median = median, q5 = q5, q95=q95)

  prior_var = fit$samples_prior$summary(c("I0_raw", "theta", "p_detect1", "p_detect2", "p_detect3" ))  %>%
    dplyr::mutate(type="Prior distribution")

  dat_g2 <- dplyr::bind_rows(post_var,prior_var) %>%
    dplyr::mutate(type=factor(type,levels=c("Posterior distribution","Chosen value","Prior distribution")),
                  variable=gsub("_raw","",variable))

  g2 =  ggplot( dat_g2 ) +
    geom_pointrange(aes(x=type,y=median,ymin=q5,ymax=q95,colour=type,shape=type)) +
    facet_wrap(~variable,scales="free",ncol=3) +
    scale_shape_manual(values=c(16,15,95),guide="none") +
    scale_colour_manual(values=c(col_post,col_simul,col_prior),guide="none") +
    coord_flip() +
    labs(x=NULL,y="Parameter values (median and 90% CrI)")

  # Special for Geneva

  post_detect = fit$samples_posterior$summary(c("p_detect1", "p_detect2", "p_detect3" )) %>%
    dplyr::mutate(type="Posterior distribution", median = median, q5 = q5, q95=q95)

  post_detect$age_group = rep(c("0-19", "20-59", "60+"), 3)
  post_detect$time_period = rep(1:3, each = 3) # c(0.8, 1, 1.2, 1.8, 2, 2.2, 2.8, 3, 3.2)

  g3 =  ggplot( post_detect ) +
    geom_pointrange(aes(x=as.character(time_period),y=median,ymin=q5,ymax=q95,colour=age_group),
                    position = position_dodge( width=0.5) ) +
    labs(x="Time period",y="Ascertainment rate (median and 90% CrI)") +
    scale_colour_manual(values=c("0-19"=col_simul,
                                 "20-59"="orange",
                                 "60+"="red" )) +
    theme(text=element_text(size=16))

  g4 =  ggplot( cum_inc ) +
    geom_pointrange(aes(x=as.character(time_period),y=median,ymin=q5,ymax=q95,colour=age_group),
                    position = position_dodge( width=0.5) ) +
    labs(x="Time period",y="Cumulative incidence (median and 90% CrI)") +
    scale_colour_manual(values=c("0-19"=col_simul,
                                 "20-59"="orange",
                                 "60+"="red" )) +
    scale_y_continuous(trans='log10') +
    theme(text=element_text(size=16))

  full_plot <- cowplot::plot_grid(g1+ theme(legend.position = "none"),g1prob, g3+ theme(legend.position = "none"),
                                  nrow=2,labels=LETTERS, byrow = FALSE, rel_widths = c(2,1))

  return(full_plot)
}

#' plot_compare
#'
#' @param fitBM
#' @param fitSpline
#' @param fitGP
#' @param data
#' @param params
#' @param transmission_prob_sim
#'
#' @return
#' @export
#'
#' @examples
plot_compare = function(fitBM, fitSpline, fitGP, data, params,transmission_prob_sim) {

  # posterior predictive check
  dat_ = data.frame(variable = "I_t_simulated", time = 1:length(data), null = "", median = data )
  pred_BM = fitBM$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  pred_spline = fitSpline$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  pred_GP = fitGP$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  g1 = ggplot() +
    geom_ribbon(data=pred_BM,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior -- BM model"),alpha=.3) +
    geom_ribbon(data=pred_spline,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior -- Spline model"),alpha=.3) +
    geom_ribbon(data=pred_GP,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior -- GP model"),alpha=.3) +
    geom_line(data=pred_BM,aes(x=as.numeric(time),y=median),colour="Blue") +
    geom_line(data=pred_spline,aes(x=as.numeric(time),y=median),colour="Orange") +
    geom_line(data=pred_GP,aes(x=as.numeric(time),y=median),colour="Red") +
    geom_point(data=dat_,aes(x=as.numeric(time),y=median,colour="Simulated data")) +
    scale_fill_manual(values=c("Posterior -- BM model"="Blue",
                               "Posterior -- Spline model"="Orange",
                               "Posterior -- GP model"="Red")) +
    scale_colour_manual(values=c("Simulated data"=col_simul)) +
    #scale_y_continuous(trans='log10') +
    labs(x="Time (weeks)",y="Incidence cases",colour=NULL,fill=NULL) +
    theme(legend.position = "bottom", text=element_text(size=16))

  dat_prob = data.frame(variable = "Transmission_probability", time = 1:length(data), null = "", median = transmission_prob_sim )
  pred_BM_prob = fitBM$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  pred_spline_prob = fitSpline$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  pred_GP_prob = fitGP$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
  g2 = ggplot() +
    geom_ribbon(data=pred_BM_prob,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior -- BM model"),alpha=.3) +
    geom_ribbon(data=pred_spline_prob,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior -- Spline model"),alpha=.3) +
    geom_ribbon(data=pred_GP_prob,aes(x=as.numeric(time),ymin=q5,ymax=q95,fill="Posterior -- GP model"),alpha=.3) +
    geom_line(data=pred_BM_prob,aes(x=as.numeric(time),y=median),colour="blue") +
    geom_line(data=pred_spline_prob,aes(x=as.numeric(time),y=median),colour="Orange") +
    geom_line(data=pred_GP_prob,aes(x=as.numeric(time),y=median),colour="Red") +
    geom_line(data=dat_prob,aes(x=as.numeric(time),y=median,colour="Simulated data")) +
    scale_fill_manual(values=c("Posterior -- BM model"="Blue",
                               "Posterior -- Spline model"="Orange",
                               "Posterior -- GP model"="Red")) +
    scale_colour_manual(values=c("Simulated data"=col_simul)) +
    labs(x="Time (weeks)",y="rho",colour=NULL,fill=NULL) +
    theme(legend.position = "bottom",text=element_text(size=16))

  full_plot <- cowplot::plot_grid(g1,g2,nrow=2,labels=LETTERS)
  full_plot
  return(full_plot)
}
