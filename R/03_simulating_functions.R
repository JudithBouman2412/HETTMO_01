# functions to simulate data using the spline-based time-implementation

#' set_parameters
#'
#' Function to create set of parameters to simulate data as presented in the
#' publication.
#'
#' @param stratified should the simulated data be stratified by age?
#' @param contact_all contact data to calculate the contact matrix using Prem et al (2021)
#'
#' @return list of parameters
#' @export
#'
#' @examples
set_parameters <- function(stratified = FALSE,
                           contact_all){

  # Set of parameters general for the various type of simulations
  params = list( order=4, # order of the simulated splines
                 knots = seq(0,48,4), # defined knot sequence to simulate spline-based transmission rate
                 I0 = 3, # number of infected individuals at t=0
                 generation_time = 5.2/7, # assumed generation time per week
                 fraction_pre = 0.5, # distribution of generation time over I and E compartment
                 theta = 5, # overdispersion, based on quasi-poisson sampling distribution
                 t_detectSwitch = 25 # time at which ascertainment rate switches
                 )

  if (stratified){ # parameters specific for stratified data by age
    strat_pop <- create_contactmatrix_GE(contact_all) # create contact matrix by age
    params$p_detect1 = 0.3 # ascertainment rate at first time period
    params$p_detect2 = 0.7 # ascertainment rate at second time period
    params$popsize <- strat_pop[[2]] # population distribution across age-classes
    params$contact <- strat_pop[[1]] *7 # contact matrix
    params$n_tested_survey <- c(1000,2000,2000) # number of serological tests per age class in seroprevalence study
    params$beta_fixed <- 0.052 # baseline probability of transmission per contact
    params$a <- rbind(c(rep(1,3), 0.8,  rep((0.5),2), rep((1),7), rep(0.6,2) ),
                         c(rep(1,3), 0.8, rep((0.6),3), rep((1),6), rep(0.6,2) ),
                         c(rep(1,2), rep(0.5, 2),  rep((0.5),4), rep((1),5), rep(0.6,2) )) # coefficients for spline to construct rho(t)
    params$num_class <- 3 # number of age-categories

  } else {
    params$beta_fixed = 0.0638 # represents the fixed probability of transmission per contact
    params$a = c(rep((0.058),3), (0.05), (0.03), rep((0.025),1), (0.04), rep((0.05),3), rep((0.03),5) )/0.058 # coefficients for spline to construct rho(t)
    params$popsize = 100000 # simulated population size
    params$contact = 77 # average number of contacts per week
    params$p_detect1 = 0.1 # ascertainment rate at first time period
    params$p_detect2 = 0.5 # ascertainment rate at second time period
    params$n_tested_survey = 0.05 * 100000 # number of serological tests in seroprevalence study
  }

  return(params)
}

#' create_contactmatrix_GE
#'
#' Function to construct a contact matrix based on Prem et al. (2021) for Switzerland and the age
#' distribution of the Canton of Geneva (2021). Age groups are fixed at 0-19, 20-59
#' and 60+.
#'
#' @param contact_all contact data to calculate the contact matrix using Prem et al (2021)
#'
#' @return List with contact matrix (1) and the simulated population (2)
#' @export
#'
#' @examples
create_contactmatrix_GE <- function(contact_all){
  # create contact matrix for three groups
  pop_2021 <- readxl::read_excel("data/su-d-01.02.03.06.xlsx", sheet = "2021")
  pop_2021 <- data.frame(age = 0:100, n = as.numeric(pop_2021[2, 3:103]))
  age_sent <- c(0, 19, 59) # Age groups for Sentinella
  pop <- numeric(length(age_sent))
  age_range <- c(age_sent, 200)

  for(i in 1:length(pop)) {
    pop[i] <- sum(pop_2021$n[pop_2021$age >= age_range[i] & pop_2021$age < age_range[i + 1]])
  }
  pop_2021 <- data.frame(lower.age.limit = age_sent, population = pop)

  sim_pop <- round(300000*pop_2021$population/sum(pop_2021$population),0)

  pop_final <- data.frame(lower.age.limit=age_sent, population=sim_pop)

  # Create hypothetical contact matrix for baseline contacts switzerland
  prem <- contact_all$CHE

  prem_adjusted = matrix(NA, 3, 3)
  group_list = c(rep(1,4),
                 rep(2,8),
                 rep(3,4))

  for (i in 1:3){
    for (j in 1:3){
      prem_adjusted[i,j] = sum(prem[group_list==i, group_list==j])/length(group_list[group_list==j])
    }
  }

  return(list(prem_adjusted, sim_pop))
}

#' Function with system of ODE equations for SEIR model
#'
#' @param time time in weeks
#' @param y number of individuals in each compartment at time
#' @param params set of parameters for SEIR model
#'
#' @return list of number of individuals per compartment
#' @export
#'
#' @examples
ODE_simple <- function(time, y, params) {
  g <- with(as.list(c(y,params)),
            c( - beta_fixed*beta*contact*(I)*S/popsize,
               (beta_fixed*beta*contact*I*(S))/popsize - tau*E,
               tau*E - gamma*I,
               gamma*I )
  )
  return(list(g))
}

#' Function with system of ODE equations for SEIR model with spline-based
#' time-dependent transmission rate
#'
#' @param time time in weeks
#' @param y number of individuals in each compartment at time
#' @param params set of parameters as generated by "set_parameters"
#'
#' @return list of number of individuals per compartment
#' @export
#'
#' @examples
ODE_spline <- function(time, y, params) {


  beta_time <- analytical_bspline( b_hat = params$b_hat , a = params$a,
                                   b_support = params$b_support, t_new = time, knots = params$knots )

  g <- with(as.list(c(y,params)),
            c( -beta_fixed*(beta_time)*contact*(I)*S/popsize,
               (beta_fixed*(beta_time)*contact*I*(S))/popsize - tau*E,
               tau*E - gamma*I,
               gamma*I )
  )
  return(list(g))
}

#' Function to determine the number of the equation corresponding to given age class
#'
#' @param comp required compartment of SEIR model (1=S, 2=E, 3=I, 4 = R)
#' @param age_sex required age-classe (1= 0-19, 2=20-59, 3=60+)
#' @param num_class total number of modelled age-classes
#'
#' @return number of the ODE function corresponding to "comp" and "age_sex"
#' @export
#'
#' @examples
ind <- function(comp, age_sex, num_class=3){

  return(age_sex+(comp-1)*num_class)
}


#' Function with system of ODE equations for age-stratified SEIR model with
#' spline-based time-dependent transmission rate.
#'
#' @param time time in weeks
#' @param y number of individuals in each compartment at time
#' @param params set of parameters as generated by "set_parameters"
#'
#' @return list of number of individuals per compartment
#' @export
#'
#' @examples
ODE_spline_strat <- function(time, y, params) {

  num_comp <- 4
  num_class = length(y)/num_comp
  beta_time <- rep(0, num_class)
  num_eq = num_class*num_comp

  for (i in 1:num_class){
    beta_time[i] <- analytical_bspline( b_hat = params$b_hat , a = params$a[i,],
                                        b_support = params$b_support, t_new = time, knots = params$knots )
  }

  dydt = rep(0,num_eq)

  for(i in 1:num_class){

    f_inf_t = 0

    for (j in 1:num_class){
      f_inf_t = f_inf_t + params$contact[i,j]*y[ind(3,j, num_class)]/params$popsize[j]
    }

    dydt[ind(1,i, num_class)] = - params$beta_fixed* (beta_time[i]) * as.numeric(f_inf_t) * y[ind(1,i, num_class)]
    dydt[ind(2,i, num_class)] =  params$beta_fixed *(beta_time[i]) *as.numeric( f_inf_t )* y[ind(1,i, num_class)] - params$tau * y[ind(2,i, num_class)]
    dydt[ind(3,i, num_class)] =  params$tau * y[ind(2,i, num_class)] - params$gamma * y[ind(3,i, num_class)]
    dydt[ind(4,i, num_class)] =  params$gamma * y[ind(3,i, num_class)]
  }

  return(list(dydt))
}

#' get_incidence
#'
#' Function to calculate the ascertained COVID-19 incidence of simulated data using the output
#' of SEIR-ode function. Incidence is defined as those moving from compartment
#' I to R.
#'
#' @param out output of SEIR-ode model
#' @param popsize population size
#' @param p_detect1 ascertainment rate of first time-period
#' @param p_detect2 ascertainment rate of second time-period
#' @param t_detectSwitch time at which ascertainment-rate is switched
#'
#' @return vector with incidence per week
#' @export
#'
#' @examples
get_incidence <- function( out, popsize, p_detect1, p_detect2, t_detectSwitch ){

  if(!is.null(dim(out))){
    num_t = dim(out)[1]
    incidence = rep( 0, num_t )
    incidence[1] = p_detect1*(out[1,"R"]+1)

    for (i in 2:num_t){
      if (i < t_detectSwitch){
        incidence[i] = p_detect1*(out[i,"R"] - out[i-1,"R"]+1)
      } else {
        incidence[i] = p_detect2*(out[i,"R"] - out[i-1,"R"]+1)
      }
    }
  } else {
    num_t = length(out)
    incidence = rep( 0, num_t )
    incidence[1] = p_detect1*(out[1]+1)

    for (i in 2:num_t){
      if (i < t_detectSwitch){
        incidence[i] = p_detect1*(out[i,"R"] - out[i-1,"R"]+1)
      } else {
        incidence[i] = p_detect2*(out[i,"R"] - out[i-1,"R"]+1)
      }
    }
  }

  return(incidence)
}


#' get_incidence_strat
#'
#' Function to calculate the ascertained COVID-19 incidence of stratified simulated
#' data using the output of SEIR-ode function. Incidence is defined as those moving from compartment
#' I to R.
#'
#' @param R array with number of individuals per age-class per week
#' @param popsize vector of population size per age-class
#' @param p_detect1 ascertainment rate of first time-period
#' @param p_detect2 ascertainment rate of second time-period
#' @param t_detectSwitch time at which ascertainment-rate is switched
#'
#' @return array with incidence per age-class per week
#' @export
#'
#' @examples
get_incidence_strat <- function( R, popsize, p_detect1, p_detect2, t_detectSwitch ){

  if(!is.null(length(R))){
    num_t = length(R)
    incidence = rep( 0, num_t )
    incidence[1] = p_detect1*(R[1]+1)

    for (i in 2:num_t){
      if (i < t_detectSwitch){
        incidence[i] = p_detect1*(R[i] - R[i-1]+1)
      } else {
        incidence[i] = p_detect2*(R[i] - R[i-1]+1)
      }
    }
  } else {
    num_t = length(R)
    incidence = rep( 0, num_t )
    incidence[1] = p_detect1*R[1]

    for (i in 2:num_t){
      if (i < t_detectSwitch){
        incidence[i] = p_detect1*(R[i] - R[i-1]+1)
      } else {
        incidence[i] = p_detect2*(R[i] - R[i-1]+1)
      }
    }
  }

  return(incidence)
}


#' rqpois
#'
#' Function to randomly sample from a quasi-poisson distribution
#'
#' @param n number of samples required
#' @param mu mean of distribution sampling from
#' @param theta overdispersion parameter of quasi-poisson distribution
#'
#' @return
#' @export
#'
#' @examples
rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

#' simulate_data
#'
#' Function to simulate unstratified data as presented in the publication
#'
#'
#' @param params list of parameters as generated by "set_parameters"
#' @param time type of time-dependence in the simulation (2=spline) (currently only one option)
#' @param ts vector with time points at which simulated data is required
#'
#' @return list with following elements: (1) vector with number of simulated positive PCR-tests
#'  (2) vector with effective reproduction number over time
#'  (3) vector with rho(t) over time
#'  (4) vector with number of individuals in compartment I
#'  (5) vector with number of individuals in compartment E
#'  (6) vector with number of individuals in compartment R
#'  (7) number of test in the seroprevalence study that were positive
#'
#' @export
#'
#' @examples
simulate_data <- function( params,
                          time = 2,
                          ts){

  # Define gamma and tau based on generation time and fraction
  params$gamma = 1/(params$generation_time*params$fraction_pre)
  params$tau = 1/(params$generation_time*(1-params$fraction_pre))

  # run ODE with correct parameters
  if (time == 2){ # use splines to simulate data

    degree = params$order-1 # degree of the spline to fit
    num_knots = length(params$knots) # number of knots
    num_basis = num_knots + degree - 1 # total number of b-splines

    int_knots = c()
    for (i in 2:num_knots){
      int_knots = c(int_knots, seq(params$knots[i-1], params$knots[i], length.out=params$order) )
    }
    int_knots <- unique(int_knots)

    ext_knots = c(rep(params$knots[1],degree), params$knots, rep(params$knots[num_knots], degree))

    # Find B splines
    B <- create_B(int_knots, params$order, params$knots )

    # Get coefficients and their support
    B_splines <- ana_B(B, ext_knots, int_knots, params$order, params$knots, num_basis )

    params$b_hat = B_splines[[1]]
    params$b_support = B_splines[[2]]

    state = c(S = params$popsize, E = 0, I = params$I0, R = 0)

    out <- ode(y = state, times = ts, func = ODE_spline, parms = params)

    # Calculate incidence
    incidence = get_incidence(out, params$popsize, params$p_detect1, params$p_detect2, params$t_detectSwitch)

    # Calculate effective reproduction number
    beta = rep(0, length(ts))
    for (i in ts){
      beta[i] = (analytical_bspline( b_hat = params$b_hat , a = params$a,
                                     b_support = params$b_support, t_new = i, knots = params$knots ))
    }

    re <- as.vector( out[,"S"]/params$popsize * beta * params$contact / (params$gamma) )

    # Add sampling distribution -- Quasi Poisson
    simulation <- rqpois(length(ts), mu = incidence, theta = params$theta)

  }

  # calculate seroprevalence estimate
  n_infected_survey = rbinom(1, params$n_tested_survey, out[length(ts),"R"]/params$popsize )

  return(list(simulation, re, beta, out[,"I"], out[,"E"], out[,"R"], n_infected_survey ) )
}

#' simulate_strat_data
#'
#' Function to simulate stratified data as presented in the publication
#'
#' @param params list of parameters as generated by "set_parameters"
#' @param timetype of time-dependence in the simulation (2=spline) (currently only one option)
#' @param ts vector with time points at which simulated data is required
#'
#' @return list with following elements: (1) array with number of simulated positive PCR-tests per age-class
#'  (2) array with effective reproduction number over time per age-class
#'  (3) number of test in the seroprevalence study that were positive per age-class
#'  (4) array with ascertained incidence per week per age-class
#'  (5) array with output of ODE model for SEIR per age-class
#' @export
#'
#' @examples
simulate_strat_data <- function( params,
                                 time = 2,
                                 ts){

  # Define gamma and tau based on generation time and fraction
  params$gamma = 1/(params$generation_time*params$fraction_pre)
  params$tau = 1/(params$generation_time*(1-params$fraction_pre))

  # run ODE with correct parameters
  if (time == 2){ # use splines to simulate data

    degree = params$order-1 # degree of the spline to fit
    num_knots = length(params$knots) # number of knots
    num_basis = num_knots + degree - 1 # total number of b-splines
    num_class = dim(params$a)[1]

    int_knots = c()
    for (i in 2:num_knots){
      int_knots = c(int_knots, seq(params$knots[i-1], params$knots[i], length.out=params$order) )
    }
    int_knots <- unique(int_knots)

    ext_knots = c(rep(params$knots[1],degree), params$knots, rep(params$knots[num_knots], degree))

    # Find B splines
    B <- create_B(int_knots, params$order, params$knots )

    # Get coefficients and their support
    B_splines <- ana_B(B, ext_knots, int_knots, params$order, params$knots, num_basis )

    params$b_hat = B_splines[[1]]
    params$b_support = B_splines[[2]]

    state = c(S1 = (params$popsize[1]),
              S2 = (params$popsize[2]),
              S3 = (params$popsize[3]),
              E1 = 0,
              E2 = 0,
              E3 = 0,
              I1 = params$I0,
              I2 = params$I0,
              I3 = params$I0,
              R1 = 0,
              R2 = 0,
              R3 = 0)

    out <- ode(y = state, times = ts, func = ODE_spline_strat, parms = params)

    # Calculate incidence
    incidence = matrix(0, nrow=length(ts), ncol = params$num_class)
    for (i in 1:params$num_class){
      incidence[,i] = get_incidence_strat(out[,ind(4,i,params$num_class)+1], params$popsize, params$p_detect1, params$p_detect2, params$t_detectSwitch)
    }

    # Calculate effective reproduction number
    beta = matrix(0, nrow = length(ts), ncol=params$num_class)
    for (i in ts){
      for (j in 1:params$num_class){
        beta[i,j] = (analytical_bspline( b_hat = params$b_hat , a = params$a[j,],
                                         b_support = params$b_support, t_new = i, knots = params$knots ))
      }
    }

    # Add sampling distribution -- Quasi Poisson
    simulation <- matrix(0, nrow = length(ts), ncol=params$num_class)
    for (i in 1:params$num_class){
      simulation[,i] <- rqpois(length(ts), mu = incidence[,i], theta = params$theta)
    }
  }

  # simulate seroprevalence data
  n_infected_survey = rep(0,length(params$n_tested_survey))
  for (i in 1:length(params$n_tested_survey)){
    n_infected_survey[i] = rbinom(1, params$n_tested_survey[i], out[45,ind(4,i,params$num_class)+1]/params$popsize[i] )
  }

  return(list(simulation, beta, n_infected_survey, incidence, out) )
}


#' reciprocal
#'
#' Function to rebalance contactmatrix based on the number of individuals in each
#' class
#'
#' @param m original contact matrix
#' @param N vector with number of individuals per age-class
#'
#' @return rebalanced contact matrix
#' @export
#'
#' @examples
reciprocal <- function(m, N) {
  m_rec <- matrix(NA, nrow = length(N), ncol = length(N))
  for(i in 1:length(N)) {
    for(j in 1:length(N)) {
      m_rec[i, j] <- (m[i, j]*N[i] + m[j, i]*N[j])/(2*N[i])
    }
  }
  return(m_rec)
}
