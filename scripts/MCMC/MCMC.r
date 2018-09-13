args = commandArgs(trailingOnly=TRUE)

libraries = c("dplyr","magrittr","tidyr","readxl","nimble","HDInterval")
for(x in libraries) { library(x,character.only=TRUE,warn.conflicts=FALSE) }

'%&%' = function(x,y)paste0(x,y)

rho = 0.19

# 1. Incubation period

## Data
#### excluding index case
onsetTimes = c(17, 19, 21, 21, 21, 22, 22, 23, 23, 23, 23, 23, 25, 26, 26, 26, 27, 31, 39, 41, 41, 42, 44)
deathTimes = c(33, 48, 31, 47, 30, 41, 37, 32, 31, 35, 54, 30, 31, 33, 30, 33, 48, 39, 51, 48, 55, 60, 62)
## Prior knowledge
exposureTimesLower = c(5,  9,  5,  9,  5,  5,  5,  11, 5,  5,  13, 11, 12, 5,  11, 11, 5,  5,  12, 25, 5,  25, 25)
exposureTimesUpper = c(15, 18, 18, 15, 18, 18, 18, 15, 18, 18, 14, 15, 18, 18, 15, 15, 18, 18, 18, 31, 18, 31, 31)

## Initialization
nimIncubationData = list(
    onsetTime = onsetTimes,
    deathTime = deathTimes,
    exposureTimeLower = exposureTimesLower,
    exposureTimeUpper = exposureTimesUpper
)
nimIncubationConsts = list(
    nCases = length(onsetTimes),
    incubation_uncertainty = 0.5,
    death_uncertainty = 0.5
)
nimIncubationInits = list(
    incubation_mean = mean(onsetTimes-.5*(exposureTimesLower+exposureTimesUpper)),
    incubation_sd = sqrt(var(onsetTimes-.5*(exposureTimesLower+exposureTimesUpper))),
    death_mean = mean(deathTimes-onsetTimes),
    death_sd = sqrt(var(deathTimes-onsetTimes)),
    exposureTime = .5*(exposureTimesLower+exposureTimesUpper),
    exposure_mu = .5*(exposureTimesLower+exposureTimesUpper),
    exposure_sigma = (exposureTimesUpper-exposureTimesLower)/1.96/2,
    incubationTime = onsetTimes-.5*(exposureTimesLower+exposureTimesUpper),
    timeToDeath = deathTimes-onsetTimes
)
## Main script
nimIncubationCode = nimbleCode({
    # priors for parameters of incubation period
    incubation_mean ~ dinvgamma(0.001, 0.001)
    incubation_sd ~ dinvgamma(0.001, 0.001)
    death_mean ~ dinvgamma(0.001, 0.001)
    death_sd ~ dinvgamma(0.001, 0.001)
    # priors for known time of exposure
    exposure_mu[1:nCases] <- (exposureTimeLower[1:nCases]+exposureTimeUpper[1:nCases])/2
    exposure_sigma[1:nCases] <- (exposureTimeUpper[1:nCases]-exposureTimeLower[1:nCases])/1.96/2
    for (k in 1:nCases) {
        exposureTime[k] ~ dnorm(exposure_mu[k], sd=exposure_sigma[k])
        incubationTime[k] ~ dgamma(mean=incubation_mean, sd=incubation_sd)
        timeToDeath[k] ~ dgamma(mean=death_mean, sd=death_sd)
        onsetTime[k] ~ dnorm(exposureTime[k] + incubationTime[k], sd=incubation_uncertainty)
        deathTime[k] ~ dnorm(exposureTime[k] + incubationTime[k] + timeToDeath[k], sd=death_uncertainty)
    }
})
## Nimble Model
nimIncubationModel = nimbleModel(nimIncubationCode,
                        data = nimIncubationData,
                        constants = nimIncubationConsts,
                        inits = nimIncubationInits)

nimIncubationConf = configureMCMC(nimIncubationModel,
                                  monitors=c("incubation_mean","incubation_sd"),
                                  monitors2=c("death_mean","death_sd"),
                                  monitors3=nimIncubationSamples %>% select(contains("exposure")) %>% colnames,
                                  thin = 50, thin2 = 50, thin3 = 10)

## Model compilation
nimIncubationMCMC = buildMCMC(nimIncubationConf)
compiledIncubationModel = compileNimble(nimIncubationModel,nimIncubationMCMC)

set.seed(0) #for reproducibility
Niter = 5e5
Nburn = 2e4
compiledIncubationModel$nimIncubationMCMC$run(niter=Niter,nburn=Nburn)
compiledIncubationModel$nimIncubationMCMC$mvSamples %>% as.matrix %>% as.data.frame -> nimIncubationSamples
compiledIncubationModel$nimIncubationMCMC$mvSamples2 %>% as.matrix %>% as.data.frame -> nimIncubationSamples2

## Credible intervals
IncubationCIs = hdi(nimIncubationSamples,.95) %>% as.data.frame %>% { tibble::rownames_to_column(.,"var") }

nimIncubationSamples %>%
    summarise_all(funs(mean)) %>%
    mutate(var="median") %>%
    rbind(IncubationCIs) %>%
    select(var,everything()) -> IncubationCIs

(IncubationCIs)

## Credible intervals
IncubationCIs2 = hdi(nimIncubationSamples2,.95) %>% as.data.frame %>% { tibble::rownames_to_column(.,"var") }

nimIncubationSamples2 %>%
    summarise_all(funs(mean)) %>%
    mutate(var="median") %>%
    rbind(IncubationCIs2) %>%
    select(var,everything()) -> IncubationCIs2

IncubationCIs %<>% left_join(IncubationCIs2)

IncubationCIs %>%
    filter(var=="median") %>%
    select(-var) -> mediansJos

# Incubation period = gamma distribution
message("incubation_mean")
mediansJos$incubation_mean
message("incubation_sd")
mediansJos$incubation_sd

incubation_shape = mediansJos$incubation_mean^2/mediansJos$incubation_sd^2
incubation_rate = mediansJos$incubation_mean/mediansJos$incubation_sd^2*7.0 #from days to weeks
message("incubation shape and scale")
print(c(incubation_shape,1/incubation_rate))

# Time from Onset to Death
message("period_to_death_mean")
mediansJos$death_mean
message("period_to_death_sd")
mediansJos$death_sd
death_shape = mediansJos$death_mean^2/mediansJos$death_sd^2
death_rate = mediansJos$death_mean/mediansJos$death_sd^2*7.0 #from days to weeks
message("period_to_death shape and scale")
print(c(death_shape,1/death_rate))

# 2. Analyzing the main dataset

yearMin = 2016

data = read_excel("../../data/Nigeria_raw.xlsx", sheet = "Incidence") %>%
    select(-one_of("Timeseries","Imputation","File in the repo"),-contains("URL")) %>%
    filter(Year>=yearMin) %>%
    group_by(Year) %>%
    mutate(Incidence_Reported = if_else(Week==1,Reported,Reported-lead(Reported)),
           Incidence_Deaths = if_else(Week==1,Deaths,Deaths-lead(Deaths))
          ) %>%
    ungroup

data %>% select(Year,Week,starts_with("Incidence")) -> Df

Df %>% head

data.frame(Year=yearMin-1,Week=1:52,Incidence_Reported=NA,Incidence_Deaths=NA) %>%
    rbind(Df %>% arrange(Year,Week)) %>%
    rowwise %>%
    mutate(Incidence_Reported_NA=ifelse(is.na(Incidence_Reported),rpois(1,30),NA),
           Incidence_Deaths_NA=ifelse(is.na(Incidence_Deaths),rpois(1,1),NA)) %>%
    ungroup -> Df

K = nrow(Df)

# Convolutions
incubation_probability = pgamma(1:K,shape=incubation_shape,rate=incubation_rate)-pgamma(1:K-1,shape=incubation_shape,rate=incubation_rate)
timeFromOnsetToDeath_probability = pgamma(1:K,shape=death_shape,rate=death_rate)-pgamma(1:K-1,shape=death_shape,rate=death_rate)
# time from Exposure event to Death is the convolution of two latter probabilities
timeFromExposureToDeath_probability = c(0)
for (x in 2:K) {
    timeFromExposureToDeath_probability = c(timeFromExposureToDeath_probability,
        sum(incubation_probability[1:(x-1)]*timeFromOnsetToDeath_probability[(x-1):1]))
}

shift = 26

nimData = list(infected = Df$Incidence_Reported,
                dead = Df$Incidence_Deaths,
                one = 1)

nimConsts = list(r = rho,
                  K = K,
                  week = (Df$Week-shift-1)%%52+1,
                  incubationPeriod = incubation_probability[1:K],
                  periodToDeath = timeFromExposureToDeath_probability[1:K])

nimInits = function(){
    list(l1 = (runif(1,44,52)-shift-1)%%52+1,
        l2 = (runif(1,5,12)-shift-1)%%52+1,
        infected = Df$Incidence_Reported_NA,
        dead = Df$Incidence_Deaths_NA,
        lambdaIncidence = rep(0,K),
        lambdaDeath = rep(0,K),
        exposure = rgamma(K,10,1),
        CFR = rgamma(K,0.08,1),
        mean_a = c(runif(1,2,12),runif(1,13,28)),
        sd_a = rgamma(2,1,1),
        mean_q = rbeta(1,1,1),
        sd_q = rbeta(1,1,1))}

nimConvolution = nimbleFunction(
     run = function(a = double(1), b = double(1)) {
         L <- dim(a)[1]
         ans <- inprod(a[L+1-1:L],b)
         return(ans)
         returnType(double(0))
     }
 )

nimCode = nimbleCode({
     ## exposure and awareness
     for (k in 1:K) {
         mean_a_realized[k] <- mean_a[1]+(mean_a[2]-mean_a[1])*equals(step(week[k]-l1),step(l2-week[k]))
         sd_a_realized[k] <- sd_a[1]+(sd_a[2]-sd_a[1])*equals(step(week[k]-l1),step(l2-week[k]))
         exposure[k] ~ dgamma(mean=mean_a_realized[k],sd=sd_a_realized[k])
         CFR[k] ~ dgamma(mean=mean_q,sd=sd_q)
     }
     ### I start from week 53 because the first 52 weeks were added with missing data
     ### in order to avoid edge effects for the convolution
     for (k in 52:(K-1)) {
         lambdaIncidence[k+1] <- nimConvolution(exposure[1:k],incubationPeriod[1:k])/(1-r)
         lambdaDeath[k+1] <- CFR[k+1]*nimConvolution(exposure[1:k],periodToDeath[1:k])/(1-r)
         infected[k+1] ~ dpois(lambdaIncidence[k+1])
         dead[k+1] ~ dpois(lambdaDeath[k+1])
     }
     ## Priors
     for(i in 1:2) {
         mean_a[i] ~ dinvgamma(0.001, 0.001)
         sd_a[i] ~ dinvgamma(0.001, 0.001)
     }
     mean_q ~ dinvgamma(0.001, 0.001)
     sd_q ~ dinvgamma(0.001, 0.001)
     l1 ~ dunif(0,53)
     l2 ~ dunif(0,53)
     # Constraints
     one ~ dbern(condition)
     condition <- step(l2-l1-2)*step(mean_a[2]-mean_a[1])
})

set.seed(as.numeric(args[1])) #so that the seed will be controllable, but would give different initial conditions
nimModel = nimbleModel(nimCode,
                       constants = nimConsts,
                       data = nimData,
                       inits = nimInits())

nimConf = configureMCMC(nimModel,setSeed=TRUE,control=list(reflective=TRUE))
# thinning parameter
nimConf$setThin(100)
## Checking that all variables are properly initialized
nimModel$initializeInfo()

## Model compilation
nimMCMC = buildMCMC(nimConf)
compiledModel = compileNimble(nimModel,nimMCMC)

Niter = 2e5
Nburn = 2e4
set.seed(as.numeric(args[1]))
compiledModel$nimMCMC$run(niter=Niter+Nburn, nburnin = Nburn)
compiledModel$nimMCMC$mvSamples %>% as.matrix %>% as.data.frame -> nimSamples

saveRDS(nimSamples, file = paste0("nimSamples-final-",args[1],".rds"))
