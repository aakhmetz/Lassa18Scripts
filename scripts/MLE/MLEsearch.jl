# load the functionality into Julia
using FileIO, DataFrames, CSV, Distributions, Optim, ShiftedArrays, StatsFuns

# Loading the data
data = CSV.read("data.csv") |> DataFrame
df = DataFrame(
    Year=[isa.(x,Number)==true ? Int64.(x) : missing for x in data[:Year]],
    Week=[isa.(x,Number)==true ? Int64.(x) : missing for x in data[:Week]],
    Reported=[isa.(x,Number)==true ? Int64.(x) : missing for x in data[:Reported]],
    Deaths=[isa.(x,Number)==true ? Int64.(x) : missing for x in data[:Deaths]]
)

minYear = 2016

# remove entities with missing values for the week
df = df[.~ismissing.(df[:Week]),:]
# consider only more or less complete data
df = df[df[:Year].>=minYear,:]

# adding one more year to the dataframe with completely missing data
# in order to regularize convolutions (so they will have periodic structure for known years)
df = [DataFrame(Year=minimum(df[:Year])-1,Week=1:52,Reported=missing,Deaths=missing);df]

# calculating incidence
df[:IncidenceReported] = ifelse.(df[:Week].!=1,df[:Reported].-lead(df[:Reported]),df[:Reported])
df[:IncidenceDeaths] = ifelse.(df[:Week].!=1,df[:Deaths].-lead(df[:Deaths]),df[:Deaths])

sort!(df, [order(:Year), order(:Week)])

### Incubation period
IncubationPeriod_mean = 12.8271965544676/7.0
IncubationPeriod_var = 21.4535562783816/7.0^2
IncubationPeriod = Gamma(IncubationPeriod_mean^2/IncubationPeriod_var,IncubationPeriod_var/IncubationPeriod_mean)

dIncubationPeriod = (x -> cdf(IncubationPeriod,x)-cdf.(IncubationPeriod,x-1)).(1:length(df[:1]))

###  From onset to death
TimeFromOnsetToDeath_mean = 13.804525254164/7.0
TimeFromOnsetToDeath_var = 60.6313724811289/7.0^2
TimeFromOnsetToDeath = Gamma(TimeFromOnsetToDeath_mean^2/TimeFromOnsetToDeath_var,TimeFromOnsetToDeath_var/TimeFromOnsetToDeath_mean)

dTimeFromOnsetToDeath = (x -> cdf(TimeFromOnsetToDeath,x)-cdf(TimeFromOnsetToDeath,x-1)).(1:length(df[:1]))

### From exposure to Deaths
dDeathPeriod = (x -> x>1 ? vecdot(dIncubationPeriod[1:(x-1)],dTimeFromOnsetToDeath[(x-1):-1:1]) : 0.0).(1:length(df[:1]))

### First try of optimization
stepFunction(x) = x >= 0 ? 1 : 0
heaviside(x) = x > 0 ? 1 : (x == 0 ? 0 : -1)
equals(x,y) = x==y ? 1 : 0

ρ = 0.19

f = open("MLEsearch.csv","w")
write(f,"i1,i2,exp1,exp2,CFR,loglk\n")
close(f)

function getNegativeLoglk(prms,indices)
    exposureLevels = exp.(prms[1:2])
    riskOfDeath = logistic(prms[3])

    indexExposureLower = indices[1]
    indexExposureUpper = indices[2]
    nrecords = nrow(df)

    # here I have a bit complicated formula, because I may assume that indexExposureUpper can be less than indexExposureLower
    exposure = exposureLevels[1]+(exposureLevels[2]-exposureLevels[1])*(equals.(stepFunction.(df[:Week]-indexExposureLower),stepFunction.(indexExposureUpper-df[:Week])))

    λReported = 1/(1-ρ).*((x -> x>1 ? vecdot(exposure[1:(x-1)],dIncubationPeriod[(x-1):-1:1]) : 0.0).(1:nrecords))
    λDeaths = riskOfDeath/(1-ρ).*((x -> x>1 ? vecdot(exposure[1:(x-1)],dDeathPeriod[(x-1):-1:1]) : 0.0).(1:nrecords))

    return(-sum(df[j,:IncidenceReported]*log(λReported[j])-λReported[j] for j in 1:nrecords if (isa.(df[j,:IncidenceReported],Number))) -
            sum(df[j,:IncidenceDeaths]*log(λDeaths[j])-λDeaths[j] for j in 1:nrecords if (isa.(df[j,:IncidenceDeaths],Number))))
end

for i1 in 1:52
    println("\t ($i1)")
    result = [0,0,1,1,1,1e10]
    for i2 in 53:-1:i1
        pars = [30,10,.1]
        pars = [log.(pars[1:2]);logit(pars[3])]
        sol = optimize((x -> getNegativeLoglk(x,[i1,i2])),pars)
        pars = sol.minimizer
        loglk = Optim.minimum(sol)
        if (loglk<result[end])
            result = [i1,i2,exp(pars[1]),exp(pars[2]),logistic(pars[3]),loglk]
            print(result)
        end
    end
    f = open("MLEsearch.csv","a")
    write(f,join(result,", "),"\n")
    close(f)
end
