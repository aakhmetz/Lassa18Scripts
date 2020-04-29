# Lassa18Scripts

Supplementary materials for Akhmetzhanov AR, Asai Y, Nishiura H: Quantifying the seasonal drivers of transmission for Lassa fever in Nigeria. Philosophical Transactions of the Royal Society B 374, 20180268, 2019

[Manuscript](http://dx.doi.org/10.1098/rstb.2018.0268) (Open access)

**Abstract:** Lassa fever (LF) is a zoonotic disease that is widespread in West Africa and involves animal-to-human and human-to-human transmission. Animal-to-human transmission occurs upon exposure to rodent excreta and secretions, i.e. urine and saliva, and human-to-human transmission occurs via the bodily fluids of an infected person. To elucidate the seasonal drivers of LF epidemics, we employed a mathematical model to analyse the datasets of human infection, rodent population dynamics and climatological variations and capture the underlying transmission dynamics. The surveillance-based incidence data of human cases in Nigeria were explored, and moreover, a mathematical model was used for describing the transmission dynamics of LF in rodent populations. While quantifying the case fatality risk and the rate of exposure of humans to animals, we explicitly estimated the corresponding contact rate of humans with infected rodents, accounting for the seasonal population dynamics of rodents. Our findings reveal that seasonal migratory dynamics of rodents play a key role in regulating the cyclical pattern of LF epidemics. The estimated timing of high exposure of humans to animals coincides with the time shortly after the start of the dry season and can be associated with the breeding season of rodents in Nigeria.

## Jupyter notebooks

### Preprocessing of the data
* [A1. Data aggregation - PlotDigitizer of the data from 2012 to the mid of 2014.ipynb](https://nbviewer.jupyter.org/github/aakhmetz/Lassa2018Scripts/blob/master/scripts/A1.%20Data%20aggregation%20-%20PlotDigitizer%20of%20the%20data%20from%202012%20to%20the%20mid%20of%202014.ipynb) (written in *R*)</br>Very preliminary check of the grabbed data for 2012–mid of 2014. To note that the data sources were quite scarce and sometime contradictive for that period. 
* [A2. Data aggregation - Final steps.ipynb](https://nbviewer.jupyter.org/github/aakhmetz/Lassa2018Scripts/blob/master/scripts/A2.%20Data%20aggregation%20-%20Final%20steps.ipynb) (written in *R*)
* [A3. Nigeria map.ipynb](https://nbviewer.jupyter.org/github/aakhmetz/Lassa2018Scripts/blob/master/scripts/A3.%20Nigeria%20map.ipynb) (written in *R*)</br>Here, we created the maps of Nigeria with incidence of LF for 2012–2018 which were showed in Supplementary Materials 

### Main analysis
* [B. Main analysis.ipynb](https://nbviewer.jupyter.org/github/aakhmetz/Lassa2018Scripts/blob/master/scripts/B.%20Main%20analysis.ipynb) (written in *R* using the package [*nimble*](https://r-nimble.org/) for Bayesian inference)</br>
Obtaining main results reported in the paper. For example, we used the data from the nosocomial outbreak in [Jos in 1970](http://dx.doi.org/10.1016/0035-9203(72)90271-4), to estimate the **incubation period** distribution of LF as the gamma distribution with a mean of 12.8 days (95% credible interval (CI): 10.7, 15.0) and a standard deviation of 4.6 days (95% CI: 2.8, 6.6)
* [C. Transmission dynamics in rodents.ipynb](https://nbviewer.jupyter.org/github/aakhmetz/Lassa2018Scripts/blob/master/scripts/C.%20Transmission%20dynamics%20in%20rodents.ipynb) (written in *Python*)</br>Modeling transmission of LF within rodent populations by solving the SIR model with seasonality and density-dependent contact rate
* [D. Precipitation.ipynb](https://nbviewer.jupyter.org/github/aakhmetz/Lassa2018Scripts/blob/master/scripts/D.%20Precipitation.ipynb) (written in *R*)</br>Showing the precipitation data for Nigeria
* [E. CCM.ipynb](https://nbviewer.jupyter.org/github/aakhmetz/Lassa2018Scripts/blob/master/scripts/E.%20CCM.ipynb) (written in *R* using the package [rEDM](https://ha0ye.github.io/rEDM/articles/rEDM.html))</br>Performing coss-map causality analysis for shared seasonality according to the method developed by George Sugihara and colleagues

---------
**Thank you for your interest to our work!** 

Few words of caution: We would like to note that our code is not supposed to work out of box, because the links used in the notebooks were user-specific, and our main intent was to show the relevance of the methods used in our paper.
