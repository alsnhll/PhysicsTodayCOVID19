# PhysicsTodayCOVID19

This article provides the code required to reproduce the figures and analysis included in the forthcoming article in Physics Today "COVID-19 and the math behind epidemics", by Alison L Hill. All code is written in Matlab. 

# Description of folders/files

__Figure1__: These files contain the data to draw the transmisison network which was described in the study "Jang S, Han SH, Rhee J-Y. _Cluster of Coronavirus Disease Associated with Fitness Dance Classes, South Korea_. Emerging Infectious Diseases Journal. 2020;26. doi:10.3201/eid2608.200633". 
I encoded the network based on the data in Figure 1. 
* korean_dance_network.xlsx - Excel file with network edges
* zumba_net.mat - Matlab data file with network edges
* zumba_names.mat - Matlab data file with node types
* plot_zumba_graph.m - Matlab script to create and plot the network
* korean_network.gephi - network file in format to be edited by network visualization software Gephi

__Box1__: These files are related to the calculations for the basic reproductive ratio R0
* dispR0.m - Matlab script to calculate the distribution of secondary infections for different values of the basic reproductive ratio R0 and the dispersion parameter, and the cumulative percent of infected individuals responsible for each percentile of total secondary infections. 
* rVsR0.m - Matlab script to calculate the R0 value consistent with observed values of the exponential growth rate and incubation and infectious periods. The large dots assume both the latent period and infectious period are exponentially distributed. The upper estimate of R0 assumes the distribution of the incubation period is narrower (gamma-distributed with shape parameter 5), then estimates of R0 will be higher, whereas if the distribution of the infectious period is narrower (gamma-distributed with shape parameter 5), then R0  is inferred to be lower (triangles). 

Further reading on this topic, which inspired these calculations:

Re dispersion in individual-R0 and superspreading:
* Lloyd-Smith JO, Schreiber SJ, Kopp PE, Getz WM. Superspreading and the effect of individual variation on disease emergence. Nature. 2005;438: 355–359. doi:10.1038/nature04153
* Bi Q, Wu Y, Mei S, Ye C, Zou X, Zhang Z, et al. Epidemiology and transmission of COVID-19 in 391 cases and 1286 of their close contacts in Shenzhen, China: a retrospective cohort study. The Lancet Infectious Diseases. 2020;0. doi:10.1016/S1473-3099(20)30287-5
*  On Kwok K, Hin Chan HH, Huang Y, Cheong Hui DS, Anantharajah Tambyah P, In Wei W, et al. Inferring super-spreading from transmission clusters of COVID-19 in Hong Kong, Japan and Singapore. J Hosp Infect. 2020. doi:10.1016/j.jhin.2020.05.027
* Endo A, Centre for the Mathematical Modelling of Infectious Diseases COVID-19 Working Group, Abbott S, Kucharski AJ, Funk S. Estimating the overdispersion in COVID-19 transmission using outbreak sizes outside China. Wellcome Open Res. 2020;5: 67. doi:10.12688/wellcomeopenres.15842.1
* Riou J, Althaus CL. Pattern of early human-to-human transmission of Wuhan 2019 novel coronavirus (2019-nCoV), December 2019 to January 2020. Eurosurveillance. 2020;25: 2000058. doi:10.2807/1560-7917.ES.2020.25.4.2000058
* Adam DC, Wu P, Wong JY, Lau EHY, Tsang TK, Cauchemez S, et al. Clustering and superspreading potential of SARS-CoV-2 infections in Hong Kong. Nature Medicine. 2020; 1–6. doi:10.1038/s41591-020-1092-0
* Althouse BM, Wenger EA, Miller JC, Scarpino SV, Allard A, Hébert-Dufresne L, et al. Stochasticity and heterogeneity in the transmission dynamics of SARS-CoV-2. arXiv:200513689 [physics, q-bio]. 2020. Available: http://arxiv.org/abs/2005.13689
* Zhang Y, Li Y, Wang L, Li M, Zhou X. Evaluating Transmission Heterogeneity and Super-Spreading Event of COVID-19 in a Metropolis of China. Int J Environ Res Public Health. 2020;17. doi:10.3390/ijerph17103705
* Lau MSY, Grenfell B, Thomas M, Bryan M, Nelson K, Lopman B. Characterizing superspreading events and age-specific infectiousness of SARS-CoV-2 transmission in Georgia, USA. PNAS. 2020;117: 22430–22435. doi:10.1073/pnas.2011802117

Re estimating R0 from the epidemic growth rate and transmission intervals:
* Wearing HJ, Rohani P, Keeling MJ. Appropriate Models for the Management of Infectious Diseases. PLOS Medicine. 2005;2: e174. doi:10.1371/journal.pmed.0020174
* Wallinga J, Lipsitch M. How generation intervals shape the relationship between growth rates and reproductive numbers. Proceedings of the Royal Society B: Biological Sciences. 2007;274: 599–604. doi:10.1098/rspb.2006.3754
* Park SW, Champredon D, Weitz JS, Dushoff J. A practical generation-interval-based approach to inferring the strength of epidemics from their speed. Epidemics. 2019;27: 12–18. doi:10.1016/j.epidem.2018.12.002
* Park SW, Bolker BM, Champredon D, Earn DJD, Li M, Weitz JS, et al. Reconciling early-outbreak estimates of the basic reproductive number and its uncertainty: framework and applications to the novel coronavirus (SARS-CoV-2) outbreak. Journal of The Royal Society Interface. 2020;17: 20200144. doi:10.1098/rsif.2020.0144
* Sanche S, Lin YT, Xu C, Romero-Severson E, Hengartner N, Ke R. High Contagiousness and Rapid Spread of Severe Acute Respiratory Syndrome Coronavirus 2. Emerging Infectious Diseases journal. 2020;26. doi:10.3201/eid2607.200282
* Pan A, Liu L, Wang C, Guo H, Hao X, Wang Q, et al. Association of Public Health Interventions With the Epidemiology of the COVID-19 Outbreak in Wuhan, China. JAMA. 2020. doi:10.1001/jama.2020.6130
* Ali ST, Wang L, Lau EHY, Xu X-K, Du Z, Wu Y, et al. Serial interval of SARS-CoV-2 was shortened over time by nonpharmaceutical interventions. Science. 2020. doi:10.1126/science.abc9004


__Box2__: These plots were generated using the same SEIRHD differential equation model as Figure 3. The total population size was N=1e5, the average duration of the incubation period was 5 days, the average duration of the infectious period was 5 days, the average duration of delay to death after infectious period  was 14 days, and 1% of all infected individuals died. For the "low r" case, we set b = 0.4/N (corresponding to r = 0.08/day or doubling time ~8.3 days). For the "high r" case, we set b = 0.8/N (corresponding to r = 0.2/day or doubling time ~3.5 days). For the "low -> high r" case we set b = 0.4/N  and eff= -1 and Tint=50 (corresponding to an intervention that led to a 2-fold increase in transmission starting 50 days after the epidemic began). Other parameters stayed at the baseline values in the code. 
* SEIRHD_COVID19.m - Matlab script to define the parameters and simulate the infection course
* SEIRHD_COVID19_eqns_v2.m - Matlab function with differential equations
* Note you need to add the toolbox folder brewer_map to your path to get the ColorBrewer color schemes

Further reading on this topic, which inspired these calculations:
* Lipsitch M, Donnelly CA, Fraser C, Blake IM, Cori A, Dorigatti I, et al. Potential Biases in Estimating Absolute and Relative Case-Fatality Risks during Outbreaks. PLoS Negl Trop Dis. 2015;9. doi:10.1371/journal.pntd.0003846
* Hauser A, Counotte MJ, Margossian CC, Konstantinoudis G, Low N, Althaus CL, et al. Estimation of SARS-CoV-2 mortality during the early stages of an epidemic: A modeling study in Hubei, China, and six regions in Europe. PLOS Medicine. 2020;17: e1003189. doi:10.1371/journal.pmed.1003189
* Verity R, Okell LC, Dorigatti I, Winskill P, Whittaker C, Imai N, et al. Estimates of the severity of coronavirus disease 2019: a model-based analysis. The Lancet Infectious Diseases. 2020;0. doi:10.1016/S1473-3099(20)30243-7
* Wu JT, Leung K, Bushman M, Kishore N, Niehus R, de Salazar PM, et al. Estimating clinical severity of COVID-19 from the transmission dynamics in Wuhan, China. Nature Medicine. 2020; 1–5. doi:10.1038/s41591-020-0822-7

__Figure3__: This “compartmental” model classifies individuals as susceptible (S), exposed (E), infectious (I), hospitalized (H), recovered (R), or dead (D). The model is simulated as a system of ordinary differential equations. The average duration of the latent period (“exposed” stage) and infectious period are 5 days. 10% of individuals progress to hospitalization, of average length 2 weeks. 1% of all individuals (10% of those hospitalized) eventually die. At baseline, R0 = 2.5. The social distancing intervention starts at day 60 and leads to a 70% reduction in transmission. For case isolation intervention, 90% of infected individuals leave the infectious stage after an average of 1 day and enter a "isolated" state where they can't transmit to others. 
* SEIRHD_COVID19_isolation.m - Matlab script to define the parameters and simulate the infection course
* SEIR_COVID19_eqns_v3.m - Matlab function with differential equations
* SEIR_COVID19_eqns_v4.m - Matlab function with differential equations including isolation
* Note you need to add the toolbox folder brewer_map to your path to get the ColorBrewer color schemes

Examples of models for COVID-19 that inspired this one:
* https://alhill.shinyapps.io/COVID19seir/
* https://covid19-scenarios.org/
* Kissler SM, Tedijanto C, Goldstein E, Grad YH, Lipsitch M. Projecting the transmission dynamics of SARS-CoV-2 through the postpandemic period. Science. 2020. doi:10.1126/science.abb5793
* Ferretti L, Wymant C, Kendall M, Zhao L, Nurtay A, Abeler-Dörner L, et al. Quantifying SARS-CoV-2 transmission suggests epidemic control with digital contact tracing. Science. 2020. doi:10.1126/science.abb6936
* Peak CM, Kahn R, Grad YH, Childs LM, Li R, Lipsitch M, et al. Individual quarantine versus active monitoring of contacts for the mitigation of COVID-19: a modelling study. The Lancet Infectious Diseases. 2020;0. doi:10.1016/S1473-3099(20)30361-3
* Ferguson NM, Laydon D, Nedjati-Gilani G, Imai N, Ainslie K, Baguelin M, et al. Impact of non-pharmaceutical interventions (NPIs) to reduce COVID- 19 mortality and healthcare demand. 2020; 20. 
* Davies NG, Kucharski AJ, Eggo RM, Gimma A, Edmunds WJ, Jombart T, et al. Effects of non-pharmaceutical interventions on COVID-19 cases, deaths, and demand for hospital services in the UK: a modelling study. The Lancet Public Health. 2020;5: e375–e385. doi:10.1016/S2468-2667(20)30133-X
* Salje H, Kiem CT, Lefrancq N, Courtejoie N, Bosetti P, Paireau J, et al. Estimating the burden of SARS-CoV-2 in France. Science. 2020;369: 208–211. doi:10.1126/science.abc3517

__Box3__: This code generates networks with different statistical properties, and then uses a discrete-time stochastic simulation of a simple SIR epidemic model. 
* make_networks.m - Matlab script that creates and plots the three different networks: i) a uniform random network where every individual has exactly 10 connectoins, ii) a small world network created with the Watt-Strogatz algorithm with mean degree 10 and rewiring probability 10%, and iii) a random netwwork with a gamma-distributed degree distribution with mean 10 and standard deviation 10. 
* networkMake.m - Matlab function that can construct many different types of networks
* stubconnect.m - Matlab function that creates random networks from an input degree distribution
* run_networkSIR.m - Matlab script that to define the parameters for a stochastic SIR model on a defined transmission network. The population size was N=1000, the transmission rate per contact per time was 0.05, and the average duration of infection was 5 days. 
* networkSIR.m - Matlab function that simulates multiple iterations of a stochastic SIR model on a defined transmission network. Infection is stimulated starting from a single infected individual, and runs until the infection is extinct. 

This simulations were inspired by:
* May RM, Lloyd AL. Infection dynamics on scale-free networks. Phys Rev E. 2001;64: 066112. doi:10.1103/PhysRevE.64.066112
* Pastor-Satorras R, Vespignani A. Epidemic spreading in scale-free networks. Physical Review Letters. 2001;86: 3200–3203. doi:10.1103/PhysRevLett.86.3200
* Watts DJ, Strogatz SH. Collective dynamics of ‘small-world’networks. Nature. 1998;393: 440–442. 
* Santos FC, Rodrigues JF, Pacheco JM. Epidemic spreading and cooperation dynamics on homogeneous small-world networks. Phys Rev E. 2005;72: 056128. doi:10.1103/PhysRevE.72.056128
* Miller JC. Spread of infectious disease through clustered populations. Journal of The Royal Society Interface. 2009;6: 1121–1134. doi:10.1098/rsif.2008.0524
* Volz EM, Miller JC, Galvani A, Ancel Meyers L. Effects of Heterogeneous and Clustered Contact Patterns on Infectious Disease Dynamics. PLoS Comput Biol. 2011;7: e1002042. doi:10.1371/journal.pcbi.1002042
* Leventhal GE, Hill AL, Nowak MA, Bonhoeffer S. Evolution and emergence of infectious diseases in theoretical and real-world networks. Nature Communications. 2015;6. doi:10.1038/ncomms7101
* Hébert-Dufresne L, Althouse BM, Scarpino SV, Allard A. Beyond $R_0$: the importance of contact tracing when predicting epidemics. arXiv:200204004 [physics, q-bio]. 2020. Available: http://arxiv.org/abs/2002.04004

