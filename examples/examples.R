popModel <- '

i=~1*se1+1*se2+1*se3+1*se4
s=~0*se1+1*se2+2*se3+3*se4

se1~0*1
se2~0*1
se3~0*1
se4~0*1

se1~~0.071*se1
se2~~0.034*se2
se3~~0.067*se3
se4~~0.025*se4

i~2.983*1
s~0.086*1

i~~0.268*i+(-0.039)*s
s~~0.023*s

'


analyzeModel <- '

i=~1*se1+1*se2+1*se3+1*se4
s=~0*se1+1*se2+2*se3+3*se4

se1~0*1
se2~0*1
se3~0*1
se4~0*1

se1~~se1
se2~~se2
se3~~se3
se4~~se4

i~1
s~1

i~~i+s
s~~s

'



wave.ex1=simPM(
  popModel=popModel,                                 #supply the population model using lavaan language
  analyzeModel=analyzeModel,                         #supply the analysis model using lavaan language
  VarNAMES=c("se1","se2","se3","se4"),               #the observed variable names, in chronological order
  Time=4,                                            #total number of waves
  Time.complete=1,                                   #number of waves completed before funding cut occurs
  k=1,                                               #number of observed variables collected at each wave
  pc=0.2,                                            #percentage of participants to provide complete data after funding cut
  pd=0,                                              #percentage of participants to provide no data after funding cut
  costmx=c(5,10,15),                                 #unit cost of each data point at the following waves
  n=323,                                             #original sample size
  nreps=1000,                                        #number of replications for simulation
  focal.param=c("i~1","s~1","i~~i","s~~s"),          #specify the focal parameters
  complete.wave=NULL,                                #specify any future wave/variables that need complete data 
  eval.budget=T,                                     #whether or not there is a budget restriction
  rm.budget=30*323*0.7,                              #the amount of remaining budget
  distal.var=NULL,                                   #specify any distal variables that are not subject to PM
  seed=12345,                                        #random seed
  engine="l",                                        #use lavaan to fit the models
  methods="wave")                                    #type of PHPM designs, "wave" indicates wave-level missing

setwd("C:/Users/yifeng94/Desktop/simPM/simPM-git/examples")
save(wave.ex1,file="linear LGM example.rda")
