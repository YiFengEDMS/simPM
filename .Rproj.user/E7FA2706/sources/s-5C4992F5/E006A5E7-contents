## Welcome to the homepage of _simPM_!

`simPM` is an R package that automates the search for optimal 'post hoc' planned missing (PHPM) designs. This R package is developed and maintained by [Yi Feng](https://terpconnect.umd.edu/~yifeng94/) & [Dr. Gregory R. Hancock](https://education.umd.edu/directory/gregory-r-hancock) from the University of Maryland.

#### The story behind _simPM_
*"What should I do when they cut my research funding after my longitudinal study is already underway?"*

`simPM` was created to help researchers survive the unexpected funding cut in the course of a longitudinal study. It can be used to find an optimal 'post hoc' planned missing design that allows the researchers to complete the study at a reduced cost, while maintaining satisfactory level of statistical power for testing the focal parameters. 

#### What does _simPM_ do?
By automizing the simulation-based power analysis for planned missing designs in longitudinal context, `simPM` can free the researchers from manually configuring the possible PM designs, determining their eligibility, setting up the simulations, and summarizing the results over replications, which can be tedious and time-consuming work especially when there is a large number of plausible PHPM designs to be evaluated.


#### How to install _simPM_? 

The source code of this R package is made public on the author's [Github page](https://github.com/YiFengEDMS/simPM). To install the R package on your local machine, please run the following code:

```r

install.packages("devtools")
library(devtools)
devtools::intall_github("YiFengEDMS/simPM")

``` 

#### How to use _simPM_? 

The main function in package `simPM` is `simPM()`. As an example, the following code can be used to search for an optimal wave-level PHPM design with 4 waves of repeated measures.       

```r
wave.ex1=simPM(
        popModel=popModel,                                 #supply the population model using lavaan language
        analyzeModel=analyzeModel,                         #supply the analysis model using lavaan language
        VarNAMES=c("se1","se2","se3","se4"),               #specify the observed variable names, in chronological order
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
``` 


To view the results, use the `summary.opt` function.
```r
summary.opt(wave.ex1)
``` 

To view the missing data patterns in the optimal PHPM design, use the `plotPM` function.
```r
plotPM(wave.ex1,Time=4,k=1)
```


To view more details of the optimal PHPM design, use the following code:
```r
summary(wave.ex1$opt.output)
```

More details are available in the [package manual](). More examples will be available on this page soon. Please stay tuned!

#### Questions or Suggestions?
Send an email to [yifeng94@umd.edu](yifeng94@umd.edu). We are happy to hear about your thoughts!

