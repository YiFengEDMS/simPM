## Welcome to the homepage of _simPM_!

`simPM` is an R package that automates the search for optimal 'post hoc' planned missing (PHPM) designs. This R package is developed and maintained by [Yi Feng](https://terpconnect.umd.edu/~yifeng94/) & [Dr. Gregory R. Hancock](https://education.umd.edu/directory/gregory-r-hancock) from the University of Maryland.

#### How to install the _simPM_ R package? 

The source code of this R package is made public on [Github](https://github.com/YiFengEDMS/simPM). To install the R package on your local machine, please use the following code:

```markdown

install.packages("devtools")
library(devtools)
devtools::intall_github("YiFengEDMS/simPM")

```

#### How to use _simPM_? 

- The main function in package `simPM` is `simPM()`. As an example, the following code can be used to search for an optimal wave-level PHPM design with 4 waves of repeated measures.

```
wave.ex1=simPM(popModel=popModel,                          #supply the population model using lavaan language
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



- To view the results, use the `summary.opt` function:
```
summary.opt(wave.ex1)
```



- To view the missing data patterns in the optimal PHPM design, use the `plotPM` function:
```
plotPM(wave.ex1,Time=4,k=1)
```



- To view more details of the optimal PHPM design, use the following code:

```
summary(wave.ex1$opt.output)
```

More details are available in the [package manual](). More examples will be available on this page soon.

#### Questions or Suggestions?
Send an email to yifeng94@umd.edu. We are happy to hear about your thoughts!


```
- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```
