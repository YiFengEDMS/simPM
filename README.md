## Welcome to the homepage of _simPM_!

'simPM' is an R package that automates the search for optimal 'post hoc' planned missing designs. This R package is developed and maintained by [Yi Feng](https://terpconnect.umd.edu/~yifeng94/) & [Dr. Gregory R. Hancock](https://education.umd.edu/directory/gregory-r-hancock) from the University of Maryland.

#### How to install the _simPM_ R package 

The source code of this R package is made public on [Github](https://github.com/YiFengEDMS/simPM). To install the R package on your local machine, please use the following code:

```install simPM
install.packages("devtools")
library(devtools)
devtools::intall_github("YiFengEDMS/simPM")
```

#### How to use _simPM_? 

The main function in package `simPM` is `simPM()`.

As an example, 

```
wave.ex1=simPM(popModel=popModel, 
        analyzeModel=analyzeModel,
        VarNAMES=c("se1","se2","se3","se4"),
        Time=4,
        Time.complete=1,
        k=1,
        pc=0.2,
        pd=0,
        costmx=c(5,10,15),
        n=323,
        nreps=1000,
        focal.param=c("i~1","s~1","i~~i","s~~s"),
        complete.wave=NULL,
        eval.budget=T,
        rm.budget=30*323*0.7,
        distal.var=NULL,
        seed=12345,
        engine="l",
        methods="wave")

summary.opt(wave.ex1)
summary(wave.ex1$opt.output)

plotPM(wave.ex1,Time=4,k=1)

```

More details are available in the [package manual]().

```


#### Questions or Suggestions?
Send an email to yifeng94@umd.edu. We are happy to hear about your thoughts!


# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).
