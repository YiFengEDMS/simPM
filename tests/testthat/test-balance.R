context("balanced item-level PM design")

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

output <- balance.miss.l(
  popModel = popModel,
  analyzeModel = analyzeModel,
  NAMES = c("se1","se2","se3","se4"),
  Time = 4,
  Time.complete = 1,
  k = 1,
  pc = 0.2,
  pd = 0,
  costmx = c(5, 10, 15),
  n = 323,
  nreps = 3,
  focal.param = c("i~1","s~1","i~~i","s~~s"),
  complete.var = NULL,
  eval.budget = TRUE,
  rm.budget = 30*323*0.7,
  distal.var = NULL,
  seed=1234
)

test_that("balanced item-level PM design returns a list", {
  expect_is(output, class = "list")
}
)