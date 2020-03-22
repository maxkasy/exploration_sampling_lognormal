##############################################
#Invoke this script from the terminal as follows:
#RScript --vanilla command_line_app.R data/test_data.csv 600
#the first argument is the input file of prior data
#the second argument is the number of units in the current wave
##############################################
library(tidyverse)
source("modified_thompson.R")

# Reading in data and arguments
command_arguments = commandArgs(trailingOnly = T)
prior_datapath = command_arguments[1]
# prior_datapath = "test/test_data.csv"
Nt = as.integer(command_arguments[2])

priordata = priordata=read_csv(prior_datapath)
k=max(priordata$D) # number of treatments

###########################
# Calculating treatment assignment
alpha = DtchoiceThompsonProbabilities(priordata$Y, priordata$D, k)
shares = DtchoiceThompson_modified(alpha)
Dt = ProportionalAssignment(shares, Nt)

filename = paste("data/", Sys.Date(), "_treatment_assignment.csv", sep = "")
write_csv(tibble(treatment = Dt), path = filename)



