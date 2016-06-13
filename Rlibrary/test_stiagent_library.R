##################################################################
######
######    MINIMAL TEST FOR 'stiagent' LIBRARY
######
######
##################################################################


library(stiagent,lib.loc = "./lib")

t0 <- Sys.time()
# path to model input files:
folder_inputs = "./test-inputs/"
folder_calib = "../calibration/"

# founder population
founder.size <- 100
founder.fem.prop <- 0.5
founder.csw.prop <- 0.02

# scenario 
# (includes intervention filenames 
#  that will be done during this simualtion)
scen.file <- "in_scenario_baseline.csv"

# simply run the model:
x <- stiagent_runsim(params = list(folder_inputs = folder_inputs,
                                   folder_calib = folder_calib,
								   founder_size = founder.size,
								   founder_femprop = founder.fem.prop,
								   founder_cswprop = founder.csw.prop,
                                   scenario_file = scen.file,
                                   MC_id = 1,  # <-- MC id (*not* number of MC iterations!!!)
								   displayProgress=1)
                     )



t1 <- Sys.time()
message(paste("time elapsed:",round(t1-t0,1),"sec"))

# check if everything went well:
msg <- ifelse(length(x)>0,
              "==> stiagent R library seems to work.",
              "THERE IS A PROBLEM WITH stiagent LIBRARY")
message(paste0(rep("=",40)))
message(msg)
