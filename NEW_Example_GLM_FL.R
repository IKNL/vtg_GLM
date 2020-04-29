remove(list=ls())

#place your own directory
source("C:/Users/mce1908.52713/Desktop/FEDERATED_GLM/GLM_FL.R")
#set your own directory to sve files
setwd("C:/Users/mce1908.52713/Desktop/FEDERATED_GLM/Example_Files/")


### create the master file. This has to be done once for all
### can be done by the server and distributed to all users
master_init(formula = num_awards ~ prog + math,family = 'poisson',maxit=25,tol= 1e-08)

repeat{
  ###this happen in node 1
  master=readRDS("master.Rds")
  data_user1=read.csv('data_user1.csv')
  node_beta(Data = data_user1,master = master,userID = 'user1')
  rm(data_user1,master)
  
  ###this happen in node 2
  master=readRDS("master.Rds")
  data_user2=read.csv('data_user2.csv')
  node_beta(Data = data_user2,master = master,userID = 'user2')
  rm(data_user2,master)
  
  ###this happen in node 3
  master=readRDS("master.Rds")
  data_user3=read.csv('data_user3.csv')
  node_beta(Data = data_user3,master = master,userID = 'user3')
  rm(data_user3,master)
  
  ###############################
  ###this happen in the master ##
  ###############################
  user_1=readRDS("user1.Rds")
  user_2=readRDS("user2.Rds")
  user_3=readRDS("user3.Rds")
  master=readRDS("master.Rds")
  master_beta(user_1,user_2,user_3,master = master)
  rm(user_1,user_2,user_3,master)
  ###############################
  ###############################
  
  ###this happen in node 1
  master=readRDS("master.Rds")
  data_user1=read.csv('data_user1.csv')
  node_deviance(Data = data_user1,master = master,userID = 'user1')
  rm(master,data_user1)
  
  ###this happen in node 2
  master=readRDS("master.Rds")
  data_user2=read.csv('data_user2.csv')
  node_deviance(Data = data_user2,master = master,userID = 'user2')
  rm(master,data_user2)
  
  ###this happen in node 3
  master=readRDS("master.Rds")
  data_user3=read.csv('data_user3.csv')
  node_deviance(Data = data_user3,master = master,userID = 'user3')
  rm(master,data_user3)
  
  
  ###############################
  ###this happen in the master ##
  ###############################
  user_1=readRDS("user1.Rds")
  user_2=readRDS("user2.Rds")
  user_3=readRDS("user3.Rds")
  master=readRDS("master.Rds")
  master_deviance(user_1,user_2,user_3,master=master)
  rm(user_1,user_2,user_3,master)
  ###############################
  ###############################
  #convergence check
  master=readRDS("master.Rds")
  if(master$converged) break
  rm(master)
}

master



