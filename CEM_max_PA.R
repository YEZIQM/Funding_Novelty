library("plyr")
library(Hmisc)
library("ggpubr")
library('cem')

# your data should contain you dep and indep vars, control vars, data type as factors for fixed effects 
load(file = '/data/home/bsw639/Expertise/R notebooks/paper0127.Rdata')

# add one comment

# add another comment
cal_delta <- function(treatment = NULL, data = NULL,vars = NULL,cutpoints = NULL,Y = 'mesh_inno'){
# you don't have to change this function;
# treatment: treatment, your main indep var;
# vars: all var names;
# cutpoints: specify number of cutpoints or cut intervals, see cem instruction;
# Y: your dep var;

time1 = Sys.time()
mat <- cem(treatment = treatment, data = data[vars],
         cutpoints = cutpoints)
    
# you can change the formula, add indep vars. However, it should be similar as using treatment only as CEM is designed for this
formula0 = as.formula(paste(Y, treatment, sep=" ~ "))
   
homo1 <- att(mat, formula = formula0 , data = data)
time3 = Sys.time()

print(time3 -time1)    
    
mean_ctl <- sum((1-data$high_max)*(data$mesh_inno*mat$w))/sum((1-data$high_max*mat$w))
    
beta = homo1$att.model['Estimate', treatment]
std = homo1$att.model['Std. Error', treatment]
p =  homo1$att.model['p-value', treatment]
CI_lower = beta - 1.96*std 
CI_higher = beta + 1.96*std 
num_C = homo1$tab['All','G0']
num_T = homo1$tab['All','G1']
num_C_matched = homo1$tab['Matched','G0']
num_T_matched = homo1$tab['Matched','G1']

    
res = c(beta,std,p,CI_lower,CI_higher,num_C,num_T,num_C_matched,num_T_matched,mean_ctl)

  
return(res)
}

# set the thresholds you want to define treatment and controlled groups
crts = c(0.01,0.05,0.1,0.25,0.5)
paras = c('beta','std','p','CI_lower','CI_higher','num_C','num_T','num_C_matched','num_T_matched','mean_ctl')
nc = length(paras)
nr = length(crts)
m <- matrix(1:nr*nc,nrow = nr ,ncol = nc)
k = 0

# your treatment var name
treatment = 'high_max'

# your dep var name
Y = 'mesh_inno'

# name of vars used to control the matching
vars = c('A' ,'B' ,'C' ,'D' ,'E', 'F', 'G', 'H' ,'I' ,'J' ,'K' ,'L' ,'M' ,'N','paper_num_MeSH','num_author',
         'year','prior_new','avg_pub_age','avg_SJR',treatment)

# 'cur' is the name of dataframe

# replace 'CV_pub_age' by the name of your treatment (4 of them)

# high_max is the dummy for treatment(1 for treatment, i.e. high in your mian var; 0 otherwise)
for (i in crts){
    k = k+1 
    tmp = cur
    low = unname(quantile(cur$CV_pub_age,i))
    high = unname(quantile(cur$CV_pub_age,1-i))
    tmp$high_max <- NA
    tmp[which(tmp$CV_pub_age>high),]$high_max<-1 
    tmp[which(tmp$CV_pub_age<low),]$high_max<-0
    tmp<-na.omit(tmp)    
    
    # the only thing you need to change is the cutpoints, for an indep var, the more important you think it is, the more cutpoints you 
    # ask for. Default is 10 for continuous vars, and every factor for factor vars.
    d= cal_delta(treatment = treatment, data = tmp,vars = vars,Y = Y,cutpoints = list(avg_pub_age=4,avg_SJR=4,paper_num_MeSH = 4,num_author=4))
#
    m[k,] = as.array(d)
}

rownames(m) = crts

colnames(m) = paras

library('foreign')

# change the out path
write.csv(m,file = '/data/home/bsw639/Expertise/notebooks/R_out/res_PA_MAX2_new.csv' )


