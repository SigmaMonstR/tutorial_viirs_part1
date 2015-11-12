dir = "/Users/jeffreychen/Google Drive/DOC/034-VIIRS"
setwd(dir)

###############
#IMPORT VIIRS##
###############

periods <- expand.grid(seq(1,9,2),c(2014,2015))
period <- paste(periods[,2],"-",periods[,1],sep="")

data <- data.frame()
for(i in period){
  temp <- read.csv(paste(dir,"/VIIRS_composites/",i,"/TNL_",i,".csv",sep=""))
  data <- rbind(data,temp)
}

#Save file
write.csv(data, paste(dir,"/VIIRS_composites/VIIRS_odd_mos.csv",sep=""),row.names=F)

##Remove empty records
data <- data[data$sum!=0,]

#Setup date var
data$date <- as.Date(paste(as.character(data$period),"-1",sep=""),"%Y-%m-%d")

#######
##BLS##
#######

#Wrangle BLS QCEW data for 2014 to 2015
years <- c(2014,2015)
for(year in years){
  
  bls_master <- data.frame()
  bls_cache <- data.frame()
  
  files <- list.files(paste(dir,"/BLS/",year,".q1-q4.by_area",sep=""))
  
  for(k in 1:length(files)){
    print(100*k/length(files))
    bls_temp<- read.csv(paste(dir,"/BLS/",year,".q1-q4.by_area/",files[k],sep=""))
    bls_temp <- bls_temp[,c("area_fips","industry_code","year","qtr",
                            "qtrly_estabs_count","month1_emplvl",
                            "month2_emplvl","month3_emplvl",  "total_qtrly_wages",
                            "taxable_qtrly_wages","qtrly_contributions","avg_wkly_wage")]
    bls_temp$area_fips <- as.character(bls_temp$area_fips) 
    bls_cache <- rbind(bls_temp,bls_cache)
    
    if(k%%50==0){
      bls_master <- rbind(bls_master,bls_cache)
      bls_cache <- data.frame()
      print("dump")
    } 
    if(k == length(files)){
      bls_master <- rbind(bls_master,bls_cache)
    }
    
  }
  
  write.csv(bls_master,paste(dir,"/bls_",year,".csv",sep=""),row.names=F)

  bls_process <- sqldf("SELECT area_fips, industry_code, year, qtr,
                               month1_emplvl, month2_emplvl, month3_emplvl
                               FROM bls_master
                               GROUP BY area_fips, industry_code, year, qtr")
  
  write.csv(bls_master,paste(dir,"/bls_processed_",year,".csv",sep=""),row.names=F)
}

#######################
###PROCESS BLS FILES###
########################
  bls_processed_2014 <- read.csv("~/Google Drive/DOC/034-VIIRS/bls_processed_2014.csv")
  bls_processed_2015 <- read.csv("~/Google Drive/DOC/034-VIIRS/bls_processed_2015.csv")
  
  bls_processed <- rbind(bls_processed_2014,bls_processed_2015)
  bls_processed$nchar <- nchar(as.character(bls_processed$industry_code))
  bls_processed <- bls_processed[bls_processed$nchar>=2 & bls_processed$nchar<=4 ,]
  
  bls_emp <- bls_processed[,c("area_fips","industry_code",
                                "year","qtr","month1_emplvl", "month2_emplvl", "month3_emplvl")]
  colnames(bls_emp) <- c("area_fips","industry_code",
                           "year","qtr","mo1",
                           "mo2","mo3")

  bls_emp2 <- reshape(bls_emp, 
             varying = c("mo1","mo2","mo3"), 
             v.names = "emp",
             timevar = "month", 
             times = c("mo1","mo2","mo3"), 
             direction = "long")
  bls_emp2$month <-  bls_emp2$qtr*3 - 3 + as.numeric(gsub("mo","",bls_emp2$month))
  bls_emp2$date <- as.Date(paste(bls_emp2$year,"-",bls_emp2$month,"-",1,sep=""),"%Y-%m-%d")
  
  bls_emp2<-bls_emp2[,c("area_fips","industry_code","date","emp")]

  bls_ready <- reshape(bls_emp2, 
               timevar = "industry_code",
               idvar = c("area_fips","date"),
               direction = "wide")
  bls_ready$GEOID <- as.numeric(as.character(bls_ready$area_fips))

#########################
##SEARCH FOR CORRELATES##
#########################

##Loop through each employment count 
##Order by Rho(Yhat_test,Y_test)

merged <- merge(data,bls_ready,by=c("GEOID","date"))
emp_list <- colnames(merged)
tester <- data.frame()
for(k in 21:ncol(merged)){
  temp <- merged[,c(k,3:18)]
  
  for(i in 1:ncol(temp)){
    temp[,i] <- log(temp[,i])
    temp <- temp[abs(temp[,i])!=Inf & !is.na(temp[,i]),]
  }
  temp$target<-temp[,1]
  temp[,1] <-NULL
  
  rho_pre_mean = cor(temp$mean, temp$target)
  rho_pre_sum = cor(temp$mean, temp$target)
  temp$rand = runif(nrow(temp))
  
  train <- temp[temp$rand<=0.7,]
  test <- temp[temp$rand>0.7,]
  
  fit = lm(target~.,data=train)
  train$yhat = predict(fit, train)
  test$yhat = predict(fit, test)
  
  rho_train = cor(train$yhat, train$target)
  rho_test = cor(test$yhat, test$target)
  
  tester <- rbind(tester,
                  data.frame(
                    var = emp_list[k],
                    rho_pre_mean = rho_pre_mean,
                    rho_pre_sum = rho_pre_sum,
                    rho_train = rho_train,
                    rho_test = rho_test,
                    n = nrow(temp)))
  
  print(paste(emp_list[k],": ",rho,"   ", nrow(temp),sep=""))
}
tester<-tester[order(-tester[,2]),]
