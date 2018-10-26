# Loading the packages          
require(Hmisc)          
require(MatchIt)          
require(fastDummies)   
require(rgenoud)   



# Loading the dataset        
dfh <- read.LIS('us16h')          
dfp <- read.LIS('us16p')          
df <- merge(dfh,dfp,by='hid')          


# selecting variables:        
df_int <- df[,c('hid','relation','emp','ptime','children', 'nchildren','age','sex','marital','disabled','educ','hours','pi','pil')]     
unique(df_int$children) 
unique(df_int$nchildren) 
dim(df_int)    


#imputing missing values of hours, using additive regression and bootstrapping     
# finding missing values     
missing_values <- function(data){     
  for (name in colnames(data)){     
    print(c(name, sum(is.na(data[name]))))     
    
  }     
}      
missing_values(data = df_int)      



impute_values <- aregImpute( ~ hours  + disabled +ptime+ age+ sex + pi+ pil,     
                             x=T,tlinear = F,type = 'regression',     
                             data = df_int,n.impute = 50)     

imputed_df <- data.frame(impute_values$x) 
df_int$imp_hours <- imputed_df$hours 
any(is.na(impute_values$x)) 
any(is.na(df_int)) 


##############################################################     
# subsetting only females     
df_mtr <- df_int[df_int$sex =='[2]female' ,]     


# creating female_caregiver variable     
df_mtr$female_caregiver <- ifelse(df_mtr$relation == '[2100]spouse'&     
                                    df_mtr$imp_hours <= 30,     
                                  df_mtr$female_caregiver <- 1, 0)     


#############################################################     



###########################################################     
# male aged 25-45 married/cohabiting is most likely to be a father.     
# so no. of treatemnt subject is going to be higher than the no. of control subject     
# susetting only father aged 25 to 45      
df_ftr <- df_int[df_int$age >= 25 & df_int$age <= 45 &     
                   df_int$sex =='[1]male' &     
                   (df_int$relation == '[1000]head'|df_int$relation == '[2100]spouse') ,]     

###########################################################     

#creating new variables     
df_ftr$father <- ifelse(df_ftr$children == '[110]living with own children aged 0-5'| 
                          df_ftr$children == '[120]living with own children 6-12'| 
                          df_ftr$children == '[130]living with own children 13-17'| 
                          df_ftr$children == '[140]living with own children 18+', 
                        df_ftr$father <- 1, 0) 

df_ftr$higher_edu <- ifelse(df_ftr$educ == '[3]high',df_ftr$higher_edu <- 1, 0) 
df_ftr$married <- ifelse(df_ftr$marital == '[110]married',df_ftr$married <- 1, 0) 
df_ftr$part_time <- ifelse(df_ftr$imp_hours <=30 ,df_ftr$part_time <- 1, 0) 
df_ftr$age_squared <- df_ftr$age^2 


summary(df_int$marital) 
summary(df_ftr$marital) 
summary(df_ftr$married) 



unique(df_int$children)   
unique(df_ftr$married)   
unique(df$marital)   


###############################################################     

# merging female_caregiver variable with father data and     
# condidering wife's income as other household income     

data <- merge(df_ftr,df_mtr[c('hid','female_caregiver','pi')],by='hid') 
colnames(data)[colnames(data)=='pi.y'] <- 'ohh_income'  
colnames(data)[colnames(data)=='pil'] <- 'total_income' 



# selecting the variable of interest   
var_int <- data[,c('father','age','age_squared','married','higher_edu','part_time','female_caregiver','disabled','ohh_income','total_income')] 
missing_values(data=var_int) 
summary(var_int$married) 



# creating dummy variables for education,part-time ,female_caregiver and disabled   
final_data <- dummy_cols(.data = var_int,select_columns = c('higher_edu','part_time','female_caregiver','married'), 
                         remove_most_frequent_dummy = TRUE) 


any(is.na(final_data)) 
missing_values(data=final_data) 
summary(final_data) 
describe(final_data)  
dim(final_data) 




#############################     
# beginnning of the project     
############################     

#################################################################     
### part-1 : treatment and control group population     
tg <- subset(final_data, father == 1) 
cg <- subset(final_data, father == 0) 
dim(tg) 
dim(cg) 
summary(tg) 
summary(cg) 


### fitting a logistic regression model to estimate the propensity scores     
# i.e. probablity of being a father     

# base model: 
base_model <- glm(formula = father  ~ age + age_squared + higher_edu_1 + married_0, 
                  family = binomial('logit'), data = final_data) 
prop_scores_base <- base_model$fitted.values 
summary(base_model) 


# full, genetic mathces all teat and control 
match_base <- matchit(formula = base_model, 
                      method = 'genetic', 
                      distance = 'logit', 
                      data = final_data) 

summary(match_base)   

# To create a dataframe containing only the matched observations, use the match.data() function:     
matched_obs_base<- match.data(match_base)     

# estimation of treatment effect     
treat_effect_base <- lm(total_income~father,data = matched_obs_base) 
summary(treat_effect_base)    



### fitting a logistic regression model to estimate the propensity scores      
# i.e. probablity of being a father      

h2a_model <- glm(formula = father  ~ age + age_squared + higher_edu_1 + married_0 + part_time_1, 
                 family = binomial('logit'), data = final_data) 
prop_scores_h2a <- h2a_model$fitted.values 
summary(h2a_model) 

match_h2a <- matchit(formula = h2a_model, 
                     method = 'genetic', 
                     distance = 'logit', 
                     data = final_data) 

summary(match_h2a) 

# To create a dataframe containing only the matched observations, use the match.data() function:     
matched_obs_h2a <- match.data(match_h2a) 



# estimation of treatment effect     
treat_effect_h2a <- lm(total_income~father,data = matched_obs_h2a) 
summary(treat_effect_h2a) 



### fitting a logistic regression model to estimate the propensity scores       
# i.e. probablity of being a father       

h2b_model <- glm(formula = father  ~ age + age_squared + higher_edu_1 + married_0 + part_time_1 + female_caregiver_1, 
                 family = binomial('logit'), data = final_data) 
prop_scores_h2b <- h2b_model$fitted.values 
summary(h2b_model) 


match_h2b <- matchit(formula = h2b_model, 
                     method = 'genetic', 
                     distance = 'logit', 
                     data = final_data) 

summary(match_h2b)  

# To create a dataframe containing only the matched observations, use the match.data() function:     
matched_obs_h2b <- match.data(match_h2b)   


# estimation of treatment effect     
treat_effect_h2b <- lm(total_income~father,data = matched_obs_h2b) 
summary(treat_effect_h2b)   




### fitting a logistic regression model to estimate the propensity scores       
# i.e. probablity of being a father       

# hypothesis H3(has a care partner):     
h3_fc1_model <- glm(formula = father  ~ age + age_squared + higher_edu_1 + married_0, 
                    family = binomial('logit'),  
                    data = final_data[final_data$female_caregiver==1,]) 
prop_scores_h3_fc1<- h3_fc1_model$fitted.values 
dim(final_data[final_data$female_caregiver==1,]) 
summary(h3_fc1_model) 

match_h3_fc1 <- matchit(formula = h3_fc1_model, 
                        method = 'genetic', 
                        distance = 'logit', 
                        data = final_data[final_data$female_caregiver==1,]) 

summary(match_h3_fc1) 

# To create a dataframe containing only the matched observations, use the match.data() function:     
matched_obs_h3_fc1 <- match.data(match_h3_fc1) 


# estimation of treatment effect     
treat_effect_h3_fc1<- lm(total_income~father,data = matched_obs_h3_fc1) 
summary(treat_effect_h3_fc1) 



### fitting a logistic regression model to estimate the propensity scores       
# i.e. probablity of being a father       

# hypothesis H3(does not have a care partner):     
h3_fc0_model <- glm(formula = father  ~ age + age_squared + higher_edu_1 + married_0, 
                    family = binomial('logit'),  
                    data = final_data[final_data$female_caregiver==0,]) 

prop_scores_h3_fc0<- h3_fc0_model$fitted.values 
summary(h3_fc0_model) 

dim(var_int[var_int$female_caregiver==0,]) 
dim(final_data[final_data$female_caregiver==0,]) 


match_h3_fc0 <- matchit(formula = h3_fc0_model, 
                        method = 'genetic', 
                        distance = 'logit', 
                        data = final_data[final_data$female_caregiver==0,]) 

summary(match_h3_fc0) 

# To create a dataframe containing only the matched observations, use the match.data() function:     
matched_obs_h3_fc0 <- match.data(match_h3_fc0)   


# estimation of treatment effect     
treat_effect_h3_fc0 <- lm(total_income~father,data = matched_obs_h3_fc0) 
summary(treat_effect_h3_fc0)   
