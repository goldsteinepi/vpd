#################
# VPDs and Proposed Bills, 2010 - 2017
# Citation: Goldstein ND, Purtle J, Suder JS. Effect of Vaccine-preventable Disease Incidence on Proposed State Vaccine Exemption Legislation. Manuscript in preparation.
# Requires: Analysis_data.csv (this dataset is more thoroughly described here: https://doi.org/10.2105/AJPH.2018.304765)
# 3/6/19 -- Neal Goldstein
#################


### FUNCTIONS ###

library("nlme") #linear mixed effects
library("lme4") #poisson mixed effects
library("tsModel") #time-series analysis

#check for overdispersion: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

library("rvest") #scrape html pages
library("tidycensus") #retrieve ACS data, note if error installing on MacOS see: https://github.com/r-quantities/units/issues/1

#takes an mmwr table and returns state case estimates
stateCases = function(mmwr_table, column, mmwr_year) {
  
  #initialize results dataframe
  results = data.frame("state"=state.name, "estimate"=NA, stringsAsFactors=F)
  
  #difference in reporting styles by year
  if (mmwr_year==2010 | mmwr_year==2011 | mmwr_year==2012 | mmwr_year==2013) {
    
    #recode missing and zero cases
    mmwr_table[which(mmwr_table[,column]=="—"), column] = 0
    mmwr_table[which(mmwr_table[,column]=="N"), column] = NA
    mmwr_table[which(mmwr_table[,column]=="U"), column] = NA
    
    #handle each state seperately
    for (i in 1:nrow(results)) {
      
      #some states have to be handled as a special case
      if (results$state[i]=="New York") {
        results$estimate[i] = as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]=="New York (Upstate)",column], fixed=T)) + as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]=="New York City",column], fixed=T))
      } else {
        results$estimate[i] = as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]==results$state[i],column], fixed=T))
      }
    }
    
  } else if (mmwr_year==2014 | mmwr_year==2015) {
    
    #recode missing and zero cases
    mmwr_table[which(mmwr_table[,column]=="—"), column] = 0
    mmwr_table[which(mmwr_table[,column]=="N"), column] = NA
    mmwr_table[which(mmwr_table[,column]=="U"), column] = NA
    
    #handle each state seperately
    for (i in 1:nrow(results)) {
      
      #some states have to be handled as a special case
      if (results$state[i]=="New York") {
        results$estimate[i] = as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]=="New York (upstate)",column], fixed=T)) + as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]=="New York City",column], fixed=T))
      } else {
        results$estimate[i] = as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]==results$state[i],column], fixed=T))
      }
    }
    
  } else if (mmwr_year==2016 | mmwr_year==2017) {
    
    #recode missing and zero cases
    mmwr_table[which(mmwr_table[,column]=="—"), column] = 0
    mmwr_table[which(mmwr_table[,column]=="N"), column] = NA
    mmwr_table[which(mmwr_table[,column]=="U"), column] = NA
    
    #handle each state seperately
    for (i in 1:nrow(results)) {
      
      #some states have to be handled as a special case
      if (results$state[i]=="New York") {
        results$estimate[i] = as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]=="New York (excluding New York City)",column], fixed=T)) + as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]=="New York City",column], fixed=T))
      } else {
        results$estimate[i] = as.numeric(gsub(",", "", mmwr_table[mmwr_table[,1]==results$state[i],column], fixed=T))
      }
    }
    
  }
  
  
  #return a dataframe
  return(results)
  
}


### READ DATA ###

data2010 = read_html("https://www.cdc.gov/mmwr/preview/mmwrhtml/mm5953a1.htm")
data2011 = read_html("https://www.cdc.gov/mmwr/preview/mmwrhtml/mm6053a1.htm")
data2012 = read_html("https://www.cdc.gov/mmwr/preview/mmwrhtml/mm6153a1.htm")
data2013 = read_html("https://www.cdc.gov/mmwr/preview/mmwrhtml/mm6253a1.htm")
data2014 = read_html("https://www.cdc.gov/mmwr/volumes/63/wr/mm6354a1.htm?s_cid=mm6354a1_w")
data2015 = read_html("https://www.cdc.gov/mmwr/volumes/64/wr/mm6453a1.htm?s_cid=mm6453a1_w")
data2016_Dip = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2e.html")
data2016_Hib = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2f.html")
data2016_HepA_HepB_HepBPeri = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2g.html")
data2016_Flu_Pne = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2h.html")
data2016_Mea = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2i.html")
data2016_MenACWY_MenB_Mum_FluA = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2j.html")
data2016_Per = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2k.html")
data2016_Rub_RubCRS = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2l.html")
data2016_Tet = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2n.html")
data2016_VarMorb_VarMort = read_html("https://wonder.cdc.gov/nndss/static/2016/annual/2016-table2o.html")
data2017_Dip = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2e.html")
data2017_Hib = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2f.html")
data2017_HepA_HepB_HepBPeri = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2g.html")
data2017_Flu_Pne = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2h.html")
data2017_Mea = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2i.html")
data2017_MenACWY_MenB_Mum_FluA = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2j.html")
data2017_Per = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2k.html")
data2017_Rub_RubCRS = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2l.html")
data2017_Tet = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2n.html")
data2017_VarMorb_VarMort = read_html("https://wonder.cdc.gov/nndss/static/2017/annual/2017-table2o.html")

#state populations
census_api_key("paste api key here") #obtain key from: http://api.census.gov/data/key_signup.html
pop2010 = get_acs(geography="state", table="B01003", survey="acs5", year=2010)
pop2011 = get_acs(geography="state", table="B01003", survey="acs5", year=2011)
pop2012 = get_acs(geography="state", table="B01003", survey="acs5", year=2012)
pop2013 = get_acs(geography="state", table="B01003", survey="acs5", year=2013)
pop2014 = get_acs(geography="state", table="B01003", survey="acs5", year=2014)
pop2015 = get_acs(geography="state", table="B01003", survey="acs5", year=2015)
pop2016 = get_acs(geography="state", table="B01003", survey="acs5", year=2016)
pop2017 = get_acs(geography="state", table="B01003", survey="acs5", year=2017)


### PARSE VPDs ###

##VPDs
#Dip: Diphtheria
#Hib: Haemophilus influenzae --> serotype b
#HepA: Hepatitis, virus, acute --> A
#HepB: Hepatitis, virus, acute --> B
#HepBPeri: Hepatitis B perinatal infection
#Flu: Influenza-associated pediatric mortality
#Pne: Invasive pneumococcal disease --> All Ages (--> Confirmed, starting 2016)
#Mea: Measles, total
#MenACWY: Meningococcal disease., invasive, all serogroups --> serogroup A,C,Y, and W-135
#MenB: Meningococcal disease., invasive, all serogroups --> serogroup B
#Mum: Mumps
#FluA: Novel influenza A virus infections
#Per: Pertussis
#Pol: Polio
#Rub: Rubella
#RubCRS: Rubella, congenital syndrome
#Tet: Tetanus
#VarMorb: Varicella (Chickenpox) --> Morbidity
#VarMort: Varicella (Chickenpox) --> Mortality
VPDs = data.frame("Year"=rep(2010:2017,each=50), "State"=rep(state.name,length(2010:2017)), "Dip"=NA, "Hib"=NA, "HepA"=NA, "HepB"=NA, "HepBPeri"=NA, "Flu"=NA, "Pne"=NA, "Mea"=NA, "MenACWY"=NA, "MenB"=NA, "Mum"=NA, "FluA"=NA, "Per"=NA, "Pol"=NA, "Rub"=NA, "RubCRS"=NA, "Tet"=NA, "VarMorb"=NA, "VarMort"=NA, "Total"=NA, stringsAsFactors=F)

#2010
#VPDs$Dip[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[9], fill=T)[[1]], 4, 2010)[2])
VPDs$Hib[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[10], fill=T)[[1]], 3, 2010)[2])
VPDs$HepA[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[11], fill=T)[[1]], 2, 2010)[2])
VPDs$HepB[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[11], fill=T)[[1]], 3, 2010)[2])
#VPDs$HepBPeri[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[11], fill=T)[[1]], 7, 2010)[2])
VPDs$Flu[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[11], fill=T)[[1]], 6, 2010)[2])
VPDs$Pne[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[16], fill=T)[[1]], 3, 2010)[2])
VPDs$Mea[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[12], fill=T)[[1]], 6, 2010)[2])
VPDs$MenACWY[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[13], fill=T)[[1]], 3, 2010)[2])
VPDs$MenB[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[13], fill=T)[[1]], 4, 2010)[2])
VPDs$Mum[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[13], fill=T)[[1]], 7, 2010)[2])
VPDs$FluA[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[13], fill=T)[[1]], 8, 2010)[2])
VPDs$Per[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[14], fill=T)[[1]], 2, 2010)[2])
VPDs$Rub[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[15], fill=T)[[1]], 2, 2010)[2])
#VPDs$RubCRS[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[15], fill=T)[[1]], 5, 2010)[2])
VPDs$Tet[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[16], fill=T)[[1]], 8, 2010)[2])
#VPDs$VarMorb[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[18], fill=T)[[1]], 4, 2010)[2])
#VPDs$VarMort[VPDs$Year==2010] = unlist(stateCases(html_table(html_nodes(data2010, "table")[18], fill=T)[[1]], 5, 2010)[2])

#2011
#VPDs$Dip[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[9], fill=T)[[1]], 4, 2011)[2])
VPDs$Hib[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[10], fill=T)[[1]], 3, 2011)[2])
VPDs$HepA[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[11], fill=T)[[1]], 2, 2011)[2])
VPDs$HepB[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[11], fill=T)[[1]], 3, 2011)[2])
#VPDs$HepBPeri[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[11], fill=T)[[1]], 7, 2011)[2])
VPDs$Flu[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[11], fill=T)[[1]], 6, 2011)[2])
VPDs$Pne[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[16], fill=T)[[1]], 3, 2011)[2])
VPDs$Mea[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[12], fill=T)[[1]], 6, 2011)[2])
VPDs$MenACWY[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[13], fill=T)[[1]], 3, 2011)[2])
VPDs$MenB[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[13], fill=T)[[1]], 4, 2011)[2])
VPDs$Mum[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[13], fill=T)[[1]], 7, 2011)[2])
VPDs$FluA[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[13], fill=T)[[1]], 8, 2011)[2])
VPDs$Per[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[14], fill=T)[[1]], 2, 2011)[2])
VPDs$Rub[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[15], fill=T)[[1]], 2, 2011)[2])
#VPDs$RubCRS[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[15], fill=T)[[1]], 5, 2011)[2])
VPDs$Tet[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[16], fill=T)[[1]], 8, 2011)[2])
#VPDs$VarMorb[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[18], fill=T)[[1]], 4, 2011)[2])
#VPDs$VarMort[VPDs$Year==2011] = unlist(stateCases(html_table(html_nodes(data2011, "table")[18], fill=T)[[1]], 5, 2011)[2])

#2012
VPDs$Dip[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[9], fill=T)[[1]], 4, 2012)[2])
VPDs$Hib[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[10], fill=T)[[1]], 5, 2012)[2])
VPDs$HepA[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[11], fill=T)[[1]], 4, 2012)[2])
VPDs$HepB[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[11], fill=T)[[1]], 5, 2012)[2])
VPDs$HepBPeri[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[11], fill=T)[[1]], 7, 2012)[2])
VPDs$Flu[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[12], fill=T)[[1]], 2, 2012)[2])
VPDs$Pne[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[12], fill=T)[[1]], 3, 2012)[2])
VPDs$Mea[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[13], fill=T)[[1]], 2, 2012)[2])
VPDs$MenACWY[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[13], fill=T)[[1]], 6, 2012)[2])
VPDs$MenB[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[13], fill=T)[[1]], 7, 2012)[2])
VPDs$Mum[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[14], fill=T)[[1]], 2, 2012)[2])
VPDs$FluA[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[14], fill=T)[[1]], 3, 2012)[2])
VPDs$Per[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[14], fill=T)[[1]], 4, 2012)[2])
VPDs$Rub[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[15], fill=T)[[1]], 4, 2012)[2])
VPDs$RubCRS[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[15], fill=T)[[1]], 5, 2012)[2])
VPDs$Tet[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[17], fill=T)[[1]], 2, 2012)[2])
VPDs$VarMorb[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[18], fill=T)[[1]], 5, 2012)[2])
VPDs$VarMort[VPDs$Year==2012] = unlist(stateCases(html_table(html_nodes(data2012, "table")[18], fill=T)[[1]], 6, 2012)[2])

#2013
#VPDs$Dip[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[9], fill=T)[[1]], 4, 2013)[2])
VPDs$Hib[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[10], fill=T)[[1]], 5, 2013)[2])
VPDs$HepA[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[11], fill=T)[[1]], 4, 2013)[2])
VPDs$HepB[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[11], fill=T)[[1]], 5, 2013)[2])
VPDs$HepBPeri[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[11], fill=T)[[1]], 7, 2013)[2])
VPDs$Flu[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[12], fill=T)[[1]], 2, 2013)[2])
VPDs$Pne[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[12], fill=T)[[1]], 3, 2013)[2])
VPDs$Mea[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[13], fill=T)[[1]], 2, 2013)[2])
VPDs$MenACWY[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[13], fill=T)[[1]], 6, 2013)[2])
VPDs$MenB[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[13], fill=T)[[1]], 7, 2013)[2])
VPDs$Mum[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[14], fill=T)[[1]], 2, 2013)[2])
VPDs$FluA[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[14], fill=T)[[1]], 3, 2013)[2])
VPDs$Per[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[14], fill=T)[[1]], 4, 2013)[2])
VPDs$Rub[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[15], fill=T)[[1]], 4, 2013)[2])
VPDs$RubCRS[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[15], fill=T)[[1]], 5, 2013)[2])
VPDs$Tet[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[17], fill=T)[[1]], 2, 2013)[2])
VPDs$VarMorb[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[18], fill=T)[[1]], 4, 2013)[2])
VPDs$VarMort[VPDs$Year==2013] = unlist(stateCases(html_table(html_nodes(data2013, "table")[18], fill=T)[[1]], 5, 2013)[2])

#2014
VPDs$Dip[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[6], fill=T)[[1]], 3, 2014)[2])
VPDs$Hib[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[7], fill=T)[[1]], 5, 2014)[2])
VPDs$HepA[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[8], fill=T)[[1]], 4, 2014)[2])
VPDs$HepB[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[8], fill=T)[[1]], 5, 2014)[2])
VPDs$HepBPeri[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[8], fill=T)[[1]], 7, 2014)[2])
VPDs$Flu[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[9], fill=T)[[1]], 3, 2014)[2])
VPDs$Pne[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[9], fill=T)[[1]], 4, 2014)[2])
VPDs$Mea[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[10], fill=T)[[1]], 6, 2014)[2])
VPDs$MenACWY[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[11], fill=T)[[1]], 3, 2014)[2])
VPDs$MenB[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[11], fill=T)[[1]], 4, 2014)[2])
VPDs$Mum[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[11], fill=T)[[1]], 7, 2014)[2])
VPDs$FluA[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[11], fill=T)[[1]], 8, 2014)[2])
VPDs$Per[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[11], fill=T)[[1]], 9, 2014)[2])
VPDs$Rub[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[12], fill=T)[[1]], 9, 2014)[2])
VPDs$RubCRS[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[13], fill=T)[[1]], 2, 2014)[2])
VPDs$Tet[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[14], fill=T)[[1]], 6, 2014)[2])
VPDs$VarMorb[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[16], fill=T)[[1]], 2, 2014)[2])
VPDs$VarMort[VPDs$Year==2014] = unlist(stateCases(html_table(html_nodes(data2014, "table")[16], fill=T)[[1]], 3, 2014)[2])

#2015
#VPDs$Dip[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[6], fill=T)[[1]], 3, 2015)[2])
VPDs$Hib[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[9], fill=T)[[1]], 3, 2015)[2])
VPDs$HepA[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[10], fill=T)[[1]], 3, 2015)[2])
VPDs$HepB[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[10], fill=T)[[1]], 4, 2015)[2])
VPDs$HepBPeri[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[10], fill=T)[[1]], 6, 2015)[2])
VPDs$Flu[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[11], fill=T)[[1]], 3, 2015)[2])
VPDs$Pne[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[15], fill=T)[[1]], 6, 2015)[2])
VPDs$Mea[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[12], fill=T)[[1]], 3, 2015)[2])
VPDs$MenACWY[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[12], fill=T)[[1]], 7, 2015)[2])
VPDs$MenB[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[12], fill=T)[[1]], 8, 2015)[2])
VPDs$Mum[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[13], fill=T)[[1]], 2, 2015)[2])
VPDs$FluA[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[13], fill=T)[[1]], 3, 2015)[2])
VPDs$Per[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[13], fill=T)[[1]], 4, 2015)[2])
VPDs$Rub[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[14], fill=T)[[1]], 4, 2015)[2])
VPDs$RubCRS[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[14], fill=T)[[1]], 5, 2015)[2])
VPDs$Tet[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[16], fill=T)[[1]], 2, 2015)[2])
VPDs$VarMorb[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[17], fill=T)[[1]], 4, 2015)[2])
VPDs$VarMort[VPDs$Year==2015] = unlist(stateCases(html_table(html_nodes(data2015, "table")[17], fill=T)[[1]], 5, 2015)[2])

#2016
VPDs$Dip[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Dip, "table")[1], fill=T)[[1]], 6, 2016)[2])
VPDs$Hib[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Hib, "table")[1], fill=T)[[1]], 5, 2016)[2])
VPDs$HepA[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_HepA_HepB_HepBPeri, "table")[1], fill=T)[[1]], 5, 2016)[2])
VPDs$HepB[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_HepA_HepB_HepBPeri, "table")[1], fill=T)[[1]], 6, 2016)[2])
VPDs$HepBPeri[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_HepA_HepB_HepBPeri, "table")[1], fill=T)[[1]], 7, 2016)[2])
VPDs$Flu[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Flu_Pne, "table")[1], fill=T)[[1]], 3, 2016)[2])
VPDs$Pne[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Flu_Pne, "table")[1], fill=T)[[1]], 4, 2016)[2])
VPDs$Mea[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Mea, "table")[1], fill=T)[[1]], 7, 2016)[2])
VPDs$MenACWY[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_MenACWY_MenB_Mum_FluA, "table")[1], fill=T)[[1]], 3, 2016)[2])
VPDs$MenB[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_MenACWY_MenB_Mum_FluA, "table")[1], fill=T)[[1]], 4, 2016)[2])
VPDs$Mum[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_MenACWY_MenB_Mum_FluA, "table")[1], fill=T)[[1]], 7, 2016)[2])
VPDs$FluA[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_MenACWY_MenB_Mum_FluA, "table")[1], fill=T)[[1]], 8, 2016)[2])
VPDs$Per[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Per, "table")[1], fill=T)[[1]], 2, 2016)[2])
VPDs$Rub[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Rub_RubCRS, "table")[1], fill=T)[[1]], 4, 2016)[2])
VPDs$RubCRS[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Rub_RubCRS, "table")[1], fill=T)[[1]], 5, 2016)[2])
VPDs$Tet[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_Tet, "table")[1], fill=T)[[1]], 5, 2016)[2])
VPDs$VarMorb[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_VarMorb_VarMort, "table")[1], fill=T)[[1]], 5, 2016)[2])
VPDs$VarMort[VPDs$Year==2016] = unlist(stateCases(html_table(html_nodes(data2016_VarMorb_VarMort, "table")[1], fill=T)[[1]], 6, 2016)[2])

#2017
VPDs$Dip[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Dip, "table")[1], fill=T)[[1]], 6, 2017)[2])
VPDs$Hib[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Hib, "table")[1], fill=T)[[1]], 5, 2017)[2])
VPDs$HepA[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_HepA_HepB_HepBPeri, "table")[1], fill=T)[[1]], 5, 2017)[2])
VPDs$HepB[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_HepA_HepB_HepBPeri, "table")[1], fill=T)[[1]], 6, 2017)[2])
VPDs$HepBPeri[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_HepA_HepB_HepBPeri, "table")[1], fill=T)[[1]], 7, 2017)[2])
VPDs$Flu[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Flu_Pne, "table")[1], fill=T)[[1]], 3, 2017)[2])
VPDs$Pne[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Flu_Pne, "table")[1], fill=T)[[1]], 4, 2017)[2])
VPDs$Mea[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Mea, "table")[1], fill=T)[[1]], 7, 2017)[2])
VPDs$MenACWY[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_MenACWY_MenB_Mum_FluA, "table")[1], fill=T)[[1]], 3, 2017)[2])
VPDs$MenB[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_MenACWY_MenB_Mum_FluA, "table")[1], fill=T)[[1]], 4, 2017)[2])
VPDs$Mum[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_MenACWY_MenB_Mum_FluA, "table")[1], fill=T)[[1]], 7, 2017)[2])
VPDs$FluA[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_MenACWY_MenB_Mum_FluA, "table")[1], fill=T)[[1]], 8, 2017)[2])
VPDs$Per[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Per, "table")[1], fill=T)[[1]], 2, 2017)[2])
VPDs$Rub[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Rub_RubCRS, "table")[1], fill=T)[[1]], 4, 2017)[2])
VPDs$RubCRS[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Rub_RubCRS, "table")[1], fill=T)[[1]], 5, 2017)[2])
VPDs$Tet[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_Tet, "table")[1], fill=T)[[1]], 5, 2017)[2])
VPDs$VarMorb[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_VarMorb_VarMort, "table")[1], fill=T)[[1]], 5, 2017)[2])
VPDs$VarMort[VPDs$Year==2017] = unlist(stateCases(html_table(html_nodes(data2017_VarMorb_VarMort, "table")[1], fill=T)[[1]], 6, 2017)[2])

#totals
VPDs$VPDs_total = rowSums(VPDs[,3:21], na.rm=T)

#join state population data
VPDs$Population = NA
for (i in 1:nrow(VPDs)) {
  
  if (VPDs$Year[i]==2010) {
    VPDs$Population[i] = pop2010$estimate[pop2010$NAME==VPDs$State[i]]
  } else if (VPDs$Year[i]==2011) {
    VPDs$Population[i] = pop2011$estimate[pop2011$NAME==VPDs$State[i]]
  } else if (VPDs$Year[i]==2012) {
    VPDs$Population[i] = pop2012$estimate[pop2012$NAME==VPDs$State[i]]
  } else if (VPDs$Year[i]==2013) {
    VPDs$Population[i] = pop2013$estimate[pop2013$NAME==VPDs$State[i]]
  } else if (VPDs$Year[i]==2014) {
    VPDs$Population[i] = pop2014$estimate[pop2014$NAME==VPDs$State[i]]
  } else if (VPDs$Year[i]==2015) {
    VPDs$Population[i] = pop2015$estimate[pop2015$NAME==VPDs$State[i]]
  } else if (VPDs$Year[i]==2016) {
    VPDs$Population[i] = pop2016$estimate[pop2016$NAME==VPDs$State[i]]
  } else if (VPDs$Year[i]==2017) {
    VPDs$Population[i] = pop2017$estimate[pop2017$NAME==VPDs$State[i]]
  }
}
rm(i)

#per capita: 100,000 population
VPDs$VPDs_percapita = VPDs$VPDs_total/VPDs$Population*100000


### SAVE DATA ###

save(VPDs, file="VPD NNDSS 2010-17.RData")


### READ DATA ###

load("VPD NNDSS 2010-17.RData")
bills = read.csv("Analysis_data.csv", as.is=T, stringsAsFactors=F, na.strings="")


# ### CREATE BEFORE AND AFTER LAW ###
# 
# Law = data.frame("State"=bills$State[bills$Law=="Yes"], "Enacted"=bills$Year[bills$Law=="Yes"], "Year"=rep(2010:2017,each=sum(bills$Law=="Yes")), "Law"=NA, "VPDs"=NA, "Population"=NA, "VPDs_percapita"=NA, stringsAsFactors=F)
# 
# for (i in 1:nrow(Law))
# {
#   #VPDs by year
#   Law$VPDs[i] = VPDs$VPDs_total[VPDs$State==state.name[which(state.abb==Law$State[i])] & VPDs$Year==Law$Year[i]]
#   
#   #VPDs per 100000
#   Law$VPDs_percapita[i] = VPDs$VPDs_percapita[VPDs$State==state.name[which(state.abb==Law$State[i])] & VPDs$Year==Law$Year[i]]
# 
#   #population by year
#   Law$Population[i] = VPDs$Population[VPDs$State==state.name[which(state.abb==Law$State[i])] & VPDs$Year==Law$Year[i]]
# 
#   #add a one year lag; that is law effects won't happen until subsequent year
#   Law$Law[i] = ifelse(Law$Year[i]<=Law$Enacted[i], 0, 1)
# }
# rm(i)
# 
# #create time indicator (2010 as base)
# Law$Time = Law$Year - 2010


### COUNT BILLS ###

#add proposed legislation
VPDs$Bills_total = NA
VPDs$Bills_anti = NA
VPDs$Bills_pro = NA
VPDs$Bills_law = NA

for (i in 1:nrow(VPDs)) {
  
  #add a one year lag; that is VPDs in a given year influence bills the following year
  proposed = bills[bills$Year==(VPDs$Year[i]+1) & bills$State==state.abb[which(state.name==VPDs$State[i])], ]
  proposed_anti = proposed[proposed$Assessment=="Anti", ]
  proposed_pro = proposed[proposed$Assessment=="Pro", ]
  proposed_law = proposed[proposed$Law=="Yes", ]
  
  #count bills
  VPDs$Bills_total[i] = nrow(proposed)
  VPDs$Bills_anti[i] = nrow(proposed_anti)
  VPDs$Bills_pro[i] = nrow(proposed_pro)
  VPDs$Bills_law[i] = nrow(proposed_law)
}
rm(i, bills, proposed, proposed_anti, proposed_pro, proposed_law)

#create dichotomous indicators
VPDs$Bills_total_YN = ifelse(VPDs$Bills_total>0, 1, 0)
VPDs$Bills_anti_YN = ifelse(VPDs$Bills_anti>0, 1, 0)
VPDs$Bills_pro_YN = ifelse(VPDs$Bills_pro>0, 1, 0)
VPDs$Bills_law_YN = ifelse(VPDs$Bills_law>0, 1, 0)

#sort by state, year
VPDs = VPDs[with(VPDs, order(State, Year)), ]

#subset to only VPD years with bills the following year
VPDs = subset(VPDs, Year<=2016)

#create time indicator (2010 as base)
VPDs$Time = VPDs$Year - 2010


### ANALYSIS: NUMBER of BILLS ###

sum(VPDs$Bills_total)
sum(VPDs$Bills_anti)
sum(VPDs$Bills_pro)
sum(VPDs$VPDs_total)

#continuous outcome 

#choose a covariance structure
summary(lme(Bills_total ~ Time*VPDs_percapita, random=~1|State,data=VPDs,method="REML")) #no correlation
summary(lme(Bills_total ~ Time*VPDs_percapita, random=~1|State,correlation=corSymm(),data=VPDs,method="REML")) #symmetric
summary(lme(Bills_total ~ Time*VPDs_percapita, random=~1|State,correlation=corAR1(),data=VPDs,method="REML")) #autoregressive
summary(lme(Bills_total ~ Time*VPDs_percapita, random=~1|State,correlation=corCompSymm(),data=VPDs,method="REML")) #compound symmetry

#MODEL 1: Total bills by VPDs in prior year
summary(lme(Bills_total ~ Time*VPDs_percapita, random=~1|State, data=VPDs, method="ML"))

#MODEL 2: Anti bills by VPDs in prior year
summary(lme(Bills_anti ~ Time*VPDs_percapita, random=~1|State, data=VPDs, method="ML"))

#MODEL 3: Pro bills by VPDs in prior year
summary(lme(Bills_pro ~ Time*VPDs_percapita, random=~1|State, data=VPDs, method="ML"))

#MODEL 4: Enacted bills by VPDs in prior year
summary(lme(Bills_law ~ Time*VPDs_percapita, random=~1|State, data=VPDs, method="ML"))

#count outcome

#MODEL 1: Total bills by VPDs in prior year
summary(glmer(Bills_total ~ (1|State) + Time + scale(VPDs_percapita), offset=log(Population), data=VPDs, family=poisson()))
summary(glmer(Bills_total ~ (1|State) + Time*scale(VPDs_percapita), offset=log(Population), data=VPDs, family=poisson()))

#MODEL 2: Anti bills by VPDs in prior year
model = glmer(Bills_anti ~ (1|State) + Time + scale(VPDs_percapita), offset=log(Population), data=VPDs, family=poisson())
#summary(glmer(Bills_anti ~ (1|State) + Time*scale(VPDs_percapita), offset=log(Population), data=VPDs, family=poisson()))
overdisp_fun(model)
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#MODEL 3: Pro bills by VPDs in prior year
model = glmer(Bills_pro ~ (1|State) + Time + scale(VPDs_percapita), offset=log(Population), data=VPDs, family=poisson())
#summary(glmer(Bills_pro ~ (1|State) + Time*scale(VPDs_percapita), offset=log(Population), data=VPDs, family=poisson()))
overdisp_fun(model)
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#MODEL 4: Enacted bills by VPDs in prior year
summary(glmer(Bills_law ~ (1|State) + Time + scale(VPDs_percapita), offset=log(Population), data=VPDs, family=poisson()))
summary(glmer(Bills_law ~ (1|State) + Time*scale(VPDs_percapita), offset=log(Population), data=VPDs, family=poisson()))

#dichotmous outcome

#MODEL 1: Total bills by VPDs in prior year
summary(glmer(Bills_total_YN ~ (1|State) + Time*VPDs_percapita, data=VPDs, family=poisson()))

#MODEL 2: Anti bills by VPDs in prior year
summary(glmer(Bills_anti_YN ~ (1|State) + Time*VPDs_percapita, data=VPDs, family=binomial()))

#MODEL 3: Pro bills by VPDs in prior year
summary(glmer(Bills_pro_YN ~ (1|State) + Time*VPDs_percapita, data=VPDs, family=binomial()))

#MODEL 4: Enacted bills by VPDs in prior year
summary(glmer(Bills_law_YN ~ (1|State) + Time*VPDs_percapita, data=VPDs, family=binomial()))


# ### ANALYSIS: BEFORE and AFTER ###
# 
# #rate outcome
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="WA" & Law$Enacted=="2011",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="CA" & Law$Enacted=="2012",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="VT" & Law$Enacted=="2012",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="OR" & Law$Enacted=="2013",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="CO" & Law$Enacted=="2014",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="CA" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="CT" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="DE" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="IL" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="MT" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="VT" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="WV" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time, data=Law[Law$State=="UT" & Law$Enacted=="2017",]))
# 
# #adjust for seasonality
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="WA" & Law$Enacted=="2011",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="CA" & Law$Enacted=="2012",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="VT" & Law$Enacted=="2012",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="OR" & Law$Enacted=="2013",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="CO" & Law$Enacted=="2014",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="CA" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="CT" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="DE" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="IL" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="MT" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="VT" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="WV" & Law$Enacted=="2015",]))
# summary(lm(VPDs_percapita ~ Law*Time + harmonic(Year,2,7), data=Law[Law$State=="UT" & Law$Enacted=="2017",]))
# 
# #count outcome
# 
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="WA" & Law$Enacted=="2011",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="CA" & Law$Enacted=="2012",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="VT" & Law$Enacted=="2012",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="OR" & Law$Enacted=="2013",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="CO" & Law$Enacted=="2014",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="CA" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="CT" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="DE" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="IL" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="MT" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="VT" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="WV" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time, offset=log(Population), data=Law[Law$State=="UT" & Law$Enacted=="2017",], family=poisson()))
# 
# #adjust for seasonality
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="WA" & Law$Enacted=="2011",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="CA" & Law$Enacted=="2012",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="VT" & Law$Enacted=="2012",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="OR" & Law$Enacted=="2013",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="CO" & Law$Enacted=="2014",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="CA" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="CT" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="DE" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="IL" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="MT" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="VT" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="WV" & Law$Enacted=="2015",], family=poisson()))
# summary(glm(VPDs ~ Law*Time + harmonic(Year,2,7), offset=log(Population), data=Law[Law$State=="UT" & Law$Enacted=="2017",], family=poisson()))


### FIGURE ###

#tiff("Figure.tif",height=6,width=10,units='in',res=1200) 

#create a table of pro/anti vax bills by year
figure_data = t(matrix(data=c(as.numeric(by(VPDs$Bills_anti, VPDs$Year, FUN=sum)), as.numeric(by(VPDs$Bills_pro, VPDs$Year, FUN=sum))), nrow=length(unique(VPDs$Year)), ncol=2, byrow=F))

#side-by-side bar plot of pro/anti
par(mar = c(5,4,4,4) + 0.1)
bar_center = barplot(figure_data, xaxt="n", beside=T, ylim=c(0,30))

#plot VPDs, shift down for visual clarity
lines(x=colMeans(bar_center), y=by((VPDs$VPDs_percapita-5), VPDs$Year, FUN=median), lwd=4)

#axis labels
axis(1, at=colMeans(bar_center), labels=2011:2017)
axis(4, at=c(0,10,20,30), labels=c(5,15,25,35))
mtext(expression(italic("Total bills proposed")), side=2, line=2)
mtext(expression(italic("Legislative year")), side=1, line=2)
mtext(expression(italic("VPDs per 100,000")), side=4, line=2)

#legend
legend(1,30, c("Expand", "Restrict"), fill=gray.colors(2))

#dev.off()
