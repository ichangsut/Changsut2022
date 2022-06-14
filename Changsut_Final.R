##STATS ON A21 DATAAA##
#######load some packages!##########
if(!require(rstatix)){
  install.packages("rstatix")
  library(rstatix)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(ggpubr)){
  install.packages("ggpubr")
  library(ggpubr)
}


##############SET YOUR WD :) setwd("yourwd")##
setwd("~/Desktop/Changsut_Biostats_Final")

#############input the data file##
data = read.delim("Changsut_Final_RawData.txt", row.names = 1)

#######we'll start with t-tests for each of our assays, so we can make some SICK boxplots#####


#############POX DATA##
data[1:10,1:5]
POX_ttest=data[,c(1,3)]
POX_ttest

##check for outliers##
POX_ttest %>% 
  group_by(Type) %>%
  identify_outliers(POX.mg.protein)

##remove the outlier##
POX_ttest
POX_ttest_noout=POX_ttest[row.names(POX_ttest) != "A21-004", , drop = FALSE]

##check for normality
POX_ttest_noout %>% 
  group_by(Type) %>%
  shapiro_test(POX.mg.protein)

##check for homogeneity of variances
POX_ttest_noout %>% levene_test(POX.mg.protein ~ Type)

##run a t-test
stat.test <- POX_ttest_noout %>% 
  t_test(POX.mg.protein ~ Type) %>%
  add_significance()
stat.test

##ok ttest is non-sig. What about a correlation to symbiont counts?
POX_corr=data[,c(2:3)]
cor.test(POX_corr$Symbiont.Density, POX_corr$POX.mg.protein)

##what happens if we transform them?
POX_corr$Sqrt.Symbiont = sqrt(POX_corr$Symbiont.Density)
cor.test(POX_corr$Sqrt.Symbiont, POX_corr$POX.mg.protein)

##without our outlier##
POX_corr_noout=POX_corr[row.names(POX_corr) != "A21-004", , drop = FALSE]
cor.test(POX_corr_noout$Symbiont.Density, POX_corr_noout$POX.mg.protein)
cor.test(POX_corr_noout$Sqrt.Symbiont, POX_corr_noout$POX.mg.protein)


############PO DATA##
data[1:10,1:5]
PPO_ttest=data[,c(1,4)]
PPO_ttest

##check for outliers##
PPO_ttest %>% 
  group_by(Type) %>%
  identify_outliers(PPO.mg.protein)

##check for normality
PPO_ttest %>% 
  group_by(Type) %>%
  shapiro_test(PPO.mg.protein)

##check for homogeneity of variances
PPO_ttest %>% levene_test(PPO.mg.protein ~ Type)

##run a t-test
stat.test <- PPO_ttest %>% 
  t_test(PPO.mg.protein ~ Type, var.equal = TRUE) %>%
  add_significance()
stat.test

##ok ttest is non-sig. What about a correlation to symbiont counts?
PPO_corr=data[,c(2:4)]
cor.test(PPO_corr$Symbiont.Density, PPO_corr$PPO.mg.protein)
##what happens if we transform them?
PPO_corr$Sqrt.Symbiont = sqrt(PPO_corr$Symbiont.Density)
cor.test(PPO_corr$Sqrt.Symbiont, PPO_corr$PPO.mg.protein)


##################CATALASE DATA##
data[1:10,1:5]
CAT_ttest=data[,c(1,5)]
CAT_ttest

##check for outliers##
CAT_ttest %>% 
  group_by(Type) %>%
  identify_outliers(CAT.mg.protein)

##remove the outlier##
CAT_ttest
CAT_ttest_noout=CAT_ttest[row.names(CAT_ttest) != "A21-008", , drop = FALSE]

##check for normality
CAT_ttest_noout %>% 
  group_by(Type) %>%
  shapiro_test(CAT.mg.protein)

##check for homogeneity of variances
CAT_ttest_noout %>% levene_test(CAT.mg.protein ~ Type)

##run a t-test
stat.test <- CAT_ttest_noout %>% 
  t_test(CAT.mg.protein ~ Type, var.equal = FALSE) %>%
  add_significance()
stat.test

##ok ttest is non-sig. What about a correlation to symbiont counts?
CAT_corr=data[,c(2,5)]
cor.test(CAT_corr$Symbiont.Density, CAT_corr$CAT.mg.protein)
##what happens if we transform them?
CAT_corr$Sqrt.Symbiont = sqrt(CAT_corr$Symbiont.Density)
cor.test(CAT_corr$Sqrt.Symbiont, CAT_corr$CAT.mg.protein)

##without our outlier##
CAT_corr_noout=CAT_corr[row.names(CAT_corr) != "A21-008", , drop = FALSE]
cor.test(CAT_corr_noout$Symbiont.Density, CAT_corr_noout$CAT.mg.protein)
cor.test(CAT_corr_noout$Sqrt.Symbiont, CAT_corr_noout$CAT.mg.protein)

####################MELANIN DATA##
data[1:10,1:6]
mel_ttest=data[,c(1,6)]
mel_ttest

##check for outliers##
mel_ttest %>% 
  group_by(Type) %>%
  identify_outliers(mg.melanin.mg.tissue)

##remove the outlier##
mel_ttest
mel_ttest_noout=mel_ttest[row.names(mel_ttest) != "A21-016", , drop = FALSE]

##check for normality
mel_ttest_noout %>% 
  group_by(Type) %>%
  shapiro_test(mg.melanin.mg.tissue)

##check for homogeneity of variances
mel_ttest_noout %>% levene_test(mg.melanin.mg.tissue ~ Type)

##run a t-test
stat.test <- mel_ttest_noout %>% 
  t_test(mg.melanin.mg.tissue ~ Type, var.equal = FALSE) %>%
  add_significance()
stat.test


#Merge correlation data with metadata into one table#
merged_data <- merge(data, CAT_corr_noout, by = "row.names")

#Filter columns out of merged data#
names(merged_data)
merged_data <- merged_data[c(1:2,6:7,12)]

################################TABLE USING STATS#########################
Assay <- c("Peroxidase", "Prophenoloxidase", "Catalase", "Antibacterial", "Melanin" )
Statistic_value <- c("-0.591", "-0.865", "2.02", "1.55", "4.40 ")
dF <- c("14.4", "17", "10.6", "17", "9.66")
p_value <- c("0.564", "0.399", "0.070", "0.139", "0.002")
Changsut_Final_Table <- data.frame(Assay, Statistic_value, dF, p_value)
names(Changsut_Final_Table)[2]<- "Statistic value" ###changing the name###
names(Changsut_Final_Table)[4]<- "p-value" ###changing the name###
Changsut_Final_Table
write.table(Changsut_Final_Table, file = "Changsut_Final_Table.txt", sep = "\t", quote = FALSE, row.names = FALSE)



#########################FIGURES#########################################
#colors first :) to match the coral
colors <- c('B' = "#895b2c",
            'W' = "white")
##Scatter plot catalase
CAT_plot <- ggplot(merged_data, aes(Sqrt.Symbiont, CAT.mg.protein.x)) +
  geom_point(pch = 21, color = "black", size = 2, alpha = 0.75, aes(fill = Type)) + 
  geom_smooth(method = "lm", color = "grey25", linetype = 5, alpha = .15) + scale_fill_manual(values = colors) +  xlab("Symbiont Density Transformed \n (cells/mL)") + ylab("Catalase Activity \n (Hydrogen peroxide scavenged/minute/mg protein)") + theme_bw() +
  theme(panel.grid = element_blank()) + theme(legend.position = "") +
  theme(strip.background = element_blank(), text = element_text(size = 8, face = "bold", color = "black")) + stat_cor(label.y=2500) + theme(strip.text = element_blank())
CAT_plot

#Scatter plot melanin
MEL_plot <- ggplot(merged_data, aes(Sqrt.Symbiont, mg.melanin.mg.tissue)) +
  geom_point(pch = 21, color = "black", size = 2, alpha = 0.75, aes(fill = Type)) + 
  geom_smooth(method = "lm", color = "grey25", linetype = 5, alpha = .15) + scale_fill_manual(values = colors) +  xlab("Symbiont Density Transformed \n (cells/mL)") + ylab("Melanin Production \n (mg Melanin/mg Tissue)") + theme_bw() +
  theme(panel.grid = element_blank()) + theme(legend.position = "") +
  theme(strip.background = element_blank(), text = element_text(size = 8, face = "bold", color = "black")) + stat_cor(label.y= 0.025) + theme(strip.text = element_blank())
MEL_plot

##Box plot for melanin
boxplotmel <- ggplot(data, aes(x=Type, y=mg.melanin.mg.tissue)) + geom_boxplot(aes(fill = Type)) + geom_point(pch = 21, size = 2, aes(fill = Type)) +
  scale_fill_manual(values = c('B' = "burlywood4", 'W' = "white")) + theme(legend.position = "none") + 
  theme_bw() + theme(legend.position = "none") + theme(panel.grid= element_blank()) + xlab("Symbiont state \n (cells/mL)") + ylab("Melanin Production \n (mg Melanin/mg Tissue)") +
  theme(strip.background = element_rect(fill = "black"), text = element_text(size = 8, face = "bold", color = "black")) + stat_cor(label.y=) + annotate(geom = "text", x= 0.75, y=0.025, label="p=0.002") + ylim(0,0.025)
boxplotmel

##Box plot for catalase
boxplotcat <- ggplot(data, aes(x=Type, y=CAT.mg.protein)) + geom_boxplot(aes(fill = Type)) + geom_point(pch = 21, size =2, aes(fill = Type)) +
  scale_fill_manual(values = c('B' = "burlywood4", 'W' = "white")) + theme(legend.position = "none") + 
  theme_bw() + theme(legend.position = "none") + theme(panel.grid= element_blank()) + xlab("Symbiotic state \n (cells/mL)") + ylab("Catalase Activity \n (Hydrogen peroxide scavenged/minute/mg protein)") +
  theme(strip.background = element_rect(fill = "black"), text = element_text(size = 8, face = "bold", color = "black")) + annotate(geom = "text", x= 0.75, y=2500, label="p=0.007") + ylim(0,2500)
boxplotcat


Multipanelfig <- ggarrange(boxplotmel, boxplotcat, MEL_plot, CAT_plot, nrow = 2, ncol = 2, align = "hv", legend = "none")
Multipanelfig

##########################SAVING THE PLOT AS A PDF##############################
ggsave("Changsut_Final_Fig.pdf", height = 7, width = 6)

