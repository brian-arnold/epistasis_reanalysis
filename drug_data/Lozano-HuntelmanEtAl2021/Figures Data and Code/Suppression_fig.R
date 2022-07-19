library(ggplot2)
library(scales)
library(forcats)
library(RColorBrewer)

setwd("~/Desktop/Lab_Stuff/Ch3_Antagonisms_3-4-5/materials/DATA/FIgures")

DA <- read.csv("~/Desktop/Lab_Stuff/Ch3_Antagonisms_3-4-5/materials/DATA/FIgures/DA.csv")
E <- read.csv("~/Desktop/Lab_Stuff/Ch3_Antagonisms_3-4-5/materials/DATA/FIgures/Emergent.csv")
#Sup_BD<- read.csv("~/Desktop/Lab_Stuff/Ch3_Antagonisms_3-4-5/materials/DATA/FIgures/BreakDown.csv")

DA$Order<- as.factor(DA$Order)
DA$Suppressive<- as.factor(DA$Suppressive)
DA$DA<- as.factor(DA$DA)

DA$og_DA<- -1
DA$og_DA[DA$DA=="Additive"]<-0
DA$og_DA[DA$DA=="Antagonistic"]<-0.5
DA$og_DA[DA$DA=="Antagonistic Suppression"]<- 1.3

DA$og_DA<- as.factor(DA$og_DA)
DA$Hid<- "yes"
DA$Hid[(DA$Assigned.Case==c("Not Suppressive"))] <- "no"
DA$Hid[(DA$Assigned.Case==c("No Hidden Suppression"))] <- "no"
DA$compare<-paste(DA$og_DA,DA$Hid)

DA$Assigned.Case[(DA$Assigned.Case==c("Suppression"))]<- "Hidden Suppression"

DA2<- DA[!(DA$Assigned.Case==c("Not Suppressive")),]

DA_1<-split(DA, DA$Order)
PerHid<-c()
Int_Hid_Sy<-c()
Int_Hid_Ad<-c()
Int_Hid_Ant<-c()
Int_Hid_AntS<-c()
Count_Hid<- c()
i=2

for (i in 1:3){
  total<- sum(DA_1[[i]]$Count)
  DA_1[[i]]$overall_per<- 100* DA_1[[i]]$Count/total
  hidden<- sum(DA_1[[i]]$Count[DA_1[[i]]$Hid==c("yes")])
  DA_1[[i]]$overall_hid<-100* DA_1[[i]]$Count/hidden
  DA_1[[i]]$overall_hid[DA_1[[i]]$Assigned.Case==c("Not Suppressive")]<-0 
  PerHid$per[i]<-hidden/total
  PerHid$Order[i]<- DA_1[[i]]$Order[1]
  Count_Hid[i]<- hidden
  Int_Hid_Sy[i]<- sum(DA_1[[i]]$Count[DA_1[[i]]$DA==c("Synergistic") & !DA_1[[i]]$Assigned.Case==c("Not Suppressive")])/sum(DA_1[[i]]$Count[DA_1[[i]]$DA==c("Synergistic")])
  Int_Hid_Ad[i]<-sum(DA_1[[i]]$Count[DA_1[[i]]$DA==c("Additive") & !DA_1[[i]]$Assigned.Case==c("Not Suppressive")])/sum(DA_1[[i]]$Count[DA_1[[i]]$DA==c("Additive")])
  Int_Hid_Ant[i]<-sum(DA_1[[i]]$Count[DA_1[[i]]$DA==c("Antagonistic") & !DA_1[[i]]$Assigned.Case==c("Not Suppressive")])/sum(DA_1[[i]]$Count[DA_1[[i]]$DA==c("Antagonistic")])
  Int_Hid_AntS[i]<-sum(DA_1[[i]]$Count[DA_1[[i]]$DA==c("Antagonistic Suppression") & !DA_1[[i]]$Assigned.Case==c("No Hidden Suppression")])/sum(DA_1[[i]]$Count[DA_1[[i]]$DA==c("Antagonistic Suppression")])
}
PerHid<- as.data.frame(PerHid)
PernotHid<- PerHid
PernotHid$per<- 1-PernotHid$per
PerHid$Hid<- c("Hid")
PernotHid$Hid<- c("notHid")
PerHid<- rbind(PerHid,PernotHid)

PerHid$Order[PerHid$Order=="1"]<-"3-Drug Combinations"
PerHid$Order[PerHid$Order=="2"]<-"4-Drug Combinations"
PerHid$Order[PerHid$Order=="3"]<-"5-Drug Combinations"

DA3<-rbind(DA_1[[1]], DA_1[[2]], DA_1[[3]])
DA3<- DA3[!(DA3$Count==0),]

#overall amount of hidden suppression
T<-sum(DA$Count)
S<- 1-sum(DA$Count[(DA$Assigned.Case==c("Not Suppressive"))])/T
ASup<- sum(DA$Count[(DA$DA==c("Antagonistic Suppression"))])/T
H<- sum(DA$Count[(DA$Hid==c("yes"))])/T
new<- 1-sum(DA$Count[(DA$DA==c("Antagonistic Suppression"))])/sum(DA$Count[(!DA$Assigned.Case==c("Not Suppressive"))])


DA3$order[DA3$Order=="3"]<-"3-Drug Combinations"
DA3$order[DA3$Order=="4"]<-"4-Drug Combinations"
DA3$order[DA3$Order=="5"]<-"5-Drug Combinations"
DA4<- DA3[!(DA3$Hid==c("no")),]

# pie chart for overall amount of hidden supression
ggplot(PerHid, aes(x="", y=per, fill= Hid))+
  geom_bar(stat="identity", width=1.2)+
  coord_polar("y", start=0)+
  facet_grid(facets=. ~ Order)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  labs(fill = NULL)+
  theme(legend.position="bottom")+
  theme(axis.text = element_blank(),axis.ticks = element_blank(), panel.grid  = element_blank())+
  scale_fill_manual(labels=c("Hidden Suppression","No Hidden Suppression"), values = c("azure4","azure3"))+
  theme(legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 18))

#bar compare 3, 4, 5
ggplot(DA3, aes(x=og_DA, y=overall_per, fill= compare))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  ylab("%")+
  xlab(NULL)+
  #theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position="bottom", legend.text=element_text(size=18))+
  scale_fill_manual(labels=c("Synergy","Synergy with HS","Additive","Additive with HS","Antagonism","Antogonism with HS","Suppression","Suppression with HS"), values = c("red", "red3", "gray80","gray60","olivedrab2","olivedrab","mediumturquoise","darkcyan"))+
  #scale_x_discrete(labels=c("Synergy","Additive","Antagomism","Suppression"))+
  theme(strip.text.x = element_text(size = 18))+
  theme(legend.title=element_blank())+
  theme( axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=18),axis.text.y=element_text(size=16),)

  


#looking at how hidden suppression is dist. based on on special cases

#by DA interaction
ggplot(DA4, aes(x=og_DA, y=overall_hid, fill= Assigned.Case))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  scale_x_discrete(labels=c("Synergy","Additive","Antagomism","Suppression"))+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("%")+
  xlab("net interaction")

# by special case
ggplot(DA4, aes(x=Assigned.Case, y=overall_hid, fill= og_DA))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size =14))+
  theme(axis.text.x=element_blank())+
  labs(fill = c("net interaction"))+
  theme(axis.text.y = element_text(size=14))+
  ylab(NULL)+
  xlab(NULL)+
  scale_fill_manual(labels=c("Synergy","Additive","Antagomism","Suppression"), values = c("red3","gray60","olivedrab","darkcyan"))+
  theme(legend.text=element_text(size=24))+
  theme(strip.text.x = element_text(size = 20))+
  theme(legend.title=element_blank())

ggplot(DA4, aes(x=Assigned.Case, y=overall_hid, fill= og_DA))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(fill = c("net interaction"))+
  ylab("%")+
  xlab("special case")+
  scale_fill_manual(labels=c("Synergy","Additive","Antagomism","Suppression"), values = c("red3","gray60","olivedrab","darkcyan"))+
  theme(legend.text=element_text(size=10))+
  theme(strip.text.x = element_text(size = 24))

  
DA5<- DA4[!(DA4$Assigned.Case==c("Fully Nested Hidden Suppression")),]
#DA5<- DA5[!(DA5$Assigned.Case==c("Hidden Suppression")),]
DA5<- DA5[!(DA5$Assigned.Case==c("Nested Hidden Suppression")),]

ggplot(DA5, aes(x=Assigned.Case, y=overall_hid, fill= og_DA))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, size=14))+
  theme(axis.text.x=element_blank())+
  labs(fill = c("net interaction"))+
  theme(axis.text.y = element_text(size=14))+
  ylab(NULL)+
  xlab(NULL)+
  scale_fill_manual(labels=c("Synergy","Additive","Antagomism","Suppression"), values = c("red3","gray60","olivedrab","darkcyan"))+
  theme(legend.text=element_text(size=24))+
  theme(strip.text.x = element_text(size = 20))+
  theme(legend.title=element_blank())


ggplot(DA5, aes(x=Assigned.Case, y=overall_hid, fill= og_DA))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(fill = c("net interaction"))+
  ylab("%")+
  xlab("special case")+
  scale_fill_manual(labels=c("Synergy","Additive","Antagomism","Suppression"), values = c("red3","gray60","olivedrab","darkcyan"))+
  theme(legend.text=element_text(size=10))+
  theme(strip.text.x = element_text(size = 12))


#Distribution Fig 6

Dist <- read.csv("~/Desktop/Lab_Stuff/Ch3_Antagonisms_3-4-5/materials/DATA/FIgures/Dist.csv")

#Dist$Amount<- as.factor(Dist$Amount)
Dist$Value<- Dist$Value*100

ggplot(Dist, aes(x=Amount, y=Value, fill=Int))+
  geom_bar(stat="identity", position = "dodge", width = .22)+
  facet_grid(facets=. ~ Order*Compare)+
  theme_minimal()+
  labs(x=NULL, y= c("%"), fill=NULL)+
  theme(legend.position="bottom")+
  scale_fill_manual(labels= c("Only Hidden Suppression", "Overall Suppressive"), values=c("gray28", "mediumturquoise")) +
  theme(strip.text.x = element_text(size = 8))+
  theme(legend.text=element_text(size=16))+
  theme(axis.title.y = element_text(size=16),axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="blank"))



#~~~~~~~~~~~~~~~~~~~~~~~~~Emergent Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


E$Order<- as.factor(E$Order)
E$Suppressive<- as.factor(E$Suppressive)
E$E<- as.factor(E$E)

E$og_E<- -1
E$og_E[E$E=="Additive"]<-0
E$og_E[E$E=="Antagonistic"]<-0.5
E$og_E[E$E=="Antagonistic Suppression"]<- 1.3

E$og_E<- as.factor(E$og_E)
E$Hid<- "yes"
E$Hid[(E$Assigned.Case==c("Not Suppressive"))] <- "no"
E$Hid[(E$Assigned.Case==c("No Hidden Suppression"))] <- "no"
E$compare<-paste(E$og_E,E$Hid)

E$Assigned.Case[(E$Assigned.Case==c("Suppression"))]<- "Hidden Suppression"

E2<- E[!(E$Assigned.Case==c("Not Suppressive")),]

E_1<-split(E, E$Order)
PerHid<-c()
Int_Hid_Sy<-c()
Int_Hid_Ad<-c()
Int_Hid_Ant<-c()
Int_Hid_AntS<-c()
Count_Hid<- c()

for (i in 1:3){
  total<- sum(E_1[[i]]$Count)
  E_1[[i]]$overall_per<- 100* E_1[[i]]$Count/total
  hidden<- sum(E_1[[i]]$Count[E_1[[i]]$Hid==c("yes")])
  E_1[[i]]$overall_hid<-100* E_1[[i]]$Count/hidden
  E_1[[i]]$overall_hid[E_1[[i]]$Assigned.Case==c("Not Suppressive")]<-0 
  PerHid$per[i]<-hidden/total
  PerHid$Order[i]<- E_1[[i]]$Order[1]
  Count_Hid[i]<- hidden
  Int_Hid_Sy[i]<- sum(E_1[[i]]$Count[E_1[[i]]$E==c("Synergistic") & !E_1[[i]]$Assigned.Case==c("Not Suppressive")])/sum(E_1[[i]]$Count[E_1[[i]]$E==c("Synergistic")])
  Int_Hid_Ad[i]<-sum(E_1[[i]]$Count[E_1[[i]]$E==c("Additive") & !E_1[[i]]$Assigned.Case==c("Not Suppressive")])/sum(E_1[[i]]$Count[E_1[[i]]$E==c("Additive")])
  Int_Hid_Ant[i]<-sum(E_1[[i]]$Count[E_1[[i]]$E==c("Antagonistic") & !E_1[[i]]$Assigned.Case==c("Not Suppressive")])/sum(E_1[[i]]$Count[E_1[[i]]$E==c("Antagonistic")])
  Int_Hid_AntS[i]<-sum(E_1[[i]]$Count[E_1[[i]]$E==c("Antagonistic Suppression") & !E_1[[i]]$Assigned.Case==c("No Hidden Suppression")])/sum(E_1[[i]]$Count[E_1[[i]]$E==c("Antagonistic Suppression")])
}
PerHid<- as.Eta.frame(PerHid)
PernotHid<- PerHid
PernotHid$per<- 1-PernotHid$per
PerHid$Hid<- c("Hid")
PernotHid$Hid<- c("notHid")
PerHid<- rbind(PerHid,PernotHid)

PerHid$Order[PerHid$Order=="1"]<-"3-Drug Combinations"
PerHid$Order[PerHid$Order=="2"]<-"4-Drug Combinations"
PerHid$Order[PerHid$Order=="3"]<-"5-Drug Combinations"

E3<-rbind(E_1[[1]], E_1[[2]], E_1[[3]])
E3<- E3[!(E3$Count==0),]

#overall amount of hidden suppression
T<-sum(E$Count)
S<- 1-sum(E$Count[(E$Assigned.Case==c("Not Suppressive"))])/T
ASup<- sum(E$Count[(E$E==c("Antagonistic Suppression"))])/T
H<- sum(E$Count[(E$Hid==c("yes"))])/T
new<- 1-sum(E$Count[(E$E==c("Antagonistic Suppression"))])/sum(E$Count[(!E$Assigned.Case==c("Not Suppressive"))])


E3$order[E3$Order=="3"]<-"3-Drug Combinations"
E3$order[E3$Order=="4"]<-"4-Drug Combinations"
E3$order[E3$Order=="5"]<-"5-Drug Combinations"
E4<- E3[!(E3$Hid==c("no")),]

#bar compare 3, 4, 5
ggplot(E3, aes(x=og_E, y=overall_per, fill= compare))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  ylab("%")+
  xlab(NULL)+
  #theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position="bottom", legend.text=element_text(size=18))+
  scale_fill_manual(labels=c("Synergy","Synergy with HS","Additive","Additive with HS","Antagonism","Antogonism with HS","Suppression","Suppression with HS"), values = c("red", "red3", "gray80","gray60","olivedrab2","olivedrab","mediumturquoise","Darkcyan"))+
  #scale_x_discrete(labels=c("Synergy","Additive","Antagonism","Suppression"))+
  theme(strip.text.x = element_text(size = 18))+
  theme(legend.title=element_blank())+
  theme( axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=18),axis.text.y=element_text(size=16),)




#looking at how hidden suppression is dist. based on on special cases

#by E interaction
ggplot(E4, aes(x=og_E, y=overall_hid, fill= Assigned.Case))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  scale_x_discrete(labels=c("Synergy","Additive","Antagomism","Suppression"))+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("%")+
  xlab("net interaction")

# by special case
ggplot(E4, aes(x=Assigned.Case, y=overall_hid, fill= og_E))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, size =14))+
  theme(axis.text.x=element_blank())+
  labs(fill = c("net interaction"))+
  theme(axis.text.y = element_text(size=14))+
  ylab(NULL)+
  xlab(NULL)+
  scale_fill_manual(labels=c("Synergy","Additive","Antagomism","Suppression"), values = c("red3","gray60","olivedrab","Darkcyan"))+
  theme(legend.text=element_text(size=24))+
  theme(strip.text.x = element_text(size = 24))+
  theme(legend.title=element_blank())



E5<- E4[!(E4$Assigned.Case==c("Fully Nested Hidden Suppression")),]
#E5<- E5[!(E5$Assigned.Case==c("Hidden Suppression")),]
E5<- E5[!(E5$Assigned.Case==c("Nested Hidden Suppression")),]

ggplot(E5, aes(x=Assigned.Case, y=overall_hid, fill= og_E))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, size=14))+
  theme(axis.text.x=element_blank())+
  labs(fill = c("net interaction"))+
  theme(axis.text.y = element_text(size=14))+
  ylab(NULL)+
  xlab(NULL)+
  scale_fill_manual(labels=c("Synergy","Additive","Antagomism","Suppression"), values = c("red3","gray60","olivedrab","Darkcyan"))+
  theme(legend.text=element_text(size=24))+
  theme(strip.text.x = element_text(size = 24))+
  theme(legend.title=element_blank())


ggplot(E5, aes(x=Assigned.Case, y=overall_hid, fill= og_E))+
  geom_bar(stat="identity", width=0.75)+
  facet_grid(facets=. ~ order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(fill = c("net interaction"))+
  ylab("%")+
  xlab("special case")+
  scale_fill_manual(labels=c("Synergy","Additive","Antagomism","Antagonistic\nSuppression"), values = c("red3","gray60","olivedrab","Darkcyan"))+
  theme(legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 24))

#~~~~~~~~~~~~~~~~~~~~~~~~~ OLD FIGURE IDEAS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#pie compare with 3,4,5
ggplot(DA3, aes(x="", y=overall_per, fill= compare))+
  geom_bar(stat="identity", width=0.75)+
  coord_polar("y", start=0)+
  facet_grid(facets=. ~ Order)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  labs(fill = NULL)+
  theme(legend.position="bottom")+
  theme(axis.text = element_blank(),axis.ticks = element_blank(), panel.grid  = element_blank())+
  scale_fill_manual(labels=c("Synergy","Synergy with HS","Additive","Additive with HS","Antagonism","Antogonism with HS","Antagonistic Suppression","Antagonistic Suppression with HS"), values = c("firebrick1", "firebrick3", "gray80","gray70","olivedrab1","olivedrab3","cyan2","cyan4"))+
  theme(strip.text.x = element_text(size = 18))

#breakdown of paths
Sup_BD1<- split(Sup_BD, Sup_BD$Order_down)

ggplot(Sup_BD1[[3]], aes(x="", y=per))+
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  facet_grid(facets=. ~ Drug)+
  theme_bw()+
  ylim(0,1)+
  xlab(NULL)+
  ylab(NULL)+
  scale_fill_manual(labels=NULL, values = c("green4"))+
  # ylab("Amount of hidden supresive interactions to the \n third lower order combiantion" )+
  theme(axis.text = element_blank(),axis.ticks = element_blank(), panel.grid  = element_blank())



