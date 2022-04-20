## From (https://github.com/junesk9/)

library(ggplot2)
library(ggthemes)
names(windowsFonts())


Colours=c("#0984ff", "#FF0000")
Colours=c("#EE87B4", "#F9C270", "#54C3F1")

tjd<-read.csv("sum_TjD.pl.stdout.A_B_CD.csv")
ggplot()+theme_classic()+
  theme(axis.title.x = element_text(size = 20, family = "Arial"), 
        axis.title.y = element_text(size = 20, family = "Arial"), 
        axis.text.x = element_text(size = 20, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 20, colour = 1, family = "Arial"))+
  geom_density(aes(x = TajimaD, fill = Pop, colour=Pop), alpha = 0.4,  data = tjd)+
  scale_fill_manual(values = Colours)+
  scale_color_manual(values = Colours)



pi<-read.csv("sum_pi.pl.stdout.A_B_CD.csv")
ggplot()+theme_classic()+
  theme(axis.title.x = element_text(size = 10, family = "Arial"), 
        axis.title.y = element_text(size = 10, family = "Arial"), 
        axis.text.x = element_text(size = 20, colour = 1, family = "Arial"), 
        axis.text.y = element_text(size = 20, colour = 1, family = "Arial"))+
  geom_density(aes(x = Pi, fill = Pop, colour=Pop), alpha = 0.4,  data = pi)+
  scale_fill_manual(values = Colours)+xlim(0,3e-5)+
  scale_color_manual(values = Colours)

  
  
 

                                    
