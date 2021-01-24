
######## ANÁLISIS COLORACIÓN ###########

#Preparar datos
BD_Color = read.table(file = "clipboard", sep = "\t", header = TRUE)
str(BD_Color)

# Separar la BD por filas por subespecie

BD_Color_split = split(BD_Color, BD_Color$Subespecie)
str(BD_Color_split)

# Escribir un archivo independiente para cada subespecie

Color_cayan = BD_Color_split$cayanensis
Color_hell = BD_Color_split$hellmayri
Color_rufi = BD_Color_split$rufipennis
Color_ery = BD_Color_split$erythropterus


# Crear una función para evaluar la moda
getmode = function(x) {
        uniqv = unique(x)
        uniqv[which.max(tabulate(match(x, uniqv)))]
}

# encontrar la moda de cada columna (colores) en los diferentes subespecies

x = sapply(Color_cayan, getmode)
writeClipboard(x)

y = sapply(Color_hell, getmode)
writeClipboard(y)

z = sapply(Color_rufi, getmode)
writeClipboard(z)

w = sapply(Color_ery, getmode)
writeClipboard(w)







######## ANÁLISIS MORFOMÉTRICO ###########

###El análisis de MANOVA re realizó tomando como referencia la información contenida en la pagina:
###https://www.datanovia.com/en/lessons/one-way-manova-in-r/#data-preparation

#paquetes necesarios
install.packages("Hmisc")
library(Hmisc)
install.packages("mice")
library(mice)
library(tidyverse)
install.packages("ggpubr")
library(ggpubr)
install.packages("rstatix")
library(rstatix)
library(car)
install.packages("broom")
library(broom)

#Directorio de trabajo
setwd(choose.dir())
BD.morpho=read.csv("BD Morpho_TOTAL.csv", header=TRUE, fill = TRUE)

#Depurar columnas
BD.morpho=BD.morpho[,c(1:5,7,8,11,12,14,17,18,19,21, 24)]


#Crear un df con datos numericos
BD.morphnumeric=BD.morpho[,c(6:15)]

#Asegurar datos as.numeric
BD.morphnumeric=apply(BD.morphnumeric, MARGIN = 2,as.double)


length(which(complete.cases(BD.morphnumeric)==TRUE))
length(which(complete.cases(BD.morphnumeric)==FALSE))


### Evalación de supuestos de MANOVA ###

## Evaluar correlación = pearson
correlation=cor(BD.morphnumeric, use = "complete.obs", method = "pearson")
corrconp=rcorr(BD.morphnumeric, type = "pearson")

write.csv(corrconp[[1]], "clipboard")
write.csv(corrconp[[3]], "clipboard")


# Evaluar NA's

p_missing <- unlist(lapply(as.data.frame(BD.morphnumeric), function(x) sum(is.na(x))))/nrow(as.data.frame(BD.morphnumeric))
sort(p_missing[p_missing > 0], decreasing = TRUE)


#creación de BD análisis
BD.analysis=cbind( BD.morpho[,c(1:5)], as.data.frame(BD.morphnumeric))

#verificar tamaño muestral
BD.analysis[,c(3,6:15)] %>%
        group_by(Subespecie) %>%
        summarise(N = n())

## Evaluar valores extremos (outliers) univariados

extremos=as.data.frame(apply(BD.analysis[,c(6:15)], MARGIN = 2, is_extreme))
listadeextremos=apply(extremos, MARGIN = 2, function(x) which(x==TRUE))

#reemplazar extremos por NA's
for (i in 6:15){BD.analysis[,i][which(extremos[,i-5]==TRUE)]=NA}

extremos2=as.data.frame(apply(BD.analysis[,c(6:15)], MARGIN = 2, is_extreme))
listadeextremos2=apply(extremos2, MARGIN = 2, function(x) which(x==TRUE))                                
listadeextremos2


## Lidiar con NA's

par(mfrow=c(5,2), mar=c(1.5,1.5,1.5,1.5),cex=0.7)
nombres=variable.names(BD.analysis)[6:15]

for (i in 6:15){hist(BD.analysis[,i], main = nombres[i-5])
        abline(v=median(BD.analysis[,i], na.rm = TRUE), col="green")
        abline(v=mean(BD.analysis[,i], na.rm = TRUE), col="red")}

#reemplazar los datos perdidos por la media
for (i in 6:15){BD.analysis[,i][is.na(BD.analysis[,i])]=mean(BD.analysis[,i], na.rm = TRUE)}

p_missing <- unlist(lapply(as.data.frame(BD.analysis), function(x) sum(is.na(x))))/nrow(as.data.frame(BD.analysis))
sort(p_missing[p_missing > 0], decreasing = TRUE)



## Evaluar Outliers Multivariados

maht=mahalanobis_distance(BD.analysis[,c(6:15)])

Spp1=subset(BD.analysis[,c(3,6:15)], Subespecie=="cayanensis")
spp2=subset(BD.analysis[,c(3,6:15)], Subespecie=="rufipennis")
spp3=subset(BD.analysis[,c(3,6:15)], Subespecie=="hellmayri")
spp4=subset(BD.analysis[,c(3,6:15)], Subespecie=="erythropterus")

MDsp1=mahalanobis_distance(Spp1[,c(2:11)])
MDsp2=mahalanobis_distance(spp2[,c(2:11)])
MDsp3=mahalanobis_distance(spp3[,c(2:11)])
MDsp4=mahalanobis_distance(spp4[,c(2:11)])

MDt=list("cayanensis"=subset(MDsp1$is.outlier, MDsp1$is.outlier==TRUE),
         "rufipennis"=subset(MDsp2$is.outlier, MDsp2$is.outlier==TRUE),
         "hellmayri"=subset(MDsp3$is.outlier, MDsp3$is.outlier==TRUE),
         "erythropterus"=subset(MDsp4$is.outlier, MDsp4$is.outlier==TRUE))

BD.analysis=BD.analysis[-190,]


## Evaluar normalidad univariada

Spp1=subset(BD.analysis[,c(3,6:15)], Subespecie=="cayanensis")
spp2=subset(BD.analysis[,c(3,6:15)], Subespecie=="rufipennis")
spp3=subset(BD.analysis[,c(3,6:15)], Subespecie=="hellmayri")
spp4=subset(BD.analysis[,c(3,6:15)], Subespecie=="erythropterus")

normalidad_uniivariada=list("cayanensis"=apply(Spp1[,-1], MARGIN = 2, shapiro.test),
                            "rufipennis"=apply(spp2[,-1], MARGIN = 2, shapiro.test),
                            "hellmayri"=apply(spp3[,-1], MARGIN = 2, shapiro.test),
                            "erythropterus"=apply(spp4[,-1], MARGIN = 2, shapiro.test))

normalidad_uniivariadadf=data.frame(row.names = c("cayanenis", "rufipennis","hellmayri","erythropterus"))
for (i in 1:4 ){for (j in 1:10) {normalidad_uniivariadadf[i,j]=normalidad_uniivariada[[i]][[j]][["p.value"]]}}
colnames(normalidad_uniivariadadf)=nombres


write.csv(normalidad_uniivariadadf, "Test_norm_univariada_p_values.csv")


par(mfrow=c(2,2))

qqnorm(spp3$Alto.pico..nivel.de.narina., main="alto del pico")
qqline(spp3$Alto.pico..nivel.de.narina., col="red")
qqnorm(spp3$Ancho.distal.del.tarso, main="Ancho distal del Tarso")
qqline(spp3$Ancho.distal.del.tarso, col="red")
hist(spp3$Alto.pico..nivel.de.narina., main="", xlab = "")
hist(spp3$Ancho.distal.del.tarso, main="", xlab="")

## Evaluar normalidad multivariada

mshapiro_test(BD.analysis[,c(6:15)])

install.packages("MVN")
library(MVN)

Mardia<-mvn(data= BD.analysis[,c(6:15)],mvnTest="mardia")
Mardia$multivariateNormality

HZ<-mvn(data= BD.analysis[,c(6:15)],mvnTest="hz")
HZ$multivariateNormality

Royston=mvn(data= BD.analysis[,c(6:15)],mvnTest="royston")
Royston$multivariateNormality

dev.off()
qqplotm=mvn(BD.analysis[,c(6:15)], multivariatePlot = "qq")

MVNtpdf=list()
for (i in 6:15){MVNtpdf[[i-5]]=lm(BD.analysis[,i]~BD.analysis$Subespecie)$residuals}

lapply(MVNtpdf, shapiro.test)


## Evaluar homogenidad de  covarianza (homocedastecidad)

box_m(BD.analysis[, c(6:15)], BD.analysis$Subespecie)

BD.analysis %>% 
        gather(key = "variable", value = "value", nombres) %>%
        group_by(variable) %>%
        levene_test(value ~ Subespecie)

GGally::ggpairs(BD.analysis[,c(6:15)])
GGally::ggpairs(Spp1[,c(2:11)])
GGally::ggpairs(spp2[,c(2:11)])
GGally::ggpairs(spp3[,c(2:11)])
GGally::ggpairs(spp4[,c(2:11)])


### Ejecutar Modelos ###

## Modelo 1 ## 

modelo1 <- lm(cbind(Culmen.desde.narina..desde.el.borde.distal.,Alto.pico..nivel.de.narina., Tarso,Ancho.distal.del.tarso,Cola,Ancho.de.la.rectriz.central,Cuerda.alar..Calibrador.,Largo.P10,Largo.S1,Largo.P1) ~ Subespecie, BD.analysis)
manovamo=Manova(modelo1, test.statistic="Pillai")
summary(manovamo)


## Modelo 3 ##

modelo3 <- lm(cbind(Culmen.desde.narina..desde.el.borde.distal., Tarso,Cola,Cuerda.alar..Calibrador.) ~ Subespecie, BD.analysis)
manova3=Manova(modelo3, test.statistic="Pillai")


BD.ana4=BD.analysis[, c(1:5,6,8,10,12)]

grouped.data <- BD.ana4[,c(3,6:9)] %>%
        gather(key = "variable", value = "value", nombres[c(1,3,5,7)]) %>%
        group_by(variable)

grouped.data %>% anova_test(value ~ Subespecie) %>% .$p  %>% p.adjust(method = "BH")  %>% data.frame("variable"=nombres[c(1,3,5,7)])

#Test post-hoc

#prueba de Welch
grouped.data %>% welch_anova_test(value ~ Subespecie) %>% .$p  %>% p.adjust(method = "BH")  %>% data.frame("variable"=nombres[c(1,3,5,7)])


pwctukey <- BD.ana4[,c(3,6:9)]  %>%
        gather(key = "variables", value = "value", nombres[c(1,3,5)]) %>%
        group_by(variables) %>%
        tukey_hsd(value ~ Subespecie) %>%
        select(-estimate, -conf.low, -conf.high) # Remove details
pwcgames_howell= BD.ana4[,c(3,6:9)]  %>%
        gather(key = "variables", value = "value", nombres[c(1,3,5)]) %>%
        group_by(variables) %>%
        games_howell_test(value ~ Subespecie) %>%
        select(-estimate, -conf.low, -conf.high) # Remove details

install.packages("DTK")
library(DTK)
BD.ana4$Subespecie=gsub("hellmayri", "h",x =BD.ana4$Subespecie )
BD.ana4$Subespecie=gsub("rufipennis", "r",x =BD.ana4$Subespecie )
BD.ana4$Subespecie=gsub("cayanensis", "c",x =BD.ana4$Subespecie )
BD.ana4$Subespecie=gsub("erythropterus", "e",x =BD.ana4$Subespecie )

PWCDTK=list()
par(mfrow=c(2,2))
for (i in 6:9){PWCDTK[[i-5]]=DTK.test(x=BD.ana4[,i], f=BD.ana4$Subespecie)
DTK.plot(PWCDTK[[i-5]])}
DTK.plot(PWCDTK[[2]])
PWCDTK=as.data.frame(do.call(cbind, PWCDTK))
PWCDTK=PWCDTK[,-c(1,4,7,10)]

# Graficos
BD.analysis2 = BD.analysis[, -218]

BD.analysis2[,c(3,4, 6, 8, 10, 12)] %>% filter ( Sexo =="M" | Sexo =="F") %>%  
        pivot_longer(c(-Subespecie, -Sexo), names_to="Variables", values_to ="Lenghtmm")%>%
        ggplot(aes(x=Subespecie, y=Lenghtmm, color=Sexo ))+geom_violin()+ylab("Length (mm)")+ 
        scale_x_discrete(labels=c("c","e","h","r"))+facet_wrap(~ Variables, scales = "free")


BD.analysis[,c(3,4, 6, 8, 10, 12)] %>% filter ( Sexo =="M" | Sexo =="F") %>%  
        pivot_longer(c(-Subespecie, -Sexo), names_to="Variables", values_to ="Lenghtmm")%>%
        ggplot(aes(x=Subespecie, y=Lenghtmm, color=Sexo ))+geom_violin()+ylab("Length (mm)")+ 
        scale_x_discrete(labels=c("c","e","h","r"))+facet_wrap(~ Variables, scales = "free")



### Dimorfismo sexual ###

ttestsexo=BD.analysis[,c(3:4,6:15)]   %>% filter ( Sexo =="M" | Sexo =="F") %>% 
        pivot_longer(-c(Subespecie, Sexo), names_to="variable", values_to="value") %>% 
        group_by(variable, Subespecie) %>% 
        t_test(value ~ Sexo) %>%
        adjust_pvalue(method = "bonferroni") %>%
        select(!c(.y., p )) %>% 
        add_significance("p.adj")


BD.analysis[,c(3:4,6:15)]   %>% filter ( Sexo =="M" | Sexo =="F") %>% 
        pivot_longer(-c(Subespecie, Sexo), names_to="variable", values_to="value") %>%
        ggplot(aes(x=Sexo, y=value))+geom_boxplot()+facet_wrap(~Subespecie + variable, scales = "free")

BD.analysis$Subespecie=factor(BD.analysis$Subespecie, levels=c("hellmayri" ,    "rufipennis"  ,  "cayanensis"   , "erythropterus"))
BD.analysis[,c(3,4,6,8,10,12)] %>% rename(Culmen=Culmen.desde.narina..desde.el.borde.distal., Tail=Cola,Wing=Cuerda.alar..Calibrador.,Tarsus=Tarso) %>%  pivot_longer(c(-Subespecie, -Sexo), names_to="Variables", values_to ="Lenghtmm")%>%
        ggplot(aes(x=Subespecie, y=Lenghtmm , col=Subespecie))+geom_violin()+ylab("Length (mm)")+xlab("Subspecies")+theme_bw()+
        theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
        geom_jitter(shape=16, position=position_jitter(0.2))+scale_x_discrete(labels=c("h", "r", "c", "e"))+
        facet_wrap(~ Variables, scales = "free")


### K means approach fro clustering ###

setwd(choose.dir())

BD.analysis=read.csv("BDanalysis.csv", header=TRUE)

#escalar las variables
BD.sc4kmeans=scale(BD.analysis[,-c(1:5)])
rownames(BD.sc4kmeans)=BD.analysis$Subespecie
BD.sc4kmeans=as.data.frame(BD.sc4kmeans)



#fijar una semilla
set.seed(123)


install.packages("factoextra")
library(factoextra)

# Elbow method
fviz_nbclust(BD.sc4kmeans, kmeans, method = "wss") +
        labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(BD.sc4kmeans, kmeans, method = "silhouette")+
        labs(subtitle = "Silhouette method")

# Gap statistic
fviz_nbclust(BD.sc4kmeans, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
        labs(subtitle = "Gap statistic method")



km=kmeans(BD.sc4kmeans, centers = 7, nstart = 50, iter.max = 1000000, algorithm = "Hartigan-Wong")
BD.sc4kmeansc=cbind(BD.sc4kmeans, cluster=km$cluster)
BD.sc4kmeansc=cbind(BD.sc4kmeansc, BD.analysis[,1:5])
write.csv(BD.sc4kmeansc, "Clusterizacion7.csv")
fviz_cluster(km, BD.sc4kmeans)


#boxplot
nombres=names(BD.analysis[,-c(1:5)])
par(mfrow=c(2,5), cex=0.5)
for (i in 1:10){boxplot(BD.analysis[,i+5]~BD.analysis$Sexo:BD.analysis$Subespecie, xlab="", ylab =nombres[i]) }

##crear un df solo de hembras y machos
BDhembras=subset(BD.analysis, Sexo=="F")
BDmachos=subset(BD.analysis, Sexo=="M")

BDhembras[,-c(1:5)]=scale(BDhembras[,-c(1:5)], center = TRUE, scale = TRUE)
row.names(BDhembras)=make.names(BDhembras$Subespecie, unique = TRUE)
BDmachos[,-c(1:5)]=scale(BDmachos[,-c(1:5)], center = TRUE, scale = TRUE)
row.names(BDmachos)=make.names(BDmachos$Subespecie, unique = TRUE)


#correr algoritmo de clusterizacion

kmhmebras=kmeans(BDhembras[, -c(1:5)], centers = 7, nstart = 50, iter.max = 1000000, algorithm = "Hartigan-Wong")
fviz_cluster(kmhmebras, BDhembras[,-c(1:5,16)])

#añadir columna de clasificacion al df original
BDhembrasclus=cbind(BDhembras, cluster=as.factor(kmhmebras$cluster))
write.csv(BDhembrasclus, "clusterizacionhembras.csv")

#machos
kmmachos=kmeans(BDmachos[,-c(1:5)],centers = 7, nstart = 50, iter.max=1000000, algorithm = "Hartigan-Wong")
fviz_cluster(kmmachos, BDmachos[,-c(1:5)])
BDmachosclus=cbind(BDmachos, cluster=as.factor(kmmachos$cluster))
write.csv(BDmachosclus, "clusterizacionmachos.csv")

tapply(BDmachosclus$cluster, BDmachosclus$Subespecie, summary)
tapply(BDhembrasclus$cluster, BDhembrasclus$Subespecie, summary)


##clusterizacion jerarquica
clusterizacion=function(x,y){plot(hclust(dist(x, method = "euclidean"),method = y), labels=row.names(x), main = paste(deparse(substitute(x)),y))}
clusterizacion(BDhembras[,-c(1:5)], "complete")

abline(h=6, col="red")

clusterizacion(BDmachos[,-c(1:5)], "complete")
abline(h=6, col="red")










######## ANÁLISIS VOCALIZACIONES ###########

#Descargar cantos de wikiaves  (www.wikiaves.com.br)

#Set directory

install.packages("remotes")
remotes::install_github("Athospd/wikiaves")
library(wikiaves)
birds_metadata <- wa_metadata("Myiozetetes cayanensis")
wa_download(birds_metadata, path = choose.dir())



##Descargar cantos de Xenocanto (www.xenocanto.org)
install.packages("warbleR")
require(warbleR)


##Set directory

## Query xeno-canto for all recordings of the species Myiozetetes cayanensis
Myio.cayan <- querxc(qword = "Myiozetetes cayanensis", download = TRUE) 
str(Myio.cayan)
write.csv(Myio.cayan, file = "BD Xenocanto.csv")


##Image map of species from data set downloaded
xcmaps(X = Myio.cayan, img = TRUE, it = "tiff")
xcmaps(X = Myio.cayan, img = FALSE)





### ANOVA ###



## cargar datos
setwd()
library(tidyverse)
BDanovaf=read.csv('BD_ANOVA_UTOS.csv', header = TRUE)
BDanovaf= BDanovaf %>% rename(Indice_div=indice_de_diversidad, Cadencia=Cadencia_de_notas_C_.cadencia_del_trino.)


BDanovaf$UTO=as.factor(BDanovaf$UTO)
violinesgrananova=BDanovaf[,-1] %>% pivot_longer(!UTO, names_to="variable", values_to="value") %>% 
        ggplot(aes(UTO, y=value, fill=UTO, col=UTO))+geom_violin()+
        facet_wrap(~variable, scales = "free")

ggsave("violinesgrananova.png", violinesgrananova, device = "png", width = 210, height = 297, units = "mm" )
ggsave("violinesgrananova.eps", violinesgrananova, device = "eps", width = 210, height = 297, units = "mm" )
ggsave("violinesgrananova.svg", violinesgrananova, device = "svg", width = 210, height = 297, units = "mm" )


#ejecición de ANOVA
variables=names(BDanovaf)[3:5]
anovas=lapply(variables, function(x) aov(BDanovaf[,x]~BDanovaf$UTO))
names(anovas)=variables
pvalues=sapply(1:3, function(x) sapply(anovas, summary)[[x]][["Pr(>F)"]][1])
pvaluesadj=p.adjust(pvalues, method = "BH")
resultados=data.frame("Variable"=names(anovas), "p-value corregido"=pvaluesadj)
resultados


## TEST post-hoc
tukeys= TukeyHSD(anovas[["Indice_div"]],  ordered=TRUE)
tukeys=as.data.frame(tukeys$`BDanovaf$UTO`)
tukeys$comp=rownames(tukeys)

tukeys$variable = "Indice_div"


tukeysplots1= tukeys %>% select(lwr, upr, comp) %>% ggplot()+geom_linerange(aes(x=comp, ymax=upr, ymin=lwr))+
        geom_point(aes(x=comp, y=lwr, color="brown"))+ geom_point(aes(x=comp, y=upr, color="brown"))+
        coord_flip()+xlab("Comparacion de UTOS")+ylab("Intervalo de confianza")+
        geom_hline(aes(yintercept = 0), lwd=0.1, color="red", alpha=0.5)+ggtitle("Tukey HSD para el indice de diversidad")+
        theme_light()+ theme(legend.position = 'n')


#Numero de Notas
tukeys2= TukeyHSD(anovas[["No_notas"]],  ordered=TRUE)
tukeys2=as.data.frame(tukeys2$`BDanovaf$UTO`)
tukeys2$comp=rownames(tukeys2)

tukeys2$variable = "No_notas"

tukeys_finales = rbind(tukeys, tukeys2)

tukeysplots2= tukeys2 %>% select(lwr, upr, comp) %>% ggplot()+geom_linerange(aes(x=comp, ymax=upr, ymin=lwr))+
        geom_point(aes(x=comp, y=lwr, color="brown"))+ geom_point(aes(x=comp, y=upr, color="brown"))+
        coord_flip()+xlab("Comparacion de UTOS")+ylab("Intervalo de confianza")+
        geom_hline(aes(yintercept = 0), lwd=0.1, color="red", alpha=0.5)+ggtitle("Tukey HSD para el Numero de Notas")+
        theme_light()+ theme(legend.position = 'n')

tukeysplots = plot_grid(tukeysplots1, tukeysplots2)

ggsave("tukeysplotsanovatotal.png", tukeysplots, device = "png", width = 210, height = 297, units = "mm" )
ggsave("tukeysplotsanovatotal.eps", tukeysplots, device = "eps", width = 210, height = 297, units = "mm" )
ggsave("tukeysplotsanovatotal.svg", tukeysplots, device = "svg", width = 210, height = 297, units = "mm" )





## Clusterización vocalizaciones ##

setwd()
BDsinutos=read.csv("BDcantosresumida.csv", header = TRUE)
BDsinutos$X.1=NULL
BDsinutos$View=NULL;BDsinutos$Begin.File=NULL
matcdf=read.csv("cordenadas y utos for match.csv", header = TRUE)


BDsinutos$file=gsub(".WAV","",BDsinutos$file);BDsinutos$file=gsub(".wav","",BDsinutos$file);BDsinutos$file=gsub(".mp3","",BDsinutos$file)
matcdf$Archivo.Mp3=gsub(".WAV","",matcdf$Archivo.Mp3);matcdf$Archivo.Mp3=gsub(".wav","",matcdf$Archivo.Mp3);matcdf$Archivo.Mp3=gsub(".mp3","",matcdf$Archivo.Mp3)
BDutos=merge(BDsinutos, matcdf, by.x = "file", by.y = "Archivo.Mp3",all = TRUE)
BDutos=BDutos[complete.cases(BDutos),]


##clusteriazcion utilizando todas las variables
library(factoextra)
library(plyr)
library(ggnewscale)
library(tidyverse)
library(cowplot)


BDutos[,-c(1,2)]=apply(BDutos[,-c(1:2)], 2, as.numeric)

row.names(BDutos)=make.unique(as.character(BDutos$ï..UTO))

BDutos=BDutos %>% filter(file != c("Myiozetetes-cayanensis-26132" ,"Myiozetetes-cayanensis-26132"))

BDutos=filter(BDutos, Note_type %in% c("A","C","D","Z"))
BDutos=split(BDutos, BDutos$Note_type)


lapply(BDutos, function(x) fviz_nbclust(x[,8:17], kmeans, method = "wss")) 

clusterkmeans4=lapply(BDutos, function(x) kmeans(x[,8:17],centers = 4,nstart = 25,iter.max = 100))

for (i in seq(length(BDutos))){BDutos[[i]]$kmeans=clusterkmeans4[[i]]$cluster}


## visualizar los clusters
kmeaclt=list()
for (i in seq(length(BDutos))){kmeaclt[[i]]=fviz_cluster(clusterkmeans4[[i]], 
                                                         BDutos[[i]][,8:17], 
                                                         ggtheme = theme_bw(), 
                                                         pointsize = 2, ellipse = TRUE, ellipse.alpha = 0, 
                                                         stand = TRUE)}

plot_grid(plotlist = kmeaclt, labels = names(BDutos), nrow = 2)


## añadir una columna al PCA de los utos originales y de los cluster hallados por el Kmean (como factor)
pca=lapply(BDutos, function(x) as.data.frame(prcomp(x[,8:17],scale. = TRUE)$x))
for(i in seq(length(pca))){pca[[i]]$UTO=BDutos[[i]]$ï..UTO
pca[[i]]$kmean=BDutos[[i]]$kmeans
pca[[i]]$Note_type=names(BDutos)[[i]]}



#unión de casco convexo
hulls=lapply(pca, function(x) x[,c(1:2,12)] %>% group_by(kmean) %>% 
                     slice(chull(PC1,PC2)))


kmeansplotslist=list()
for ( i in seq(length(pca))){kmeansplotslist[[i]]=ggplot(mapping=aes(PC1,PC2))+
        geom_polygon(data=hulls[[i]], aes(color=as.factor(kmean),  fill=as.factor(kmean)), alpha=0.1)+
        scale_colour_grey()+
        new_scale_color()+
        geom_point(data = pca[[i]], aes(color=as.factor(UTO), shape=as.factor(UTO) ), size=1.5)+
        scale_colour_hue(l = 60, c = 150)+theme_minimal()}

kmeanscantos=plot_grid(plotlist = kmeansplotslist, labels = names(BDutos) , nrow=2)
ggsave("kmeanscantos.png", kmeanscantos, device = "png", width = 210, height = 297, units = "mm" )
ggsave("kmeanscantos.eps", kmeanscantos, device = "eps", width = 210, height = 297, units = "mm" )
ggsave("kmeanscantos.svg", kmeanscantos, device = "svg", width = 210, height = 297, units = "mm" )


lapply(pca, function(x) table(paste("UTO", x$UTO, sep=""), paste("k",x$kmean, sep = "")))




##Cluster jerarquico (clados) http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
library(dendextend)
distancias=lapply(BDutos, function(x) dist(x[,8:17], method = "euclidean"))
hc=lapply(distancias, function(x)hclust(x, method = "complete"))
arbol=list()
for (i in 1:4) {arbol[[i]]=hc[[i]] %>% as.dendrogram %>%
        set("branches_k_color", k=4) %>% set("branches_lwd", 1.2) %>%
        set("labels_colors") %>% set("labels_cex", c(.9,1.2)) %>% 
        set("leaves_pch", 19) %>% set("leaves_col", c(as.factor(pca[[i]]$UTO)))}
arbol=lapply(arbol, as.ggdend)

##visualizar
arbol[[1]];arbol[[2]];arbol[[3]];arbol[[4]]



hclusters=lapply(hc, function(x) cutree(x, k=4))

for (i in 1:4){BDutos[[i]]$hclus=as.factor(hclusters[[i]])
pca[[i]]$hclus=as.factor(hclusters[[i]])}
hullshc=lapply(pca, function(x) x[,c(1:2,14)] %>% group_by(hclus) %>% 
                       slice(chull(PC1,PC2)))
hclustcantoslist=list()
for ( i in seq(length(pca))){hclustcantoslist[[i]]=ggplot(mapping=aes(PC1,PC2))+
        geom_polygon(data=hullshc[[i]], aes(color=as.factor(hclus),  fill=as.factor(hclus)), alpha=0.1)+
        scale_colour_grey()+
        new_scale_color()+
        geom_point(data = pca[[i]], aes(color=as.factor(UTO), shape=as.factor(UTO) ), size=1.5)+
        scale_colour_hue(l = 60, c = 150)+theme_minimal()}

hclustcantos=plot_grid(plotlist = hclustcantoslist, labels = names(BDutos) , nrow=2)
ggsave("hclustcantos.png", hclustcantos, device = "png", width = 210, height = 297, units = "mm" )
ggsave("hclustcantos.eps", hclustcantos, device = "eps", width = 210, height = 297, units = "mm" )
ggsave("hclustcantos.svg", hclustcantos, device = "svg", width = 210, height = 297, units = "mm" )



varaov=which( names(BDutos[[1]]) %in% c("Low.Freq..Hz."	,"High.Freq..Hz.",	"Peak.Freq..Hz.",	
                                        "PFC.Max.Freq..Hz.",	"PFC.Min.Freq..Hz.",	"Dur.90...s.",
                                        "BW.90...Hz."	,"Delta.Freq..Hz."	,"Delta.Time..s."	,"Center.Freq..Hz."	,
                                        "Max.Freq..Hz."))
for (i in 1:4) {BDutos[[i]]$ï..UTO=as.factor(BDutos[[i]]$ï..UTO)}


violineslist=list()

for (i in 1:4){violineslist[[i]]=BDutos[[i]][,c(varaov, 19)] %>% rename(UTO=ï..UTO) %>% 
        pivot_longer(cols = !UTO, names_to="variables", values_to="value") %>% 
        ggplot(aes(x=UTO, y=value, col=UTO, fill=UTO))+geom_violin()+facet_wrap(~variables, scales = "free")+
        theme_light()+ylab("")+theme(axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank(),
                                     panel.grid.major = element_blank(),
                                     strip.background = element_rect(fill="white", color = "black"),
                                     strip.text = element_text(colour = "black"))}


violines=plot_grid(plotlist = violineslist, labels = names(BDutos) , nrow=2)
ggsave("violines.png", violines, device = "png", width = 210, height = 297, units = "mm" )
ggsave("violines.eps", violines, device = "eps", width = 210, height = 297, units = "mm" )
ggsave("violines.svg", violines, device = "svg", width = 210, height = 297, units = "mm" )




anovas=lapply(BDutos, function(x) lapply(varaov, function(y) aov(x[,y]~x$ï..UTO)))
for ( i in 1:4) {names(anovas[[i]])=names(BDutos[[i]])[varaov]}


pvalues=list()
for ( i in 1:4) {pvalues[[i]]=sapply(1:10, function(x) sapply(anovas[[i]], summary)[[x]][["Pr(>F)"]][1])}


pvalajustados=lapply(pvalues,function(x) p.adjust(x, method = "bonferroni"))


resultados=list()
for (i in 1:4) {resultados[[i]]=data.frame("Variable"=names(anovas[[i]]), "p-value corregido"=pvalajustados[[i]])}



signif=lapply(resultados, function(x) x %>% filter(p.value.corregido <= 0.05 ) %>% select(Variable))

signif=unlist(signif[[3]])
signif=anovas[["D"]][signif]

#test pos-hoc
tukeys=sapply(signif, TukeyHSD,  ordered=TRUE)
tukeys=lapply(tukeys, as.data.frame)
names(tukeys)=names(signif)
for (i in seq(length(tukeys))){tukeys[[i]]$variable=names(tukeys)[i]
tukeys[[i]]$comp=rownames(tukeys[[i]])}
tukeys=do.call(rbind, tukeys)
tukeysplotsTOTAL=tukeys %>% select(lwr, upr, variable, comp) %>% ggplot()+geom_linerange(aes(x=comp, ymax=upr, ymin=lwr))+
        geom_point(aes(x=comp, y=lwr, color="brown"))+ geom_point(aes(x=comp, y=upr, color="brown"))+
        coord_flip()+facet_wrap(~variable, scales = "free_x")+xlab("UTO comparison")+geom_hline(aes(yintercept = 0), lwd=0.1, color="red", alpha=0.5)+
        theme_light()+
        theme(legend.position = 'n',
              axis.title.x = element_blank())

ggsave("tukeysplotstotal.png", tukeysplotsTOTAL, device = "png", width = 210, height = 297, units = "mm" )
ggsave("tukeysplotstotal.eps", tukeysplotsTOTAL, device = "eps", width = 210, height = 297, units = "mm" )
ggsave("tukeysplotstotal.svg", tukeysplotsTOTAL, device = "svg", width = 210, height = 297, units = "mm" )

