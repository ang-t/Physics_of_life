library(tidyverse)
library(ggplot2)
library(rgl)
library(scatterplot3d)
library(plotly)
library(simecol)
library(ggpubr)
library(patchwork)
library(plot3D)

## Hodkin-Huxley model

HH <- odeModel(
  main = function(time, init, parms) {
    with(as.list(c(init, parms)),{
      
      am <- function(v) 0.1*(v+40)/(1-exp(-(v+40)/10))
      bm <- function(v) 4*exp(-(v+65)/18)
      ah <- function(v) 0.07*exp(-(v+65)/20)
      bh <- function(v) 1/(1+exp(-(v+35)/10))
      an <- function(v) 0.01*(v+55)/(1-exp(-(v+55)/10))
      bn <- function(v) 0.125*exp(-(v+65)/80)
      
      dv <- (I - gna*h*(v-Ena)*m^3-gk*(v-Ek)*n^4-gl*(v-El))/C
      dm <- am(v)*(1-m)-bm(v)*m
      dh <- ah(v)*(1-h)-bh(v)*h
      dn <- an(v)*(1-n)-bn(v)*n
      
      return(list(c(dv, dm, dh, dn)))
    })
  },
  ## Set parameters
  parms = c(Ena=50, Ek=-77, El=-54.4, gna=120, gk=36, gl=0.3, C=1, I=0),
  ## Set integrations times
  times = c(from=0, to=50, by = 0.25),
  ## Set initial state
  init = c(v=-65, m=0.052, h=0.596, n=0.317),
  solver = "lsoda"
)

HH <- sim(HH)
plot(HH)
mtext("Initial H-H model",                
      side = 3,
      line = - 2,
      outer = TRUE)


times(HH)["to"] <- 100

## Stimuli
Ilow<- c(2,2,2,2,2,2,2,2)
Imed <- c(5,5,5,5,5,5,5,5)
Ihigh<- c(6.5,6.5,6.5,6.5,6.5,6.5,6.5,6.5)
Iirreg <- c(6.5,2,6.2,5,5,2,6.5,2)

sims_low <- do.call("rbind",
                lapply(Ilow, function(i){
                  parms(HH)["I"] <- i
                  HH <- sim(HH)
                  cbind(I=paste("I =", i), out(HH))
                }))

sims_med <- do.call("rbind",
                   lapply(Imed, function(i){
                     parms(HH)["I"] <- i
                     HH <- sim(HH)
                     cbind(I=paste("I =", i), out(HH))
                   }))

sims_high <- do.call("rbind",
                   lapply(Ihigh, function(i){
                     parms(HH)["I"] <- i
                     HH <- sim(HH)
                     cbind(I=paste("I =", i), out(HH))
                   }))


sims_irreg <- do.call("rbind",
                   lapply(Iirreg, function(i){
                     parms(HH)["I"] <- i
                     HH <- sim(HH)
                     cbind(I=paste("I =", i), out(HH))
                   }))


simslow <- do.call("rbind",
                   lapply(I, function(i){
                     parms(HH)["I"] <- i
                     HH <- sim(HH)
                     cbind(I=paste("I =", i), out(HH))
                   }))

tot_time <- seq(0,(100*8+1.75),0.25)
df_stimuli<- cbind.data.frame("time"=tot_time,"low" = sims_low$v,"med"=sims_med$v,"high"=sims_high$v,"irreg"=sims_irreg$v)

## Plot stimuli

low <- ggplot(df_stimuli)+geom_line(mapping = aes(x=time,y=low))+
  ylim(-90,50)+ xlab(NULL) + ylab(NULL)+ theme_light()+ggtitle("I = 2")+
  theme(plot.title=element_text(hjust=0.5))

med <- ggplot(df_stimuli)+geom_line(mapping = aes(x=time,y=med))+
  ylim(-90,50)+ xlab(NULL) + ylab(NULL)+ theme_light()+ ggtitle("I = 5")+
  theme(plot.title=element_text(hjust=0.5))


high <- ggplot(df_stimuli)+geom_line(mapping = aes(x=time,y=high))+
  ylim(-90,50)+ xlab(NULL) + ylab(NULL)+ theme_light()+ ggtitle("I = 6.5")+
  theme(plot.title=element_text(hjust=0.5))


irreg <- ggplot(df_stimuli)+geom_line(mapping = aes(x=time,y=irreg))+
  ylim(-90,50)+ xlab(NULL) + ylab(NULL)+ theme_light() + ggtitle("I = 2 to 6.5")+
  theme(plot.title=element_text(hjust=0.5))

p = ggarrange(low,med,high,irreg)


annotate_figure(p,left="mV", bottom= "Time (ms)")

## attrtactor

f_low <- plot_ly(sims_low, x = ~v, y = ~m, z = ~(n+h), type = 'scatter3d', mode = 'lines',
                 scene="scene",name = 'I = 2')%>%
  layout(title= "I = 2",
         scene=(list(
           camera = list(
             eye = list(
               x = -1.2,
               y = 2,
               z = 0.3
             )
           ),
           xaxis = list(title = 'x'),
           yaxis = list(title = 'y'),
           zaxis = list(title = 'z'),aspectmode='cube')))


f_med <- plot_ly(sims_med, x = ~v, y = ~m, z = ~(n+h), type = 'scatter3d', mode = 'lines',scene="scene2",
                 name = 'I = 5')%>%
  layout(title= "I = 5",
         scene2=(list(
           camera = list(
             eye = list(
               x = -1.2,
               y = 2,
               z = 0.3
             )
           ),
           xaxis = list(title = 'x'),
           yaxis = list(title = 'y'),
           zaxis = list(title = 'z'),aspectmode='cube')))


f_high <- plot_ly(sims_high, x = ~v, y = ~m, z = ~(n+h), type = 'scatter3d', mode = 'lines',
                  scene="scene3",name = 'I = 6.5')%>%
  layout(title= "I = 6.5",
         scene3=(list(
           camera = list(
             eye = list(
               x = -1.2,
               y = 2,
               z = 0.3
             )
           ),
           xaxis = list(title = 'x'),
           yaxis = list(title = 'y'),
           zaxis = list(title = 'z'),aspectmode='cube')))


f_irreg <- plot_ly(sims_irreg, x = ~v, y = ~m, z = ~(n+h), type = 'scatter3d', mode = 'lines',
                   scene="scene4",name = 'I = 2 to 6.5')%>%
  layout(title= "I = 2 to 6.5",
         scene4=(list(
           camera = list(
             eye = list(
               x = -1.2,
               y = 2,
               z = 0.3
             )
           ),
           xaxis = list(title = 'x'),
           yaxis = list(title = 'y'),
           zaxis = list(title = 'z'),aspectmode='cube')))



fig_tot <- subplot(f_low,f_med,f_high,f_irreg)
fig_tot <- fig_tot %>% layout(title = "Dynamical System of the Hodgkin and Huxley model",
                      scene = list(domain=list(x=c(0,0.5),y=c(0.5,1))),
                      scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1))),
                      scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5))),
                      scene4 = list(domain=list(x=c(0.5,1),y=c(0,0.5))))

fig_tot

