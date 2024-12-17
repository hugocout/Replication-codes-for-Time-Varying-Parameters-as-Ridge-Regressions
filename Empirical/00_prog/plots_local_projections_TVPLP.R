# =============================================================================================================
# INITIALIZATION
# =============================================================================================================

rm(list=ls()) # To clear the current environment
yymmdd = paste(substr(Sys.Date(),3,4),substr(Sys.Date(),6,7),substr(Sys.Date(),9,10),sep='')

wd = 'INSERT YOUR PATH HERE/Empirical/'
setwd(wd)
paths <- list(pro = "00_prog",
              dat = "10_data",
              too = "20_tools",
              fun = "20_tools/functions",
              out = "30_output",
              rst = "40_results")

# Functions & Libraries
# Array of packages to be checked and installed if not already
myPKGs <- c('plotly','pracma','extrafont','ggplot2','htmlwidgets','reshape2') 

# check to see if each package is installed, and if not add it to a list to install
InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]

# install any needed packages
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

library(plotly)
library(pracma)
library(extrafont)
library(ggplot2)
library(htmlwidgets)
library(reshape2)

# Plot custom functions
theme_bluewhite <- function (base_size = 11, base_family = "") {
  theme_linedraw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "white"),
      panel.grid.minor  = element_line(color = "white"),
      panel.background = element_rect(fill = "grey90"),
      panel.border = element_rect(color = "grey40", fill = NA)
    )
}


theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.3), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = 'black'),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.6, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.margin=margin(-20,5,5,-20),
            legend.box.margin=margin(-5,-5,-5,-5),
            #legend.text=element_text(size=12),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale(name = "","fill","Publication",manual_pal(values = c("#386cb0","#fdb462",brewer.pal('Dark2',n=3)[1],"#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462",brewer.pal('Dark2',n=3)[1],"#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# =============================================================================================================
# LOAD RESULTS
# =============================================================================================================

# Load initial results data
load(paste(paths$out, '/figure3/cs18results12.Rdata', sep='/'))
cBETAS_temp <- 0 * cBETAS
cBETAS_temp2 <- 0 * cBETAS

# Loop through various versions and models
for (v in 1:3) {
  for (mod in 2) {
    
    # Load TVP results
    load(paste(paths$out, '/figure3/cs18results', v, mod, '.Rdata', sep=''))
    cBETAS_temp[, mod, v, ] <- cBETAS[, mod, v, ]
    
    # BETA^2srr - BETA^ols
    cBETAS_temp2[,mod,v,]=cBETAS[,mod,v,]-t(repmat(cumsum(BOLS[50,]),dim(cBETAS_temp)[4],1))
    
    
  }
}
cBETAS <- cBETAS_temp
cBETAS2 <- cBETAS_temp2

# =============================================================================================================
# RESPONSE FUNCTIONS : Figure 3 (f, g, and h)
# =============================================================================================================

# Set up the graphics layout
par(mfrow=c(3, 2))

# Set parameters for the plot
mod <- 2 # Model index
graph_names <- c("GDP","Unemp","Inflation")

for(v in 1:3) {
  
  # Prepare the data matrix
  irfmat <- cBETAS[, mod, v, ]
  if (is.element(v, c(1, 3))) {
    irfmat <- irfmat * 100
  }
  
  # Create a 3D surface plot
  fig <- plot_ly(z = ~irfmat, colorbar = list(title = '')) %>% add_surface(
    contours = list(
      z = list(
        show = TRUE, usecolormap = TRUE, highlightcolor = "#ff0000",
        project = list(z = TRUE)
      )))
  
  # Camera and font settings
  f <- list(family = "Arial", size = 18, color = "#7f7f7f")
  x <- list(title = "Time", titlefont = f)
  y <- list(title = "Horizon", titlefont = f)
  
  # X-axis tick settings
  xticks <- as.Date(seq.Date(as.Date('1976-12-01'), as.Date('2020-12-01'), by = "month")[1:419], format = "%y")
  xticks <- format(xticks, "%y")
  
  # Layout adjustments for the plot
  fig <- fig %>% layout(
    title = "",
    scene = list(
      xaxis = list(
        title = "Years",
        tickmode = 'array',
        tickvals = c(seq(2, 419, by = 24)),
        ticktext = xticks[c(seq(2, 419, by = 24))]
      ),
      yaxis = list(title = "Horizon"),
      zaxis = list(title = "")
    )) 
  
  # Display the figure
  fig
  
  # Save the 3D figure into 2D figure
  saveWidget(fig, paste(paths$rst, '/figure3_TVP-LP_',graph_names[v],'.html', sep=''), selfcontained = TRUE)
}

# =============================================================================================================
# BETA^2srr  − BETA^ols for the cumulative effect of MP shocks (Figure 6)
# =============================================================================================================

xticks=as.Date(seq.Date(as.Date('1976-12-01'), as.Date('2020-12-01'), by = "month")[1:419])[1:419]
toplot = array(NA,dim=c(6,3,419))
hors=c(1,6,12,24,36,48)
mult=c(100,1,100)
reorderthis = c(2,3,1)
for(v in 1:3){
  for(mod in 2){
    toplot[,v,]=mult[reorderthis[v]]*(cBETAS2[hors,mod,reorderthis[v],])
  }}

dimnames(toplot)[[1]] = paste('h=',hors,sep='')
dimnames(toplot)[[2]] = c('Unemployment','CPI Inflation','GDP')
dimnames(toplot)[[3]] = xticks 

tp.long = reshape2::melt(toplot[,,])
colnames(tp.long)=c('Horizon','Variable','Years','value')
tp.long$Years=as.Date(tp.long$Years,origin='1970-01-01')

p=ggplot() + 
  geom_line(data = tp.long[,], mapping=aes(x=Years,y=value,color=Horizon), size = 1.5) + 
  facet_wrap( ~ Variable,scales="free")+
  scale_color_brewer(palette = "Dark2")+
  theme_bluewhite()+
  geom_vline(xintercept = as.Date('1992-01-01',origin='1970-01-01'), linetype="dashed", 
             color = "black", size=1)+
  geom_hline(yintercept = 0, #linetype="dashed", 
             color = "black", size=0.5)+
  labs(color='')+
  theme(axis.text.x=element_text(size=10,hjust=1,angle=0,family="Arial"), #angle=
        axis.text.y=element_text(size=10,family="Arial"),
        strip.text = element_text(face="bold", size=15,family="Arial"),
        legend.text=element_text(size=15,family="Arial"),
        legend.position="top")+
  theme_Publication()+
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0), #,family="Palatino Linotype"
    strip.text = element_text(face="bold", colour = "white",size=15,family="Palatino Linotype"), #
    strip.background=element_rect(colour="black",fill="black"))+
  xlab('')+ylab('')+theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=15),strip.text = element_text(face="bold", size=15,family="helvetica"))

plot(p)

# save
ggsave(paste(paths$rst, '/figure6_TVP-LP.png', sep=''), width = 300, height = 150,
       units = "mm", dev='png')

