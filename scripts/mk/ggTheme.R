ggTheme<-function(theme=1,size_title=12){
  if(theme==1){
    ggTheme<-theme_minimal() +
             theme(plot.title = element_text(size = size_title, face = "bold",hjust = 0.5),
                   legend.title = element_text(color = "Black", size = 12, face = "bold"),
                   legend.text=element_text(color = "Black", size = 12, face = "bold"),
                   legend.position = "bottom",
                   axis.title=element_text(size = 12,face="bold", colour = "black"),
                   axis.text = element_text(size = 12),
                   axis.ticks = element_line(colour='black'),
                   panel.background = element_blank(),
                   panel.grid.major =  element_line(colour = "grey90", size = 0.2),
                   panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
                   panel.border = element_rect(color='black',fill=NA))
    return(ggTheme)
  }
}
