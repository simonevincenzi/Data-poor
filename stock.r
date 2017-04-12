##### themes and parameter values for plots
library(tidyverse)

size.title = 15
line.lwd = 1.2
size.label.x = 18
size.text.x = 14
size.point = 4
size.label.y = 18
size.text.y = 14
size.legend.text = 15
size.legend.title = 20
unit.legend.h = 1.8
unit.legend.w = 1.8
size.ann = 10
colour.axis = "gray20"
colour.theme = "black"
colour.axis.line = "gray20"
colour.line = "gray50"
max_size_dot = 8

## Theme to be used for all plots

theme.stock =  theme(plot.title = element_text(lineheight=.8, face="bold", size = size.title,hjust = 0.5), 
                   plot.background = element_blank()
                   ,panel.grid.major = element_blank()
                   ,panel.grid.minor = element_blank()
                   ,panel.border = element_blank()
                   ,panel.background = element_blank(),
                   axis.line = element_line(color = 'black'),
                   plot.margin = unit(c(1,2,1,1), "cm"),
                   axis.title.x = element_text(size=size.label.x,vjust=-80),
                   axis.text.x  = element_text(size=size.text.x, vjust = 0.5),
                   axis.title.y = element_text(size=size.label.x, vjust = 50),
                   axis.text.y  = element_text(size=size.text.x),
                   legend.title = element_blank(),
                   legend.text = element_text(size = size.legend.text),
                   legend.position = c(0.8, 0.9),
                   legend.key = element_rect(fill = "white")) 




stock = seq(from = 0, to = 5000, by = 50)
sons = rep(0, length(stock))
a = 20
Rp = 20000  

for (i in 1:length(stock)) {
sons[i] = a*stock[i]*exp(-a*(stock[i]/(Rp*exp(1))))
}




## create the data.frame

stock_df = data.frame(stock = stock, sons = sons)


stock_gg = ggplot(stock_df, aes(x = stock, y = sons)) +
  geom_point() +
  theme.stock +
  guides(size = guide_legend(override.aes = list(alpha = 0.2))) +
  scale_y_continuous(limits = c(0,(max(stock_df$sons) + max(stock_df$sons)/10))) +
  scale_x_continuous(limits = c(0,(max(stock_df$stock) + max(stock_df$stock)/10))) +
  labs(y = "Recruits") +
  labs(x = "Stock") 

plot(stock_gg)
