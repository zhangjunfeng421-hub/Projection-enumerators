D_1 <- c(1,2,3,0,7,4,5,6,
         1,3,7,5,2,0,4,6,
         0,6,2,4,3,5,1,7)
D_1 <- matrix(D_1,nrow = 3,byrow = T)
D_1 <- t(D_1)
D_2 <- matrix(c(0.5,1.5,2.5,3.5,-3.5,-2.5,-1.5,-0.5,
                #-1.5,0.5,-3.5,2.5,-2.5,3.5,-0.5,1.5,
                3.5,2.5,-1.5,-0.5,0.5,1.5,-2.5,-3.5,
                2.5,-3.5,-0.5,1.5,-1.5,0.5,3.5,-2.5),nrow = 3,byrow = TRUE)
D_2 <- t(D_2)+3.5

points <- D_1
colnames(points) <- c("x", "y", "z")
library(ggplot2)
library(gridExtra)

df <- as.data.frame(points) 

vlines <- 0:7   
hlines <- 0:7  


my_theme <- theme_minimal() +
  theme(
    # plot.background  = element_rect(fill = "white", color = "black", size = 1),  
    panel.background  = element_rect(fill = "white", color = "gray90"),        
    panel.grid.major  = element_line(color = "white", linetype = "dotted"),    
    panel.grid.minor  = element_blank(),                                       
    plot.title  = element_text(hjust = 0.5, face = "bold")                      
  )


p1 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = c(-0.5,3.5,7.5), linetype = 2, color = "black") +
  geom_hline(yintercept = c(-0.5,1.5,3.5,5.5,7.5), linetype = 2, color = "black") +
  labs(x = "Column 1", y = "Column 2") +
  xlim(-0.5, 7.5) + ylim(-0.5, 7.5) +
  coord_fixed(ratio = 1) +
  my_theme 



points <- D_2
colnames(points) <- c("x", "y", "z")
library(ggplot2)
library(gridExtra)

df <- as.data.frame(points) 


my_theme <- theme_minimal() +
  theme(
    # plot.background  = element_rect(fill = "white", color = "black", size = 1),  
    panel.background  = element_rect(fill = "white", color = "gray90"),        
    panel.grid.major  = element_line(color = "white", linetype = "dotted"),    
    panel.grid.minor  = element_blank(),                                       
    plot.title  = element_text(hjust = 0.5, face = "bold")                      
  )


p2 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = c(-0.5,3.5,7.5), linetype = 2, color = "black") +
  geom_hline(yintercept =  c(-0.5,1.5,3.5,5.5,7.5), linetype = 2, color = "black") +
  labs(x = "Column 1", y = "Column 2") +
  xlim(-0.5, 7.5) + ylim(-0.5, 7.5) +
  coord_fixed(ratio = 1) +
  my_theme 
  


P <- grid.arrange( 
  p1, p2,
  ncol = 2,
  padding = unit(0.5, "line") 
)

ggsave("bac_GSOA.png", P,
       width = 8, height = 4, dpi=300, bg="white")
