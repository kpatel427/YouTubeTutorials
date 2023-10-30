
.libPaths()
# [1] "/Users/kr/Library/R/x86_64/4.1/library"                        
# [2] "/Library/Frameworks/R.framework/Versions/4.1/Resources/library"

library(tidyverse)

df <- iris

df %>% 
  ggplot(., aes(Sepal.Length, Sepal.Width, color = Species)) +
  geom_point()

sessionInfo()
