library(renv)
renv::activate()

# Restore the project environment
renv::restore()

.libPaths()
# [1] "/Users/kr/Desktop/demo/conda_renv/test2/renv/library/R-4.1/x86_64-apple-darwin17.0"            
# [2] "/Users/kr/Library/Caches/org.R-project.R/R/renv/sandbox/R-4.1/x86_64-apple-darwin17.0/2af913d7"

library(tidyverse)

renv::install('tidyverse')

library(tidyverse)

df <- iris

df %>% 
  ggplot(., aes(Sepal.Length, Sepal.Width, color = Species)) +
  geom_point()

sessionInfo()

renv::snapshot()






