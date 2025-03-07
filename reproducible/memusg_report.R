library(tidyverse);  theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))
library(magrittr)

dd <- readRDS("memusg.rds")
print(length(dd))

set_run <- . %>% mutate(across(run, ~ forcats::fct_inorder(factor(.))))
metadata <- read.csv("reproducible/testvals.csv") |>
    mutate(run = seq(n())) |>
    set_run()

res <- map_dfr(dd, ~.[["memory_use"]], .id = "run") |>
    set_run() |>
    right_join(x = metadata, by = "run")

ggplot(res, aes(time, rss, colour = factor(ntax))) + geom_line() +
    facet_wrap(~nsubj, labeller = label_both) + zmargin
## colorspace::scale_colour_discrete_sequential()
 
ggsave("memusg_report1.pdf", width = 10, height = 10)
