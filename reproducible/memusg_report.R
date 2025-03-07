library(tidyverse);  theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))
library(magrittr)

pdf("memusg.pdf")

dd <- readRDS("memusg.rds")
print(length(dd))

set_run <- . %>% mutate(across(run, ~ forcats::fct_inorder(factor(.))))
metadata <- read.csv("reproducible/testvals.csv") |>
    mutate(run = seq(n())) |>
    set_run()

res <- map_dfr(dd, ~.[["memory_use"]], .id = "run") |>
    set_run() |>
    right_join(x = metadata, by = "run") |>
    mutate(rss_gb = rss/1e9)

print(ggplot(res, aes(time, rss_gb, colour = factor(ntax))) + geom_line() +
      facet_wrap(~nsubj, labeller = label_both) + zmargin +
      labs(title= "memory use trace (Gb)")
      )
## colorspace::scale_colour_discrete_sequential()



## find max (? quantile ?)
sum_res <- (res
    |> summarise(rss_max_gb = max(rss_gb), .by = c("ntax", "nsubj"))
)

print(
    ggplot(sum_res, aes(ntax, rss_max_gb, colour = factor(nsubj))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x) +
    labs(title= "peak memory (Gb) vs ntaxa")
)

print(
    ggplot(sum_res, aes(nsubj, rss_max_gb, colour = factor(ntax))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x) +
    labs(title= "peak memory (Gb) vs nsubj")
)

print(
    lm(log(rss_max_gb) ~ log(ntax)+ factor(nsubj), data = sum_res)
)
## scaling 1.83

print(
    lm(log(rss_max_gb) ~ factor(ntax)+ log(nsubj), data = sum_res)
)

## scaling 0.92
print(
    lm(log(rss_max_gb) ~ factor(ntax) + factor(ntax):log(nsubj), data = sum_res)
)
## slightly larger: 0.79, 0.94, 0.96, ..., 0.92

full <- lm(log(rss_max_gb) ~ log(ntax)*log(nsubj), data = sum_res)
print(exp(predict(full, newdata = list(ntax=969, nsubj = 8))))  ## 588 Gb

## could do another scaling to predict time ...

res_time <- map_dfr(dd, ~.[["result"]], .id = "run") |>
    set_run()

print(
    ## nsubj == 5 is a little wonky (too short to be reliable for model fit)
    ggplot(filter(res_time, nsubj > 5),
           aes(nsubj, time_sec, colour = factor(ntax))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x) +
    facet_wrap(~task) +
    labs(title= "time (sec) vs nsubj")
)

print(
    ggplot(filter(res_time, nsubj>5),
           aes(ntax, time_sec, colour = factor(nsubj))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x) +
    facet_wrap(~task) +
    labs(title= "time (sec) vs ntax")
)

full_time <- lm(log(time_sec) ~ log(ntax)*log(nsubj), data = res_time,
                subset = (task == "leverage"))
print(exp(predict(full_time, newdata = list(ntax=969, nsubj = 8))))  ## 11000 seconds
## = 3+hours

dev.off()
