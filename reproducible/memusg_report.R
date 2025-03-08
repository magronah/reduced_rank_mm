library(tidyverse);  theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))
library(magrittr)

pdf("memusg.pdf")

get_data <- function(metadata_fn = "testvals.csv", data_fn = "memusg.rds", dir = "reproducible",
                     by_vars = c("ntax", "nsubj")) {
    
    set_run <- . %>% mutate(across(run, ~ forcats::fct_inorder(factor(.))))

    dd <- readRDS(file.path(dir, data_fn))
    metadata <- read.csv(file.path(dir, metadata_fn)) |>
        mutate(run = seq(n())) |>
        set_run()

    res <- map_dfr(dd, ~.[["memory_use"]], .id = "run") |>
        set_run() |>
        right_join(x = metadata, by = "run") |>
        mutate(rss_gb = rss/1e9)

    res_time <- map_dfr(dd, ~.[["result"]], .id = "run") |>
        set_run()

    sum_res <- (res
        |> summarise(rss_max_gb = max(rss_gb), .by = all_of(by_vars))
    )
    
    return(list(mem_trace = res, mem = sum_res, time = res_time))
}

res <- get_data()

res2 <- get_data("testvals2.csv", data_fn = "memusg2.rds",
                 by_vars = c("ntax", "nsubj", "d", "include_ttt")) |>
print(nrow(res2$mem))

print(ggplot(res2$mem_trace, aes(time, rss_gb, colour = factor(ntax), linetype = include_ttt)) +
      geom_line() +
      facet_grid(d~nsubj, labeller = label_both) + zmargin +
      labs(title= "memory use trace (Gb)") +
      scale_x_continuous(limits = c(NA, 300), oob = scales::squish)
      )

## did I guess right on increasing order? no ...
plot(1:nrow(res2$mem), res2$mem$rss_max_gb)

print(
    ggplot(res2$mem, aes(ntax, rss_max_gb, colour = factor(nsubj), linetype = factor(d))) +
    facet_wrap(~include_ttt, labeller = label_both) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x, alpha = 0.1) +
    labs(title= "peak memory (Gb) vs ntaxa")
)



print(ggplot(res$mem_trace, aes(time, rss_gb, colour = factor(ntax))) + geom_line() +
      facet_wrap(~nsubj, labeller = label_both) + zmargin +
      labs(title= "memory use trace (Gb)")
      )
## colorspace::scale_colour_discrete_sequential()


print(
    ggplot(res$mem, aes(ntax, rss_max_gb, colour = factor(nsubj))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x) +
    labs(title= "peak memory (Gb) vs ntaxa")
)


print(
    ggplot(res$mem, aes(nsubj, rss_max_gb, colour = factor(ntax))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x) +
    labs(title= "peak memory (Gb) vs nsubj")
)

print(
    lm(log(rss_max_gb) ~ log(ntax)+ factor(nsubj), data = res$mem)
)
## scaling 1.83

print(
    lm(log(rss_max_gb) ~ factor(ntax)+ log(nsubj), data = res$mem)
)

## scaling 0.92
print(
    lm(log(rss_max_gb) ~ factor(ntax) + factor(ntax):log(nsubj), data = res$mem)
)
## slightly larger: 0.79, 0.94, 0.96, ..., 0.92

full <- lm(log(rss_max_gb) ~ log(ntax)*log(nsubj), data = res$mem)
print(exp(predict(full, newdata = list(ntax=969, nsubj = 8))))  ## 588 Gb

## could do another scaling to predict time ...


print(
    ## nsubj == 5 is a little wonky (too short to be reliable for model fit)
    ggplot(filter(res$time, nsubj > 5),
           aes(nsubj, time_sec, colour = factor(ntax))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x) +
    facet_wrap(~task) +
    labs(title= "time (sec) vs nsubj")
)

print(
    ggplot(filter(res$time, nsubj>5),
           aes(ntax, time_sec, colour = factor(nsubj))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", formula = y ~ x) +
    facet_wrap(~task) +
    labs(title= "time (sec) vs ntax")
)

full_time <- lm(log(time_sec) ~ log(ntax)*log(nsubj), data = res$time,
                subset = (task == "leverage"))
print(exp(predict(full_time, newdata = list(ntax=969, nsubj = 8))))  ## 11000 seconds
## = 3+hours

dev.off()
