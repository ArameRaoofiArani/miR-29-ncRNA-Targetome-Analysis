library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

cor_hacat <- readRDS("cor_hacat.rds")
cor_tigk <- readRDS("cor_tigk.rds")

make_df <- function(cor_mat){
  as.data.frame(cor_mat) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to="lncRNA", values_to="cor") %>%
    mutate(abs_cor=abs(cor))
}

df_hacat <- make_df(cor_hacat)
df_tigk  <- make_df(cor_tigk)

plot_fun <- function(df, title){
  ggplot(df, aes(lncRNA, gene, size=abs_cor, color=cor)) +
    geom_point() +
    scale_color_gradient2(low="blue", mid="white", high="red") +
    theme_classic() +
    labs(title=title)
}

p1 <- plot_fun(df_hacat,"HaCaT")
p2 <- plot_fun(df_tigk,"TIGK")

ggsave("bubble_hacat.pdf", p1)
ggsave("bubble_tigk.pdf", p2)
