meta.table <- read.csv(file = "metamorphosis_table.csv", header = TRUE, row.names = 1)
class(meta.table)
dim(meta.table)
colnames(meta.table) <- c("0", "1", "2", "3")

meta.table.percent <- apply(meta.table, 1, function(x) x / sum(x) * 100)
meta.table.percent <- t(meta.table.percent)

data <- meta.table.percent

library(ggplot2)

## Treatment group

df <- tibble::tribble(
  ~specimen,    ~`0`, ~`1`, ~`2`, ~`3`,
  "treat_MAC",     44,   44,   12,    0,
  "treat_IVA",  10.53, 10.53, 47.37, 31.58,
  "treat_HMA",   6.45, 12.90, 61.29, 19.35,
  "treat_HIV",  13.64, 31.82, 36.36, 18.18
)

df_treat <- tidyr::pivot_longer(df, cols = c(`0`, `1`, `2`, `3`), names_to = "category", values_to = "percentage")

df_treat$specimen <- factor(df_treat$specimen, levels = c("treat_MAC", "treat_IVA", "treat_HMA", "treat_HIV"))

x <- ggplot(df_treat) +
     geom_bar(aes(x = category, y = percentage, fill = category),
           position = "dodge", stat = "identity", width = 0.8) +
     labs(x = "Category", y = "Percentage %", fill = "Category") +
     facet_grid(~specimen)


## Control group

df.cont <- tibble::tribble(
  ~specimen,    ~`0`, ~`1`, ~`2`, ~`3`,
  "cont_MAC",    80,   20,    0,    0,
  "cont_IVA", 36.84, 26.32, 36.84,  0,
  "cont_HMA",    40,   40,   20,    0,
  "cont_HIV", 54.54, 22.73, 22.73,  0
)

df_cont <- tidyr::pivot_longer(df.cont, cols = c(`0`, `1`, `2`, `3`), names_to = "category", values_to = "percentage")

df_cont$specimen <- factor(df_cont$specimen, levels = c("cont_MAC", "cont_IVA", "cont_HMA", "cont_HIV"))

y <- ggplot(df_cont) +
     geom_bar(aes(x = category, y = percentage, fill = category),
           position = "dodge", stat = "identity", width = 0.8) +
     labs(x = "Category", y = "Percentage %", fill = "Category") +
     facet_grid(~specimen) 

gridExtra::grid.arrange(y, x, nrow = 2)

f.test <- fisher.test(meta.table, hybrid = TRUE, simulate.p.value = TRUE, B = 10000)

f.test.treatment <- fisher.test(meta.table[1:4,], hybrid = TRUE, simulate.p.value = TRUE, B = 10000)
f.test.control <- fisher.test(meta.table[5:8,], hybrid = TRUE, simulate.p.value = TRUE, B = 10000)
f.test.mac <- fisher.test(meta.table[c(1,5),], hybrid = TRUE, simulate.p.value = TRUE, B = 10000)
f.test.iva <- fisher.test(meta.table[c(2,6),], hybrid = TRUE, simulate.p.value = TRUE, B = 10000)
f.test.hma <- fisher.test(meta.table[c(3,7),], hybrid = TRUE, simulate.p.value = TRUE, B = 10000)
f.test.hiv <- fisher.test(meta.table[c(4,8),], hybrid = TRUE, simulate.p.value = TRUE, B = 10000)