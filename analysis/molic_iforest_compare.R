## CT: molic vs iforest

library(tidyverse)
library(patchwork)

theme_set(theme_bw(base_size = 16))

iforest <- readRDS("objects/pm_iforest_exact.Rds")
molic <- readRDS("objects/pm_odmgm_exact.Rds")

iforest_df <- iforest %>% as_tibble(rownames = "true_class") %>%
  pivot_longer(-true_class, names_to = "tested_class", values_to = "iForest")

molic_df <- molic %>% as_tibble(rownames = "true_class") %>%
  pivot_longer(-true_class, names_to = "tested_class", values_to = "ODMGM")

compare_df <- molic_df %>% full_join(iforest_df, by = c("true_class", "tested_class"))

plot_xy <- compare_df %>% 
  mutate(sh = ifelse(true_class == tested_class, "True case", "Outlier")) %>% 
  ggplot(aes(x = ODMGM, y = iForest, shape = sh, colour = true_class)) +
  geom_abline() + 
  geom_point(size = 2) + 
  facet_wrap(~tested_class, labeller = label_both, nrow = 2) + 
#  theme(legend.position = c(0.8, 0.1)) + 
  guides(shape = FALSE)
plot_xy

plot_li <- compare_df %>% 
  gather(method, rejection_fraction, ODMGM:iForest) %>% 
  mutate(sh = ifelse(true_class == tested_class, "True case", "Outlier")) %>% 
  ggplot(aes(x = true_class, y = rejection_fraction, colour = method, shape = sh)) +
  geom_hline(yintercept = 0.05, lty = 2) +
  geom_line(aes(group = method), alpha = 0.2) +
  geom_point(size = 2) + facet_wrap(~tested_class, labeller = label_both, nrow = 2)
plot_li

CT_delta <- 0.3
bands <- tibble(CT = 1:7, xm = CT-CT_delta, xM = CT+CT_delta, ym = -Inf, yM = Inf) %>% 
  mutate(`Cover type` = paste(CT))


plot_il <- compare_df %>% 
  # gather(method, rejection_fraction, ODMGM:iForest) %>%
  pivot_longer(cols = -c("true_class", "tested_class"), names_to = "method", values_to = "rejection_fraction") %>%
  mutate(sh = ifelse(true_class == tested_class, "True case", "Outlier")) %>% 
  rename(`Cover type` = true_class) %>% 
  ggplot(aes(x = tested_class, y = rejection_fraction, fill = method)) +
  geom_rect(data = bands, show.legend = FALSE, fill = "#999999", alpha = 0.4,
            aes(y = NULL, fill = NULL, x = NULL, ymin = ym, ymax = yM, xmin = xm, xmax = xM)) + 
  geom_hline(yintercept = 0.05, lty = 2) +
  geom_line(aes(group = method), alpha = 0.2) +
  geom_point(size = 2, pch = 21) + facet_wrap(~`Cover type`, labeller = label_both, nrow = 2) + 
  scale_fill_grey() + labs(y = "Rejection fraction", x = "Tested class", fill = "Method:") + 
  theme(legend.position = "top")
plot_il

pdf("../../paper/incl/fig/molic_iforest.pdf", width = 16, height = 8)
plot_il
# (plot_xy | plot_li + labs(shape = "type")) + plot_layout(guides = 'collect')
dev.off()
