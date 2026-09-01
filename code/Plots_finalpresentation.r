# ============================================================
# MOTILITY GRAPHS
# ============================================================

library(ggplot2)

# ============================================================
# GRAPH 1: CONTROL vs CUTANEOUS MODEL
# ============================================================

data_cutaneous <- data.frame(
  Model = c("Control", "Cutaneous"),
  Motility = c(0, 3),
  No_Motility = c(7, 12),
  Total = c(7, 15)
)

# Calculate percentage
data_cutaneous$Percent_Motility <- 
  (data_cutaneous$Motility / data_cutaneous$Total) * 100

# Order groups
data_cutaneous$Model <- factor(
  data_cutaneous$Model,
  levels = c("Control", "Cutaneous")
)

# Fisher's exact test
table_cutaneous <- matrix(
  c(0, 7,
    3, 12),
  nrow = 2,
  byrow = TRUE
)

fisher_cutaneous <- fisher.test(table_cutaneous)

# Print Fisher result
print(fisher_cutaneous)


# Create Graph 1
graph_cutaneous <- ggplot(
  data_cutaneous,
  aes(x = Model, y = Percent_Motility)
) +
  
  geom_col(
    width = 0.65
  ) +
  
  geom_text(
    aes(
      label = paste0(
        round(Percent_Motility, 1),
        "%\n(",
        Motility,
        "/",
        Total,
        ")"
      )
    ),
    vjust = -0.5,
    size = 4.5
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Model",
    y = "Larval zebrafish showing motility (%)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(15, 15, 15, 15)
  )


# Show Graph 1
print(graph_cutaneous)



# ============================================================
# GRAPH 2: CONTROL vs ACRAL vs CRKL
# ============================================================

data_acral <- data.frame(
  Model = c("Control", "Acral", "CRKL"),
  Motility = c(0, 18, 8),
  No_Motility = c(17, 10, 5),
  Total = c(17, 28, 13)
)

# Calculate percentage
data_acral$Percent_Motility <- 
  (data_acral$Motility / data_acral$Total) * 100

# Order groups
data_acral$Model <- factor(
  data_acral$Model,
  levels = c("Control", "Acral", "CRKL")
)


# ------------------------------------------------------------
# Fisher's exact test: Control vs Acral
# ------------------------------------------------------------

table_acral <- matrix(
  c(0, 17,
    18, 10),
  nrow = 2,
  byrow = TRUE
)

fisher_acral <- fisher.test(table_acral)

print(fisher_acral)


# ------------------------------------------------------------
# Fisher's exact test: Control vs CRKL
# ------------------------------------------------------------

table_crkl <- matrix(
  c(0, 17,
    8, 5),
  nrow = 2,
  byrow = TRUE
)

fisher_crkl <- fisher.test(table_crkl)

print(fisher_crkl)


# Create Graph 2
graph_acral <- ggplot(
  data_acral,
  aes(x = Model, y = Percent_Motility)
) +
  
  geom_col(
    width = 0.65
  ) +
  
  geom_text(
    aes(
      label = paste0(
        round(Percent_Motility, 1),
        "%\n(",
        Motility,
        "/",
        Total,
        ")"
      )
    ),
    vjust = -0.5,
    size = 4.5
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Model",
    y = "Larval zebrafish showing motility (%)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(15, 15, 15, 15)
  )


# Show Graph 2
print(graph_acral)



# ============================================================
# SAVE BOTH GRAPHS
# ============================================================

ggsave(
  "Cutaneous_Motility.png",
  graph_cutaneous,
  width = 5,
  height = 5,
  dpi = 600,
  bg = "transparent"
)

ggsave(
  "Acral_CRKL_Motility.png",
  graph_acral,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "transparent"
)
