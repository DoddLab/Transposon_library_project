# 必要なパッケージを読み込む
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(ggpubr)

# データの読み込み
df_2ndscreening_B1 <- read.csv('/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/Second_screening/20250207_Final Analysis/20250207_TL2ndscreening_B1.csv')
df_2ndscreening_B2 <- read.csv('/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/Second_screening/20250207_Final Analysis/20250207_TL2ndscreening_B2.csv')
df_2ndscreening_B3 <- read.csv('/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/Second_screening/20250207_Final Analysis/20250207_TL2ndscreening_B3.csv')
df_2ndscreening_B4 <- read.csv('/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/Second_screening/20250207_Final Analysis/20250207_TL2ndscreening_B4.csv')
df_2ndscreening_B5 <- read.csv('/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/Second_screening/20250206_Additional/20250206_additional 2nd screening.csv')

df_2ndscreening_B1_long <- reshape2::melt(df_2ndscreening_B1, 
                                                  id.vars = c("well", "type", "TLnumber","batch", "gene", "insertion.point"), 
                                                  variable.name = "compound", value.name = "value") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%　# NAを0に変換
  filter(!grepl("Carnitine|t.Cinnamic.Acid|Asn|IAA|PAA|Phenylethylamine|PPA|PPA\\.D6|Propionic\\.acid|Isocaproic\\.acid", compound))

df_2ndscreening_B2_long <- reshape2::melt(df_2ndscreening_B2, 
                                          id.vars = c("well", "type", "TLnumber","batch", "gene", "insertion.point"), 
                                          variable.name = "compound", value.name = "value")%>%
  mutate_all(~replace(., is.na(.), 0)) %>%　# NAを0に変換
  filter(!grepl("Carnitine|t.Cinnamic.Acid|Asn|IAA|PAA|Phenylethylamine|PPA|PPA\\.D6|Propionic\\.acid|Isocaproic\\.acid", compound))

df_2ndscreening_B3_long <- reshape2::melt(df_2ndscreening_B3, 
                                          id.vars = c("well", "type", "TLnumber","batch", "gene", "insertion.point"), 
                                          variable.name = "compound", value.name = "value") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%　# NAを0に変換
  filter(!grepl("Carnitine|t.Cinnamic.Acid|Asn|IAA|PAA|Phenylethylamine|PPA|PPA\\.D6|Propionic\\.acid|Isocaproic\\.acid", compound))

df_2ndscreening_B4_long <- reshape2::melt(df_2ndscreening_B4, 
                                          id.vars = c("well", "type", "TLnumber","batch", "gene", "insertion.point"), 
                                          variable.name = "compound", value.name = "value") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%　# NAを0に変換
  filter(!grepl("Carnitine|t.Cinnamic.Acid|Asn|IAA|PAA|Phenylethylamine|PPA|PPA\\.D6|Propionic\\.acid|Isocaproic\\.acid", compound))

df_2ndscreening_B5_long <- reshape2::melt(df_2ndscreening_B5, 
                                          id.vars = c("well", "type", "TLnumber","batch", "gene", "insertion.point"), 
                                          variable.name = "compound", value.name = "value") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%　# NAを0に変換
  filter(!grepl("Carnitine|t.Cinnamic.Acid|Asn|IAA|PAA|Phenylethylamine|PPA|PPA\\.D6|Propionic\\.acid|Isocaproic\\.acid", compound))

substrates <- c("Choline", "Asn", "Arg", "Citruline", "Ornithine",
                "X4.OH.Proline", "Pro", "Ser","Thr", "Met", 
                "Tyr", "Ile","Leu","Uric.acid", "Trp", "Val", "Phe.D8")

products <- c("X5.Aminovaleric.acid","IPA", "X2.aminobutyric.acid",
              "Isovaleric.acid", 
              "Isobutyric.acid", "Tryptamine", 
              "IAA", "PAA","Isocaproic.acid", "PPA","PPA.D6",
              "Propionic.acid", "Phenylethylamine")

calculate_volcano_data <- function(df, compound_name) {
  df <- df %>% mutate(compound = as.character(compound))  # Ensure compound is a character
  df_compound <- df %>% filter(compound == compound_name)  # Now filtering works
  
  # WT の値を取得
  WT_values <- df_compound$value[df_compound$gene == "WT"]
  WT_mean <- mean(WT_values, na.rm = TRUE)
  
  SACC_values <- df_compound$value[df_compound$gene == "SACC"]
  SACC_mean <- mean(SACC_values, na.rm = TRUE)
  
  # Determine scaling factor: 1 for substrates, -1 for products
  scaling_factor <- ifelse(compound_name %in% substrates, 1, -1)
  
  # 各変異体との log2FC と P 値を計算 (keeping TLnumber & insertion.point)
  df_results <- df_compound %>%
    filter(gene != "WT") %>%
    group_by(gene, TLnumber, insertion.point) %>%  # Include insertion.point
    summarise(
      Compound = compound_name,  
      Mutant_mean = mean(value, na.rm = TRUE),
      log2FC = log2(Mutant_mean / WT_mean),
      AbsoluteChange = (Mutant_mean - WT_mean),
      RelativeChange = scaling_factor * (Mutant_mean - WT_mean) / (SACC_mean - WT_mean),  # Apply scaling factor
      p_value = if (length(WT_values) > 1 & length(value) > 1) {
        t.test(WT_values, value)$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      neg_log10_p = -log10(p_value),
      Significance = case_when(
        !is.na(p_value) & p_value < 0.05 & abs(log2FC) > 1 ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )
  
  return(df_results)  # Return data including TLnumber and insertion.point
}

# すべての `compound` に対して処理し、結果をリストに格納
unique_compounds <- unique(df_2ndscreening_B1_long$compound)

df_results_list_B1<- lapply(unique_compounds, function(compound_name) {
  calculate_volcano_data(df_2ndscreening_B1_long, compound_name)
})

df_results_list_B2<- lapply(unique_compounds, function(compound_name) {
  calculate_volcano_data(df_2ndscreening_B2_long, compound_name)
})

df_results_list_B3<- lapply(unique_compounds, function(compound_name) {
  calculate_volcano_data(df_2ndscreening_B3_long, compound_name)
})

df_results_list_B4<- lapply(unique_compounds, function(compound_name) {
  calculate_volcano_data(df_2ndscreening_B4_long, compound_name)
})

df_results_list_B5<- lapply(unique_compounds, function(compound_name) {
  calculate_volcano_data(df_2ndscreening_B5_long, compound_name)
})

# すべての `df_results` を1つのデータセットに統合
df_results_all <- bind_rows(df_results_list_B1, df_results_list_B2, 
                            df_results_list_B3, df_results_list_B4, df_results_list_B5)

df_results_all$Compound <- factor(df_results_all$Compound,
                                levels = c("Choline", "Asn", "Arg", "Citruline", "Ornithine",
                                           "X4.OH.Proline", "Pro", "Ser","Thr", "Met", 
                                           "Tyr", "Ile","Leu","Uric.acid", "Trp", "Val", "Phe.D8",
                                           "X5.Aminovaleric.acid","IPA", "X2.aminobutyric.acid",
                                           "Isovaleric.acid", 
                                           "Isobutyric.acid", "Tryptamine", 
                                           "IAA", "PAA","Isocaproic.acid", "PPA","PPA.D6",
                                           "Propionic.acid", "Phenylethylamine"))

# NaNのp_valueを1に変換
df_results_all <- df_results_all %>%
  mutate(
    p_value = ifelse(is.na(p_value) | is.nan(p_value), 1, p_value),  # Replace NaN p-values with 1
    log2FC = case_when(
      is.infinite(log2FC) & log2FC > 0 ~ 10,  # Replace Inf with 10
      is.infinite(log2FC) & log2FC < 0 ~ -10, # Replace -Inf with -10
      TRUE ~ log2FC  # Keep other values unchanged
    )
  ) %>%
  mutate_all(~replace(., is.na(.), 0)) 

# TLnumber list to filter
valid_TLnumbers <- c(
  "TL046-A04", "TLExp031-C06", "TL022-G06", "TL034-E04", "TL054-E10",
  "TLExp003-E08", "TLExp001-F12", "TLExp010-G04", "TLExp034-F12", "TL041-E11",
  "TLExp035-H10", "TLExp038-H12", "TLExp013-H10", "TLExp036-A10", "TLExp038-D11",
  "TL026-D07", "TL067-A03", "TLExp041-A11", "TL061-H11", "TLExp017-C07",
  "TL014-D04", "TLExp040-C02", "TL054-C05", "TL022-H02", "TLExp037-D01",
  "TL019-A02", "TL001-E02", "TL022-A12", "TL038-E01", "TLExp041-B12",
  "TL039-H06", "TL044-B09", "TLExp040-D01", "TL019-F09", "TL021-A06",
  "TLExp008-H12", "TL017-D04", "TL069-D10", "TL068-G08", "TL014-D11",
  "TL031-B09", "TLExp041-E08", "TL004-B01", "TL005-A09", "TL041-C06",
  "TL012-B04", "TLExp008-A04", "TL018-F07", "TL037-F09", "TL063-H10",
  "TLExp042-E05", "TLExp011-H10", "TL058-C04", "TLExp003-A03", "TLExp001-G01",
  "TL034-F12", "TLExp003-F08", "TL029-B11", "TL069-D06", "TL067-E06",
  "TL001-A04", "TLExp033-B02", "TL014-D12", "TLExp040-E09", "TL044-D09",
  "TL023-F04", "TL001-A12", "TL001-H01", "TL002-A01", "TL002-A02",
  "TL002-A03", "TL002-A05", "TL002-A07", "TL002-A12", "TL004-A12",
  "TL005-E05", "TL010-B03", "TL014-C02", "TL016-F09", "TL019-G07",
  "TL019-G08", "TL019-H12", "TL021-A12", "TL022-H06", "TL024-E12"
)

# List of TLnumbers to remove
remove_TLnumbers <- c(
  "TLExp001-F12", "TLExp010-G04", "TLExp034-F12", "TLExp035-H10", "TLExp013-H10",
  "TLExp036-A10", "TLExp038-D11", "TLExp017-C07", "TL022-A12", "TL012-B04",
  "TL018-F07", "TLExp011-H10", "TL023-F04", "TL001-A12", "TL001-H01",
  "TL002-A01", "TL002-A02", "TL002-A03", "TL002-A05", "TL002-A07",
  "TL002-A12", "TL004-A12", "TL001-E02", "TL001-A04", "TL019-G08", "TL019-G07",
  "TL010-B03", "TL019-H12", "TL021-A12", "TL022-H06", "TL014-C02", "TLExp031-C06",
  "TL069-D10", "TL046-A04"
)

df_filtered <- df_results_all %>%
  filter(TLnumber %in% valid_TLnumbers) %>%
  filter(!TLnumber %in% remove_TLnumbers) %>%
  mutate(
    insertion.point = as.character(insertion.point),  # Ensure it's a character
    insertion.point = gsub(",", "", insertion.point),  # Remove commas
    insertion.point = as.numeric(insertion.point)  # Convert to numeric
  ) %>%
  arrange(desc(insertion.point)) %>%
  mutate(TL_label = paste0(TLnumber, " (", gene, "_", insertion.point, ")"))  # Format: TLnumber (gene_insertion.point)

df_filtered$TL_label <- factor(df_filtered$TL_label, levels = unique(df_filtered$TL_label))

write.csv(df_filtered, "/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/Second_screening/20250206_Additional/df_filtered 2nd screening.csv", 
          row.names = FALSE)

# Create the plot with combined labels on Y-axis
ggplot(df_filtered, aes(x = TL_label, y = Compound, color = RelativeChange, 
                        size = ifelse(-log10(p_value) >= 1.3, -log10(p_value), NA))) +  
  geom_point() +  
  scale_color_gradientn(
    colors = c("#0000FF", "#007FFF", "#FFFFFF", "#F67C6D", "#D54046"),  # Ensure 0 is white
    values = scales::rescale(c(min(df_filtered$RelativeChange, na.rm = TRUE), 
                               -2.5, -1, 0, 1, 2, 
                               max(df_filtered$RelativeChange, na.rm = TRUE)))  # Ensures proper positioning
  ) +  
  scale_size_continuous(
    range = c(1, 6),  
    breaks = c(1.3, 2, 3, 4),  
    labels = c("p = 0.05", "p = 0.01", "p = 0.001", "p = 0.0001")  
  ) +
  labs(title = "Filtered Heatmap (2nd screening)",
       x = "TLnumber (Gene_insertion point)",
       y = "",
       size = "-log10(p-value)",  
       color = "Adjusted change") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, face = "bold"),  # Move X labels to the top
    axis.text.y = element_text(hjust = 0, vjust = 0, face = "bold"),  # Move Y labels to the right
    axis.ticks.y.right = element_line(),  # Ensure ticks appear on the right
    axis.text.y.right = element_text(margin = margin(l = 10)),  # Add spacing for right labels
    axis.title.x = element_text(vjust = 1, hjust = 0.5),  # Move X-axis title up
    axis.title.y = element_blank(),  # Remove default Y-axis title
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center title
  ) +
  scale_y_discrete(position = "right") +  # Move Y labels to the right
  scale_x_discrete(position = "top", limits = rev(unique(df_filtered$TL_label)))  # Reverse X-axis order and keep labels on top

ggplot(df_filtered, aes(x = insertion.point, y = 0)) +
  geom_segment(aes(x = insertion.point, xend = insertion.point, y = -0.1, yend = 0.1), 
               color = "black") +  # Short vertical lines for insertion points
  geom_text(aes(label = TL_label), vjust = -1, angle = 90, size = 3) +  # Labels above points
  scale_x_continuous(limits = c(1, 4141984), breaks = seq(0, 4141984, by = 500000)) +  # Genome range
  scale_y_continuous(limits = c(-0.2, 0.2)) +  # Keep y-axis minimal
  labs(title = "Genome Map of Insertions",
       x = "Genomic Position",
       y = "") +  # No y-axis label (purely positional)
  theme_pubr() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  # Hide y-axis ticks & labels


df_unique <- df_filtered %>%
  distinct(insertion.point, .keep_all = TRUE) %>%# Keep only unique insertion points
  mutate(
    insertion.point = as.numeric(insertion.point)  # Convert to numeric
  ) %>%
  arrange(insertion.point) %>%  # Sort in ascending order
  mutate(heat = seq(0, by = 83000, length.out = n()))  # Generate heat values
  
ggplot(df_unique) +
  # Draw lines connecting (insertion.point, 0) to (heat, 1)
  geom_segment(aes(x = insertion.point, y = 0, xend = heat, yend = 1), color = "gray50") +
  # Plot insertion points at y = 0
  geom_point(aes(x = insertion.point, y = 0), pch = 6, size = 2, color = "black") +
  # Plot heat points at y = 1
  geom_point(aes(x = heat, y = 1), shape = 16, size = 2, color = "blue") +
  # X-axis settings
  scale_x_continuous(limits = c(0, 4150000), 
                     breaks = c(seq(0, 4200000, by = 500000), 4141984)
                     ) +# Add 4141984 explicitly
  # Remove y-axis labels for clarity
  theme_pubr() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  # Title and labels
  labs(title = "Genome Map of Unique Insertions with Heat Mapping",
       x = "Genomic Position",
       y = "")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
