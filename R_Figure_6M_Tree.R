library(dplyr)
library(tidyr)
library(ggtree)
library(ggplot2)
library(ggnewscale)

# 系統樹を読み込む
tree <- read.tree("/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For Submission/RpoB Tree for Dodd Strain Library/194 protein sequences alignment consensus tree.newick")
# シングルクォートを除去する
tree$tip.label <- gsub("'", "", tree$tip.label)

data1 <- read.csv("/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/NCAIR/BlastP/20240823_DoddLabLibrary_larABCE_NPNOCD_proteinblast.csv")
  
# データをワイド形式に変換
wide_data <- data1 %>%
  pivot_wider(names_from = gene, values_from = AA.identity)

wide_data <- wide_data %>%
  rename(tip.label = strain) 
  
missing_labels <- setdiff(tree$tip.label, wide_data$tip.label)

# missing_labelsをデータフレームに変換
missing_df <- data.frame(tip.label = missing_labels)

# wide_data2にmissing_dfを追加
# 新しい列を作成してNAで埋める
wide_data2 <- bind_rows(wide_data, missing_df)

data_phylum <- read.csv("/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/NCAIR/BlastP/Phylum_doddlibrary.csv")

wide_data2 <- wide_data2 %>%
  replace(is.na(.), 0) %>%  
  rowwise() %>%
  mutate(strain = ifelse(NPN.OCD >= 0.1, tip.label, NA)) %>%
  ungroup()%>%
  mutate(in_tree = ifelse(tip.label %in% tree$tip.label, "In Tree", "Not in Tree"))

data <- subset(wide_data2, in_tree == "In Tree") 
data <- merge(data, data_phylum, by = "tip.label") 
label_df <- data.frame(tip.label = tree$tip.label, label_number = 1:length(tree$tip.label))

# データフレームに番号をマージ（結合）
data <- data %>%
  left_join(label_df, by = "tip.label")

write.csv(data, "/Users/kazuma/Library/CloudStorage/OneDrive-Stanford/Stanford_Transposon Muntant Library Project/For R/R_studio/NCAIR/BlastP/For_paper.csv")
#### Data frame is now ready.###


p <- ggtree(tree, layout = "circular") %<+% data +
  xlim(NA, 1.5) +
  geom_tippoint(aes(color = Phylum), size = 0.5) +
  geom_tiplab(aes(label = label_number, color = Phylum), 
              size = 1.5, offset = 0.001, align = TRUE) +
  scale_color_manual(
    name = "Phylum",                       
    values = c("springgreen3", "slateblue1", "#B22222","lightsalmon","#CDAD00","dodgerblue"),
    breaks = c("Firmicutes","Verrucomicrobia", "Bacteroidetes", 
               "Actinobacteria", "Fusobacteria","Proteobacteria"),
    labels = c("Firmicutes","Verrucomicrobia", "Bacteroidetes", 
               "Actinobacteria", "Fusobacteria","Proteobacteria")) +  
  geom_tiplab(aes(label=strain) ,
              color = "red",
              offset = 0.22,
              size = 2,
              linetype = "blank",
              geom = "text",
              align = TRUE)+
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      size = 8,
      face = "bold",
      hjust = 0.5,
      vjust = -15)) 
p

larB <- data.fPhylumlarB <- data.frame("larB" = data[,c("larB")])
rownames(larB) <- data$tip.label

h1 <- p + new_scale_fill()
p2 <- gheatmap(h1, larB,  
               offset = 0.01, 
               width = 0.05,
               colnames = FALSE)+
  scale_fill_continuous(name = "AA identity",  # ここでは、MIC の連続変数のグラデーションカラースキームを定義する
                        low = "#EAF5FA", high = "#273871",
                        breaks = c(0.3, 0.65, 1.0),
                        limits = c(0.3, 1.0),
                        na.value = "white") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())

p2

larC <- data.frame("larC" = data[,c("larC")])
rownames(larC) <- data$tip.label

h2 <- p2 + new_scale_fill()
p3 <- gheatmap(h2, larC,  
               offset = 0.05, 
               width = 0.05,
               colnames = FALSE)+
  scale_fill_continuous(name = "AA identity",  # ここでは、MIC の連続変数のグラデーションカラースキームを定義する
                        low = "#EAF5FA", high = "#273871",
                        breaks = c(0.3, 0.65, 1.0),
                        limits = c(0.3, 1.0),
                        na.value = "white") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())
p3

larE <- data.frame("larE" = data[,c("larE")])
rownames(larE) <- data$tip.label

h3 <- p3 + new_scale_fill()
p4 <- gheatmap(h3, larE,  
               offset = 0.09, 
               width = 0.05,
               colnames = FALSE)+
  scale_fill_continuous(name = "AA identity",  # ここでは、MIC の連続変数のグラデーションカラースキームを定義する
                        low = "#EAF5FA", high = "#273871",
                        breaks = c(0.3, 0.65, 1.0),
                        limits = c(0.3, 1.0),
                        na.value = "white") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())

p4

larA <- data.frame("larA" = data[,c("larA")])
rownames(larA) <- data$tip.label

h4 <- p4 + new_scale_fill()
p5 <- gheatmap(h4, larA,  
               offset = 0.13, 
               width = 0.05,
               colnames = FALSE)+
  scale_fill_continuous(name = "AA identity [larA]",  # ここでは、MIC の連続変数のグラデーションカラースキームを定義する
                        low = "#EAF5FA", high = "#273871",
                        breaks = c(0.3, 0.65, 1.0),
                        limits = c(0.3, 1.0),
                        na.value = "white") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin()) 
p5

NPN <- data.frame("NPN" = data[,c("NPN.OCD")])
rownames(NPN) <- data$tip.label

h5 <- p5 + new_scale_fill()
p6 <- gheatmap(h5, NPN,  
               offset = 0.17, 
               width = 0.05,
               colnames = FALSE)+
  scale_fill_continuous(name = "AA identity [NPN-OCD]",  # ここでは、MIC の連続変数のグラデーションカラースキームを定義する
                        low = "#FFF8E1", high = "#FF6F00",
                        breaks = c(0.3, 0.65, 1.0),
                        limits = c(0.3, 1.0),
                        na.value = "white") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin()) 
p6

#####ベン図
# VennDiagram パッケージをインストール
install.packages("VennDiagram")

# 必要なパッケージを読み込む
library(VennDiagram)
# 各群を定義する
group1 <- data %>% 
  filter(larB >= 0.3 & larC >= 0.3 & larE >= 0.3) %>% 
  pull(tip.label)

group2 <- data %>% 
  filter(larA >= 0.3) %>% 
  pull(tip.label)

group3 <- data %>% 
  filter(NPN.OCD >= 0.3) %>% 
  pull(tip.label)

# ベン図を描く
venn.plot <- venn.diagram(
  x = list(Group1 = group1, Group2 = group2, Group3 = group3),
  category.names = c("larB/larC/larE >= 0.3", "larA >= 0.3", "NPN.OCD >= 0.3"),
  filename = NULL,  # ファイル保存せずにR内で表示
  fill = c("blue", "yellow", "red"),
  alpha = 0.3,
  cex = 1,
  cat.cex = 1,
  cat.col = c("blue", "yellow", "red")
)

# ベン図を表示
grid.draw(venn.plot)
