library("data.table")
library("ggplot2")
library("dplyr")
library("tidyr")
library("forcats")

#### Read in data ####

# Picrust IMG genome annotations
picrust <- fread("ko_13_5_precalculated_nometadata.tab",
                      header = TRUE, sep = "\t",
                      colClasses = list(factor = 1, integer = 2:6910, numeric = 6911))
picrust <- picrust %>%
  filter(metadata_NSTI <= 0.10789) # trying to only select the genomes that were in IMG

picrust_taxonomy <- fread("97_otu_taxonomy.txt",header = FALSE, sep = "\t")
# KO and Module Descriptions
KO_desc <- fread("ko_13_5_precalculated_kegg_descriptors_colcorr.tab",
                      header = TRUE, sep = "\t")
M_desc <- fread("ko_13_5_precalculated_keggM_descriptors_colcorr.tab",
                      header = TRUE, sep = "\t")

# My Bin annotations
bins_ko <- fread("ko_untargets_all_bins_r.txt",
                      header = TRUE, sep = "\t")

#### Clean and prepare data ####

# Prepare KO product annotations
KO_desc_t <- data.frame(t(KO_desc[,-1]))
names(KO_desc_t) <- c("Product")
KO_desc_t <- KO_desc_t %>%
   add_rownames(var = "KO")

# Prepare KO module annotations
M_desc_t <- data.frame(t(M_desc[,-1]))
names(M_desc_t) <- c("Module")
M_desc_t <- M_desc_t %>%
   separate(col = Module, into = paste("Module",1:19,sep="."), sep = "\\|") %>%
   add_rownames(var = "KO") %>%
   mutate_each(funs(gsub("Unclassified;", "", .)), -KO) %>%
   gather(Colname, Module, -KO) %>%
   separate(col = Module, into = paste("Level",1:3,sep="."), sep = ";") %>%
   filter(!is.na(Level.1))# %>%
   #select(-Colname)
   
# Join product and module annotations
KO_annot <- M_desc_t %>%
   inner_join(KO_desc_t, by = "KO")

# Prepare picrust annotations
picrust_pa <- picrust %>%
   select(-metadata_NSTI) %>%
   mutate_each(funs(replace(., . > 0, 1)), -OTU_IDs) %>%
   mutate(bin_gg = "GG")
  

picrust_pa_melt <- picrust_pa %>%
   gather(KO, presence, -OTU_IDs, -bin_gg) %>%
   mutate(KO = as.factor(KO), presence = as.integer(presence),
          OTU_IDs = as.factor(OTU_IDs))

# Prepare picrust taxonomy
picrust_tax <- picrust_taxonomy %>%
   rename(OTU_IDs = V1, Taxonomy = V2) %>%
   mutate(OTU_IDs = as.factor(OTU_IDs))

picrust_pa_melt_tax <- picrust_pa_melt %>%
   inner_join(picrust_tax)

# Prepare Bins annotations
bins_ko_compact <- bins_ko %>% group_by(bin, KO) %>%
   tally() %>%
   mutate(presence = replace(n, n > 0, 1)) %>%
   select(-n) %>%
   mutate(bin_gg = "Bin")

# Join Bins and picrust annotations
merge()
### There are 376 KO assignments in my genoems that are not in the database. 
### Need to figure out what to do about it...
bins_picrust_ko <- bins_ko_compact %>%
   full_join(picrust_pa_melt, by = (c("bin" = "OTU_IDs", "KO",
                                      "presence", "bin_gg"))) %>%
   #mutate(bin_gg = ifelse(grepl("bin", bin), "Bin", "GG"))

#### Calculate Rarity of KOs ####
picrust_pa_melt_tally <- picrust_pa_melt %>% 
   group_by(KO, presence) %>%
   summarise(n = n()) %>%
   mutate(freq = n / sum(n)) %>%
   filter(presence == 1)
picrust_rare <- picrust_pa_melt_tally %>%
   filter(freq < 0.05) %>%
   mutate(Rare = "Rare")

#### Combine Bin and rare picrust data ####
# Joining picrust counts and bin counts
bins_picrust_ko_otu <-  picrust_pa_melt_tax %>%
   filter(grepl("g__Nostoc", Taxonomy)) %>%
   inner_join(bins_picrust_ko)


rare_bins_ko <- bins_picrust_ko %>% 
   full_join(picrust_rare, by = c("KO",
                                  "presence")) %>%
   filter(Rare == "Rare") %>%
   filter(bin_gg == "Bin") %>%
   select(bin, KO, bin_gg, freq, Rare)

rare_by_bins <- rare_bins_ko %>%
   mutate(presence = 1) %>%
   spread(bin, presence, fill = 0)

#### Combine KO annotations and Rare picrust data ####
rare_ko_annot <- picrust_rare %>%
   select(KO, freq) %>%
   inner_join(KO_annot, by = c("KO"))

# Join KO annotations to rare bins
# I've calcuated that all of the rare kos are contained in teh set of KO_annot
# so an inner_join is appropriate
rare_by_bins_annot_nolevels <- rare_by_bins %>%
   inner_join(rare_ko_annot, by = c("KO", "freq")) %>%
   unite(KeggModules, Level.1, Level.2, Level.3, sep = ";") %>%
   spread(Colname, KeggModules)

rare_by_bins_annot <- rare_by_bins %>%
   inner_join(rare_ko_annot, by = c("KO", "freq"))

#### Plot Results ####

# Ordered Heatmap
plot_rare_df <- as.matrix(rare_by_bins[c(7,12,9)])
rownames(plot_rare_df) <- rare_by_bins$KO
row.order <- hclust(dist(plot_rare_df))$order
col.order <- hclust(dist(t(plot_rare_df)))$order
plot_rare_df_ord <- data.frame(plot_rare_df[row.order, col.order])
plot_rare_df_ord$KO <- rownames(plot_rare_df_ord)

rare_ord_melt <- plot_rare_df_ord %>%
   gather(bin, value, -KO) %>%
   inner_join(KO_annot, by = "KO") %>%
   mutate(Level.1 = replace(Level.1, value == 0, NA),
          Level.2 = replace(Level.2, value == 0, NA),
          Level.3 = replace(Level.3, value == 0, NA))  %>%
   mutate(Level.1 = factor(Level.1), Level.2 = factor(Level.2), 
          Level.3 = factor(Level.3)) %>%
   filter(Level.2 == "Energy Metabolism")

ggplot(data = rare_ord_melt,
       aes(x = bin, y = fct_reorder(KO, Level.3), fill = Level.3)) + 
   geom_raster() + geom_label(aes(label = Level.3))
   scale_fill_discrete(na.value = "black") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         axis.text.y = element_blank()) +
   ggtitle("Clustered Rare KOs")

### Ordered by KEGG orthology
plot_rare_df_ord <- rare_by_bins_annot %>%
   select(KO, bin.029, bin.033, bin.044, bin.061, Level.1, Level.2, Level.3, Product) %>%
   arrange(Level.3, Product)

porph_chlorophyll <- plot_rare_df_ord %>% 
   filter(Level.3 == "Porphyrin and chlorophyll metabolism") %>%
   select(KO) %>%
   inner_join(plot_rare_df_ord)


rare_ord_melt <- plot_rare_df_ord %>%
   gather(bin, value, -KO, -Level.1, -Level.2, -Level.3, -Product) %>%
   mutate(Level.1 = replace(Level.1, value == 0, NA),
          Level.2 = replace(Level.2, value == 0, NA),
          Level.3 = replace(Level.3, value == 0, NA))  %>%
   mutate(Level.1 = factor(Level.1), Level.2 = factor(Level.2), 
          Level.3 = factor(Level.3)) %>%
   arrange(Level.3, KO) %>%
   add_rownames(var = "order") %>%
   mutate(order = as.numeric(order)) %>%
   filter(Level.2 == "Metabolism of Cofactors and Vitamins")

ggplot(data = rare_ord_melt,
       aes(x = bin, y = fct_reorder(KO, order), fill = Level.3)) + 
   geom_raster() + 
   geom_text(aes(label = Level.3)) +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         legend.position = "none") +
   ggtitle("Clustered Rare KOs")

