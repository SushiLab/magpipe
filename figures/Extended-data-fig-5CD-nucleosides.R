# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 5CD - Nucleosides BGCs ==============================================

go_micro_antismash_file = "data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-antismash_summary.tsv.gz"
scaffolds_length = read_tsv("data/raw/go_microbiomics/summaries/genomes_bgcs-scaffolds_length.tsv")

genomic_bgc_table = read_tsv(go_micro_antismash_file) %>%
  left_join(scaffolds_length %>% dplyr::rename(scaffold_length = length)) %>%
  filter(scaffold_length >= 5000) %>%
  filter(!genome %in% blacklisted_genomes) %>%
  select(-scaffold_length) %>%
  select(-contains("CDS"))

genomic = read_tsv("https://microbiomics.io/ocean/db/1.0/MGEs-predictions/go_microbiomics-genomic_bgc_scaffolds-with-MGEs.tsv") %>%
  filter(scaffold_name %in% genomic_bgc_table$scaffold) %>%
  left_join(genomic_bgc_table, by = c("scaffold_name" = "scaffold")) %>%
  select(scaffold = scaffold_name, products, region, eukrep, plasflow_prediction, cbar_prediction, plasmidfinder_plasmid, deepvirfinder_score, deepvirfinder_pvalue, virsorter_score, virsorter_category) %>%
  mutate(Genomic = TRUE)

metagenomic = read_tsv("https://microbiomics.io/ocean/db/1.0/MGEs-predictions/go_microbiomics-metagenomic_bgcs-with-MGEs.tsv") %>%
  select(scaffold, region, products, eukrep, plasflow_prediction, cbar_prediction, plasmidfinder_plasmid, deepvirfinder_score, deepvirfinder_pvalue, virsorter_score, virsorter_category) %>%
  mutate(Genomic = FALSE)

table_summary = metagenomic %>%
  rbind(genomic) %>%
  mutate(Chromosome = 
           (eukrep == "Prokarya" & plasflow_prediction == "chromosome" & cbar_prediction == "Chromosome" & is.na(plasmidfinder_plasmid) & deepvirfinder_pvalue > 0.05 & is.na(virsorter_score)),
         Plasmid = 
           ((!is.na(plasmidfinder_plasmid) | (plasflow_prediction == 'plasmid' & cbar_prediction == "Plasmid")) & eukrep == "Prokarya" & (!virsorter_score %in% c(1, 2)) & deepvirfinder_pvalue > 0.05),
         Phage = 
           ((virsorter_category == "phage" & !is.na(virsorter_score) & deepvirfinder_pvalue < 0.01) & cbar_prediction != "Plasmid" & plasflow_prediction != "plasmid" & eukrep == "Prokarya"),
         `Putative phage` = 
           !Phage & (((virsorter_category == "phage" & !is.na(virsorter_score)) | deepvirfinder_pvalue < 0.05) & cbar_prediction != "Plasmid" & plasflow_prediction != "plasmid" & eukrep != "Eukarya"),
         Prophage = 
           (virsorter_category == "prophage" & !is.na(virsorter_score) & cbar_prediction != "Plasmid" & plasflow_prediction != "plasmid" & eukrep == "Prokarya")) %>%
  select(scaffold, region, products, Genomic, Chromosome, Plasmid, Phage, `Putative phage`, Prophage) %>%
  replace(is.na(.), FALSE)

table_to_plot1 = table_summary %>%
  separate_rows(products, sep = ";") %>%
  group_by(products) %>%
  summarize(n = n(),
            Genomic = sum(Genomic)/n(),
            Chromosome = sum(Chromosome)/n(),
            Plasmid = sum(Plasmid)/n(),
            Phage = sum(Phage)/n(),
            `Putative phage` = sum(`Putative phage`)/n(),
            Prophage = sum(Prophage)/n()) %>%
  filter(n > 100)

p1 = ggplot() +
  geom_point(aes(x = n, y = Genomic), size = 1, data = table_to_plot1 %>% filter(products != 'nucleoside'), alpha = .8) +
  geom_point(aes(x = n, y = Genomic), size = 1, color = "#17BECF", data = table_to_plot1 %>% filter(products == 'nucleoside'), alpha = .8) +
  ggrepel::geom_label_repel(aes(x = n, y = Genomic, label = products), size = 2, label.padding = unit(.5, 'mm'), data = table_to_plot1 %>% filter(products != 'nucleoside')) +
  ggrepel::geom_label_repel(aes(x = n, y = Genomic, label = products), size = 2, label.padding = unit(.5, 'mm'), color = "#17BECF", data = table_to_plot1 %>% filter(products == 'nucleoside')) +
  ylim(0, 1) +
  scale_x_log10() +
  ylab("Proportion of genomic BGCs") +
  xlab("Number of BGCs predicted to encode the product type (prdoct types with > 100 BGCs)") +
  theme_bw() +
  theme(text = element_text(size = 6),
        line = element_line(size = unit(.3, 'mm')))
p1

table_to_plot2 = table_summary %>%
  separate_rows(products, sep = ";") %>%
  filter(!Genomic) %>%
  group_by(products) %>%
  summarize(n = n(),
            Chromosome = sum(Chromosome)/n(),
            Plasmid = sum(Plasmid)/n(),
            Phage = sum(Phage)/n(),
            `Putative phage` = sum(`Putative phage`)/n(),
            Prophage = sum(Prophage)/n()) %>%
  filter(n > 100) %>%
  select(-n) %>%
  gather(key = type, value = freq, -products) %>%
  mutate(type = factor(type, levels = c("Chromosome", "Plasmid", "Prophage", "Phage", "Putative phage")))

p2 = table_to_plot2 %>%
  ggplot() +
  geom_boxplot(aes(x = type, y = freq), outlier.shape = NA) +
  geom_point(aes(x = type, y = freq), data = table_to_plot2 %>% filter(products != 'nucleoside'), alpha = .8) +
  geom_point(aes(x = type, y = freq), color = "#17BECF", data = table_to_plot2 %>% filter(products == 'nucleoside'), alpha = .8) +
  ylab("Proportion of BGCs") +
  xlab("Remaining metagenomic fragments (excl. unannotated)") +
  theme_bw() +
  theme(text = element_text(size = 6),
        line = element_line(size = unit(.3, 'mm')))
p2

p = p1 + p2 + plot_layout(widths = c(.6, .4))
p

ggsave("data/processed/figures/Figure-S11-nucleosides/Figure-S11-nucleosides.raw.pdf", p, width = 183, height = 90, unit = 'mm')

# How many gcfs?
table_gcf = read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_metadata.tsv") %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_to_gcf_gcc.tsv")) %>%
  mutate(region = paste0(scaffold, "-biosynth_", gsub(".*.region0+|.gbk", "", file))) %>%
  arrange(dataset) %>%
  filter(!duplicated(region)) %>%
  select(region, gcf, gcc)

table_summary %>% 
  left_join(table_gcf) %>%
  filter(products == "nucleoside") %>%
  filter(!Genomic) %>% #View()
  pull(gcf) %>%
  unique %>%
  length

table_summary %>% 
  left_join(table_gcf) %>%
  filter(products == "nucleoside") %>%
  filter(!Genomic) %>% #View()
  pull(gcf) %>%
  table %>%
  sort

table_summary %>% 
  left_join(table_gcf) %>%
  filter(products == "nucleoside") %>%
  filter(!Genomic) %>% #View()
  pull(gcc) %>%
  table

table_summary %>% 
  left_join(table_gcf) %>%
  filter(gcc == 58) %>% View()
