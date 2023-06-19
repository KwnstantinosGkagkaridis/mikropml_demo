library(tidyverse)
library(broom)
library(ggtext)

set.seed(19760620)

#load the shared file, count info for each taxon in each of the 490 samples
shared <- read_tsv("raw_data/baxter.subsample.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double())) %>%
  rename_all(tolower) %>%
  select(group, starts_with("otu")) %>%
  pivot_longer(-group, names_to = "otu", values_to = "count")


#info about the taxonomic name for all of the genera in the analysis
taxonomy <- read_tsv("raw_data/baxter.cons.taxonomy") %>%
  rename_all(tolower) %>%
  select(-size) %>%
  mutate(otu = tolower(otu),
         taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""), #\\( --> match an actual parenthesis, 
         #\\d+ --> matches one or more digit characters like a number, 
         #\\) --> close the parenthesis, 
         #"" --> replace it with nothing
         #remove the numbers from cells of the taxonomy columns
         taxonomy = str_replace(taxonomy, ";unclassified", "_unclassified"),
         taxonomy = str_replace(taxonomy, ";unclassified", ""),
         taxonomy = str_replace(taxonomy, ";$", ""),
         taxonomy = str_replace(taxonomy, ".*;", "")) 


#info about the disease status
#create two columns: srn & lesion
metadata <-read_tsv("raw_data/baxter.metadata.tsv",
                    col_types = cols(sample = col_character(), 
                                     Hx_Prev = col_logical())) %>%
  rename_all(tolower) %>%
  rename(group = sample) %>%
  mutate(srn = dx_bin == "Adv Adenoma" | dx_bin == "Cancer",
         lesion = dx_bin == "Adv Adenoma" | dx_bin == "Cancer" | dx_bin == "Adenoma")  # srn --> screen relavent neoplasia


#join them all together
composite <- inner_join(shared, taxonomy, by = "otu") %>%
  group_by(group, taxonomy) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  group_by(group) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by= "group")
