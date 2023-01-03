# download gd_holygrail_filtered_final.rds from here: https://drive.google.com/file/d/16SeXuyOkRm3y33TqxwIcWewmBY_Nbpfw/view?usp=share_link

multimodal_gd_gex <- readRDS("gd_holygrail_filtered_final.rds")

#figure 4a
Idents(multimodal_gd_gex) <- "celltypes"
factor(Idents(multimodal_gd_gex), levels= c("others" ,"Areg+","Rorc+","Proliferating","Klrg1+ effector","Gzmb+ IELs","Cd160+","Sell+Ly6c+","Sell+Ly6c-"))
Idents(multimodal_gd_gex) <- factor(Idents(multimodal_gd_gex), levels= c("others" ,"Areg+","Rorc+","Proliferating","Klrg1+ effector","Gzmb+ IELs","Cd160+","Sell+Ly6c+","Sell+Ly6c-"))

DotPlot(multimodal_gd_gex, features = c("Sell","Ly6c2","Cd160","Gzmb","Klrg1","Mki67","Rorc","Areg"),cols = c("RdYlBu")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + RotatedAxis() +theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  text = element_text(size = 10),
  panel.grid.major.x = element_line(color = "grey80"),
  panel.grid.major.y = element_line(color = "grey80") 
  
)

#figure 4b

table_organ_by_celltypes <- multimodal_gd_gex@meta.data %>%
  group_by(organ, celltypes) %>%
  summarize(count = n()) %>%
  spread(celltypes, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('organ', 'total_cell_count', everything())) %>%
  arrange(factor(organ, levels = levels(multimodal_gd_gex@meta.data$organ)))

library(dplyr)
x <- c("LN","spleen","liver","lung","skin","LI","SI","pooled")

table_organ_by_celltypes <- table_organ_by_celltypes %>% mutate(organ =  factor(organ, levels = x)) %>% arrange(organ) 
temp_labels <- multimodal_gd_gex@meta.data %>% group_by(organ) %>% tally()
temp_labels <-  temp_labels %>% mutate(organ =  factor(organ, levels = x)) %>% arrange(organ)

table_organ_by_celltypes %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'organ') %>%
  ggplot(aes(organ, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = organ, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'celltypes', values = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471')) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )
