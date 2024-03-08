filter2_BH_spear <- df3[, c(1:3)]
leg_connects <- subset(filter2_BH_spear, origin == "Legionella" | destination == "Legionella") 

acan_count1 <- subset(filter2_BH_spear, origin == "Legionella" | origin == "Acanthamoeba") 
acan_count2 <- subset(acan_count1, destination == "Legionella" | destination == "Acanthamoeba") 
acan_count3 <- acan_count2 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(origin == destination,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)
acan_count4 <- subset(acan_count3, value >0.5) 
acan <- nrow(acan_count4)

echi_count1 <- subset(filter2_BH_spear, origin == "Legionella" | origin == "Echinamoeba") 
echi_count2 <- subset(echi_count1, destination == "Legionella" | destination == "Echinamoeba") 
echi_count3 <- echi_count2 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(origin == destination,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)
echi_count4 <- subset(echi_count3, value >0.5) 
echi <- nrow(echi_count4)

koro_count1 <- subset(filter2_BH_spear, origin == "Legionella" | origin == "Korotnevella") 
koro_count2 <- subset(koro_count1, destination == "Legionella" | destination == "Korotnevella") 

koro_count3 <- koro_count2 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(origin == destination,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)
koro_count4 <- subset(koro_count3, value >0.5) 
koro <- nrow(koro_count4)

naeg_count1 <- subset(filter2_BH_spear, origin == "Legionella" | origin == "Naegleria") 
naeg_count2 <- subset(naeg_count1, destination == "Legionella" | destination == "Naegleria") 
naeg_count3 <- naeg_count2 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(origin == destination,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)
naeg_count4 <- subset(naeg_count3, value >0.5) 
naeg <- nrow(naeg_count4)


tetra_count1 <- subset(filter2_BH_spear, origin == "Legionella" | origin == "Tetrahymena") 
tetra_count2 <- subset(tetra_count1, destination == "Legionella" | destination == "Tetrahymena") 

tetra_count3 <- tetra_count2 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(origin == destination,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)
tetra <- nrow(tetra_count3)

vann_count1 <- subset(filter2_BH_spear, origin == "Legionella" | origin == "Vannella") 
vann_count2 <- subset(vann_count1, destination == "Legionella" | destination == "Vannella") 
vann_count3 <- vann_count2 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(origin == destination,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)
vann_count4 <- subset(vann_count3, value >0.5) 
vann <- nrow(vann_count4)


verma_count1 <- subset(filter2_BH_spear, origin == "Legionella" | origin == "Vermamoeba") 
verma_count2 <- subset(verma_count1, destination == "Legionella" | destination == "Vermamoeba") 

verma_count3 <- verma_count2 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(origin == destination,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)
verman <- nrow(verma_count3)



leg_count1 <- subset(filter2_BH_spear, origin == "Legionella" | destination == "Legionella") 
leg_count2 <- leg_count1 %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(origin == destination,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)
leg_count3 <- subset(leg_count2, value >0.5) 
leg <- nrow(leg_count3)


first_column <- c(acan, echi, koro, naeg, tetra, vann, verman)
second_column <- c(acan/leg, echi/leg, koro/leg, naeg/leg, tetra/leg, vann/leg, verman/leg)

edges <- data.frame(first_column, second_column)

