# load tidyverse package
library(tidyverse)

pca_input="all_samples_QC_MAC_LD.eigenvec"
eigenval_input="all_samples_QC_MAC_LD.eigenval"

# read in data 
pca <- read_table(pca_input, col_names = FALSE)
eigenval <- scan(eigenval_input)

# remove nuisance column
pca <- pca[,-1]

# Rename col1 to ind and following columns PC1, PC2, etc.
names(pca)[1] <- "ind" 
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Add location information
loc <- rep(NA, length(pca$ind))
loc[grep("BP", pca$ind)] <- "DRC"
loc[grep("AG", pca$ind)] <- "Gambia"
loc[grep("AK", pca$ind)] <- "Kenya"
loc[grep("BL", pca$ind)] <- "Tanzania"


# remake data.frame to tibble and add location information
pca <- as_tibble(data.frame(pca, loc))

# convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# plot amount of variance explained by PCs
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()


# plot PC1 and PC2
b <- ggplot(pca, aes(PC1, PC2, col =loc)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("#1B3F7D", "#BC413A", "#2E6B34", "#EF8632"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# plot PC1 and PC3
b <- ggplot(pca, aes(PC1, PC3, col =loc)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("#1B3F7D", "#BC413A", "#2E6B34", "#EF8632"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

# plot PC2 and PC3
b <- ggplot(pca, aes(PC2, PC3, col =loc)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("#1B3F7D", "#BC413A", "#2E6B34", "#EF8632"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))




