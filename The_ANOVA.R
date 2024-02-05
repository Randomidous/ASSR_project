# Load tidyverse package
library(tidyverse)
library(readxl)
library(broom)
library(knitr)
library(xlsx)

# Load Excel file
anova_left_VR <- read.xlsx('C:/Users/ericw/Desktop/Master_Thesis/Data/CURRENT_DATA/10_STATS/anova_left_VR.xlsx', 1, header=FALSE)
anova_right_VR <- read.xlsx('C:/Users/ericw/Desktop/Master_Thesis/Data/CURRENT_DATA/10_STATS/anova_right_VR.xlsx', 1, header=FALSE)
anova_left_noVR <- read.xlsx('C:/Users/ericw/Desktop/Master_Thesis/Data/CURRENT_DATA/10_STATS/anova_left_noVR.xlsx', 1, header=FALSE)
anova_right_noVR <- read.xlsx('C:/Users/ericw/Desktop/Master_Thesis/Data/CURRENT_DATA/10_STATS/anova_right_noVR.xlsx', 1, header=FALSE)

# Assuming you already have the tables loaded into anova_left_VR, anova_right_VR, anova_left_noVR, and anova_right_noVR

# Function to split even and odd rows
split_rows <- function(data) {
  even_rows <- data[seq(2, nrow(data), by = 2), ]
  odd_rows <- data[seq(1, nrow(data), by = 2), ]
  return(list(even_rows, odd_rows))
}

# Split rows for each table
split_anova_left_VR <- split_rows(anova_left_VR)
split_anova_right_VR <- split_rows(anova_right_VR)
split_anova_left_noVR <- split_rows(anova_left_noVR)
split_anova_right_noVR <- split_rows(anova_right_noVR)

# Create new variables with appended suffixes
anova_left_VR_39 <- split_anova_left_VR[[2]]
anova_left_VR_41 <- split_anova_left_VR[[1]]

anova_right_VR_39 <- split_anova_right_VR[[2]]
anova_right_VR_41 <- split_anova_right_VR[[1]]

anova_left_noVR_39 <- split_anova_left_noVR[[2]]
anova_left_noVR_41 <- split_anova_left_noVR[[1]]

anova_right_noVR_39 <- split_anova_right_noVR[[2]]
anova_right_noVR_41 <- split_anova_right_noVR[[1]]


# Now we can create data frames
df39 <- data.frame(group=rep(c('VR', 'noVR'), each=8),
                 direction=rep(c('left', 'right'), each=4),
                 freqs=c(anova_left_VR_39[,1], anova_right_VR_39[,1], anova_left_noVR_39[,1], anova_right_noVR_39[,1]),
                 power=c(anova_left_VR_39[,2], anova_right_VR_39[,2], anova_left_noVR_39[,2], anova_right_noVR_39[,2]))


df41 <- data.frame(group=rep(c('VR', 'noVR'), each=8),
                   direction=rep(c('left', 'right'), each=4),
                   freqs=c(anova_left_VR_41[,1], anova_right_VR_41[,1], anova_left_noVR_41[,1], anova_right_noVR_41[,1]),
                   power=c(anova_left_VR_41[,2], anova_right_VR_41[,2], anova_left_noVR_41[,2], anova_right_noVR_41[,2]))


# Convert group and direction to factors
df39$group <- as.factor(df39$group)
df39$direction <- as.factor(df39$direction)

df41$group <- as.factor(df41$group)
df41$direction <- as.factor(df41$direction)


# Perform two-way ANOVA
output_directory <- "C:/Users/ericw/Desktop/Master_Thesis/Data/CURRENT_DATA/10_STATS/"

anova_result_39 <- aov(power ~ group * direction, data = df39)
anova_tidy_39 <- tidy(anova_result_39)
kable(anova_tidy_39, format = "markdown", digits = 3)
write.table(anova_tidy_39, file = paste0(output_directory, "anova_table_39.txt"), sep = "\t", quote = FALSE)

anova_result_41 <- aov(power ~ group * direction, data = df41)
anova_tidy_41 <- tidy(anova_result_41)
kable(anova_tidy_41, format = "markdown", digits = 3)
write.table(anova_tidy_41, file = paste0(output_directory, "anova_table_41.txt"), sep = "\t", quote = FALSE)

#### NICE #####

library(apaTables)

lm_output_39 <- lm(power ~ group * direction, data = df39)
apa.aov.table(lm_output_39, filename = "C:/Users/ericw/Desktop/Master_Thesis/Data/CURRENT_DATA/10_STATS/anova_table_39.doc")

lm_output_41 <- lm(power ~ group * direction, data = df41)
apa.aov.table(lm_output_41, filename = "C:/Users/ericw/Desktop/Master_Thesis/Data/CURRENT_DATA/10_STATS/anova_table_41.doc")


#### Lets do only 1 and 10

# Now we can create data frames
df39 <- data.frame(group=rep(c('VR', 'noVR'), each=4),
                   direction=rep(c('left', 'right'), each=1),
                   freqs=c(anova_left_VR_39[1,1], anova_left_VR_39[3,1], anova_right_VR_39[1,1], anova_right_VR_39[3,1], anova_left_noVR_39[2,1], anova_left_noVR_39[3,1], anova_right_noVR_39[2,1], anova_right_noVR_39[3,1]),
                   power=c(anova_left_VR_39[1,2], anova_left_VR_39[3,2], anova_right_VR_39[1,2], anova_right_VR_39[3,2], anova_left_noVR_39[2,2], anova_left_noVR_39[3,2], anova_right_noVR_39[2,2], anova_right_noVR_39[3,2]))


df41 <- data.frame(group=rep(c('VR', 'noVR'), each=4),
                   direction=rep(c('left', 'right'), each=1),
                   freqs=c(anova_left_VR_41[1,1], anova_left_VR_41[3,1], anova_right_VR_41[1,1], anova_right_VR_41[3,1], anova_left_noVR_41[2,1], anova_left_noVR_41[3,1], anova_right_noVR_41[2,1], anova_right_noVR_41[3,1]),
                   power=c(anova_left_VR_41[1,2], anova_left_VR_41[3,2], anova_right_VR_41[1,2], anova_right_VR_41[3,2], anova_left_noVR_41[2,2], anova_left_noVR_41[3,2], anova_right_noVR_41[2,2], anova_right_noVR_41[3,2]))

anova_result_39 <- aov(power ~ group * direction, data = df39)
anova_result_41 <- aov(power ~ group * direction, data = df41)

anova_tidy_39 <- tidy(anova_result_39)
anova_tidy_41 <- tidy(anova_result_41)

kable(anova_tidy_39, format = "markdown", digits = 3)
kable(anova_tidy_41, format = "markdown", digits = 3)


### t-test 39 and 41

# Convert the columns to numeric
data39[,5] <- as.numeric(as.character(data39[,5]))
data41[,5] <- as.numeric(as.character(data41[,5]))

# Now perform the paired t-test
t_test_result <- t.test(data39[,5], data41[,5], paired = TRUE)

# View the result
print(t_test_result)

# Mean and SD for the 5th column of data39
mean_data39 <- mean(data39[,5], na.rm = TRUE)
sd_data39 <- sd(data39[,5], na.rm = TRUE)

# Print the results for data39
print(paste("Mean for data39: ", mean_data39))
print(paste("SD for data39: ", sd_data39))

# Mean and SD for the 5th column of data41
mean_data41 <- mean(data41[,5], na.rm = TRUE)
sd_data41 <- sd(data41[,5], na.rm = TRUE)

# Print the results for data41
print(paste("Mean for data41: ", mean_data41))
print(paste("SD for data41: ", sd_data41))
