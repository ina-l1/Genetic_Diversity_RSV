# Data visualization
# Sequence count (GenBank) and case count (WHO FluNet) 

library(tidyr)
library(ggplot2)
library(lubridate)
library(MMWRweek)

#################################################################################

# Base directory

base_dir <- "~/Yale_Projects/Genetic_Diversity_RSV/"

##################################################

# Sequence count: Read metadata

path_rsvAB <- file.path(base_dir, "Germany", "rsvAB_metadata_GER.csv")
rsvAB <- read.csv(path_rsvAB, row.names = 1)

rsvAB$Collection_Date <- as_date(rsvAB$Collection_Date)
rsvAB$MMWR_start_date <- MMWRweek2Date(MMWRyear = rsvAB$MMWRyear, MMWRweek = rsvAB$MMWRweek)
rsvAB$MMWR_month_year <- as_date(paste(rsvAB$Collection_Year, rsvAB$Collection_Month, 1, sep = "-"))

path_rsvAB_EU <- file.path(base_dir, "Europe", "rsvAB_metadata_EU.csv")
rsvAB_EU <- read.csv(path_rsvAB_EU, row.names = 1)

rsvAB_EU$Collection_Date <- as_date(rsvAB_EU$Collection_Date)
rsvAB_EU$MMWR_start_date <- MMWRweek2Date(MMWRyear = rsvAB_EU$MMWRyear, MMWRweek = rsvAB_EU$MMWRweek)

rsvAB_EU$'EU/GER' <- ifelse(rsvAB_EU$Country == "Germany", "GER", "EU")

rsvAB_EU_noGer <- subset(rsvAB_EU, Country != "Germany")

rsvAB_ESP <- subset(rsvAB_EU, Country == "Spain")

seq_count_DEU <- as.data.frame(table(rsvAB$MMWR_start_date))
colnames(seq_count_DEU) <- c("MMWR_start_date", "count")
seq_count_DEU$MMWR_start_date <- as_date(seq_count_DEU$MMWR_start_date)

seq_count_ESP <- as.data.frame(table(rsvAB_ESP$MMWR_start_date))
colnames(seq_count_ESP) <- c("MMWR_start_date", "count")
seq_count_ESP$MMWR_start_date <- as_date(seq_count_ESP$MMWR_start_date)

#################################################################################

# Case count: Read metadata

start_date <- MMWRweek2Date(2014, 40)
end_date <- MMWRweek2Date(2023, 39) #Most recent sequence data available from 2022/2023 season

path_case_flunet <- file.path(base_dir, "Europe", "VIW_FNT.csv")
case_flunet <- read.csv(path_case_flunet)

case_flunet$ISO_WEEKSTARTDATE <- as_date(case_flunet$ISO_WEEKSTARTDATE)
case_flunet$MMWR_WEEKSTARTDATE <- as_date(case_flunet$MMWR_WEEKSTARTDATE)

case_DEU <- case_flunet %>% subset(COUNTRY_CODE == "DEU") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date)

case_NLD <- case_flunet %>% subset(COUNTRY_CODE == "NLD") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date)

case_ESP <- case_flunet %>% subset(COUNTRY_CODE == "ESP") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date)

case_EU <- case_flunet %>% subset(COUNTRY_CODE == "AUT" | COUNTRY_CODE == "FIN" | COUNTRY_CODE == "FRA" | COUNTRY_CODE == "ITA" |
                                    COUNTRY_CODE == "NLD" | COUNTRY_CODE == "RUS" | COUNTRY_CODE == "ESP" | COUNTRY_CODE == "X09") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date) #GenBank sequences only available for these countries

#################################################################################

# Sequence count (GenBank)

tab_rsvAB_season <- as.data.frame.matrix(table(rsvAB$Collection_Season, rsvAB$Type))
tab_rsvAB_EU_season <- as.data.frame.matrix(table(rsvAB_EU$Collection_Season, rsvAB_EU$Type))

# Time series

timeseries_plot <- ggplot(rsvAB, aes(x = Collection_Date, color = Type, fill = Type)) +
  geom_histogram(binwidth = 7)+
  scale_x_date(
    date_breaks = "3 months", 
    date_minor_breaks = "1 week", 
    date_labels = "%W"
  ) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  facet_grid(Type ~ .)
timeseries_plot

timeseries_stack_plot <- ggplot() +
  geom_bar(data = rsvAB, aes(x = MMWR_start_date, fill = Type)) +
  scale_fill_manual(values = c("red", "blue")) +
  xlab(label = "collection year") +
  ylab(label = "sequence count") +
  theme_minimal()
timeseries_stack_plot

ggplot(rsvAB, aes(x = Collection_Date, color = Type, fill = Type)) +
  geom_histogram(binwidth = 7)+
  scale_x_date(
    date_breaks = "1 year", 
    date_minor_breaks = "1 month", 
    date_labels = "%Y"
  ) +
  scale_y_continuous() +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal()

#EU + GER
custom_date_labels <- function(dates) {
  ifelse(is.na(dates), NA, ifelse(format(dates, "%m") == "01", format(dates, "%b %Y"), format(dates, "%b")))
}

ggplot() +
#  geom_bar(data = rsvAB_EU, aes(x = MMWR_start_date)) +
  geom_bar(data = rsvAB, aes(x = MMWR_month_year, fill = Type)) +
  scale_fill_manual(values = c("red", "blue"), name = "RSV subtype") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_date(
    date_breaks = "3 months",
    labels = function(x) custom_date_labels(x), # Apply robust custom labels
    guide = guide_axis(angle = -90),
    limits = c(as_date("2014-9-01"), max(rsvAB$MMWR_month_year))
  ) +
  xlab(label = "collection date") +
  ylab(label = "sequence count") +
  theme_minimal() +
  theme(legend.position = "top",
        text = element_text(size = 30),
        axis.title.x = element_text(margin = margin(t = 10)))

#################################################################################

# Case count (FluNet)

case_count_DEU <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_DEU, sum)
case_count_NLD <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_NLD, sum)
case_count_ESP <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_ESP, sum)
case_count_EU <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_EU, sum)

ggplot(case_count_DEU, aes(x = MMWR_WEEKSTARTDATE, y = RSV)) +
  geom_bar(stat = "identity") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  xlab("year") +
  ylab("RSV case count") +
  theme_minimal()

ggplot(case_count_DEU, aes(x = MMWR_WEEKSTARTDATE, y = RSV)) +
  geom_line() +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  xlab("year") +
  ylab("RSV case count") +
  theme_minimal()

#################################################################################

# Seq vs case count

ggplot() +
  geom_line(data = case_count_DEU, aes(x = MMWR_WEEKSTARTDATE, y = RSV), colour = "orange") +
  geom_line(data = seq_count_DEU, aes(x = MMWR_start_date, y = count), colour = "darkgreen") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  xlab("year") +
  ylab("RSV sequence count")

ggplot() +
  geom_bar(data = case_count_DEU, aes(x = MMWR_WEEKSTARTDATE, y = RSV), fill = "darkgreen", colour = "darkgreen", stat = "identity") +
  geom_bar(data = seq_count_DEU, aes(x = MMWR_start_date, y = count), fill = "orange", colour = "orange", stat = "identity") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  xlab("year") +
  ylab("count") +
  theme_minimal()
  

ggplot() +
  geom_line(data = case_count_ESP, aes(x = MMWR_WEEKSTARTDATE, y = RSV), colour = "orange") +
  geom_line(data = seq_count_ESP, aes(x = MMWR_start_date, y = count), colour = "darkgreen") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  xlab("year") +
  ylab("RSV sequence count")