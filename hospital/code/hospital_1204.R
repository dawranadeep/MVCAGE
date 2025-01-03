library(tidycensus)
library(tidyverse)
library(sf)

setwd("/Users/rdaw/Documents/MVCAGE/hospital")
options(tigris_use_cache = TRUE)


# Read the original data
df = read.csv("data/county_hedis_6575ffs.csv")

colnames(df)

summary(df$adjusted_rate)

df$adjusted_rate[df$adjusted_rate == -99999] = NA
df$adjusted_rate[df$adjusted_rate == -88888] = NA

rate = df$adjusted_rate


idx1 = df$Year == 2015
idx2 = !is.na(df$adjusted_rate)
idx = idx1 & idx2

#plot(density(rate[idx]))


ratedf = df[idx, ]
ratedf = ratedf[order(ratedf$Geo_code), ]
ratedf = ratedf[ratedf$Race == "Overall", ]
ratedf = select(ratedf, c(Geo_code, adjusted_rate))
ratedf = ratedf %>% rename(GEOID = Geo_code)

readRenviron("~/.Renviron")
Sys.getenv("CENSUS_API_KEY")

income <- get_acs(
  geography = "county", 
  variables = "B19013_001", 
  year = 2015, 
  geometry = T
)

income2 = income
crs = st_crs(income)

income = as.data.frame(income)
income$NAME = NULL
income$variable = NULL
income$moe = NULL
income$GEOID = as.numeric(income$GEOID)
income = income %>% rename(med_income = estimate)
centroids = st_centroid(income2) %>% st_coordinates()
income$x = centroids[,1]
income$y = centroids[,2] 
income$areas = st_area(income2)
jdf = inner_join(income, ratedf, by="GEOID")


outdf = read.csv("output/hosp.csv", header = F)
names(outdf) = c("GEOID", "Predicted_Income", "Predicted_Hospital_Quality", "e1hat", "e2hat", "e1hat_A", "e2hat_A", "clusters")


# Values are from the Matlab code
mY1 = 9.8693 
MY1 = 11.7236
mY2 = 22.2045
MY2 = 93.6169
outdf$Predicted_Income = mY1 + (MY1 - mY1)*outdf$Predicted_Income
outdf$Predicted_Hospital_Quality = mY2 + (MY2 - mY2)*outdf$Predicted_Hospital_Quality  

income3 = inner_join(jdf, outdf, by = "GEOID")
income3$med_income = log(income3$med_income)


## all Fips
maps::county.fips %>%
  as_tibble %>% 
  extract(polyname, c("region", "subregion"), "^([^,]+),([^,]+)$") -> dfips

income3$fips = income3$GEOID

## some county maps & left_join

map_data("county") %>% 
  inner_join(dfips) -> data


## more left_join
data<-inner_join(data, income3, by="fips")


library(cowplot)

library(ggpubr)

## map with random, fictive cows
p1 = data %>%   ggplot( aes(long, lat, group = group)) +
  geom_polygon(aes(fill=med_income)) +
  coord_map() +
  theme_void() +
  theme(legend.key.size = unit(1.5, 'cm'))+
  scale_fill_gradientn(
    colours = hcl.colors(10, "YlGnBu"),
    limits = c(mY1, MY1), name="")+
  ggtitle("Log Median Income")


p2 = data %>%   ggplot( aes(long, lat, group = group)) +
  geom_polygon(aes(fill=Predicted_Income)) +
  coord_map() +
  theme_void() +
  theme(legend.key.size = unit(1.5, 'cm'))+
  scale_fill_gradientn(
    colours = hcl.colors(10, "YlGnBu"),
    limits = c(mY1, MY1), name="")+
  ggtitle("Predicted Log Median Income")




## map with random, fictive cows
p3 = data %>%   ggplot( aes(long, lat, group = group)) +
  geom_polygon(aes(fill=adjusted_rate)) +
  coord_map() +
  theme_void() +
  theme(legend.key.size = unit(1.5, 'cm'))+
  scale_fill_gradientn(
    colours = hcl.colors(10, "YlGnBu"),
    limits = c(mY2, MY2), name="")+
  ggtitle("Hospital Ratings")


p4 = data %>%   ggplot( aes(long, lat, group = group)) +
  geom_polygon(aes(fill=Predicted_Hospital_Quality)) +
  coord_map() +
  theme_void()+
  theme(legend.key.size = unit(1.5, 'cm')) +
  scale_fill_gradientn(
    colours = hcl.colors(10, "YlGnBu"),
    limits = c(mY2, MY2), name="")+
  ggtitle("Predicted Hospital Ratings")


png("output/income1204_2.png", width = 480, height = 960)
ggpubr::ggarrange(p1, p2, ncol=1,
                  common.legend = T, # COMMON LEGEND
                  legend = "bottom", # legend position
                  align = "v", # Align them both, horizontal and vertical
                  nrow = 2)  # number of rows
dev.off()



png("output/hospital1204_2.png", width = 480, height = 960)
ggpubr::ggarrange(p3, p4, ncol=1,
                  common.legend = T, # COMMON LEGEND
                  legend = "bottom", # legend position
                  align = "v", # Align them both, horizontal and vertical
                  nrow = 2)  # number of rows
dev.off()




######## Plot Regions


p5 = data %>%   ggplot( aes(long, lat, group = group)) +
  geom_polygon(aes(fill=clusters)) +
  coord_map() +
  theme_void() + theme(legend.position = "none") + 
  scale_fill_gradientn(
    colours = hcl.colors(10, "YlGnBu"), name="")+
  ggtitle("Joint Regionalization of Log Median Income and Hospital Ratings")


png("output/hosp_regions1204.png", width = 480, height = 480)
p5
dev.off()


