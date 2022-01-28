# 2017_tomato_visits

## DESCRIPTION

Analysis for Aguirre, LA, and Adler, LS. 2022. Interacting Antagonims: Parasite Infection Alters Bombus Impaties (Hymenoptera: Apidae) Responses to herbivory on Tomato Plants. In Review.\

Code Written by LA Aguirre\
Data Cleaning, Analysis and Visualizations\

---
## METADATA 
### counts.csv variables\
bid: Unique bee identifier\
date: Date of observation\
count: Crithidia cell counts per 0.02uL\
ctime: Time of day when Crithidia was counte\
dtime: Time of day when bee was dissected\
wing_raw: Length of marginal cell (right wing) in ocular units at 400x magnification\

#### Variables created in script:\
wing_mm: Ocular units converted from wing_raw (13x/20 convertion rate)\
jdate: Julian date from date\

---
### visits.csv variables
bid: Unique bee identifier\
start: Plant treatment of first observation\
plant: Plant treatment\
seconds: Seconds on flower\
frames: Frames spent on flower (GarageBand application measured recording time in seconds and frames. 25 frames per second)\

#### Variables created in script:
hund: hundredths of second spend on flower from frames\
count: from count in `counts.csv`\
jdate: Julian date\
wing_mm: ocular units converted from wing_raw (13x/20 convertion rate)\
status: Infection status (infected/uninfected) from count\
total_visits: number of visits observed for individual bee\
fjdate: Julian date as factor\
fbid: Unique bee identifier as factor\
