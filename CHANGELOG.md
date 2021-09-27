# Changelog

## 2021-09-27

**Enhancements:**

- We've added contact matrices for all national and subnational locations aggregated to 18 age brackets
- An [example script](https://github.com/mobs-lab/mixing-patterns/blob/matrices_18/analysis_results/scripts/make_matrices_diff_aggregations.py) shows how to create different age mappings and aggregate the matrices and age count distributions to other age brackets.
- Our original age distributions represented samples of the population which we generated as synthetic populations. We have now added [age distributions](https://github.com/mobs-lab/mixing-patterns/tree/matrices_18/data/population_rescaled_age_distributions) scaled up to the real world population sizes for each location as well to aid anyone interested in that data. Note that our study focused on using census data from each country around the period of 2010-2011. While many countries have had additional censuses or population estimates since then, some of the countries we focus on here (for example China, India, Russia, South Africa) have yet to release
new estimates on age distributions for subnational locations. For this reason, we focus on providing estimates for the 2010-2011 period. When new age estimates are released for the majority of locations in our study, we will document those data here as well so keep an eye on this page.
- In our [Nature Communications article](https://www.nature.com/articles/s41467-020-20544-y), we provide a table documenting the sources for all of our data. To add further visibility, we have also added an online [table](https://github.com/mobs-lab/mixing-patterns/blob/matrices_18/data/population_size.csv) reiterating the source of age distribution and population size data.

**Fixed:**
- One of the contact matrices (Canada_country_level_F_school_setting_85.csv) was delimited using spaces instead of commas as with the rest of the contact matrices in the repository. Fixed this file to be consistent. This addresses issue [#12](https://github.com/mobs-lab/mixing-patterns/issues/12).

**Closed issues:**
- The online table for age distribution data addresses issue [#4](https://github.com/mobs-lab/mixing-patterns/issues/4) asking about data sources. For those interested in specific sources, please reference this table. If further clarification is needed, please reach out to us.


## 2021-01-19

**Fixed:**

- Updated the link for the article to point to the official Nature Communications publication doi.

## 2020-11-23

**Enhancements:**
- We're sharing our data on age mixing contact matrices in 277 subnational locations for social settings important to understanding the transmission of infectious diseases, in particular airborne pathogens and other diseases communicable through shared space and contact. Also included are age distributions, national estimates for the countries of focus in this study, and national contact matrices for many countries in Europe from a [previous study](https://doi.org/10.1371/journal.pcbi.1002673). 