# README
Different locations might include only a single country level matrix or additional subnational matrices as indicated below.

For each location, 5 different matrices are reported:

1. a household setting contact matrix (filenames: *_F_household_setting_85.csv) 
2. a school setting contact matrix (filenames: *_F_school_setting_85.csv) 
3. a work setting contact matrix (filenames: *_F_work_setting_85.csv) 
4. a community setting contact matrix (filenames: *_F_community_setting_85.csv)
5. the overall contact matrix obtained as a linear combination of the previous matrices (filenames: *_M_contact_matrix_flu_85.csv). 

6. the age distributions used to build the matrices (filenames: *_age_distribution_85.csv). 

NOTE: the weights used to obtain the overall contact matrix are: 4.11 contacts for the household setting, 11.41 contacts for the school setting, 8.07 contacts for the work setting, and 2.79 contacts for the general community setting.

Each contact matrix has dimensions 85x85. For the 4 settings (households, schools, workplaces, and the community) contact matrices, the matrix element F<sup>k</sup><sub>ij</sub> represents the per capita probability of contact for an individual of age i with individuals of age j in that setting k. Thus the element F<sup>k</sup><sub>00</sub> is the per capita probability of contact for an individual of age [0, 1) with individuals of age [0, 1), the element F<sup>k</sup><sub>01</sub> is the per capita probability of contact for an individual of age [0, 1) with individuals of age [1, 2), etc. Each row and column refers to a single year of age, except for the last row and column which both refer to ages 84 years old and over. 

The elements M<sub>ij</sub> of the overall contact matrix are obtained as a linear combination of the 4 setting (per capita) contact matrices. The weight of each setting contact matrix represent the average number of effective contacts made in the setting that can lead to disease transmission. The overall contact matrix element M<sub>ij</sub> can then be thought of as representing the per capita number of contacts an individual of age i has with individuals of age j.

For more details, please refer to [Mistry, D., Litvinova, M., Pastore y Piontti, A. et al. "Inferring high-resolution human mixing patterns for disease modeling." Nat Commun 12, 323 (2021)](https://www.nature.com/articles/s41467-020-20544-y).

[![DOI](https://zenodo.org/badge/239899144.svg)](https://zenodo.org/badge/latestdoi/239899144)


## List of matrices included in this dataset
1. Australia
    * Australian Capital Territory
    * New South Wales
    * Northern Territory
    * Queensland
    * South Australia
    * Tasmania
    * Victoria
    * Western Australia
2. Austria
3. Bulgaria
4. Canada
    * Alberta
    * British Columbia
    * Manitoba
    * New Brunswick
    * Newfoundland and Labrador
    * Northwest Territories
    * Nova Scotia
    * Nunavut
    * Ontario
    * Prince Edward Island
    * Quebec
    * Saskatchewan
    * Yukon
5. China
    * Chongqing
    * Fujian
    * Gansu
    * Guangdong
    * Guangxi
    * Guizhou
    * Hainan
    * Hebei
    * Heilongjiang
    * Henan
    * Hubei
    * Hunan
    * Inner Mongolia
    * Jiangsu
    * Jiangxi
    * Jilin
    * Liaoning
    * Ningxia
    * Qinghai
    * Shaanxi
    * Shandong
    * Shanghai
    * Shanxi
    * Sichuan
    * Tianjin
    * Tibet
    * Xinjiang
    * Yunnan
    * Zhejiang
6. Cyprus
7. Czech
8. Denmark
9. Estonia
10. Finland
11. France
12. Germany
13. Greece
14. Hungary
15. India
    * Andaman and Nicobar Islands
    * Andhra Pradesh (split into Andhra Pradesh and Telegana as of 2014)
    * Arunachal Pradesh
    * Assam
    * Bihar
    * Chhattisgarh
    * Daman and Diu
    * Goa
    * Gujarat
    * Haryana
    * Himachal Pradesh
    * Jammu and Kashmir
    * Jharkhand
    * Karnataka
    * Kerala
    * Madhya Pradesh
    * Maharashtra
    * Manipur
    * Meghalaya
    * Mizoram
    * Nagaland
    * Nct of Delhi
    * Odisha
    * Puducherry
    * Punjab
    * Rajasthan
    * Sikkim
    * Tamil Nadu
    * Tripura
    * Uttar Pradesh
    * Uttarakhand
    * West Bengal
16. Ireland
17. Israel
18. Italy
19. Japan
    * Aichi-ken
    * Akita-ken
    * Aomori-ken
    * Chiba-ken
    * Ehime-ken
    * Fukui-ken
    * Fukuoka-ken
    * Fukushima-ken
    * Gifu-ken
    * Gumma-ken
    * Hiroshima-ken
    * Hokkaido
    * Hyogo-ken
    * Ibaraki-ken
    * Ishikawa-ken
    * Iwate-ken
    * Kagawa-ken
    * Kagoshima-ken
    * Kanagawa-ken
    * Kochi-ken
    * Kumamoto-ken
    * Kyoto-fu
    * Mie-ken
    * Miyagi-ken
    * Miyazaki-ken
    * Nagano-ken
    * Nagasaki-ken
    * Nara-ken
    * Niigata-ken
    * Oita-ken
    * Okayama-ken
    * Okinawa-ken
    * Osaka-fu
    * Saga-ken
    * Saitama-ken
    * Shiga-ken
    * Shimane-ken
    * Shizuoka-ken
    * Tochigi-ken
    * Tokushima-ken
    * Tokyo-to
    * Tottori-ken
    * Toyama-ken
    * Wakayama-ken
    * Yamagata-ken
    * Yamaguchi-ken
    * Yamanashi-ken
20. Latvia
21. Lithuania
22. Luxembourg
23. Netherlands
24. Norway
25. Portugal
26. Romania
27. Russia
    * Adygea
    * Altai Krai
    * Altai Republic
    * Amur Oblast
    * Arkhangelsk Oblast
    * Astrakhan Oblast
    * Bashkortostan
    * Belgorod Oblast
    * Bryansk Oblast
    * Buryatia
    * Chechnya
    * Chelyabinsk Oblast
    * Chukotka
    * Chuvashia
    * Dagestan
    * Ingushetia
    * Irkutsk Oblast
    * Ivanovo Oblast
    * Jewish Autonomous Oblast
    * Kabardino Balkaria
    * Kaliningrad Oblast
    * Kalmykia
    * Kaluga Oblast
    * Kamchatka Krai
    * Karachay Cherkessia
    * Karelia
    * Kemerovo Oblast
    * Khabarovsk Krai
    * Khakassia
    * Khanty Mansi Autonomous Okrug
    * Kirov Oblast
    * Komi Republic
    * Kostroma Oblast
    * Krasnodar Krai
    * Krasnoyarsk Krai
    * Kurgan Oblast
    * Kursk Oblast
    * Leningrad Oblast
    * Lipetsk Oblast
    * Magadan Oblast
    * Mari El
    * Mordovia
    * Moscow
    * Moscow Oblast
    * Murmansk Oblast
    * Nenets Autonomous Okrug
    * Nizhny Novgorod Oblast
    * North Ossetia Alania
    * Novgorod Oblast
    * Novosibirsk Oblast
    * Omsk Oblast
    * Orenburg Oblast
    * Oryol Oblast
    * Penza Oblast
    * Perm Krai
    * Primorsky Krai
    * Pskov Oblast
    * Rostov Oblast
    * Russian Federation
    * Ryazan Oblast
    * Sakha
    * Sakhalin Oblast
    * Samara Oblast
    * Saratov Oblast
    * Smolensk Oblast
    * St. Petersburg
    * Stavropol Krai
    * Sverdlovsk Oblast
    * Tambov Oblast
    * Tatarstan
    * Tomsk Oblast
    * Tula Oblast
    * Tuva
    * Tver Oblast
    * Tyumen Oblast
    * Udmurtia
    * Ulyanovsk Oblast
    * Vladimir Oblast
    * Volgograd Oblast
    * Vologda Oblast
    * Voronezh Oblast
    * Yamalo Nenets Autonomous Okrug
    * Yaroslavl Oblast
    * Zabaykalsky Krai
28. Slovakia
29. Slovenia
30. South Africa
    * Eastern Cape
    * Free State
    * Gauteng
    * KwaZulu-Natal
    * Limpopo
    * Mpumalanga
    * North West
    * Northern Cape
    * Western Cape
31. Spain
32. Sweden
33. Switzerland
34. United-Kingdom
35. United States
    * Alabama
    * Alaska
    * Arizona
    * Arkansas
    * California
    * Colorado
    * Connecticut
    * Delaware
    * District of Columbia
    * Florida
    * Georgia
    * Hawaii
    * Idaho
    * Illinois
    * Indiana
    * Iowa
    * Kansas
    * Kentucky
    * Louisiana
    * Maine
    * Maryland
    * Massachusetts
    * Michigan
    * Minnesota
    * Mississippi
    * Missouri
    * Montana
    * Nebraska
    * Nevada
    * New Hampshire
    * New Jersey
    * New Mexico
    * New York
    * North Carolina
    * North Dakota
    * Ohio
    * Oklahoma
    * Oregon
    * Pennsylvania
    * Puerto Rico
    * Rhode Island
    * South Carolina
    * South Dakota
    * Tennessee
    * Texas
    * Utah
    * Vermont
    * Virginia
    * Washington
    * West Virginia
    * Wisconsin
    * Wyoming

For some of the figures we utilize shapefiles from the Database of Global Administrative Areas. Please check out their website (https://gadm.org/) to download your own copy of the necessary shapefiles.
