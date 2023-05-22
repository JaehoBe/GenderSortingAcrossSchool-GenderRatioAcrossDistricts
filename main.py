# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/


##################################################
# import modules

import os
import pandas as pd
import requests
# from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from pyproj import CRS
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import seaborn as sns
from collections import Counter
import re

from matplotlib import font_manager, rc
plt.rc('font', family='NanumGothic')
print(plt.rcParams['font.family'])

# Set the max_columns option to None
pd.set_option('display.max_columns', None)

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Malgun Gothic'


##################################################
# set working directory

cwd = os.getcwd()
# print(cwd)
'''
schools that have closed
서울화양초등학교 -> 서울장안초서울성수초공동통학구역
'''
'''
schools that don't have own school district
서울삼각산초등학교-> 서울미양초서울삼각산초공동통학구역
'''
'''
schools that are missing in the school info data
서울반포초등학교, 서울반포초등학교학교군
'''
'''
missing in all data: 초등학교	학구(통학구역) 지정 예외 학교 : 대학 부설초등학교와 사립초등학교
서울대학교사범대학부설초등학교
서울교육대학교부설초등학교
경기초등학교
경복초등학교
경희초등학교
계성초등학교
광운초등학교
금성초등학교
대광초등학교
동광초등학교
동북초등학교
동산초등학교
리라초등학교
매원초등학교
명지초등학교
상명대학교사범대학부속초등학교
상명초등학교
선일초등학교
성동초등학교
성신초등학교
세종초등학교
숭의초등학교
신광초등학교
영훈초등학교
예일초등학교
우촌초등학교
운현초등학교
유석초등학교
은석초등학교
이화여자대학교사범대학부속초등학교
중앙대학교사범대학부속초등학교
청원초등학교
추계초등학교
충암초등학교
태강삼육초등학교
한신초등학교
한양초등학교
홍익대학교사범대학부속초등학교
화랑초등학교
'''

'''
schools that are missing in the processed elementary school data: need to know reasons of missing
경복초등학교
서울교육대학교부설초등학교
서울대학교사범대학부설초등학교
성신초등학교
태강삼육초등학교
'''
'''
Elementary school districts that for more than one school choices: 55 districts
서울강솔초서울강일초공동통학구역
서울경동초서울경수초공동통학구역
서울경동초서울경일초공동통학구역
서울공연초서울용원초공동통학구역
서울광장초서울양진초공동통학구역
서울광희초서울창신초공동통학구역
서울군자초서울답십리초공동통학구역
서울금옥초서울금호초공동통학구역
서울금옥초서울동호초공동통학구역
서울당산초서울선유초공동통학구역
서울당서초서울선유초공동통학구역
서울당서초서울영동초공동통학구역
서울대도초서울도곡초공동통학구역
서울대영초서울우신초공동통학구역
서울동자초서울신양초공동통학구역
서울두산초서울안천초공동통학구역
서울마장초서울사근초공동통학구역
서울매헌초서울언남초공동통학구역
서울면남초서울면동초서울면중초공동통학구역
서울면남초서울중곡초공동통학구역
서울묵현초서울원묵초공동통학구역
서울미양초서울삼각산초공동통학구역
서울반원초서울원촌초공동통학구역
서울방일초서울방현초공동통학구역
서울삼일초서울남성초공동통학구역
서울삼전초서울송전초공동통학구역
서울상계초서울중계초공동통학구역
서울서이초서울서일초공동통학구역
서울서초초서울원명초공동통학구역
서울석계초서울석관초공동통학구역
서울석촌초서울가락초공동통학구역
서울석촌초서울송전초공동통학구역
서울석촌초서울해누리초공동통학구역
서울송례초서울거원초공동통학구역
서울송례초서울위례별초공동통학구역
서울수암초서울을지초공동통학구역
서울숭신초서울무학초공동통학구역
서울신답초서울전농초공동통학구역
서울신학초서울창경초공동통학구역
서울신학초서울초당초공동통학구역
서울양남초서울성자초공동통학구역
서울염경초서울염동초공동통학구역
서울염리초서울용강초공동통학구역
서울원묵초서울공릉초공동통학구역
서울위례별초서울거원초공동통학구역
서울위례솔초서울거원초공동통학구역
서울은로초서울흑석초공동통학구역
서울장안초서울성수초공동통학구역
서울전곡초서울전농초공동통학구역
서울중광초서울중마초공동통학구역
서울중동초서울창서초공동통학구역
서울천왕초서울하늘숲초공동통학구역
서울태랑초서울태릉초공동통학구역
서울토성초서울풍성초공동통학구역
서울화곡초서울화일초공동통학구역
'''


base_path = "/Users/USER/PycharmProjects/genderSortingAcrossElementarySchoolInKorea"
path_engineered_data = os.path.join(base_path, r'engineered_data')

if not os.path.exists(path_engineered_data):
   os.makedirs(path_engineered_data)


##################################################
# read middle school district shp file and check distribution
file_name = "data/middleSchoolDistrict/중학교학교군.shp"
file_path = os.path.join(base_path, file_name)
shapefileMiddle = gpd.read_file(file_path)
# print(shapefileMiddle.columns)

columns_to_drop = ['CRE_DT', 'UPD_DT', 'BASE_DT']
shapefileMiddle_subset = shapefileMiddle.drop(columns_to_drop, axis=1)
# print(shapefileMiddle_subset.columns)
subsetMiddle = shapefileMiddle_subset[shapefileMiddle_subset['EDU_UP_NM'] == "서울특별시교육청"]
# print(subsetMiddle.shape)
# print(len(subsetMiddle['HAKGUDO_NM']))


##################################################
# merge shapefile with point data set(merged_df, or "merged_df_with_coordinates.csv")
data_types = {
    '위도': float,
    '경도': float,
    'longitude': float,
    'latitude': float,
    'geometry': object,
    'OBJECTID': int,
    'HAKGUDO_ID': int
}

# 2022년도 초등학생 정보 자료
file_path = os.path.join(base_path + '/engineered_data', 'merged_df_with_coordinates.csv')
merged_df = pd.read_csv(file_path, encoding='utf-8', dtype = data_types)
merged_df_with_coordinates = merged_df

# convert merged_df_with_coordinates to a GeoDataFrame
points = gpd.GeoDataFrame(merged_df_with_coordinates,
                          geometry=gpd.points_from_xy(merged_df_with_coordinates.longitude,
                                                      merged_df_with_coordinates.latitude))

# infer CRS from latitude and longitude columns using pyproj
points.crs = CRS.from_user_input('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs').to_wkt()

# convert subset to a GeoDataFrame
polygons = gpd.GeoDataFrame(subsetMiddle)

# reproject points to match the CRS of polygons
points = points.to_crs(polygons.crs)

# perform spatial join
joined = gpd.sjoin(points, polygons, op='within')
# print(joined)
# joined_columns = list(joined.columns)
# for col in joined_columns:
#     print(col)


# Group the DataFrame by 'HAKGUDO_NM' and calculate the sum of '계(남)' and '계(여)'
sum_by_district = joined.groupby('HAKGUDO_NM')['계(남)', '계(여)'].sum()

# Print the resulting DataFrame
# print(sum_by_district)

# Merge the two DataFrames based on 'HAKGUDO_NM'
merged = subsetMiddle.merge(sum_by_district, on='HAKGUDO_NM', how='left')

merged['total'] = merged['계(남)']+merged['계(여)']

merged['boys_ratio'] = merged['계(남)'] / merged['total'] * 100
# print(merged)


# Visualization: Choropleth map
fig, ax = plt.subplots(figsize=(12, 8))
merged.plot(column='boys_ratio', cmap='coolwarm', linewidth=0.8, ax=ax, edgecolor='0.8', legend=True)

# Set plot title and axis labels
ax.set_title("Boys' Ratio Across Districts", fontsize=16)
ax.set_xlabel("Longitude", fontsize=12)
ax.set_ylabel("Latitude", fontsize=12)

# Show the plot
# plt.show()

# ##################################################
# number of birth, boy/girl ratio of korea
data = {
    'year': ['2010', '2011', '2012', '2013', '2014',
            '2015', '2016', '2017', '2018', '2019', '2020'],
    'sKoreaBirthTotal': [470171, 471265, 484550, 436455, 435435,
                   438420, 406243, 357771, 326822, 302676, 272337],
    'sKoreaBoyGirlRatio': [106.9, 105.7, 105.7, 105.3, 105.3,
                   105.3, 105.0, 106.3, 105.4, 105.5, 104.8]
}

korean_births_df = pd.DataFrame(data)
print(korean_births_df)

# korean_births_df['sKoreaBoy'] = korean_births_df['sKoreaBoyGirlRatio']/100 * korean_births_df['sKoreaBirthTotal']
# print(korean_births_df['sKoreaBoy'])

korean_births_df['sKoreaBoyGirlRatio'] = korean_births_df['sKoreaBoyGirlRatio'] / (100 + korean_births_df['sKoreaBoyGirlRatio']) * 100
korean_births_df['sKoreaBirthBoy'] = (korean_births_df['sKoreaBirthTotal'] / (1 + 100/ korean_births_df['sKoreaBoyGirlRatio'])).astype(int)
korean_births_df['sKoreaBirthGirl'] = korean_births_df['sKoreaBirthTotal'] - korean_births_df['sKoreaBirthBoy']
print(korean_births_df['sKoreaBirthBoy'])
print(korean_births_df['sKoreaBirthGirl'])
# print(korean_births_df['sKoreaBoyGirlRatio'])

# Filter the DataFrame for the years 2010 to 2015
filtered_df = korean_births_df[(korean_births_df['year'] >= '2010') & (korean_births_df['year'] <= '2015')]

# Calculate the sum of 'sKoreaBirthTotal'
birth_total_sum = filtered_df['sKoreaBirthTotal'].sum()

print("Sum of sKoreaBirthTotal from 2010 to 2015:", birth_total_sum)


##################################################
# open files: seoul birth data
file_name = "data/seoulBirthData/출산순위별+출생_20230508155221.csv"
file_path = os.path.join(base_path, file_name)
df = pd.read_csv(file_path, encoding='utf-8', header=[0,1,2])

# print(df.head())

# Split the data into df1 and df2
df1 = df.iloc[:, :2]
df2 = df.iloc[:, 2:]

df2.columns = df2.columns.droplevel(level=1)

# Translate lower level headers
df2 = df2.rename(columns={'계': 'total', '남자': 'male', '여자': 'female'})

df2 = df2.drop(index=0)

# Get the upper header
upper_header = df2.columns.get_level_values(0)

# Get the lower header
lower_header = df2.columns.get_level_values(1)

# Get the year information
year_info = upper_header.str.split(" ").str[0]

# Add the year information to the lower header
new_lower_header = year_info + lower_header

# Assign the new header to the DataFrame
df2.columns = [upper_header, new_lower_header]

# drop upper header
df2.columns = df2.columns.droplevel(0)

# drop the first row
df1 = df1.iloc[1:]

# drop the first column
df1 = df1.iloc[:, 1:]
df1.columns = ['bySeoulDistrict']
df2 = df2.set_index(df1['bySeoulDistrict'])
df2.loc['SeoulTotal'] = df2.sum(axis=0)


total_cols = df2.filter(like='total').columns
df2_total = df2[total_cols]
seoul_total = df2.loc['SeoulTotal', total_cols]

male_cols = df2.filter(regex='^(?!.*female).*male$').columns
df2_male = df2[male_cols]
seoul_male = df2.loc['SeoulTotal', male_cols]

female_cols = df2.filter(like='female').columns
df2_female = df2[female_cols]
seoul_female = df2.loc['SeoulTotal', female_cols]

# reset index of each series
total = seoul_total.reset_index(drop=True)
male = seoul_male.reset_index(drop=True)
female = seoul_female.reset_index(drop=True)

# concatenate horizontally
dfSeoulBirth = pd.concat([total, male, female], axis=1)
dfSeoulBirth['year'] = range(2000, 2022)

dfSeoulBirth.columns = ['SeoulTotal', 'SeoulMaleTotal', 'SeoulFemaleTotal', 'year']

dfSeoulBirth = dfSeoulBirth.reindex(columns=['year', 'SeoulBirthTotal', 'SeoulBirthMale', 'SeoulBirthFemale'])

# convert the "year" column to a string type
dfSeoulBirth['year'] = dfSeoulBirth['year'].astype(str)

# remove the "년" string from the "year" column
dfSeoulBirth['year'] = dfSeoulBirth['year'].str.replace('년', '')

print(dfSeoulBirth)


# one-sample t-test:
from scipy.stats import ttest_1samp, kruskal

# Assuming you have the merged DataFrame 'merged' containing the required data
# and 'boys_ratio' column representing the boys' ratio for each polygon

# Perform one-sample t-test
t_statistic, p_value = ttest_1samp(merged['boys_ratio'], popmean=desired_mean)

# Print the t-test result
print("One-sample t-test p-value:", p_value)

# Kruskal-Wallis test
statistic, p_value = kruskal(*boys_ratio_by_polygon)

# Print the Kruskal-Wallis test result
print("Kruskal-Wallis p-value:", p_value)
