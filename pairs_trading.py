import csv

google_dictionary = {}
meta_dictionary = {}

with open('GOOGL stocks.csv', 'r') as csvfile1:
    google_stocks_reader = csv.reader(csvfile1)

    for index, row in enumerate(google_stocks_reader):
        if index == 0:
            continue

        google_dictionary_date = row[0]
        google_dictionary_close = row[4]
        google_dictionary[google_dictionary_date] = google_dictionary_close


with open('META stocks.csv', newline='') as csvfile2:
    meta_stocks_reader = csv.reader(csvfile2)

    for index, row in enumerate(meta_stocks_reader):
        if index == 0:
            continue

        meta_dictionary_date = row[0]
        meta_dictionary_close = row[4]
        meta_dictionary[meta_dictionary_date] = meta_dictionary_close


#merge data
#check between data sets for the same date
#we only compute the spread dates that are in both data sets/matching days
#calc spread by finding difference
spread_list = []

google_grouped_by_year = {}
meta_grouped_by_year = {}


for date, close in google_dictionary.items():
    year = date.split('-')[0]

    if year not in google_grouped_by_year:
        google_grouped_by_year[year] = {}
    
    google_grouped_by_year[year][date] = close


for date, close in meta_dictionary.items():
    year = date.split('-')[0]

    if year not in meta_grouped_by_year:
        meta_grouped_by_year[year] = {}
    
    meta_grouped_by_year[year][date] = close



for date in google_dictionary:
    if date in meta_dictionary:
        google_price = google_dictionary[date]
        meta_price = meta_dictionary[date]
        spread = float(meta_price) - float(google_price)
        spread_list.append(spread)



#calc range
range = max(spread_list) - min(spread_list)
print("Range:", range)

#calc mean
total = 0
for x in spread_list:
    total += x
    mean = total / len(spread_list)
print("Mean:", mean)

#calc standard deviation:
#-get mean
#-square differences from the mean
#-get variance (mean of the squared differences minus 1 for length)
#-square root the variance

squared_diffs = [(x - mean) ** 2 for x in spread_list]

squared_diffs_total = 0
n = len(squared_diffs)

for x in squared_diffs:
    squared_diffs_total += x
    variance = squared_diffs_total / (n - 1)
    
standard_dev = variance ** 0.5

print("Standard Deviation:", standard_dev)
print("Variance:", variance)