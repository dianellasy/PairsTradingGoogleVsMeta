import csv

google_dictionary = {}
meta_dictionary = {}

with open('google_stocks.csv', 'r') as csvfile1:
    google_stocks_reader = csv.reader(csvfile1)

    for index, row in enumerate(google_stocks_reader):
        if index == 0:
            continue

        google_dictionary_date = row[0]
        google_dictionary_close = row[4]
        google_dictionary[google_dictionary_date] = google_dictionary_close


with open('meta_stocks.csv', newline='') as csvfile2:
    meta_stocks_reader = csv.reader(csvfile2)

    for index, row in enumerate(meta_stocks_reader):
        if index == 0:
            continue

        meta_dictionary_date = row[0]
        meta_dictionary_close = row[4]
        meta_dictionary[meta_dictionary_date] = meta_dictionary_close

print(google_dictionary)


#merge data
#check between data sets for the same date
#we only compute the spread dates that are in both data sets/matching days
#calc spread/range by finding difference
spread_list = []

for date in google_dictionary:
    if date in meta_dictionary:
        google_price = google_dictionary[date]
        meta_price = meta_dictionary[date]
        spread = float(google_price) - float(meta_price)
        spread_list.append(spread)
        
print(spread_list)