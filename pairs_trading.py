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


def calculate_iqr(data):
    sorted_data = sorted(data)
    n = len(sorted_data)

    def percentile(p):
        pos = p * (n-1)
        lower = int(pos)
        upper = lower + 1

        if upper >= n:
            return sorted_data[lower]

        fractional = pos - lower
        return sorted_data[lower] * (1 - fractional) + sorted_data[upper] * fractional
    
    quartile_one = percentile(0.25)
    quartile_three = percentile(0.75)
    return quartile_three - quartile_one


#calc the spread stats for each year
#ex: we can call calculate_year_stats(2022) (any year) whenever we want
def calculate_year_stats(year):
    google_year_data = google_grouped_by_year.get(year, {})
    meta_year_data = meta_grouped_by_year.get(year, {})
    spread_list = []

#only loop through the dates to get the specific year
#before it didnt work since it looped through the dictionaries
    for date in google_year_data:
        if date in meta_year_data:
            google_price = google_dictionary[date]
            meta_price = meta_dictionary[date]
            spread = float(meta_price) - float(google_price)
            spread_list.append(spread)
            
    if not spread_list:
        return None

#calculations:

#calc range
    range = max(spread_list) - min(spread_list)


#calc mean
    total = 0
    for x in spread_list:
        total += x
        mean = total / len(spread_list)


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

#call calculate_iqr to get iqr
    iqr = calculate_iqr(spread_list)

#access certain stats/store data from the dict later
    return {
        'trading_days': len(spread_list),
        'range': range,
        'mean': mean,
        'variance': variance,
        'standard_dev': standard_dev,
        'min_spread': min(spread_list),
        'max_spread': max(spread_list),
        'IQR': iqr
    }

#print results
#make sure the years match across both data sets
print("Yearly Stock Spread Analysis: \n")

all_years = set(google_grouped_by_year.keys()) & set(meta_grouped_by_year.keys())

for year in sorted(all_years):
    print(f"Year {year}:")
    print("-" * 30)
    
    stats = calculate_year_stats(year)
    
    if stats:
        print(f"Trading Days: {stats['trading_days']}")
        print(f"Range: ${stats['range']:.2f}")
        print(f"Mean Spread: ${stats['mean']:.2f}")
        print(f"Standard Deviation: ${stats['standard_dev']:.2f}")
        print(f"Variance: ${stats['variance']:.2f}")
        print(f"Min Spread: ${stats['min_spread']:.2f}")
        print(f"Max Spread: ${stats['max_spread']:.2f}")
        print(f"IQR: ${stats['IQR']:.2f}")
    else:
        print("No matching data for this year exists")
    
    print()

#prints the overall stats/what we had before
print("Stock Spread Across All Years:")
print("-" * 40)

spread_list = []
for date in google_dictionary:
    if date in meta_dictionary:
        google_price = google_dictionary[date]
        meta_price = meta_dictionary[date]
        spread = float(meta_price) - float(google_price)
        spread_list.append(spread)

range_val = max(spread_list) - min(spread_list)
print("Range:", range_val)

mean = sum(spread_list) / len(spread_list)
print("Mean:", mean)

squared_diffs = [(x - mean) ** 2 for x in spread_list]
variance = sum(squared_diffs) / (len(squared_diffs) - 1)
standard_dev = variance ** 0.5
iqr = calculate_iqr(spread_list)

print("Standard Deviation:", standard_dev)
print("Variance:", variance)
print("IQR:", iqr)