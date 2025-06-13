import csv
import math

# Initialize dictionaries for storing closing price data for each stock
google_dictionary = {}
meta_dictionary = {}


# Read Google and Meta stock data from CSV file
# Each row (except the header) is parsed to extract the date (column 0) and closing price (column 4), which are then stored in a dictionary
with open('GOOGL stocks.csv', 'r') as csv_file_one:
    google_stocks_reader = csv.reader(csv_file_one)

    for index, row in enumerate(google_stocks_reader):
        if index == 0:
            continue

        google_dictionary_date = row[0]
        google_dictionary_close = row[4]
        google_dictionary[google_dictionary_date] = google_dictionary_close


with open('META stocks.csv', 'r') as csv_file_two:
    meta_stocks_reader = csv.reader(csv_file_two)

    for index, row in enumerate(meta_stocks_reader):
        if index == 0:
            continue

        meta_dictionary_date = row[0]
        meta_dictionary_close = row[4]
        meta_dictionary[meta_dictionary_date] = meta_dictionary_close



# Merge data by grouping closing prices by year, which is done separately for Google and META
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



# Define helper functions for statistical calculations
def calculate_iqr(data):
    """
    Calculate the interquartile range of the input data
    Sort the data and interpolates the 25th and 75th percentiles.
    """
    
    sorted_data = sorted(data)
    n = len(sorted_data)

    def percentile(p):
        pos = p * (n - 1)
        lower = int(pos)
        upper = lower + 1

        if upper >= n:
            return sorted_data[lower]

        fractional = pos - lower
        return sorted_data[lower] * (1 - fractional) + sorted_data[upper] * fractional
    
    quartile_one = percentile(0.25)
    quartile_three = percentile(0.75)
    return quartile_three - quartile_one



def calculate_correlation(data_x, data_y):
    """
    Calculate Pearson's correlation coefficient between two equally-sized datasets.
    """

    n = len(data_x)
    mean_x = sum(data_x) / n
    mean_y = sum(data_y) / n

    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(data_x, data_y)) / (n - 1)
    var_x = sum((x - mean_x) ** 2 for x in data_x) / (n - 1)
    var_y = sum((y - mean_y) ** 2 for y in data_y) / (n - 1)

    if var_x == 0 or var_y == 0:
        return 0
    
    return cov / (math.sqrt(var_x) * math.sqrt(var_y))



def ols_regression(x, y):
    """
    Perform ordinary least squares (OLS) regression to estimate the relationship between x and y
    Return the intercept (alpha) and slope (beta).
    """
    n = len(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    numerator = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    denominator = sum((xi - mean_x) ** 2 for xi in x)
    beta = numerator / denominator if denominator != 0 else 0
    alpha = mean_y - beta * mean_x
    return alpha, beta



def compute_cointegration(x, y):
    """
    Conduct a simplified Engle-Granger cointegration test between two time series
    1. Computes the residuals from the OLS regression of y on x
    2. Computes the first differences of the residuals (de) and lags (e_lag)
    3. Performs a regression of de on e_lag to estimate gamma
    4. Computes the t-statistic for the estimated gamma

    Return the t-statistic, gamma, and residual values.
    """
    n = len(x)

    alpha, beta = ols_regression(x, y)
    residuals = [yi - (alpha + beta * xi) for xi, yi in zip(x, y)]

    if n < 2: 
        return None, None, residuals
    
    de = []
    e_lag = []

    for i in range(1, n):
        de.append(residuals[i] - residuals[i - 1])
        e_lag.append(residuals[i - 1])
    
    m_e = sum(e_lag) / len(e_lag)
    m_de = sum(de) / len(de)

    num = sum((e - m_e) * (d - m_de) for e, d in zip(e_lag, de))
    den = sum((e - m_e) ** 2 for e in e_lag)

    if den == 0:
        gamma = 0

    else:
        gamma = num / den

    reg_errors = [d - gamma * e for d, e in zip(de, e_lag)]
    sse = sum(err ** 2 for err in reg_errors)

    df = len(e_lag) - 1

    if df <= 0:
        stderr = float('inf')
    else:
        variance_err = sse / df
        stderr = math.sqrt(variance_err) / math.sqrt(den) if den != 0 else float('inf')
    
    if stderr == 0:
        t_stat = float('inf')
    else:
        t_stat = gamma / stderr
    
    return t_stat, gamma, residuals



# Calculate yearly statistics for stock spread between META and Google
# The spread is defined as META closing price minus Google closing price for matching dates
# Returns a dictionary of spread statistics
def calculate_year_stats(year):
    """
    For a specific year, calculate:
    - Spread metrics (range, mean, standard deviation, variance, IQR)
    - Correlation between closing prices
    - Cointegration test (t-statistic and gamma)

    Only dates that are common to both stock datasets are analyzed.
    """

    google_year_data = google_grouped_by_year.get(year, {})
    meta_year_data = meta_grouped_by_year.get(year, {})
    
    spread_list = []
    google_prices = []
    meta_prices = []

    # Only consider dates present in both datasets
    for date in google_year_data:
        if date in meta_year_data:
            g_price = float(google_dictionary[date])
            m_price = float(meta_dictionary[date])
            spread = m_price - g_price
            spread_list.append(spread)
            google_prices.append(g_price)
            meta_prices.append(m_price)
            
    if not spread_list:
        return None


    # Calculate range of spread
    range = max(spread_list) - min(spread_list)


    # Calculate mean spread
    total = 0
    for x in spread_list:
        total += x
        mean = total / len(spread_list)


    # Calculate standard deviation and variance from the spread data
    squared_diffs = [(x - mean) ** 2 for x in spread_list]

    squared_diffs_total = 0
    n = len(squared_diffs)

    for x in squared_diffs:
        squared_diffs_total += x
        variance = squared_diffs_total / (n - 1)
        
    standard_dev = variance ** 0.5


    # Get IQR using the helper function
    iqr = calculate_iqr(spread_list)


    # Calculate correlation between Google and META closing prices
    correlation = calculate_correlation(google_prices, meta_prices)


    # Compute cointegration using the Engle-Granger test
    t_stat, gamma, _ = compute_cointegration(google_prices, meta_prices)


    return {
        'trading_days': len(spread_list),
        'range': range,
        'mean': mean,
        'variance': variance,
        'standard_dev': standard_dev,
        'min_spread': min(spread_list),
        'max_spread': max(spread_list),
        'IQR': iqr,
        'correlation': correlation,
        'cointegration_t_stat': t_stat,
        'cointegration_gamma': gamma
    }



# Calculate overall correlation and cointegration for all matching dates
overall_google_prices = []
overall_meta_prices = []
dates_common = sorted(set(google_dictionary.keys()) & set(meta_dictionary.keys()))

for date in dates_common:
    overall_google_prices.append(float(google_dictionary[date]))
    overall_meta_prices.append(float(meta_dictionary[date]))

if overall_google_prices:
    overall_corr = calculate_correlation(overall_google_prices, overall_meta_prices)
    print("Overall Pearson Correlation between Google and META closing prices: {:.4f}".format(overall_corr))

    coint_t_stat, coint_gamma, _ = compute_cointegration(overall_google_prices, overall_meta_prices)
    print("\nManual Engle-Granger Cointegration Test (Simplified):")
    print("   t-statistic (for gamma): {:.4f}".format(coint_t_stat))
    print("   gamma: {:.4f}".format(coint_gamma))

    if coint_t_stat is not None:
        if coint_t_stat <= -3.45:
            print("   -> Evidence suggests the series are cointegrated.")
        else:
            print("   -> Evidence suggests the series are not cointegrated (at 5% level, roughly).")
else:
    print("No matching data found for overall correlation and cointegration tests.")



# Print Yearly Stock Spread Analysis
# For each common year between both datasets, print the stats
print("\nYearly Stock Spread Analysis:\n")

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
        print(f"Correlation (Closing Prices): {stats['correlation']:.4f}")
        print("Cointegration Test:")
        print("   t-statistic: {:.4f}".format(stats['cointegration_t_stat']))
        print("   gamma: {:.4f}".format(stats['cointegration_gamma']))

    else:
        print("No matching data for this year exists")
    
    print()



# Print overall stock spread statistics across all years
# Aggregates the spread between META and Google over all matching dates
print("Stock Spread Across All Years:")
print("-" * 40)

spread_list = []
for date in google_dictionary:
    if date in meta_dictionary:
        google_price = google_dictionary[date]
        meta_price = meta_dictionary[date]
        spread = float(meta_price) - float(google_price)
        spread_list.append(spread)

if spread_list:
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