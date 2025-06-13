import csv
import math

# Initialize dictionaries for storing closing price data for each stock
google_closing_prices = {}
meta_closing_prices = {}


# Read Google and Meta stock data from CSV file
# Each row (except the header) is parsed to extract the date (column 0) and closing price (column 4), which are then stored in a dictionary
with open('GOOGL stocks.csv', 'r') as google_csv_file:
    google_stocks_reader = csv.reader(google_csv_file)

    for row, column in enumerate(google_stocks_reader):
        if row == 0:
            continue

        date = column[0]
        closing_price = column[4]
        google_closing_prices[date] = closing_price


with open('META stocks.csv', 'r') as meta_csv_file:
    meta_stocks_reader = csv.reader(meta_csv_file)

    for row, column in enumerate(meta_stocks_reader):
        if row == 0:
            continue

        date = column[0]
        closing_price = column[4]
        meta_closing_prices[date] = closing_price



# Merge data by grouping closing prices by year, which is done separately for Google and META
google_prices_by_year = {}
meta_prices_by_year = {}

for date, close in google_closing_prices.items():
    year = date.split('-')[0]

    if year not in google_prices_by_year:
        google_prices_by_year[year] = {}
    
    google_prices_by_year[year][date] = close


for date, close in meta_closing_prices.items():
    year = date.split('-')[0]

    if year not in meta_prices_by_year:
        meta_prices_by_year[year] = {}
    
    meta_prices_by_year[year][date] = close



# Helper functions for statistical calculations
def calculate_iqr(data):
    """
    Calculates the interquartile range (IQR) of the input data
    Sorts the data and interpolates the 25th and 75th percentiles.
    """

    # Sort the data in ascending order
    sorted_data = sorted(data)
    n = len(sorted_data)


    # A helper function to compute a percentile value
    # Calculates the position in the sorted data based on the desired percentile and interpolates if necessary
    def percentile(p):
        pos = p * (n - 1)
        lower_index = int(pos)
        upper_index = lower_index + 1

        # If computed upper index exceeds list bounds, return the last available value
        if upper_index >= n:
            return sorted_data[lower_index]

        fractional = pos - lower_index
        
        # Return interpolated value based on lower and upper indices
        return sorted_data[lower_index] * (1 - fractional) + sorted_data[upper_index] * fractional
    
    # Calculate the 25th and 75th percentiles (Q1 and Q3)
    quartile_one = percentile(0.25)
    quartile_three = percentile(0.75)

    # Compute and return the interquartile range (IQR) as the difference between Q3 and Q1
    return quartile_three - quartile_one



def calculate_pearson_correlation(series_x, series_y):
    """ Calculate Pearson's correlation coefficient between two equally-sized datasets. """

    n = len(series_x)
    mean_x = sum(series_x) / n
    mean_y = sum(series_y) / n


    # Calculate covariance
    covariance_sum = 0

    for x, y in zip(series_x, series_y):
        covariance_sum += (x - mean_x) * (y - mean_y)
    
    covariance = covariance_sum / (n - 1)

    # Calculate variance for series_x
    variance_x_sum = 0

    for x in series_x:
        variance_x_sum += (x - mean_x) ** 2
    
    variance_x = variance_x_sum / (n - 1)


    # Calculate variance for series_y
    variance_y_sum = 0

    for y in series_y:
        variance_y_sum += (y - mean_y) ** 2
    
    variance_y = variance_y_sum / (n - 1)


    # Calculate Pearson's correlation coefficient
    if variance_x == 0 or variance_y == 0:
        return 0
    
    return covariance / (math.sqrt(variance_x) * math.sqrt(variance_y))



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

    google_year_data = google_prices_by_year.get(year, {})
    meta_year_data = meta_prices_by_year.get(year, {})
    
    spread_list = []
    google_prices = []
    meta_prices = []

    # Only consider dates present in both datasets
    for date in google_year_data:
        if date in meta_year_data:
            g_price = float(google_closing_prices[date])
            m_price = float(meta_closing_prices[date])
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
    correlation = calculate_pearson_correlation(google_prices, meta_prices)


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
dates_common = sorted(set(google_closing_prices.keys()) & set(meta_closing_prices.keys()))

for date in dates_common:
    overall_google_prices.append(float(google_closing_prices[date]))
    overall_meta_prices.append(float(meta_closing_prices[date]))

if overall_google_prices:
    overall_corr = calculate_pearson_correlation(overall_google_prices, overall_meta_prices)
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

all_years = set(google_prices_by_year.keys()) & set(meta_prices_by_year.keys())

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
for date in google_closing_prices:
    if date in meta_closing_prices:
        google_price = google_closing_prices[date]
        meta_price = meta_closing_prices[date]
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