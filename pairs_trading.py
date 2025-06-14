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
    """ Calculates Pearson's correlation coefficient between two equally-sized datasets. """

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
    Performs ordinary least squares (OLS) regression to estimate the relationship between x and y
    Returns the intercept (alpha) and slope (beta).
    """

    n = len(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n


    # Calculate the numerator for beta, which is the sum of (xi - mean_x) * (yi - mean_y) for all data points
    numerator = 0

    for xi, yi in zip(x, y):
        numerator += (xi - mean_x) * (yi - mean_y)


    # Calculate the denominator for beta, which is the sum of squared differences (xi - mean_x)^2 for all x values
    denominator = 0

    for xi in x:
        denominator += (xi - mean_x) ** 2


    # Calculate beta (slope)
    # If the denominator is zero, set beta to 0 to avoid division by zero
    beta = numerator / denominator if denominator != 0 else 0


    # Calculate alpha (intercept) using the means of x and y
    alpha = mean_y - beta * mean_x

    return alpha, beta



def compute_cointegration(series_x, series_y):
    """
    Conducts a simplified Engle-Granger cointegration test between two time series
    1. Computes the residuals from the OLS regression of y on x
    2. Computes the first differences of the residuals (residuals_first_difference) and lags (lagged_residuals)
    3. Performs a regression of residuals_first_difference on lagged_residuals to estimate gamma
    4. Computes the t-statistic for the estimated gamma

    Return the t-statistic, gamma, and residual values.
    """

    # Long-run regression: y = α + βx
    n = len(series_x)
    alpha, beta = ols_regression(series_x, series_y)


    # Compute residuals: e_t = y_t − (α + βx_t)
    residuals = []

    for x_value, y_value in zip(series_x, series_y):
        residual = y_value - (alpha + beta * x_value)
        residuals.append(residual)

    # If there are not enough observations, return early
    if n < 2: 
        return None, None, residuals
    

    # Build ADF(1) design: Δe_t and e_{t-1}
    # ADF(1) means one lag of the level term and no extra lagged differences (p = 1)
    residuals_first_difference = [] # Δe_t
    lagged_residuals = []   # e_{t-1}

    for i in range(1, n):
        difference_residual = residuals[i] - residuals[i - 1]
        lagged_residual = residuals[i - 1]
        residuals_first_difference.append(difference_residual)
        lagged_residuals.append(lagged_residual)
    

    # Estimate γ in Δe_t = γ e_{t-1} + ε_t
    # a. Means
    sum_lagged = 0.0

    for lag in lagged_residuals:
        sum_lagged += lag
    
    mean_lag = sum_lagged / len(lagged_residuals)


    sum_first_difference = 0.00

    for difference in residuals_first_difference:
        sum_first_difference += difference
    
    mean_difference = sum_first_difference / len(residuals_first_difference)

    # b. Covariance numerator and variance denominator
    covariance_numerator = 0.0
    variance_denominator = 0.0

    for lag, first_difference in zip(lagged_residuals, residuals_first_difference):
        covariance_numerator += (lag - mean_lag) * (first_difference - mean_difference)
        variance_denominator += (lag - mean_lag) ** 2
    
    gamma = covariance_numerator / variance_denominator if variance_denominator != 0 else 0.0


    # Regression residuals ε_t and their variance
    regression_errors = []
    
    for lag, first_difference in zip(lagged_residuals, residuals_first_difference):
        regression_errors.append(first_difference - gamma * lag)

    sse = 0.0   # Σ ε_t²

    for error in regression_errors:
        sse += error ** 2

    degrees_of_freedom = len(lagged_residuals) - 1

    if degrees_of_freedom <= 0 or variance_denominator == 0:
        standard_error_of_gamma = float('inf')
    else:
        variance_error = sse / degrees_of_freedom
        standard_error_of_gamma = math.sqrt(variance_error) / math.sqrt(variance_denominator)


    # t-statistic for γ
    t_statistic = gamma / standard_error_of_gamma if standard_error_of_gamma != 0 else float('inf')
    
    return t_statistic, gamma, residuals



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