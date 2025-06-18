import csv

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
    
    return covariance / ((variance_x ** 0.5) * (variance_y ** 0.5))



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
        standard_error_of_gamma = (variance_error ** 0.5) / (variance_denominator ** 0.5)


    # t-statistic for γ
    t_statistic = gamma / standard_error_of_gamma if standard_error_of_gamma != 0 else float('inf')
    
    return t_statistic, gamma, residuals



# Calculate yearly statistics for stock spread between META and Google
# The spread is defined as META closing price minus Google closing price for matching dates
# Returns a dictionary of spread statistics
def calculate_year_statistics(year):
    """
    For a specific year, calculate:
    - Spread metrics (range, mean, standard deviation, variance, IQR)
    - Correlation between closing prices
    - Cointegration test (t-statistic and gamma)

    Only dates that are common to both stock datasets are analyzed.
    """

    google_year_data = google_prices_by_year.get(year, {})
    meta_year_data = meta_prices_by_year.get(year, {})
    
    spread_values = []
    google_prices = []
    meta_prices = []


    # Consider only dates that exist in both Google ansd META data for that year
    for date in google_year_data:
        if date in meta_year_data:
            google_price = float(google_closing_prices[date])
            meta_price = float(meta_closing_prices[date])
            spread = abs(meta_price - google_price)
            spread_values.append(spread)
            google_prices.append(google_price)
            meta_prices.append(meta_price)
            
    if not spread_values:
        return None


    # Explicitly find min and max spread values before calculating the range
    min_spread = min(spread_values)
    max_spread = max(spread_values)

    # Calculate range of spread
    range = max_spread - min_spread


    # Calculate mean spread
    total = 0

    for x in spread_values:
        total += x
        mean = total / len(spread_values)


    # Calculate standard deviation and variance 
    squared_differences = []

    for value in spread_values:
        squared_difference = (value - mean) ** 2
        squared_differences.append(squared_difference)

    squared_differences_total = 0
    n = len(squared_differences)

    for squared_difference in squared_differences:
        squared_differences_total += squared_difference
        variance = squared_differences_total / (n - 1)
        
    standard_deviation = variance ** 0.5


    # Calculate IQR
    iqr = calculate_iqr(spread_values)


    # Calculate correlation between Google and META closing prices
    closing_price_correlation = calculate_pearson_correlation(google_prices, meta_prices)


    # Compute cointegration using the Engle-Granger test
    t_statistic, gamma, _ = compute_cointegration(google_prices, meta_prices)


    return {
        'trading_days': len(spread_values),
        'range': range,
        'mean': mean,
        'variance': variance,
        'standard_deviation': standard_deviation,
        'min_spread': min(spread_values),
        'max_spread': max(spread_values),
        'IQR': iqr,
        'correlation': closing_price_correlation,
        'cointegration_t_statistic': t_statistic,
        'cointegration_gamma': gamma
    }



# Calculate overall correlation and cointegration for all matching dates
overall_google_prices = []
overall_meta_prices = []
common_dates = sorted(set(google_closing_prices.keys()) & set(meta_closing_prices.keys()))

for date in common_dates:
    overall_google_prices.append(float(google_closing_prices[date]))
    overall_meta_prices.append(float(meta_closing_prices[date]))

if overall_google_prices:
    overall_correlation = calculate_pearson_correlation(overall_google_prices, overall_meta_prices)
    print("-----------------------------------------------------------------------------")
    print("Overall Pearson Correlation between Google and META Closing Prices: {:.4f}".format(overall_correlation))

    overall_t_statistic, overall_gamma, _ = compute_cointegration(overall_google_prices, overall_meta_prices)
    print("\nEngle-Granger Cointegration Test:")
    print("• t-statistic (for gamma): {:.4f}".format(overall_t_statistic))
    print("• gamma: {:.4f}".format(overall_gamma))

    if overall_t_statistic is not None:
        if overall_t_statistic <= -3.45:
            print("\nEvidence suggests the series are cointegrated.")
        else:
            print("\nEvidence suggests the series are not cointegrated (at 5% level, roughly).")
else:
    print("\nNo matching data found for overall correlation and cointegration tests.")

print("-----------------------------------------------------------------------------")
print("\n\n-----------------------------------------------------------------------------")



# Print yearly stock spread analysis for each common year
print("Yearly Stock Spread Analysis:\n")
common_years = set(google_prices_by_year.keys()) & set(meta_prices_by_year.keys())

for index, year in enumerate(sorted(common_years)):
    print(f"☆ Year {year}:")
    print("---------------------------------------")
    
    statistics = calculate_year_statistics(year)
    
    if statistics:
        print(f"Trading Days: {statistics['trading_days']}")
        print(f"Range: ${statistics['range']:.2f}")
        print(f"Mean Spread: ${statistics['mean']:.2f}")
        print(f"Standard Deviation: ${statistics['standard_deviation']:.2f}")
        print(f"Variance: ${statistics['variance']:.2f}")
        print(f"Min Spread: ${statistics['min_spread']:.2f}")
        print(f"Max Spread: ${statistics['max_spread']:.2f}")
        print(f"IQR: ${statistics['IQR']:.2f}")
        print(f"Correlation (Closing Prices): {statistics['correlation']:.4f}")
        print("Cointegration Test:")
        print("• t-statistic: {:.4f}".format(statistics['cointegration_t_statistic']))
        print("• gamma: {:.4f}".format(statistics['cointegration_gamma']))
    else:
        print("No matching data for this year exists")
    
    if index < len(common_years) - 1:
        print()
    
print("-----------------------------------------------------------------------------")



# Print overall stock spread statistics across all years
# Aggregates the spread differences between META and Google over all matching dates
print("\n\n-----------------------------------------------------------------------------")
print("☆ Stock Spread Across All Years:")
print("-------------------------------------")

all_spreads = []

for date in google_closing_prices:
    if date in meta_closing_prices:
        google_price = google_closing_prices[date]
        meta_price = meta_closing_prices[date]
        spread = abs(float(meta_price) - float(google_price))
        all_spreads.append(spread)

if all_spreads:
    all_spreads_range = max(all_spreads) - min(all_spreads)
    print("Range:", all_spreads_range)

    all_spreads_total = 0

    for x in all_spreads:
        all_spreads_total += x
        all_spreads_mean = all_spreads_total / len(all_spreads)

    print("Mean:", all_spreads_mean)


    all_spreads_squared_differences = []

    for x in all_spreads:
        squared_difference = (x - all_spreads_mean) ** 2
        all_spreads_squared_differences.append(squared_difference)
 

    all_spreads_variance = sum(all_spreads_squared_differences) / (len(all_spreads_squared_differences) - 1)
    print("Variance:", all_spreads_variance)
    

    all_spreads_standard_deviation = all_spreads_variance ** 0.5
    print("Standard Deviation:", all_spreads_standard_deviation)
    

    overall_iqr = calculate_iqr(all_spreads)
    print("IQR:", overall_iqr)

print("-----------------------------------------------------------------------------")