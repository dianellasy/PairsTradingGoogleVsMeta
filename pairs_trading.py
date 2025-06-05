import csv

with open('google_stocks.csv', newline='') as csvfile1:
    google_stocks_reader = csv.reader(csvfile1)

with open('meta_stocks.csv', newline='') as csvfile2:
    meta_stocks_reader = csv.reader(csvfile2)