import csv

google_dictionary = {}
date = "Date"
close = "Close"

with open('google_stocks.csv', newline='') as csvfile1:
    google_stocks_reader = csv.reader(csvfile1)

    for row in google_stocks_reader:
        if row and row[0] == date:
            google_dictionary_date = date
            print(f"Row found where the first column is {google_dictionary_date}: {close}")
        
        if row and row[4] == close:
            google_dictionary_close = close

        google_dictionary[google_dictionary_date] = google_dictionary_close

with open('meta_stocks.csv', newline='') as csvfile2:
    meta_stocks_reader = csv.reader(csvfile2)