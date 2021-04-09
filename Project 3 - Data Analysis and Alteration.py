#-----------------------------------------------------------------------#
# SCRIPT: projectc_template.py						#
# PURPOSE: Template for project C					#
# S. Miller, March 2, 2021						#
#									#
#-----------------------------------------------------------------------#

#----------#
# Notes:   #
#----------#


#----------------------#
# Required libraries   #
#----------------------#

import csv
import math


#------------------#
# Main function    #
#------------------#

def main():

    filename = 'NOAA_data.csv'

    # Run the read_data function
    data = read_data(filename=filename)

    # Find the list of unique site codes
    usites = unique_sites(sitecode=data[0])
    print("Unique site codes")
    print(usites)

    # Calculate the data counts
    counts = data_counts(filename)
    print(counts)
    # Calculate data statistics
    summarystats = data_stats(filename)

    # Add a new data point to the data file
    add_new_data(filename)


def transpose(A): #switch the rows and columns of the input matrix
    transpose = []

    for i in range(len(A[0])):
        transpose.append([])

    for i in range(len(A[0])):
        for x in range(len(A)):
            transpose[i].append(A[x][i])
            
    return transpose

#----------#
# PART A   #
#----------#

# Purpose: read in the data and clean the data

def read_data(filename):
    sitecodes = {'wisconsin':'LEF', 'wisc':'LEF', 'waco':'WKT', 'maryland':'TMD', 'sutro':'STR', 'oklahoma':'SGP'} #dictionary of sitecodes
    
    with open('NOAA_data.csv', 'r') as csvfile: #open NOAA data to be read
        NOAA_data = csv.reader(csvfile, delimiter=',') #set up reader object
        NOAA = [] 
        for row in NOAA_data: #transfer csv file data to the list NOAA
            NOAA.append(row)
        
        for row in range(len(NOAA)): #check through each item in the lists and replace any descriptive sitecodes
            for char in range(len(NOAA[0])):
                if NOAA[row][char] in sitecodes:
                    NOAA[row][char] = sitecodes[NOAA[row][char]]
        
        NOAAT = transpose(NOAA)
        
        for row in range(len(NOAAT)): #delete the original column headers
            del NOAAT[row][0]
        
        for char in range(1, 3): #convert year and doy values to int
            for x in range(len(NOAAT[0])):
                NOAAT[char][x] = int(NOAAT[char][x])
        
        for char in range(len(NOAAT[3])): #convert obs values to float
            NOAAT[3][char] = float(NOAAT[3][char])
            
        return NOAAT


#----------#
# PART B   #
#----------#

# Purpose: Create a unique list of sites

def unique_sites(sitecode):
    uniques = []
    for site in sitecode:
        if site not in uniques: #append the sitecodes to uniques if and only if it is not already in to avoid duplicates
            uniques.append(site)
    uniques.sort() #sort in alphabetical order
    return uniques


#----------#
# PART C   #
#----------#

def data_counts(filename):
    data = read_data(filename)
    usites = unique_sites(data[0])
    counts = []
    for i in range(len(usites)): #set up list of 4 int elements, each corresponding to a unique sitecode
        counts.append(0)
    for x in data[0]:
        if x in usites:
            counts[usites.index(x)] += 1 #only increase the int that corresponds to the sitecode of the current iteration
    largest = max(counts)
    smallest = min(counts)
    print('Site {} has the largest number of observations with {} observations'.format(usites[counts.index(largest)], largest))
    print('Site {} has the smallest number of observations with {} observations'.format(usites[counts.index(smallest)], smallest))
    return counts


#----------#
# PART D   #
#----------#

# Calculate statistics on each site and write the data to a file as a visually-pleasing table.

def data_stats(filename):
    data = read_data(filename)
    usites = unique_sites(data[0])
    
    counts = []
    for i in range(len(usites)): #find the number of data points at each site
        counts.append(0)
    for x in data[0]:
        if x in usites:
            counts[usites.index(x)] += 1
    
    raw_stats = []
    for i in range(len(usites)):
        raw_stats.append([]) #give raw_stats as many empty sets as there are unique sites
    
    for x in range(len(data[0])): #add data values to their list elements according to sitecode
        if data[0][x] == usites[0]:
            raw_stats[0].append(data[3][x])
        elif data[0][x] == usites[1]:
            raw_stats[1].append(data[3][x])
        elif data[0][x] == usites[2]:
            raw_stats[2].append(data[3][x])
        elif data[0][x] == usites[3]:
            raw_stats[3].append(data[3][x])
        
    for row in range(len(raw_stats)): #convert all values to float
        for char in range(len(raw_stats[row])):
            raw_stats[row][char] = float(raw_stats[row][char])
    
    stats = [[], [], []] #min, max, mean
    for x in raw_stats:
        stats[0].append(min(x))
        stats[1].append(max(x))
        stats[2].append((sum(x) / len(x)))
 
    format_string = '|{site:^8}|{mean:^8}|{mini:^8}|{maxi:^8}|\n' #set up field
    output_list = [format_string.format(site = 'site', mean = 'mean', mini = 'min', maxi = 'max')]
    for x in range(4): #printing out the rows of data
        output_list.append(format_string.format(site = usites[x], mean = round(stats[2][x], 2), mini = stats[0][x], maxi = stats[1][x]))
    output = ''
    with open('data_summary.txt','w') as textFile:
        textFile.write(output.join(output_list))
    return stats
        

#----------#
# PART E   #
#----------#

# Add a new data point to the data file

def add_new_data(filename):
    
    try:
        site = input('Site code of new data:')
        yr = int(input('Year of new data:')) #will raise ValueError if a non int value is input
    except ValueError:
        print('Error: year must be an integer.')
        return
        
    try:
        doty = int(input('Day of new data:'))#will raise ValueError if a non int value is input
    except ValueError:
        print('Error: day of year must be an integer.')
        return
    
    try:
        obs = float(input('Observation value of new data:'))#will raise ValueError if a non float value is input
    except ValueError:
        print('Error: observation entered is not numeric.')
        return
    data = read_data(filename)
    #add back the column and new data values    
    data[0].append(site)
    data[0].insert(0, 'sitecode')
    data[1].append(str(yr))
    data[1].insert(0, 'yr')
    data[2].append(str(doty))
    data[2].insert(0, 'doy')
    data[3].append(str(obs))
    data[3].insert(0, 'obs')
        
    output_data = transpose(data)
    for x in range(len(output_data)):
        for i in range(len(output_data[0])):
            output_data[x][i] = str(output_data[x][i])
    
    with open('NOAA_data_new.csv', mode = 'w') as textFile: #write new csv file
        dataFile = csv.writer(textFile)
        for x in range(len(output_data)):
            dataFile.writerow(output_data[x])


#----------------#
# Run the code   #
#----------------#

if __name__ == "__main__": main()


#-----------------------------------------------------------------------#
# END OF SCRIPT  							#
#-----------------------------------------------------------------------#


