#Date of data from ERA5 file

#IMPORTS
#data format
import datetime

def date_data(File_name):
    #Inputs
        #File_name: File path of data from ERA5
        
    #calculus for data date 
    initial_date= datetime.datetime(1970, 1, 1)
    delta = datetime.timedelta(seconds=int(File_name['valid_time'][-1]))
    date=initial_date+delta
    return date
