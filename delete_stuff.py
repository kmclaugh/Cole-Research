from math import *
from pylab import *
import datetime
import os
import subprocess
import sys
sys.path.append('/home/kmclaugh/Desktop/Scripts')
from textsender import textsender as textsender
import pickle

current_directory = os.getcwd()
plots_directory = current_directory+"/Plots"
animation_directory = plots_directory +"/Animations"
current_data_directory = current_directory+"/Current_Data"

orbit_filename = "/home/kmclaugh/Desktop/Cole Research/Current_Data/orbit(n=2).dat"
force_x_filename = "/home/kmclaugh/Desktop/Cole Research/Current_Data/Forces_x(n=2).dat"
force_y_filename = "/home/kmclaugh/Desktop/Cole Research/Current_Data/Forces_y(n=2).dat"
semimajor_filename = "/home/kmclaugh/Desktop/Cole Research/Current_Data/smajor(n=2).dat"
eccentricity_filename = "/home/kmclaugh/Desktop/Cole Research/Current_Data/eccen(n=2).dat"

files = [semimajor_filename,eccentricity_filename]
def find_time_value(desired_time,time_list):
    
    min_index = 0
    max_index = int(len(time_list))
    current_index = int((max_index - min_index)/2+min_index)
    condition = False
    print(min_index)
    if time_list[min_index] >= desired_time:
       current_index = min_index
    else:
        while condition == False:
            if min_index+1 == max_index:
                print("find time fail")
                return(False)
            
            if time_list[current_index] >= desired_time and time_list[current_index-1] <= desired_time:
                condition = True
            else:
                if time_list[current_index] > desired_time:
                    max_index = current_index
                    current_index = int((max_index-min_index)/2)+min_index
                elif time_list[current_index] < desired_time:
                    min_index = current_index
                    current_index = int((max_index-min_index)/2)+min_index
    return(current_index)

def scan_file_values(filename,max_count=None):
    
    file = open(filename,"r")
    
    file_values = []
    time_values = []
    max_counter = 0

    for line in file:
        if max_counter % 100 ==0:
            print(max_counter)
        split_values = line.split(",")
        if len(split_values) == 2:
            time_values.append(double(split_values[0]))
            file_values.append(double(split_values[1]))
            max_counter += 1
        if max_counter == max_count:
            break
    file.close()
    print(time_values)
    return([time_values,file_values])

def scan_position_values(orbit_filename,max_count=None):
    
    orbit_file = open(orbit_filename,"r")
    position_x_values = []
    position_y_values = []
    time_values = []
    counter = 0
    for line in orbit_file:
        split_values = line.split(",")
        x_val = double(split_values[1].replace(" ",""))
        position_x_values.append(x_val)
        y_val = double(split_values[2].replace(" ",""))
        position_y_values.append(y_val)
        t_val = double(split_values[0].replace(" ",""))
        time_values.append(t_val)
        counter += 1
        if counter == max_count:
            break
    values = [position_x_values,position_y_values,time_values]
    return(values)

def scan_vector_values(vector_filename,start_index,max_count=None,end_count=None):
    
    vector_file = open(vector_filename,"r")
    
    columbic_values = []
    radiation_values = []
    applied_wave_values = []
    velocity_values = []
    my_count = 0
    max_counter = 0
    end_counter = 0

    for line in vector_file:
        if my_count == start_index:
            split_values = line.split(",")
            if len(split_values) == 4:
                ##columbic
                double_string = split_values[0]
                double_value = double(double_string.replace(" ",""))
                columbic_values.append(double_value)
            
                ##radiation
                double_string = split_values[1]
                double_value = double(double_string.replace(" ",""))
                radiation_values.append(double_value)
            
                ##applied wave
                double_string = split_values[2]
                double_value = double(double_string.replace(" ",""))
                applied_wave_values.append(double_value)
            
                ##velocity
                double_string = split_values[3]
                double_value = double(double_string.replace(" ",""))
                velocity_values.append(double_value)
                max_counter += 1
            if max_counter == max_count:
                break
        else:
            my_count += 1
        end_counter += 1
        if end_counter == end_count:
            break
    vector_file.close()
    return((columbic_values,radiation_values,applied_wave_values,velocity_values))


start_time=(6.48081*pow(10,-11))
end_time=(6.48436*pow(10,-11))

for a_file in files:
    print("scanning", a_file)
    return_values = scan_file_values(a_file)
    time = return_values[0]
    values = return_values[1]
    print("finding time")
    start_index = find_time_value(start_time,time)
    end_index = find_time_value(end_time,time)
    print("making lists")
    time = time[start_index:end_index]
    values = values[start_index:end_index]
    new_file_name = a_file[:-4] + "_test.dat"
    print("making file")
    new_file = open(new_file_name,"w")
    for t,v in zip(time,values):
          print_string = "%f, %f"%(t,v)
          print(print_string,new_file)

print("orbit")
return_values = scan_position_values(orbit_filename)
time = return_values[0]
values1 = return_values[1]
values2 = return_values[2]
start_index = find_time_value(start_time,time)
end_index = find_time_value(end_time,time)
time = time[start_index:end_index]
values1 = values1[start_index:end_index]
values2 = values2[start_index:end_index]
new_file_name = a_file[:-4] + "_test.dat"
new_file = open(new_file_name,"w")
for t,v1,v2 in zip(time,values1,values2):
      print_string = "%f, %f, %f"%(t,v1, v2)
      print(print_string,new_file)
