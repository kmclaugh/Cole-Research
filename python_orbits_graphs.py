from math import *
from numpy import *
from matplotlib import *
from pylab import *
import mpl_toolkits.mplot3d.axes3d
import datetime
import os

## Convets a file with lines formated as "xvalue, yvalue\n" to a list of x_values and a list of y_values. Returns the tuple (x_values,y_values)
def comma_seperated_file_to_xy_list(filename):
    x_vals = []
    y_vals = []
    the_file = open(filename,"r")
    for line in the_file:
        split = line.split(",")
        try:
            x_value = float(split[0].replace(" ",""))
            y_value = float(split[1].replace(" ",""))
            x_vals.append(x_value)
            y_vals.append(y_value)
        except:
            pass
    return((x_vals,y_vals))

## Converts a file with lines format as "value1, value2, value3..., valuen\n" to a plot with the formating particulars given in the function. Saves the plot in the /Plots directory as "PLOT_TITLE-CURRENT_TIME.pdf" to ensure graphs are not accidently overwritten. NOTE still does not work for more than two values.
def comma_seperated_file_to_graph(filename,the_xlabel="x",the_ylabel="y",the_zlabel=None, the_title="title",marker=None,linestyle="-",color="b"):
    current_date_and_time = datetime.datetime.now()
    the_file = open(filename,"r")
    for line in the_file:
        split_values = line.split(",")
        break
    values = []
    for one_values in split_values:
        values.append([])
    for line in the_file:
        counter = 0
        split_values = line.split(",")
        for string_value in split_values:
            float_value = float(string_value.replace(" ",""))
            values[counter].append(float_value)
            counter += 1
    if len(values) == 2:
        plot(values[0], values[1], marker=marker, linestyle=linestyle, color=color)
        ylabel(the_ylabel)
        xlabel(the_xlabel)
        title(the_title)
        todays_date = datetime.datetime.now()
        savefig("/Users/kevin/Desktop/School/Spring 2013/Cole Research/Plots/{}-{}.pdf".format(the_title,current_date_and_time), bbox_inches=0)
    if len(values) == 3:
        print("more than 3 values doesn't work yet")
#        fig=figure()
#        ax = mpl_toolkits.mplot3d.axes3d.Axes3D(fig)
#        ax.plot_surface(values[0],values[1],values[2])
#        ax.set_xlabel('X')
#        ax.set_ylabel('Y')
#        ax.set_zlabel('Z')
#        show()


### Eccentricity
eccentricity_filename = "/Users/kevin/Desktop/School/Spring 2013/Cole Research/eccen.dat"
comma_seperated_file_to_graph(eccentricity_filename,the_xlabel="x",the_ylabel="eccentricity",the_title="Eccentricity")

### Semimajor
#semimajor_filename = "/Users/kevin/Desktop/School/Spring 2013/Cole Research/smajor.dat"
#comma_seperated_file_to_graph(semimajor_filename,the_xlabel="x",the_ylabel="semimajor axis",the_title="Semimajor Axis")
#
### Theta0
#theta0_filename = "/Users/kevin/Desktop/School/Spring 2013/Cole Research/theta0.dat"
#comma_seperated_file_to_graph(theta0_filename,the_xlabel="x",the_ylabel="theta0",the_title="Theta0")

## Orbit
#orbit_filename = "/Users/kevin/Desktop/School/Spring 2013/Cole Research/orbit.dat"
#comma_seperated_file_to_graph(orbit_filename,the_title="orbit")




