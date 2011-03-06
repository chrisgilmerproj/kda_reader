#! /opt/local/bin/python

"""
    KDA File Reader (Python)
    
    Author: Chris Gilmer, Megan Marsh
    Last Modified: June 26th, 2010

    Copyright (c) 2010 All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

#--- Import the necessary libraries for this file
import glob     # used for loading files in a directory
import numpy    # used for loading text and the trapezoid rule
import optparse # used for parsing options
import os       # used for file operations
import pylab    # used for plotting and math functions
import re       # used for regular expression matching
import sys      # used for exiting program

#--- List of file extensions
FILE_EXT_LIST = ['*.KDA','*.kda']

#--- Set up the conversion factors for each column type
#    These are specific the the equipment being used
#    The labeling is (Axis)(conversion factor)(Reading), so for
#    Xc1 read "X-axis conversion factor 1"
#    Units: N/mV
Xc1 = 524.4
Xc2 = 529.9

Yc1 = 521.6
Yc2 = 523.3

Zc1 = 1037.9
Zc2 = 1027.7
Zc3 = 1029.9
Zc4 = 1039.0

#--- Default settings
GRAVITY = 9.8 # m/s/s

#--- List of files to reverse
REVERSE_LIST = ['HIT_009.KDA','GR1_094.KDA']

#--- A function to do a continuous sigma edit
def sigmaEdit(x,sigmaThresh = None):
    if sigmaThresh <= sqrt(3):
        print 'Error in sigmaEdit'
        print 'Sigma Threshold must be greater than sqrt(3)'
        print 'Sigma Threshold given: %s' % sigmaThresh
    
    nPointsRemoved = 1
    
    while nPointsRemoved > 0:
        # Get stats
        # UpperBound = mean + sigmaThresh*std
        # LowerBound = mean - sigmaThresh*std
        # Apply bounds and determine points removed
        nPointsRemoved = 0
    
    return x

#--- A function to look at the file header
def readHeader(filename): # Name of the file to read
    
    #--- Set a number for the maximum size of the header in the file
    max_header_lines = 10
    
    #--- Set up regular expressions to match contents of header
    #    Regular expressions are nice ways of specifying exactly what you
    #    are looking for in a line of text and pulling out the data you
    #    want, in this case the values between the < and > characters
    rexp_header = re.compile(r'(?P<header>[A-Z]{2}[0-9]{1})')
    rexp_datetime = re.compile(r'(?P<month>\d{1,2})/(?P<day>\d{1,2})/(?P<year>\d{4})\s+[-]{2}\s+(?P<hour>\d{1,2}):(?P<minute>\d{1,2})\s+(?P<meridiem>[A-Z]{2})')
    rexp_scanrate = re.compile(r'Scan Rate = (?P<scan_rate>[0-9]{1,4})')
    rexp_totalscans = re.compile(r'Total Scans = (?P<total_scans>[0-9]{1,4})')
    rexp_triggerstate = re.compile(r'Digital Trigger (?P<trigger_state>\w{2,3})')
    
    #--- Open the file in a read-only state
    file = open(filename, 'r')
    
    #--- Set the following default variables
    header = 'GRF'
    timestamp = '2000-01-01T12:00'
    scan_rate = 1200
    total_scans = 4800
    trigger_state = 'Off'
    
    #--- Read the header information from the file by looking at each line
    count = 0
    for line in file:
        #--- Get the test header
        rmatch_header = rexp_header.match(line)
        if rmatch_header:
            header = rmatch_header.group('header')
        
        #--- Get the date and time of the test to create a timestamp
        rmatch_datetime = rexp_datetime.match(line)
        if rmatch_datetime:
            year =   int(rmatch_datetime.group('year'))
            month =  int(rmatch_datetime.group('month'))
            day =    int(rmatch_datetime.group('day'))
            hour =   int(rmatch_datetime.group('hour'))
            minute = int(rmatch_datetime.group('minute'))
            meridiem = rmatch_datetime.group('meridiem')
            if meridiem == 'PM' and hour != 12:
                hour += 12
            elif meridiem == 'AM' and hour == 12:
                hour -= 12
            timestamp = '%04d_%02d_%02dT%02d_%02d' % \
                        (year,month,day,hour,minute)
        
        #--- Record the scan rate in Hz
        rmatch_scanrate = rexp_scanrate.match(line)
        if rmatch_scanrate:
            scan_rate = int(rmatch_scanrate.group('scan_rate'))
        
        #--- Record the total number of expected scans
        rmatch_totalscans = rexp_totalscans.match(line)
        if rmatch_totalscans:
            total_scans = int(rmatch_totalscans.group('total_scans'))
        
        #--- Determine the trigger state (perhaps not used)
        rmatch_triggerstate = rexp_triggerstate.match(line)
        if rmatch_triggerstate:
            trigger_state = rmatch_triggerstate.group('trigger_state')
        
        #--- Increment the counter until the max number of header lines
        #    has been read.  Then exit the loop.
        count += 1
        if count > max_header_lines:
            break
    
    #--- Return the important information
    return header, timestamp, scan_rate, total_scans, trigger_state

def parseFile(filename,             # Name of the file to parse
              range = (-1, -1)):    # Frame range
    
    #--- We will save the data to a dictionary
    data_dict = {}
    
    #--- Read the header using another function
    #    This returns useful information to verify the file and name it
    header, timestamp, scan_rate, total_scans, trigger_state = readHeader(filename)
    data_dict['filename'] = filename
    data_dict['header'] = header
    data_dict['timestamp'] = timestamp
    data_dict['scan_rate'] = scan_rate
    data_dict['total_scans'] = total_scans
    data_dict['trigger_state'] = trigger_state
    
    #--- Create a test name identifier for naming and saving plots
    #    The convention used here is "FILENAME_TIMESTAMP"
    title = os.path.basename(filename).split('.')[0]
    data_dict['title'] = title
    identifier = '%s_%s' % (title,timestamp)
    data_dict['identifier'] = identifier
    
    #--- Load the data from the file name that was given
    #    Skip the first six rows of header data
    #    Use a comma as the delimiter
    data = numpy.loadtxt(filename, skiprows=6, delimiter=',')
    
    #--- Do some internal checking on the data
    #    compare the data length against the file's header information
    if len(data) != total_scans:
        print "\nA problem exists with your file:"
        print "\tNumber of records expected: %d" % (total_scans)
        print "\tNumber of records found: %d" % (len(data))
        print "\tCheck that file is not missing data or corrupted"
        print "\tThe easiest fix is to adjust the header if everything looks ok"
        sys.exit()
    
    #--- Set Frame Defaults to encompass all the data
    start_frame = 0
    end_frame = len(data)-1

    #--- Determine the range of the data to use for the plot
    #    Do nothing if either frame value is set to the default
    if range[0] != -1 and range[1] != -1:
        #--- Ensure start_frame not less than zero
        if int(range[0]) < 0:
            print "\nYou cannot specify an initial frame less than zero"
            sys.exit()
        else:
            start_frame = int(range[0])
        
        #--- Ensure end_frame not outside data boundary
        if int(range[1]) >= len(data):
            print "\nYou cannot specify an ending frame larger than the data set"
            print "\tNumber of records found: %d" % (len(data))
            print "\tGiven End Frame: %s" % (range[1]) 
            sys.exit()
        else:
            end_frame = int(range[1])
        
        #--- Make sure that start_frame is less than end_frame
        if start_frame >= end_frame:
            print "\nYou must specify the start frame before the end frame"
            print "\tGiven Start Frame: %s" % (start_frame)
            print "\tGiven End Frame: %s" % (end_frame)
            sys.exit()
        
        #--- Slice the data to the desired size
        data = data[start_frame:end_frame]

    #--- Now put the GRF data for each column into a variable
    #    p1_ and p2_ are prefixes for plates 1 and 2 respectively
    #    X1, X2 - Raw force data for left-right direction of plate
    #    Y1, Y2 - Raw force data for front-back direction of plate
    #    Z1, Z2, Z3, Z4 - Raw force data for up-down direction of plate
    p1_X1,p1_X2,p1_Y1,p1_Y2 =  data[:,0], data[:,1], data[:,2], data[:,3]
    p1_Z1,p1_Z2,p1_Z3,p1_Z4 =  data[:,4], data[:,5], data[:,6], data[:,7]

    p2_X1,p2_X2,p2_Y1,p2_Y2 =  data[:,8], data[:,9], data[:,10], data[:,11]
    p2_Z1,p2_Z2,p2_Z3,p2_Z4 =  data[:,12], data[:,13], data[:,14], data[:,15]
    
    #--- Apply the conversion factors and add the data
    #    Units: N = mV * N/mV
    p1_X = p1_X1*Xc1 + p1_X2*Xc2
    p1_Y = p1_Y1*Yc1 + p1_Y2*Yc2
    p1_Z = p1_Z1*Zc1 + p1_Z2*Zc2 + p1_Z3*Zc3 + p1_Z4*Zc4
    
    p2_X = p2_X1*Xc1 + p2_X2*Xc2
    p2_Y = p2_Y1*Yc1 + p2_Y2*Yc2
    p2_Z = p2_Z1*Zc1 + p2_Z2*Zc2 + p2_Z3*Zc3 + p2_Z4*Zc4
    
    #--- Reverse files that need to be reversed
    if os.path.basename(filename) in REVERSE_LIST:
        p1_X = -1 * p1_X
        p2_X = -1 * p2_X
    
    #--- Put this data in the dictionary
    #    Units: N
    data_dict['p1_X'] = p1_X
    data_dict['p1_Y'] = p1_Y
    data_dict['p1_Z'] = p1_Z
    
    data_dict['p2_X'] = p2_X
    data_dict['p2_Y'] = p2_Y
    data_dict['p2_Z'] = p2_Z
    
    
    
    #--- Get the magnitude of the force data in each direction
    #    Units: N
    data_dict['p1_XY_mag']  = (p1_X**2 + p1_Y**2)**0.5
    data_dict['p1_XZ_mag']  = (p1_X**2 + p1_Z**2)**0.5
    data_dict['p1_YZ_mag']  = (p1_Y**2 + p1_Z**2)**0.5
    data_dict['p1_XYZ_mag'] = (p1_X**2 + p1_Y**2 + p1_Z**2)**0.5
    
    data_dict['p2_XY_mag']  = (p2_X**2 + p2_Y**2)**0.5
    data_dict['p2_XZ_mag']  = (p2_X**2 + p2_Z**2)**0.5
    data_dict['p2_YZ_mag']  = (p2_Y**2 + p2_Z**2)**0.5
    data_dict['p2_XYZ_mag'] = (p2_X**2 + p2_Y**2 + p2_Z**2)**0.5

    #--- Set the frequency and find the time delta between frames
    freq = scan_rate   # Hz
    delta_t = 1.0/freq # s
    data_dict['delta_t'] = delta_t
    
    #--- From the number of frames generate a time series for plotting
    #    If the user is inspecting the data then do not use the delta_t,
    #    this will cause the plots to use the frame number on the x-axis
    #
    #    Don't forget to add the original start_frame to the beginning
    #    to ensure the graph looks the same as before
    frame = pylab.arange(0,len(data)) # This is an array-range data structure
    frame += start_frame
    data_dict['frame'] = frame
    
    #--- Put the time information in the dictionary
    time = frame*delta_t
    data_dict['time'] = time
    
    #--- Calculate the total time for the given data range
    #    Units: s
    data_dict['total_time'] = len(frame)*delta_t
    
    #--- Use the Trapezoid Rule to calculate the net impulse for each plate
    #    NOTE: Do not set 'x' unless using variable sampling rate
    #    Units: N*s
    data_dict['p1_X_imp_net'] = numpy.trapz(p1_X, x = None, dx = delta_t, axis = -1)
    data_dict['p1_Y_imp_net'] = numpy.trapz(p1_Y, x = None, dx = delta_t, axis = -1)
    data_dict['p1_Z_imp_net'] = numpy.trapz(p1_Z, x = None, dx = delta_t, axis = -1)
    
    data_dict['p2_X_imp_net'] = numpy.trapz(p2_X, x = None, dx = delta_t, axis = -1)
    data_dict['p2_Y_imp_net'] = numpy.trapz(p2_Y, x = None, dx = delta_t, axis = -1)
    data_dict['p2_Z_imp_net'] = numpy.trapz(p2_Z, x = None, dx = delta_t, axis = -1)
    
    
    #--- First Derivative
    #
    #data_dict['p1_X_dx'] = numpy.diff(p1_X)
    #data_dict['p1_X_dx'] = numpy.append(data_dict['p1_X_dx'],data_dict['p1_X'][-1])
    #
    #print data_dict['time'].shape, data_dict['p1_X_dx'].shape
    #pylab.plot(data_dict['time'],data_dict['p1_X'],'-r')
    #pylab.grid(True)
    #pylab.twinx()
    #pylab.plot(data_dict['time'],data_dict['p1_X_dx'],'-k')
    #pylab.grid(True)
    #pylab.show()
    #sys.exit()
    
    return data_dict

#--- Define headers that will be printed for weight or impulse data
def weight_header():
    """This is the header for weight information"""
    print "#%s%s%s%s%s" % (str("File"     ).rjust(14,' '),
                          str("P1 (N)"    ).rjust(15,' '),
                          str("P2 (N)"    ).rjust(15,' '),
                          str("Total (N)" ).rjust(15,' '),
                          str("Total (kg)").rjust(15,' '))

def impulse_header():
    """This is the header for impulse information"""
    print "#%s%s%s%s%s%s%s%s%s" % (str('File'    ).rjust(14,' '),
                                  str('Timestamp').rjust(20,' '),
                                  str('Time (s)' ).rjust(10,' '),
                                  str('P1 X (Ns)').rjust(10,' '),
                                  str('P1 Y (Ns)').rjust(10,' '),
                                  str('P1 Z (Ns)').rjust(10,' '),
                                  str('P2 X (Ns)').rjust(10,' '),
                                  str('P2 Y (Ns)').rjust(10,' '),
                                  str('P2 Z (Ns)').rjust(10,' ')) 


#--- Determine the mass of the subject on each plate
#    After determining the weight exit the program
def weight(data_dict):
    
    filename = data_dict['filename']
    p1_Z = data_dict['p1_Z']
    p2_Z = data_dict['p2_Z']
    
    p1_weight = pylab.mean(p1_Z) # N
    p2_weight = pylab.mean(p2_Z) # N
    total_weight = p1_weight + p2_weight # N
    total_mass = total_weight / GRAVITY  # kg
    
    print "%s%s%s%s%s" % (str('%s' % os.path.basename(filename)).rjust(15,' '),
                          str('%.3f' % p1_weight   ).rjust(15,' '),
                          str('%.3f' % p2_weight   ).rjust(15,' '),
                          str('%.3f' % total_weight).rjust(15,' '),
                          str('%.3f' % total_mass  ).rjust(15,' '))
    return 0

def impulse(data_dict):
    """A method to calculate the total impulse"""
    
    filename = data_dict['filename']
    timestamp = data_dict['timestamp']
    total_time = data_dict['total_time']
    p1_X_imp_net = data_dict['p1_X_imp_net']
    p1_Y_imp_net = data_dict['p1_Y_imp_net']
    p1_Z_imp_net = data_dict['p1_Z_imp_net']
    p2_X_imp_net = data_dict['p2_X_imp_net']
    p2_Y_imp_net = data_dict['p2_Y_imp_net']
    p2_Z_imp_net = data_dict['p2_Z_imp_net']
    
    #--- Print the results of the total impulse
    print "%s%s%s%s%s%s%s%s%s" % (str('%s' % os.path.basename(filename)).rjust(15,' '),
                                  str('%s'   % timestamp   ).rjust(20,' '),
                                  str('%.3f' % total_time  ).rjust(10,' '),
                                  str('%.3f' % p1_X_imp_net).rjust(10,' '),
                                  str('%.3f' % p1_Y_imp_net).rjust(10,' '),
                                  str('%.3f' % p1_Z_imp_net).rjust(10,' '),
                                  str('%.3f' % p2_X_imp_net).rjust(10,' '),
                                  str('%.3f' % p2_Y_imp_net).rjust(10,' '),
                                  str('%.3f' % p2_Z_imp_net).rjust(10,' '))
    
def plot_plates(data_dict,
                inspect = False,      # Inspect the graphs by frame number
                t_range = None,       # The time range to set for plots
                plate_1 = False,      # Plot plate 1 forces
                plate_2 = False,      # Plot plate 2 forces
                x_plot = False,       # Plot forces in x-axis
                y_plot = False,       # Plot forces in y-axis
                z_plot = False,       # Plot forces in z-axis
                save_plot = False):   # Save plots at *.png files to working dir
    """
    This method is useful for plotting all the XYZ-axis data for a single
    plate on one graph.  This method will also print all the plate data,
    Plate 1 and Plate 2, for a single axis (X, Y or Z) on a single graph.
    """
    
    #--- If inspecting the frames set the time field to the frame number
    time = data_dict['time']
    if inspect:
        time = data_dict['frame']
    
    title = data_dict['title']
    timestamp = data_dict['timestamp']
    identifier = data_dict['identifier']
    p1_X = data_dict['p1_X']
    p1_Y = data_dict['p1_Y']
    p1_Z = data_dict['p1_Z']
    p2_X = data_dict['p2_X']
    p2_Y = data_dict['p2_Y']
    p2_Z = data_dict['p2_Z']
    
    #--- Make plots of the force on each plate
    if plate_1:
        pylab.figure()
        pylab.plot(time, p1_X, '.-b', label='X-axis')
        pylab.plot(time, p1_Y, '.-g', label='Y-axis')
        pylab.plot(time, p1_Z, '.-r', label='Z-axis')
        if t_range:
            pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
        pylab.legend(loc='best')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Force (N)')
        pylab.title('%s, %s\nPlate 1 Force Plot' % (title, timestamp))
        pylab.grid(True)
        
        #--- Save the plot if requested
        if save_plot:
            figure_name = '%s_plot_p1' % (identifier)
            pylab.savefig(figure_name)
            
    if plate_2:
        pylab.figure()
        pylab.plot(time, p2_X, '.-b', label='X-axis')
        pylab.plot(time, p2_Y, '.-g', label='Y-axis')
        pylab.plot(time, p2_Z, '.-r', label='Z-axis')
        if t_range:
            pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
        pylab.legend(loc='best')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Force (N)')
        pylab.title('%s, %s\nPlate 2 Force Plot' % (title, timestamp))
        pylab.grid(True)
        
        #--- Save the plot if requested
        if save_plot:
            figure_name = '%s_plot_p2' % (identifier)
            pylab.savefig(figure_name)
    
    if x_plot:
        pylab.figure()
        pylab.plot(time, p1_X, '.-b', label='Plate 1')
        pylab.plot(time, p2_X, '.-r', label='Plate 2')
        if t_range:
            pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
        pylab.legend(loc='best')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Force (N)')
        pylab.title('%s, %s\nX-axis Force Plot' % (title, timestamp))
        pylab.grid(True)
        
        #--- Save the plot if requested
        if save_plot:
            figure_name = '%s_plot_x' % (identifier)
            pylab.savefig(figure_name)
        
    if y_plot:
        pylab.figure()
        pylab.plot(time, p1_Y, '.-b', label='Plate 1')
        pylab.plot(time, p2_Y, '.-r', label='Plate 2')
        if t_range:
            pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
        pylab.legend(loc='best')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Force (N)')
        pylab.title('%s, %s\nY-axis Force Plot' % (title, timestamp))
        pylab.grid(True)
        
        #--- Save the plot if requested
        if save_plot:
            figure_name = '%s_plot_y' % (identifier)
            pylab.savefig(figure_name)
        
    if z_plot:
        pylab.figure()
        pylab.plot(time, p1_Z, '.-b', label='Plate 1')
        pylab.plot(time, p2_Z, '.-r', label='Plate 2')
        if t_range:
            pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
        pylab.legend(loc='best')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Force (N)')
        pylab.title('%s, %s\nZ-axis Force Plot' % (title, timestamp))
        pylab.grid(True)
        
        #--- Save the plot if requested
        if save_plot:
            figure_name = '%s_plot_z' % (identifier)
            pylab.savefig(figure_name)

#--- Plot the collection of plots
def plot_collection(file_dict = {},       # The dictionary of file names and data
                    inspect = False,      # Inspect the graphs by frame number
                    t_range = None,       # The time range to set for plots
                    plate_1 = False,      # Plot plate 1 forces
                    plate_2 = False,      # Plot plate 2 forces
                    mag_plot = False,     # Plot the magnitude of the forces
                    x_plot = False,       # Plot forces in x-axis
                    y_plot = False,       # Plot forces in y-axis
                    z_plot = False,       # Plot forces in z-axis)
                    save_plot = False):   # Save plots at *.png files to working dir
    """
    This method is used for plotting multiple collections of data separated
    by both axis and plate.  This means up to six graphs will be printed
    depending on options given with data from every file on each graph.
    """

    #--- Get a list of all the file names
    file_names = file_dict.keys()
    file_names.sort()
    
    #--- These are the available colors we can use for the plot lines
    #    The program will cycle through them to make reading the graphs easier
    color_list = ['b','g','r','c','m','y','k']
    
    #--- Create a list of plates to graph
    plate_list = []
    if plate_1:
        plate_list.append('1')
    if plate_2:
        plate_list.append('2')
        
    #--- Loop through each plate in the list
    for plate in plate_list:
        
        if x_plot:
            pylab.figure()
            count = 0
            for file in file_names:
                
                color = color_list[count]
                
                #--- If inspecting the frames set the time field to the frame number
                time = file_dict[file]['time'] - file_dict[file]['delta_t'] * file_dict[file]['align']
                if inspect:
                    time = file_dict[file]['frame'] - file_dict[file]['align']
                
                label = file_dict[file]['title']
                pylab.plot(time, file_dict[file]['p%s_X' % plate], '-%s' % (color) , label=label)
                count += 1
                if count >= len(color_list):
                    count = 0
            
            if t_range:
                pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
            pylab.legend(loc='best')
            pylab.xlabel('Time (s)')
            pylab.ylabel('Force (N)')
            pylab.title('X-axis Force Plot Plate %s' % (plate))
            pylab.grid(True)
            
            #--- Save the plot if requested
            if save_plot:
                figure_name = 'plate_%s_plot_x' % (plate)
                pylab.savefig(figure_name)
            
        if y_plot:
            pylab.figure()
            count = 0
            for file in file_names:
                
                color = color_list[count]
                
                #--- If inspecting the frames set the time field to the frame number
                time = file_dict[file]['time'] - file_dict[file]['delta_t'] * file_dict[file]['align']
                if inspect:
                    time = file_dict[file]['frame'] - file_dict[file]['align']
                
                label = file_dict[file]['title']
                pylab.plot(time, file_dict[file]['p%s_Y' % plate], '-%s' % (color) , label=label)
                count += 1
                if count >= len(color_list):
                    count = 0
            
            if t_range:
                pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
            pylab.legend(loc='best')
            pylab.xlabel('Time (s)')
            pylab.ylabel('Force (N)')
            pylab.title('Y-axis Force Plot Plate %s' % (plate))
            pylab.grid(True)
            
            #--- Save the plot if requested
            if save_plot:
                figure_name = 'plate_%s_plot_y' % (plate)
                pylab.savefig(figure_name)
            
        if z_plot:
            pylab.figure()
            count = 0
            for file in file_names:
                
                color = color_list[count]
                
                #--- If inspecting the frames set the time field to the frame number
                time = file_dict[file]['time'] - file_dict[file]['delta_t'] * file_dict[file]['align']
                if inspect:
                    time = file_dict[file]['frame'] - file_dict[file]['align']
                    
                label = file_dict[file]['title']
                pylab.plot(time, file_dict[file]['p%s_Z' % plate], '-%s' % (color) , label=label)
                count += 1
                if count >= len(color_list):
                    count = 0
            
            if t_range:
                pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
            pylab.legend(loc='best')
            pylab.xlabel('Time (s)')
            pylab.ylabel('Force (N)')
            pylab.title('Z-axis Force Plot Plate %s' % (plate))
            pylab.grid(True)
            
            #--- Save the plot if requested
            if save_plot:
                figure_name = 'plate_%s_plot_z' % (plate)
                pylab.savefig(figure_name)
        
        #--- Check that the required set of axis are given
        do_plot_mag = (options.x_plot and options.y_plot)
        do_plot_mag = do_plot_mag or (options.x_plot and options.z_plot)
        do_plot_mag = do_plot_mag or (options.y_plot and options.z_plot)
        if mag_plot and do_plot_mag:
            pylab.figure()
            count = 0
            for file in file_names:
                
                #--- The options given determine the dataset to choose from
                if all([options.x_plot, options.y_plot, options.z_plot]):
                    axis = 'XYZ'
                elif all([options.x_plot, options.y_plot]):
                    axis = 'XY'
                elif all([options.x_plot, options.z_plot]):
                    axis = 'XZ'
                elif all([options.y_plot, options.z_plot]):
                    axis = 'YZ'
                
                color = color_list[count]
                
                #--- If inspecting the frames set the time field to the frame number
                time = file_dict[file]['time'] - file_dict[file]['delta_t'] * file_dict[file]['align']
                if inspect:
                    time = file_dict[file]['frame'] - file_dict[file]['align']
                
                label = file_dict[file]['title']
                pylab.plot(time, file_dict[file]['p%s_%s_mag' % (plate,axis)], '-%s' % (color) , label=label)
                count += 1
                if count >= len(color_list):
                    count = 0
            
            if t_range:
                pylab.axis([t_range[0],t_range[1],t_range[2],t_range[3]])
            pylab.legend(loc='best')
            pylab.xlabel('Time (s)')
            pylab.ylabel('Force (N)')
            pylab.title('%s-axis Magnitude Force Plot Plate %s' % (axis,plate))
            pylab.grid(True)
            
            #--- Save the plot if requested
            if save_plot:
                figure_name = 'plate_%s_%s_mag_plot' % (plate,axis)
                pylab.savefig(figure_name)

#--- Declare the program that will run
if __name__ == '__main__':
    
    #--- set up the command line arguments
    description = 'Program to integrate force data from a *.kda file'
    usage = "%prog [-a] [-c] [-d dirName] [-f filename] [-i] [-n ARG1 ARG2] \
[-p 0|1|2] [-r ARG1 ARG2] [-s] [-t ARG1 ARG2 ARG3 ARG4] [-w] [-x] [-y] [-z]"
    
    p = optparse.OptionParser(usage,description=description)
    
    p.add_option('-a', action='store', type='string', dest='align',
                 help='Frames numbers to use when aligning the data')
    p.add_option('-c', action="store_true", dest="collect", default=False,
                 help='Plot forces in a collection')
    p.add_option('-d', '--dir', action='store', type='string', dest='dirName',
                 default = '.', help='Directory containing KDA files with csv data to batch process')
    p.add_option('-f','--file', action='store', type='string', dest='filename',
                 help='KDA File containing csv data')
    p.add_option('-i', action="store_true", dest="inspect", default=False,
                 help='Inspect the data by viewing plots with frame count instead of time on x-axis')
    p.add_option('-m', action="store_true", dest="mag_plot", default=False,
                 help='Plot force magnitude for both plates')
    p.add_option('-p', action="store", type='int', dest="plate",
                 help='Plot force in x and z for specified plate number (0 [both], 1 [plate 1] or 2 [plate 2])')
    p.add_option('-r', action="store", dest="range", nargs=2, default=(-1,-1),
                 help='Parse a range of data in the file between given frame numbers (start frame -> end frame)')
    p.add_option('-s', action="store_true", dest="save_plot", default=False,
                 help='Save each requested plot to a file.  Plots save automatically when batch processing.')
    p.add_option('-t', action="store", dest="t_range", nargs=4, default=(None,None,None,None),type='float',
                 help='Set the time range for the file using the start and end time and also include the force range')
    p.add_option('-w', action="store_true", dest="weight", default = False,
                 help='Determine the weight of a player from the given range of data')
    p.add_option('-x', action="store_true", dest="x_plot", default=False,
                 help='Plot force in x direction for both plates')
    p.add_option('-y', action="store_true", dest="y_plot", default=False,
                 help='Plot force in y direction for both plates')
    p.add_option('-z', action="store_true", dest="z_plot", default=False,
                 help='Plot force in z direction for both plates')
    
    options,arguments = p.parse_args()
    
    #--- If any extra arguments exist print the help and exit
    if len(arguments):
        p.print_help()
        sys.exit()
    
    #--- Determine which plates will be used to produced plots
    #    By default both plates will produce plots
    plate_1 = False
    plate_2 = False
    if options.plate >= 0:
        plate_1 = True
        plate_2 = True
        if options.plate == 1:
            plate_2 = False
        elif options.plate == 2:
            plate_1 = False
        elif options.plate != 0:
            p.print_help()
            sys.exit() 
    
    #--- If a user specifies a filename then only parse that file
    file_list = []
    if options.filename and os.path.isfile(options.filename):
        file_list.append(options.filename)
    
    #--- If no filename is specified then look through the given directory
    #    The default is to look in the current directory if none is given
    elif os.path.isdir(options.dirName):
        
        #--- Add any files in the list that match the File Extension List
        for ext in FILE_EXT_LIST:
            files = glob.glob(os.path.join(options.dirName,ext))
            file_list.extend(files)

    #--- Ensure unique values for the list (important if os is not case sensitive)
    #    Sort the list to ensure alignment values correspond to files
    file_list = list(set(file_list))
    file_list.sort()
    
    #--- Set up a regular expression for the file name
    #    (ie three letters and number, an underscore, three numbers, file ext)
    #    Example: GR1_020.KDA 
    rexp_filename = re.compile(r'(?P<filename>[A-Z0-9]{3}_[0-9]{3})\.[KDA|kda]')
    
    #--- Check each item in the list to see if it matches the naming convention
    for file in file_list:
        rmatch_filename = rexp_filename.match(os.path.basename(file))
        
        #--- If there is no match then remove file
        #if rmatch_filename:
        #    print rmatch_filename.group('filename')
        if not rmatch_filename:
            file_list.remove(file)
    
    #--- Grab the frame numbers
    align_list = []
    if options.align:
        align_list = options.align.split(',')
    
    #--- Make sure that there is one value for every file
    while len(align_list) < len(file_list):
        align_list.append('0')

    #--- Do not save any plots if looking at weight
    if options.weight:
        options.save_plot = False
    
    #--- Create a dictionary that holds the file name and time data
    file_dict = {}
    
    #--- Keep an index of the files and cycle through the list
    count = 0
    for file in file_list:
        
        #--- Parse the data into a dictionary
        file_dict[file] = parseFile(file, range = options.range)
        file_dict[file]['align'] = float(align_list[count])
        
        #--- Print weight or impulse data for the file
        if options.weight:
            #--- Print header on first time through
            if count == 0:
                weight_header()
            weight(file_dict[file])
        else:
            #--- Print header on first time through
            if count == 0:
                impulse_header()
            impulse(file_dict[file])
        
        #--- Increase the counter to get the next alignment value
        count += 1
    
    #--- If only looking at one file use this
    if len(file_list) == 1 and not options.collect:
        
        plot_plates(file_dict[file_list[0]],
                    inspect = options.inspect,
                    t_range = options.t_range,
                    plate_1 = plate_1,
                    plate_2 = plate_2,
                    x_plot = options.x_plot,
                    y_plot = options.y_plot,
                    z_plot = options.z_plot,
                    save_plot = options.save_plot)
    
    #--- Plot the data
    elif options.collect:
        
        plot_collection(file_dict = file_dict,
                        inspect = options.inspect,
                        t_range = options.t_range,
                        plate_1 = plate_1,
                        plate_2 = plate_2,
                        mag_plot = options.mag_plot,
                        x_plot = options.x_plot,
                        y_plot = options.y_plot,
                        z_plot = options.z_plot,
                        save_plot = options.save_plot)
    
    #--- Finally, show the plots unless calculating weight
    #    Plot 1 or 2 (or both) must be chose
    #    You must not be doing weight calculations
    #    Finally, if any axis is chosen then plot
    do_plot = (plate_1 or plate_2) and not options.weight
    do_plot = do_plot and (options.x_plot or options.y_plot or options.z_plot)
    if do_plot and (options.collect or len(file_list) == 1):
        pylab.show()
    