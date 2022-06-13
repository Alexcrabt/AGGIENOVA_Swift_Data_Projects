'''
@author: Alexander Crabtree
            acrabtree15@gmail.com
'''
import plotly.graph_objects as go

import plotly.express as px
from plotly.validators.scatter.marker import SymbolValidator

from plotly.subplots import make_subplots

import pandas as pd
import numpy as np

from itertools import cycle
import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

#Must have plotly version 5.1.0 higher installed.
#To do so use pip install plotly==5.1.0

#To install all dash components use:
# pip install dash
# pip install dash_core_components
# pip install dash_html_components

def sn_Data(data):
    '''
    Grabs all data on the supernovae and returns it as a dataframe.
    '''
    sndata= pd.DataFrame(columns= ['filter', 'mjd', 'mag', 'magerr'])

    for sublist in data:
        if not sublist[0] == '#':
            sndata= sndata.append({'filter': sublist[0], 'mjd': sublist[1], 'mag': sublist[2], 'magerr': sublist[3]}, ignore_index=True)

    sndata=sndata.sort_values(by="mjd", ascending= True, na_position='last')
    sndata.reset_index(drop=True, inplace=True)
    return sndata

def color_Data(sn_filter,filter1, filter2, mjd, mag, magerr):
    '''
    creates temp lists and filters all the supernovae data
    based on the 2 given filters and returns all data just for thhose filters.
    '''

    mjd1_temp=[]
    mjd2_temp=[]
    mag1_temp= []
    mag2_temp= []
    magerr1_temp=[]
    magerr2_temp=[]
    for i in range(len(sn_filter)):
        if not sn_filter[i] == filter1 and not sn_filter[i] == filter2:
            continue
        elif sn_filter[i] == filter1:
            mjd1_temp.append(mjd[i])
            mag1_temp.append(mag[i])
            magerr1_temp.append(magerr[i])
        else:
            mjd2_temp.append(mjd[i])
            mag2_temp.append(mag[i])
            magerr2_temp.append(magerr[i])

    if not len(mjd1_temp) or not len(mag1_temp):
        raise
    elif not len(mjd2_temp) or not len(mag2_temp):
        raise
    else:
        return mjd1_temp, mag1_temp, magerr1_temp, mjd2_temp, mag2_temp, magerr2_temp

def filt_Data(sn_filter,filter1, mjd, mag, magerr):
    '''
    creates temp lists and filters all the supernovae data
    based on the given filter and returns all data just for that filter.
    '''
    mjd_temp=[]
    mag_temp= []
    magerr_temp=[]


    for i in range(len(sn_filter)):
        if not sn_filter[i] == filter1:
            continue
        else:
            mjd_temp.append(mjd[i])
            mag_temp.append(mag[i])
            magerr_temp.append(magerr[i])

    if not len(mjd_temp) or not len(mag_temp):
        raise
    else:
        return mjd_temp, mag_temp, magerr_temp

def null_Color_Filter(mjd1_temp, mag1_temp, magerr1_temp, mjd2_temp, mag2_temp, magerr2_temp, j):
    '''
    Takes mjd, mag, and magerr and filters out NULL values.
    '''
    mjd1=[]
    mjd2=[]
    mag1= []
    mag2= []
    magerr1=[]
    magerr2=[]

    for i in range(len(mjd1_temp)):
        if mjd1_temp[i]== "NULL" or mag1_temp[i]== "NULL" :

            continue
        else:
            mjd1.append(float(mjd1_temp[i]))
            mag1.append(float(mag1_temp[i]))
            magerr1.append(float(magerr1_temp[i]))

    for i in range(len(mjd2_temp)):
        if mjd2_temp[i]== "NULL" or mag2_temp[i]== "NULL" :

            continue
        else:
            mjd2.append(float(mjd2_temp[i]))
            mag2.append(float(mag2_temp[i]))
            magerr2.append(float(magerr2_temp[i]))

    '''
    If mjd or mag lists are empty after filtering or all NULL values return not enough data.
    '''
    if not len(mjd1) or not len(mag1):
        raise #print("Not enough data on " + snname +'_'+ colors[j][0]+ ' for ' + colors[j][0]+ '-'+ colors[j][1])
    elif not len(mjd2) or not len(mag2):
        raise #print("Not enough data on " + snname +'_'+ colors[j][1]+ ' for ' + colors[j][0]+ '-'+ colors[j][1])
    else:
        return mjd1, mag1, magerr1, mjd2, mag2, magerr2

def null_Filter(mjd_temp, mag_temp, magerr_temp, j):
    '''
    Takes mjd, mag, and magerr and filters out NULL values.
    '''
    mjd= []
    mag= []
    magerr= []

    for i in range(len(mjd_temp)):
        if mjd_temp[i]== "NULL" or mag_temp[i]== "NULL" :

            continue
        else:
            mjd.append(float(mjd_temp[i]))
            mag.append(float(mag_temp[i]))
            magerr.append(float(magerr_temp[i]))

    '''
    If mjd or mag lists are empty after filtering or all NULL values return not enough data.
    '''
    if not len(mjd) or not len(mag):
        raise #print("Not enough data on " + snname + '_'+ filt1[j])
    else:
        return mjd, mag, magerr

def mjd_Check(mjd1, mag1, magerr1, mjd2, mag2, magerr2, dt):
    '''
    Checks if the mjd of both filters are in a certain ranage of each other, if they are
    it saves all data (mjd, mag, magerr) for them and returns them.
    '''

    if not dt=="":
        mjd1_s= [mjd-float(dt) for mjd in mjd1]
        mjd1_a= [mjd+float(dt) for mjd in mjd1]


        #_c stands for checked
        mjd_c=[]
        mag1_c= []
        magerr1_c= []

        mag2_c= []
        magerr2_c=[]

        for i in range(len(mjd1)):
            mjd_temp=[]
            mag1_temp= []
            magerr1_temp= []
            mag2_temp= []
            magerr2_temp= []
            for j in range(len(mjd2)):
                if min(mjd1_s[i], mjd1_a[i]) < mjd2[j] < max(mjd1_s[i], mjd1_a[i]):
                    mjd_temp.append(mjd1[i])
                    mag1_temp.append(mag1[i])
                    magerr1_temp.append(magerr1[i])

                    mag2_temp.append(mag2[j])
                    magerr2_temp.append(magerr2[j])
                else:
                    continue
            '''
            This if-else section checks for mag2_temp and magerr2_temp that has multiple values associated with
            the corresponding mjd_temp and averages the mag2_temp and magerr2_temp.
            '''
            if not len(mjd_temp):
                continue
            else:

                mjd_c.append(mjd_temp[0])
                mag1_c.append(mag1_temp[0])
                magerr1_c.append(magerr1_temp[0])

                mag2_mean= sum(mag2_temp)/len(mag2_temp)
                magerr2_mean= sum(magerr2_temp)/len(magerr2_temp)

                mag2_c.append(mag2_mean)
                magerr2_c.append(magerr2_mean)

        '''
        This section checks for multiple mjd_c that share the same repeated mag2 values and averages those mjd_c values together.
        '''
        comb_data= pd.DataFrame({'mjd': mjd_c, 'mag1': mag1_c, 'magerr1': magerr1_c, 'mag2': mag2_c, 'magerr2': magerr2_c})
        comb_data= comb_data.groupby('mag2').mean().reset_index()
        comb_data= comb_data.sort_values(by= 'mjd', ascending= True).reset_index(drop=True)

        mjd_c= comb_data['mjd'].to_list()
        mag1_c= comb_data['mag1'].to_list()
        magerr1_c= comb_data['magerr1'].to_list()
        mag2_c= comb_data['mag2'].to_list()
        magerr2_c= comb_data['magerr2'].to_list()

        return mjd_c, mag1_c, magerr1_c, mag2_c, magerr2_c

    else:
        dt_def= 0.1500
        mjd1_s= [mjd-dt_def for mjd in mjd1]
        mjd1_a= [mjd+dt_def for mjd in mjd1]


        #_c stands for checked
        mjd_c=[]
        mag1_c= []
        magerr1_c= []

        mag2_c= []
        magerr2_c=[]

        for i in range(len(mjd1)):
            mjd_temp=[]
            mag1_temp= []
            magerr1_temp= []
            mag2_temp= []
            magerr2_temp= []
            for j in range(len(mjd2)):
                if min(mjd1_s[i], mjd1_a[i]) < mjd2[j] < max(mjd1_s[i], mjd1_a[i]):
                    mjd_temp.append(mjd1[i])
                    mag1_temp.append(mag1[i])
                    magerr1_temp.append(magerr1[i])

                    mag2_temp.append(mag2[j])
                    magerr2_temp.append(magerr2[j])
                else:
                    continue

            '''
            This if-else section checks for mag2_temp and magerr2_temp that has multiple values associated with
            the corresponding mjd_temp and averages the mag2_temp and magerr2_temp.
            '''
            if not len(mjd_temp):
                continue
            else:
                mjd_c.append(mjd_temp[0])
                mag1_c.append(mag1_temp[0])
                magerr1_c.append(magerr1_temp[0])

                mag2_mean= sum(mag2_temp)/len(mag2_temp)
                magerr2_mean= sum(magerr2_temp)/len(magerr2_temp)

                mag2_c.append(mag2_mean)
                magerr2_c.append(magerr2_mean)

        '''
        This section checks for multiple mjd_c that share the same repeated mag2 values and averages those mjd_c values together.
        '''
        comb_data= pd.DataFrame({'mjd': mjd_c, 'mag1': mag1_c, 'magerr1': magerr1_c, 'mag2': mag2_c, 'magerr2': magerr2_c})
        comb_data= comb_data.groupby('mag2').mean().reset_index()
        comb_data= comb_data.sort_values(by= 'mjd', ascending= True).reset_index(drop=True)

        mjd_c= comb_data['mjd'].to_list()
        mag1_c= comb_data['mag1'].to_list()
        magerr1_c= comb_data['magerr1'].to_list()
        mag2_c= comb_data['mag2'].to_list()
        magerr2_c= comb_data['magerr2'].to_list()

        return mjd_c, mag1_c, magerr1_c, mag2_c, magerr2_c

def color_Val_ConverterCCD(mjd, mag1, magerr1, mag2, magerr2):
    '''
    Converts mjd to days_since_detection, both mags to color, and both magerrs to a final error
    I commented out prints that will tell you what SNe that there isn't enough data on, uncomment to
    turn it back on.
    '''
    if not len(mjd) or not len(mag1) or not len(mag2):
        raise #print("Not enough data on  " + snname)
    else:
        color= [mag_a - mag_b for mag_a, mag_b in zip(mag1, mag2)]
        error= [np.sqrt((error1)**2 + (error2)**2) for error1, error2 in zip(magerr1, magerr2)]

        '''
        If days_since_detection or color are nan lists it returns not enough data.
        I commented out prints that will tell you what SNe that there isn't enough data on, uncomment to
        turn it back on.
        '''
        if np.isnan(mjd).all() == True or np.isnan(color).all() == True:
            raise #print("Not enough data on  " + snname)
        else:
            return mjd, color, error

def color_Val_ConverterCL(mjd, mag1, magerr1, mag2, magerr2):
    '''
    Converts mjd to days_since_detection, both mags to color, and both magerrs to a final error
    '''
    if not len(mjd) or not len(mag1) or not len(mag2):
        raise #print("Not enough data on  " + snname)
    else:
        mjd_0= mjd[0]
        days_since_detection= [mjd- mjd_0 for mjd in mjd]
        color= [mag_a - mag_b for mag_a, mag_b in zip(mag1, mag2)]
        error= [np.sqrt((error1)**2 + (error2)**2) for error1, error2 in zip(magerr1, magerr2)]

        '''
        If days_since_detection or color are nan lists it returns not enough data
        '''
        if np.isnan(days_since_detection).all() == True or np.isnan(color).all() == True:
            raise #print("Not enough data on  " + snname)
        else:
            return days_since_detection, color, error

def filt_Value_Converter(mjd, mag, magerr, dist_mod):
    '''
    Converts mjd to days_since_detection, mag to absolute_mag
    '''
    if not len(mjd) or not len(mag):
        raise #print("Not enough data on  " + snname)
    else:
        mjd_0= mjd[0]
        days_since_detection= [mjd- mjd_0 for mjd in mjd]
        absolute_mag= [mag - dist_mod for mag in mag]
        final_magerr= [abs(error) for error in magerr]


        '''
        If days_since_detection or absolute_mag are nan lists it returns not enough data
        '''
        if np.isnan(days_since_detection).all() == True or np.isnan(absolute_mag).all() == True:
            raise #print("Not enough data on  " + snname)
        else:
            return days_since_detection, absolute_mag, final_magerr

def filter_Graphs_Comb(filt_dat):
    '''
    Takes the 6 individual filter data and reoganizes them to match with
    the color data
    '''
    filt_dat2= filt_dat[1:]
    i=0
    filt= []
    while i <= len(filt_dat):
        if not len(filt_dat2):
            break
        else:
            for j in range(len(filt_dat2)):
                if filt_dat[i] =="NULL":
                    filt.append("NULL")
                elif filt_dat2[j]== "NULL":
                    filt.append("NULL")
                else:
                    filt.append(filt_dat[i] + filt_dat2[j])
            filt_dat2.pop(0)
            i= i+1

    return filt
def FL_Data(snname, dist_mod, type, SNtype):
    global filt1

    #Grabs data file of given supernovae
    # /home/Alexcrabt/mysite/
    try:
        data= [i.strip().split() for i in open("data/"+snname+ "_uvotB15.1.dat")]
    except:
        #print(str(snname) + ' is not in our data files')
        return "NULL"

    try:
        sndata= sn_Data(data)

        #Turns specified columns in dataframe to lists
        sn_filter= sndata['filter'].to_list()
        mjd= sndata['mjd'].to_list()
        mag= sndata['mag'].to_list()
        magerr= sndata['magerr'].to_list()
    except:
        #prints out which supernovae we lack data on
        #print(str(snname) + ' is not in our data files')
        return "NULL"

    filt1= ['UVW2', 'UVM2', 'UVW1', 'U', 'B', 'V']
    filt_dat=[]

    for i in range(len(filt1)):
        try:
            filt_dat.append(filt_Data(sn_filter, filt1[i], mjd, mag, magerr))
        except:
            filt_dat.append("NULL")
            #prints out SNe filter we lack data on
            #print("Not enough data on " + snname+"_"+filt1[i])

    #Removes NAN/Null values for the filter data
    for i in range(len(filt_dat)):
        try:
            filt_dat[i]= null_Filter(filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], i)

        except:
            #returns Null if null_Filter returns empty list
            filt_dat[i]= "NULL"

    #Coverts values for filter light curves
    for i in range(len(filt_dat)):
        #if list has Null in it skips it
        if filt_dat[i] == "NULL":
            continue
        else:
            try:
                filt_dat[i]= filt_Value_Converter(filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], dist_mod)
            except:
                filt_dat[i]= "NULL"

    plots_dat=[]

    #Takes filter and color data and makes data to graph in plotly
    for i in range(len(filt_dat)):
        #If filter list has a null value skip graphing that as well as corresponding color data
        if filt_dat[i] == "NULL" or (len(filt_dat[i])==0):
            plots_dat.append(np.nan)
        else:
            if type.find(SNtype) == 0:
                if SNtype =="":
                    plots_dat.append((snname, filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], filt1[i]))
                elif type == SNtype:
                    plots_dat.append((snname, filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], filt1[i]))
            else:
            #If no SNe type match append Null
                plots_dat.append(np.nan)

    return plots_dat

def CL_Data(snname, dist_mod, type, SNtype, dt):
    global colors, filt1

    #Grabs data file of given supernovae
    # /home/Alexcrabt/mysite/
    try:
        data= [i.strip().split() for i in open("data/"+snname+ "_uvotB15.1.dat")]
    except:
        #print(str(snname) + ' is not in our data files')
        return "NULL"

    try:
        sndata= sn_Data(data)

        #Turns specified columns in dataframe to lists
        sn_filter= sndata['filter'].to_list()
        mjd= sndata['mjd'].to_list()
        mag= sndata['mag'].to_list()
        magerr= sndata['magerr'].to_list()
    except:
        #prints out which supernovae we lack data on
        #print(str(snname) + ' is not in our data files')
        return "NULL"

    filt1= ['UVW2', 'UVM2', 'UVW1', 'U', 'B', 'V']
    filt2= ['UVM2', 'UVW1', 'U', 'B', 'V']

    i=0
    color_dat=[]
    colors=[]

    #filters data based on given filter names to use as colors to graph later
    while i <= len(filt1):
        if not len(filt2):
            break
        else:
            for j in range(len(filt2)):
                try:
                    color_dat.append((color_Data(sn_filter, filt1[i], filt2[j], mjd, mag, magerr)))
                    colors.append((filt1[i], filt2[j]))
                except:
                    color_dat.append("NULL")
                    colors.append((filt1[i], filt2[j]))
                    #prints out exact SNe color we lack enough data on
                    #print("Not enough data on " + snname+"_"+filt1[i] +"-"+filt2[j])
            filt2.pop(0)
            i= i+1

    filt_dat=[]
    for i in range(len(filt1)):
        try:
            filt_dat.append(filt_Data(sn_filter, filt1[i], mjd, mag, magerr))
        except:
            filt_dat.append("NULL")
            #prints out SNe filter we lack data on
            #print("Not enough data on " + snname+"_"+filt1[i])

    #Removes NAN/Null values for the color data
    for i in range(len(color_dat)):
        try:
            color_dat[i]= null_Color_Filter(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], color_dat[i][5], i)

        except:
            #returns Null if null_Color_Filter returns empty list
            color_dat[i]= "NULL"

    #Removes NAN/Null values for the filter data
    for i in range(len(filt_dat)):
        try:
            filt_dat[i]= null_Filter(filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], i)

        except:
            #returns Null if null_Filter returns empty list
            filt_dat[i]= "NULL"

    #Checks if both filter mjds are in a certain range of each other for color data
    for i in range(len(color_dat)):
        #if list has Null in it skips it
        if color_dat[i] == "NULL":
            continue
        else:
            color_dat[i]= mjd_Check(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], color_dat[i][5], dt)

    #Converts values for color light curve
    for i in range(len(color_dat)):
        #if list has Null in it skips it
        if color_dat[i] == "NULL":
            continue
        else:
            try:
                color_dat[i]= color_Val_ConverterCL(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4])
            except:
                color_dat[i]= "NULL"

    #Coverts values for filter light curves
    for i in range(len(filt_dat)):
        #if list has Null in it skips it
        if filt_dat[i] == "NULL":
            continue
        else:
            try:
                filt_dat[i]= filt_Value_Converter(filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], dist_mod)
            except:
                filt_dat[i]= "NULL"

    #Pairs the filter data with corresponding color data

    filt_dat= filter_Graphs_Comb(filt_dat)
    plots_dat=[]

    #Takes filter and color data and makes data to graph in plotly
    for i in range(len(color_dat)):
        #If color list has a null value skip graphing that as well as corresponding filter data
        if color_dat[i] == "NULL":
            plots_dat.append(np.nan)
        #If filter list has a null value skip graphing that as well as corresponding color data
        elif filt_dat[i] == "NULL":
            plots_dat.append(np.nan)
        else:
            if type.find(SNtype) == 0:
                if SNtype =="":
                    plots_dat.append((snname, filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], filt_dat[i][3], filt_dat[i][4], filt_dat[i][5], color_dat[i][0], color_dat[i][1], color_dat[i][2], colors[i]))
                elif type == SNtype:
                    plots_dat.append((snname, filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], filt_dat[i][3], filt_dat[i][4], filt_dat[i][5], color_dat[i][0], color_dat[i][1], color_dat[i][2], colors[i]))
            else:
            #If no SNe type match append Null
                plots_dat.append(np.nan)

    return plots_dat

def CCD_Data(snname, dist_mod, type, SNtype, dt):
    global colors, filt1

    #Grabs data file of given supernovae
    # /home/Alexcrabt/mysite/
    try:
        data= [i.strip().split() for i in open("data/"+snname+ "_uvotB15.1.dat")]
    except:
        #print(str(snname) + ' is not in our data files')
        return "NULL"


    #Grabs specified data from datafile and saves as dataframe
    try:
        sndata= sn_Data(data)

        #Turns specified columns in dataframe to lists
        sn_filter= sndata['filter'].to_list()
        mjd= sndata['mjd'].to_list()
        mag= sndata['mag'].to_list()
        magerr= sndata['magerr'].to_list()
    except:
        #prints out which supernovae we lack data on
        #print(str(snname) + ' is not in our data files')
        return "NULL"


    filt1= ['UVW2', 'UVM2', 'UVW1', 'U', 'B', 'V']
    filt2= ['UVM2', 'UVW1', 'U', 'B', 'V']

    i=0
    color_dat=[]
    colors=[]

    #filters data based on given filter names to use as colors to graph later
    while i <= len(filt1):
        if not len(filt2):
            break
        else:
            for j in range(len(filt2)):
                try:
                    color_dat.append((color_Data(sn_filter, filt1[i], filt2[j], mjd, mag, magerr)))
                    colors.append((filt1[i], filt2[j]))
                except:
                    color_dat.append("NULL")
                    colors.append((filt1[i], filt2[j]))
                    #prints out exact SNe color we lack enough data on
                    #print("Not enough data on " + snname+"_"+filt1[i] +"-"+filt2[j])
            filt2.pop(0)
            i= i+1

    #Removes NAN/Null values for the color data
    for i in range(len(color_dat)):
        try:
            color_dat[i]= null_Color_Filter(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], color_dat[i][5], i)

        except:
            #returns Null if null_Color_Filter returns empty list
            color_dat[i]= "NULL"

    #Checks if both filter mjds are in a certain range of each other for color data

    for i in range(len(color_dat)):
        #if list has Null in it skips it
        if color_dat[i] == "NULL":
            continue
        else:
            color_dat[i]= mjd_Check(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], color_dat[i][5], dt)

    #Converts values for color light curve
    for i in range(len(color_dat)):
        #if list has Null in it skips it
        if color_dat[i] == "NULL":
            continue
        else:
            try:
                color_dat[i]= color_Val_ConverterCCD(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4])
            except:
                color_dat[i]= "NULL"

    plots_dat= []

    for i in range(len(color_dat)):

        #If color list has a null value skip graphing that as well as corresponding filter data
        if color_dat[i] == "NULL":
            plots_dat.append(np.nan)
        #If filter list has a null value skip graphing that as well as corresponding color data
        else:
            if type.find(SNtype) == 0:
                if SNtype =="":
                    plots_dat.append((snname, color_dat[i][0], color_dat[i][1], color_dat[i][2]))
                elif type == SNtype:
                    plots_dat.append((snname, color_dat[i][0], color_dat[i][1], color_dat[i][2]))
            else:
            #If no SNe type match append Null
                plots_dat.append(np.nan)

    return plots_dat

'''
This section runs the dash app by first setting up a basic plot layout before runing the various
inputs through the last 2 functions. It will repeat this for one or both functions depending on
the inputs given in the dash app.

I didn't encase this section in a main() because I was having issues getting it to run properly,
so think of everything past this section as main().

In making this dash app I used an asset folder to style it a bit. This asset folder is in the
repo, if you want to remove it delete external_stylesheets variable and the same varaible in
the app variable.
'''
external_stylesheets = ['assets/ColorColor.css']
app = dash.Dash(__name__, external_stylesheets= external_stylesheets, assets_folder='assets')

server = app.server

fig1=go.Figure()
fig1.update_layout(height=625,template= "plotly_dark")

fig2=go.Figure()
fig2.update_layout(height=625,template= "plotly_dark")


fig3=go.Figure()
fig3.update_layout(height=625,template= "plotly_dark")

app.layout= html.Div([

    html.Div([html.H1(children=[
            html.P(['These graphs were created using data from the Swift Optical Ultraviolet Supernova Archive ', html.A('(SOUSA)', href='https://pbrown801.github.io/SOUSA/'), '.']),
            html.P('For all graphs type in a preferred SNe Type and date range value (dt), then click both submit buttons.'),
            html.P('If you leave the dt value field empty the default value will be used, you must still click both buttons.'), 
            html.P ('After graph appears, if you want to change SNe or dt \
            values only click the submit button of value you changed!'),
            html.P('Enter a SNe type and dt value for one graph at a time, once one is done you can do it for the next one.'),
            html.P('Note: The average run time for setting up a graph is about 90 seconds.')], style={'color':'white', 'font-size': '1.8vw', 'text-align':'center'})]),

        html.Br(),
        html.Br(),

        #this section of the app layout is for the Filter Light Curve Graph
        html.H1('Supernovae Filter Light Curve Graph', style={'color': 'white'}),
        html.Div([html.Label(['SNe Type:'], style={'font-weight': 'bold', 'color':'white', 'background':'black'}), dcc.Input(id='SNe-FLtype', type='text', placeholder= 'Give SNe Type'), html.Button('Submit SNe Type', id='SNe-FLsub')], style={'width': '25%', 'display':'inline-block'}),

        html.Br(),
        html.Br(),

        #dcc loading gives the loading effect seen. Includes the filter dropdown lists inside
        dcc.Loading(id='load_FLdat', children=[
        html.Div(id='filter', style={'width': '25%', 'display':'inline-block'}),
        html.Div(id='fl-space', style={'width':'4px', 'height':'auto', 'display':'inline-block'}),
        html.Div(id='FLtemplate', style={'width':'25%', 'display':'inline-block'})], type="default"),

        html.Div(dcc.Graph(id="SNe-FL", figure=fig1)),

        html.Br(),
        html.Br(),

        #This section of the app layout is for the Color Light Curve Graph
        html.H1('Supernovae Color Light Curve Graph', style={'color':'white'}),
        html.Div([html.Label(['SNe Type:'], style={'font-weight': 'bold', 'color':'white', 'background':'black'}), dcc.Input(id='SNe-CLtype', type= 'text', placeholder= 'Give SNe Type'), html.Button('Submit SNe Type', id='SNe-CLsub')], style={'width': '25%', 'display': 'inline-block'}),

        html.Div([html.Label(['dt value:'], style={'font-weight': 'bold', 'color':'white'}), dcc.Input(id='dtCL', type='text', placeholder= 'empty = 0.15 value', value= ''), html.Button('Submit dt Value', id='dt-CLsub')], style={'right':'auto','width': '25%', 'display': 'inline-block', 'position':'absolute'}),

        html.Br(),
        html.Br(),

        #dcc loading gives the loading effect seen. Includes the color dropdown lists inside
        dcc.Loading(id='load_CLdat', children=[
        html.Div(id= 'color', style={'width': '25%', 'display': 'inline-block'}),
        html.Div(id='cl-space', style={'width':'4px', 'height':'auto', 'display':'inline-block'}),
        html.Div(id='CLtemplate', style={'width': '25%', 'display': 'inline-block'})], type="default"),

        #Makes the color light curve graph
        html.Div(dcc.Graph(id="SNe-CL", figure=fig2)),

        html.Br(),
        html.Br(),

        #This section of the app layout is for the Color Color Diagram:

        #The text field where you type in SNe name
        html.H1('Supernovae Color-Color Diagram', style={'color':'white'}),
        html.Div([html.Label(['SNe Type:'], style={'font-weight': 'bold', 'color':'white', 'background':'black'}), dcc.Input(id='SNe-CCDtype', type='text', placeholder= 'Give SNe Type'), html.Button('Submit SNe Type', id='SNe-sub')], style={'width': '25%', 'display': 'inline-block'}),
        #The text field where you give a dt range
        html.Div([html.Label(['dt value:'], style={'font-weight': 'bold', 'color':'white'}), dcc.Input(id='dt', type='text', placeholder= 'empty = 0.15 value', value= ''), html.Button('Submit dt Value', id='dt-sub')], style={'right':'auto','width': '25%', 'display': 'inline-block', 'position':'absolute'}),

        html.Br(),
        html.Br(),
        #dcc loading gives the loading effect seen. Includes the color dropdown lists inside
        dcc.Loading(id='load_dat', children=[
        html.Div(id= 'y-axis', style={'width': '25%', 'display': 'inline-block'}),
        html.Div(id='ccd-space1', style={'width':'4px', 'height':'auto', 'display':'inline-block'}),
        html.Div(id= 'x-axis', style={'width': '25%', 'display': 'inline-block'}),
        html.Div(id='ccd-space2', style={'width':'4px', 'height':'auto', 'display':'inline-block'}),
        html.Div(id='CCDtemplate', style={'width': '25%', 'display': 'inline-block'})], type="default"),

        #Makes the graph
        html.Div(dcc.Graph(id="SNe-scatter", figure=fig3))
    ])


'''
@app.callback is very important and is what makes dash run! They take inputs that are put into the function directly below it and
returns the outputs of that function to the specific dash html.Div referenced.
'''
@app.callback(
    [Output('filter', 'children'), Output('FLtemplate', 'children')],
    [Input('SNe-FLsub', 'n_clicks'), State('SNe-FLtype', 'value')])
def plotly_FLGraph(n_clicks, SNtype):
    global filter_dat, all_FLdat

    if n_clicks is None:
        raise PreventUpdate
    else:
        #pd.set_option('display.max_rows',1000)

        #reads in swift data to use and replaces Nans with white space ''.
        #swift= pd.read_csv('/home/Alexcrabt/mysite/NewSwiftSNweblist.csv')
        swift= pd.read_csv('NewSwiftSNweblist.csv')
        swift= swift.replace(np.nan, '', regex=True)

        #Itterates SNname, Dist_mod_cor, and type from NewSwiftSNweblist into graph_Data function and returns the data needed to graph.
        #plots_data= swift.apply(lambda row: CL_Data(row['SNname'], row['Distance_best'], row['type'], SNtype), axis=1)

        #Itterates SNname, Dist_mod_cor, and type from NewSwiftSNweblist into graph_Data function and returns the data needed to graph.
        plots_data= swift.apply(lambda row: FL_Data(row['SNname'], row['Distance_best'], row['type'], SNtype), axis=1)

        filter_dat= pd.DataFrame(columns=['UVW2', 'UVM2', 'UVW1', 'U', 'B', 'V'])

        all_dat=[]
        #This is the function that matches the data too the right filter column in the above dataframe.
        for i in range(len(plots_data)):
            if plots_data[i] == "NULL" or (len(plots_data[i]) < 6):
                continue
            elif all(i != i for i in plots_data):
                continue
            else:
                filter_dat.loc[len(filter_dat)]= plots_data[i]
                for j in range(len(plots_data[i])):
                    all_dat.append(plots_data[i][j])

        all_FLdat= pd.DataFrame({'All Filters': all_dat})

        '''
        Here I set up the label button for the color dropdown in the dash app and use a for loop to put all the color columns with data into it in the options section.
        value at the ends just specifies which color the program will start on.
        '''
    filter_FLoptions = pd.DataFrame(columns=['All Filters', 'UVW2', 'UVM2','UVW1','U','B', 'V'])

    themes=[["Basic Plotly","plotly"], ["Plotly White","plotly_white"], ["Plotly Dark","plotly_dark"], ["GGPlot2","ggplot2"], ["Seaborn", "seaborn"], ["Simple White", "simple_white"], ["No Theme","none"]]

    label1= html.Label(['Filters:'], style={'font-weight': 'bold', 'color':'white'}), dcc.Dropdown(id= 'filter-dropdown', options=[{'label': i, 'value': i} for i in filter_FLoptions],
            value='All Filters',
            searchable=False)

    label2= html.Label(['Plot Background Theme:'], style={'font-weight':'bold', 'color':'white'}), dcc.Dropdown(id= 'fl-template', options=[{'label': i[0], 'value': i[1]} for i in themes],value= 'plotly_dark',searchable=False)

    return label1, label2



@app.callback(
    [Output('color', 'children'), Output('CLtemplate', 'children')],
    [Input('SNe-CLsub', 'n_clicks'), Input('dt-CLsub', 'n_clicks'), State('SNe-CLtype', 'value'), State('dtCL', 'value')])
def plotly_CLGraph(n_clicks1, n_clicks2, SNtype, dt):
    global color_CLdat, all_CLdat

    if n_clicks1 is None or n_clicks2 is None:
        raise PreventUpdate
    else:
        #pd.set_option('display.max_rows',1000)

        #reads in swift data to use and replaces Nans with white space ''.
        #swift= pd.read_csv('/home/Alexcrabt/mysite/NewSwiftSNweblist.csv')
        swift= pd.read_csv('NewSwiftSNweblist.csv')
        swift= swift.replace(np.nan, '', regex=True)

        #Itterates SNname, Dist_mod_cor, and type from NewSwiftSNweblist into graph_Data function and returns the data needed to graph.
        plots_data= swift.apply(lambda row: CL_Data(row['SNname'], row['Distance_best'], row['type'], SNtype, dt), axis=1)

        #Sets up a pandas dataframe to put the data into their respective color columns.
        color_CLdat=pd.DataFrame(columns=['UVW2-UVM2', 'UVW2-UVW1', 'UVW2-U','UVW2-B','UVW2-V','UVM2-UVW1','UVM2-U','UVM2-B','UVM2-V','UVW1-U','UVW1-B','UVW1-V','U-B','U-V','B-V'])
        all_dat=[]
        #This is the function that matches the data too the right color column in the above dataframe.
        for i in range(len(plots_data)):
            if plots_data[i] == "NULL":
                continue
            elif all(i != i for i in plots_data[i]):
                continue
            else:
                color_CLdat.loc[len(color_CLdat)]= plots_data[i]
                for j in range(len(plots_data[i])):
                    all_dat.append(plots_data[i][j])

        all_CLdat= pd.DataFrame({'All Colors': all_dat})
        '''
        Here I set up the label button for the color dropdown in the dash app and use a for loop to put all the color columns with data into it in the options section.
        value at the ends just specifies which color the program will start on.
        '''
        color_CLoptions = pd.DataFrame(columns=['All Colors', 'UVW2-UVM2', 'UVW2-UVW1', 'UVW2-U','UVW2-B','UVW2-V','UVM2-UVW1','UVM2-U','UVM2-B','UVM2-V','UVW1-U','UVW1-B','UVW1-V','U-B','U-V','B-V'])

        themes=[["Basic Plotly","plotly"], ["Plotly White","plotly_white"], ["Plotly Dark","plotly_dark"], ["GGPlot2","ggplot2"], ["Seaborn", "seaborn"], ["Simple White", "simple_white"], ["No Theme","none"]]

        label1= html.Label(['Colors:'], style={'font-weight': 'bold', 'color':'white'}), dcc.Dropdown(id= 'color-dropdown', options=[{'label': i, 'value': i} for i in color_CLoptions],
            value='All Colors',
            searchable=False)

        label2= html.Label(['Plot Background Theme:'], style={'font-weight':'bold', 'color':'white'}), dcc.Dropdown(id= 'cl-template', options=[{'label': i[0], 'value': i[1]} for i in themes],value= 'plotly_dark',searchable=False)

        return label1, label2

@app.callback(
    [Output("y-axis", "children"), Output("x-axis", "children"), Output("CCDtemplate", "children")],
    [Input('SNe-sub', 'n_clicks'), Input('dt-sub', 'n_clicks'), State('SNe-CCDtype', 'value'), State('dt', 'value')])
def plotly_Graph(n_clicks1,n_clicks2, SNtype, dt):
    global color_CCDdat

    if n_clicks1 is None or n_clicks2 is None:
        raise PreventUpdate
    else:
        #pd.set_option('display.max_rows',1000)

        #reads in swift data to use and replaces Nans with white space ''.
        #swift= pd.read_csv('/home/Alexcrabt/mysite/NewSwiftSNweblist.csv')
        swift= pd.read_csv('NewSwiftSNweblist.csv')
        swift= swift.replace(np.nan, '', regex=True)

        #Itterates SNname, Dist_mod_cor, and type from NewSwiftSNweblist into graph_Data function and returns the data needed to graph.
        plots_data= swift.apply(lambda row: CCD_Data(row['SNname'], row['Distance_best'], row['type'], SNtype, dt), axis=1)

        #Sets up a pandas dataframe to put the data into their respective color columns.
        color_CCDdat=pd.DataFrame(columns=['UVW2-UVM2', 'UVW2-UVW1', 'UVW2-U','UVW2-B','UVW2-V','UVM2-UVW1','UVM2-U','UVM2-B','UVM2-V','UVW1-U','UVW1-B','UVW1-V','U-B','U-V','B-V'])

        #This is the function that matches the data too the right color column in the above dataframe.
        for i in range(len(plots_data)):
            if plots_data[i] == "NULL":
                continue
            elif all(i != i for i in plots_data[i]):
                continue
            else:
                color_CCDdat.loc[len(color_CCDdat)]= plots_data[i]

        '''
        Here I set up the label buttons for the color dropdowns in the dash app and use a for loop to put all the color columns with data into it in the options section.
        value at the ends just specifies which color the program will start on.
        '''
        themes=[["Basic Plotly","plotly"], ["Plotly White","plotly_white"], ["Plotly Dark","plotly_dark"], ["GGPlot2","ggplot2"], ["Seaborn", "seaborn"], ["Simple White", "simple_white"], ["No Theme","none"]]
        
        label1= html.Label(['Y-Axis Color:'], style={'font-weight': 'bold', 'color':'white'}),dcc.Dropdown(id= 'y-dropdown', options=[{'label': i, 'value': i} for i in color_CCDdat], value='UVW2-UVM2', searchable=False)
        
        label2= html.Label(['X-Axis Color:'], style={'font-weight': 'bold', 'color':'white'}),dcc.Dropdown(id= 'x-dropdown', options=[{'label': i, 'value': i} for i in color_CCDdat], value='UVW2-UVW1', searchable=False)
        
        label3= html.Label(['Plot Background Theme:'], style={'font-weight':'bold', 'color':'white'}), dcc.Dropdown(id= 'ccd-template', options=[{'label': i[0], 'value': i[1]} for i in themes],value= 'plotly_dark',searchable=False)
        
        return label1, label2, label3

#List of line styles to cycle through in update_Graph.
line_styles_names = ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']

#List of marker styles to cycle through in update_Graph.
symbols = [x for x in SymbolValidator().values[::-2] if not isinstance(x, int)]
symbols= [i for i in symbols if not i.isdigit()]
del symbols[8:12]
symbols_names = list(set([i.replace("-thin", "") for i in symbols]))
trace_color= px.colors.qualitative.Plotly

@app.callback(
    Output('SNe-FL', 'figure'),
    [Input('filter-dropdown', 'value'), Input('fl-template', 'value'), State('SNe-FLtype', 'value')]
)
def update_FLGraph(y, templ, SNtype):
    #Cycles through line styles
    line_style= cycle(line_styles_names)

    #Cycles through marker styles
    marker= cycle(symbols_names)

    c= cycle(trace_color)

    #Makes copy of color_dat dataframe so the original is not tampered with and can be reused.
    filter_data= filter_dat
    all_data= all_FLdat
    
    #sets up trace list.
    traces=[]

    if y != 'All Filters':
        #For loop makes and appends traces to traces list.
        for i in range(len(filter_data[y])):
            if pd.isna(filter_data[y][i]):
                continue
            else:
                mark= next(marker)
                line= next(line_style)
                color= next(c)

                traces.append(go.Scatter(x= filter_data[y][i][1], y= filter_data[y][i][2],  error_y= dict(array= filter_data[y][i][3], type= 'data', visible= False), marker= {'symbol': mark, 'color': color }, line={'dash': line, 'color': color}, mode= "lines+markers", name= filter_data[y][i][0]))
        
        layout= go.Layout(title= y + " Supernovae " + SNtype + " Light Curves", yaxis=dict(title="Absolute Magnitude", autorange="reversed"), xaxis=dict(title="Days since first detection"), height=625,template= templ, legend_title= y + " Supernovae", showlegend=True, updatemenus=[
        dict(
                type="buttons", showactive=False, xanchor="left", yanchor="top", x=0, y=1.08, buttons=list(
                    [   dict(
                            label="Error Bar",
                            method="update",
                            args=[{"error_y.visible": False}],
                            args2= [{"error_y.visible": True}],

                            )
                ]))])
    else:
        for i in range(len(all_data)):
            if pd.isna(all_data[y][i]):
                continue
            else:
                mark= next(marker)
                line= next(line_style)
                color= next(c)
                traces.append(go.Scatter(x= all_data[y][i][1], y= all_data[y][i][2],  error_y= dict(array= all_data[y][i][3], type= 'data', visible= False), marker= {'symbol': mark, 'color': color }, line={'dash': line, 'color': color}, mode= "lines+markers", name= all_data[y][i][0]+'_'+ all_data[y][i][4]))

        layout= go.Layout(title= "All Supernovae " + SNtype + " Filter Light Curves", yaxis=dict(title="Absolute Magnitude", autorange="reversed"), xaxis=dict(title="Days since first detection"), height=625,template= templ, legend_title= "Supernovae Type " + SNtype + "'s", showlegend=True, updatemenus=[
        dict(
                type="buttons", showactive=False, xanchor="left", yanchor="top", x=0, y=1.08, buttons=list(
                    [   dict(
                            label="Error Bar",
                            method="update",
                            args=[{"error_y.visible": False}],
                            args2= [{"error_y.visible": True}],

                            )
                ]))])

    return {'data': traces, 'layout':layout}

@app.callback(
    Output('SNe-CL', 'figure'),
    [Input('color-dropdown', 'value'),Input('cl-template', 'value'), State('dtCL', 'value'), State('SNe-CLtype', 'value')])
def update_CLGraph(y,templ, dt, SNtype):
    fig= make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0)

    #Edit axis orientation and set yaxis visibility and domain size
    fig['layout']['yaxis']['autorange'] = "reversed"
    fig['layout']['yaxis2']['autorange'] = "reversed"
    fig['layout']['yaxis3']['autorange'] = "reversed"

    #Edit axis labels
    fig['layout']['xaxis3']['title']='Days since first detection'

    #Cycles through line styles
    line_style= cycle(line_styles_names)

    #Cycles through marker styles
    marker= cycle(symbols_names)

    c= cycle(trace_color)
    #Makes copy of color_dat dataframe so the original is not tampered with and can be reused.
    color_data= color_CLdat
    all_data= all_CLdat

    #color_comb= [('UVW2-UVM2', 'UVW2', 'UVM2'), ('UVW2-UVW1', 'UVW2', 'UVW1'), ('UVW2-U', 'UVW2', 'U') ,('UVW2-B', 'UVW2', 'B') ,('UVW2-V', 'UVW2', 'V'), ('UVM2-UVW1', 'UVM2', 'UVW1') ,('UVM2-U', 'UVM2', 'U') ,('UVM2-B', 'UVM2', 'B'), ('UVM2-V', 'UVM2', 'V'), ('UVW1-U', 'UVW1', 'U') ,('UVW1-B', 'UVW1', 'B'), ('UVW1-V', 'UVW1', 'V'), ('U-B', 'U', 'B') , ('U-V', 'U', 'V'), ('B-V', 'B', 'V')]

    #sets up trace list.
    #traces=[]

    if y != 'All Colors':

        for i in range(len(color_data[y])):
            if pd.isna(color_data[y][i]):
                continue
            else:
                mark= next(marker)
                line= next(line_style)
                color= next(c)
                fig.add_trace(go.Scatter(x= color_data[y][i][1], y= color_data[y][i][2], error_y= dict(array= color_data[y][i][3], type= 'data', visible= False), marker= {'symbol': mark , 'color':color}, line={'dash': line, 'color':color}, mode= "lines+markers", visible=False, name= color_data[y][i][0]+'_'+color_data[y][i][10][0], legendgroup=color_data[y][i][0]+'_'+color_data[y][i][10][0] +'-'+color_data[y][i][10][1]), row=1, col=1)

                fig.add_trace(go.Scatter(x= color_data[y][i][4], y= color_data[y][i][5], error_y= dict(array= color_data[y][i][6], type= 'data', visible= False), marker= {'symbol': mark , 'color':color}, line={'dash': line, 'color':color}, mode= "lines+markers", visible=False, name= color_data[y][i][0]+'_'+color_data[y][i][10][1], legendgroup=color_data[y][i][0]+'_'+color_data[y][i][10][0] +'-'+color_data[y][i][10][1]), row=2, col=1)

                fig.add_trace(go.Scatter(x= color_data[y][i][7], y= color_data[y][i][8], error_y= dict(array= color_data[y][i][9], type= 'data', visible= False), marker= {'symbol': mark , 'color':color}, line={'dash': line, 'color':color}, mode= "lines+markers", visible=True,name= color_data[y][i][0]+'_'+color_data[y][i][10][0] +'-'+color_data[y][i][10][1], legendgroup=color_data[y][i][0]+'_'+color_data[y][i][10][0] +'-'+color_data[y][i][10][1]), row=3, col=1)

                filt= (color_data[y][i][10][0], color_data[y][i][10][1])

        #Layout sets up layout of the graph such as the error buttons.
        fig.update_layout(title= y+" Supernovae "+SNtype+ " Color Light Curves", xaxis3=dict(title="Days since first detection"), yaxis=dict(visible=False), yaxis2=dict(visible=False), yaxis3=dict(title= y, domain= [0,1]), legend_title= "Supernovae" + " Type "+ SNtype+"'s", height=625,template= templ, showlegend=True, updatemenus=[
            dict(
                    type="buttons", showactive=False, xanchor="left", yanchor="top", x=0.0, y=1.08, buttons=list(
                        [   dict(
                                label="Colors + Filters",
                                method="update",
                                args=[{"visible": [False, False, True]},{'title': y +" Supernovae "+ SNtype+" Color Light Curves", 'yaxis3.title': y, "yaxis.visible": False, "yaxis2.visible": False,  "yaxis3.domain":[0,1]}],
                                args2= [{"visible": [True, True, True]},{'title': y+" Supernovae "+SNtype+ " Light Curves", 'yaxis.title': filt[0] + ' Absolute Mag', 'yaxis2.title': filt[1] + ' Absolute Mag', 'yaxis3.title': y,"yaxis.visible": True, "yaxis2.visible": True,"yaxis.domain":[.65, 1], "yaxis2.domain": [.35,.6],"yaxis3.domain": [0,.3]}]

                                        )
                        ])),

            dict(
                    type="buttons", showactive=False, xanchor="left", yanchor="top", x=0.10, y=1.08, buttons=list(
                        [   dict(
                                label="Error Bar",
                                method="update",
                                args=[{"error_y.visible": False}],
                                args2= [{"error_y.visible": True}],

                                )
                        ]))])
    else:
        fig['layout']['yaxis']['visible']= False
        fig['layout']['yaxis2']['visible']= False
        fig['layout']['yaxis3']['domain']= [0,1]

        for i in range(len(all_data)):
            if pd.isna(all_data[y][i]):
                continue
            else:
                mark= next(marker)
                line= next(line_style)
                color= next(c)
                fig.add_trace(go.Scatter(x= all_data[y][i][1], y= all_data[y][i][2], error_y= dict(array= all_data[y][i][3], type= 'data', visible= False),visible=False, marker= {'symbol': mark , 'color':color}, line={'dash': line, 'color':color}, mode= "lines+markers", name= all_data[y][i][0]), row=1, col=1)

                fig.add_trace(go.Scatter(x= all_data[y][i][4], y= all_data[y][i][5], error_y= dict(array= all_CLdat[y][i][6], type= 'data', visible= False),visible=False, marker= {'symbol': mark , 'color':color}, line={'dash': line, 'color':color}, mode= "lines+markers", name= all_data[y][i][0]), row=2, col=1)

                fig.add_trace(go.Scatter(x= all_data[y][i][7], y= all_data[y][i][8], error_y= dict(array= all_data[y][i][9], type= 'data', visible= False), marker= {'symbol': mark , 'color':color}, line={'dash': line, 'color':color}, mode= "lines+markers", name= all_data[y][i][0]+'_'+all_data[y][i][10][0] +'-'+all_data[y][i][10][1]), row=3, col=1)

        #Layout sets up layout of the graph such as the error buttons.
        fig.update_layout(title= "All Supernovae " + SNtype+ " Color Light Curves", xaxis3=dict(title="Days since first detection"),yaxis3=dict(title="All Type "+ SNtype+ " Colors"), legend_title= "Supernovae" + " Type "+ SNtype+"'s", height=625,template= templ, showlegend=True, updatemenus=[

        dict(
                type="buttons", showactive=False, xanchor="left", yanchor="top", x=0.0, y=1.08, buttons=list(
                    [   dict(
                            label="Error Bar",
                            method="update",
                            args=[{"error_y.visible": False}],
                            args2= [{"error_y.visible": True}],

                            )
                    ]))])

    #returns traces and layout to update graph.
    return fig

@app.callback(
    Output("SNe-scatter", "figure"),
    [Input("y-dropdown", "value"), Input('x-dropdown', 'value'), Input('ccd-template', 'value'), State('dt', 'value'), State('SNe-CCDtype', 'value')])
def update_CCDGraph(y, x, templ, dt, SNtype):
    #Cycles through line styles
    line_style= cycle(line_styles_names)

    #Cycles through marker styles
    marker= cycle(symbols_names)

    #Makes copy of color_dat dataframe so the original is not tampered with and can be reused.
    color_data= color_CCDdat

    #sets up trace list.
    traces=[]

    #For loop makes and appends traces to traces list.
    for i in range(len(color_data[x])):
        if pd.isna(color_data[x][i]) or pd.isna(color_data[y][i]):
            continue
        else:
            color_c= mjd_Check(color_data[x][i][1], color_data[x][i][2], color_data[x][i][3], color_data[y][i][1], color_data[y][i][2], color_data[y][i][3], dt)

            traces.append(go.Scatter(x= color_c[1], y= color_c[3],  error_x= dict(array= color_c[2], type= 'data', visible= False), error_y= dict(array= color_c[4], type= 'data', visible= False), marker= {'symbol': next(marker) }, line={'dash': next(line_style)}, mode= "lines+markers", name= color_data[y][i][0]+'_'+SNtype))

    #Layout sets up layout of the graph such as the error buttons.
    layout= go.Layout(title= '('+y+')'+ '-' + '('+x+')' +' Diagram', yaxis=dict(title='('+y+')'), xaxis=dict(title='('+x+')'), height=625,template= templ, legend_title= "Supernovae", showlegend=True, updatemenus=[
        dict(
            type="dropdown", showactive=False, xanchor="left", yanchor="top", x=0, y=1.08, buttons=list(
                [dict(
                        label="Hide Error Bars",
                        method="update",
                        args=[{"error_x.visible": False,"error_y.visible": False}]
                            ),
                dict(
                        label= y + " & " + x + " Error Bars",
                        method="update",
                        args=[{"error_x.visible": True,"error_y.visible": True}]
                            ),
                dict(
                        label= y+ " Error Bars",
                        method="update",
                        args=[{"error_x.visible": False,"error_y.visible": True}]
                            ),
                dict(
                        label= x+ " Error Bars",
                        method="update",
                        args=[{"error_x.visible": True,"error_y.visible": False}]
                            )
                ]))])

    #returns traces and layout to update graph.
    return {'data': traces, 'layout':layout}

#To turn off debug feature just remove debug variable below.
if __name__ == "__main__": app.run_server(debug=True)

