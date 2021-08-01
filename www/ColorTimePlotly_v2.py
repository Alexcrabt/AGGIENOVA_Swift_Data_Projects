'''
@author: Alexander Crabtree
            acrabtree15@gmail.com
'''
import plotly.graph_objects as go
import plotly.io as pio
from plotly.validators.scatter.marker import SymbolValidator
from plotly.validators.scatter.line import DashValidator
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import itertools
from itertools import cycle

#Must have plotly version 5.1.0 or higher installed.
#To do so use pip install plotly==5.1.0

def data_Grapher(snname, dist_mod):
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

        return mjd1_temp, mag1_temp, magerr1_temp, mjd2_temp, mag2_temp, magerr2_temp

    def filt_Data(sn_filter,filter1, mjd, mag, magerr):
        '''
        creates temp lists and filters all the supernovae data
        based on the 2 given filters and returns all data just for thhose filters.
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
            raise print("Not enough data on " + snname +'_'+ colors[j][0]+ ' for ' + colors[j][0]+ '-'+ colors[j][1])
        elif not len(mjd2) or not len(mag2):
            raise print("Not enough data on " + snname +'_'+ colors[j][1]+ ' for ' + colors[j][0]+ '-'+ colors[j][1])
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
            raise print("Not enough data on " + snname + '_'+ filt1[j])
        else:
            return mjd, mag, magerr

    def mjd_Check(mjd1, mag1, magerr1, mjd2, mag2, magerr2):
        '''
        Checks if the mjd of bothe filters are in a certain ranage of each other, if they are
        it saves all data (mjd, mag, magerr) for them and returns them. 
        '''

        if not dt=="":
            mjd1_s= [mjd-float(dt) for mjd in mjd1]
            mjd1_p= [mjd+float(dt) for mjd in mjd1]
                
                #mjd_check= [min(mjd1_s, mjd1_p) < mjd2 < max(mjd1_s, mjd1_p) for mjd1_s, mjd1_p, mjd2 in zip(mjd1_s, mjd1_p, mjd2)]
                #_c stands for checked
            mjd_c=[]
            mag1_c= []
            magerr1_c= []

            mag2_c= []
            magerr2_c=[]
                
            for i in range(len(mjd1)):
                for j in range(len(mjd2)):
                    if min(mjd1_s[i], mjd1_p[i]) < mjd2[j] < max(mjd1_s[i], mjd1_p[i]):
                        mjd_c.append(mjd1[i])
                        mag1_c.append(mag1[i])
                        magerr1_c.append(magerr1[i])

                        mag2_c.append(mag2[j])
                        magerr2_c.append(magerr2[j])
                    else:
                        continue
            return mjd_c, mag1_c, magerr1_c, mag2_c, magerr2_c

        else:
            dt_def= 0.1500
            mjd1_s= [mjd-dt_def for mjd in mjd1]
            mjd1_p= [mjd+dt_def for mjd in mjd1]
                
                #mjd_check= [min(mjd1_s, mjd1_p) < mjd2 < max(mjd1_s, mjd1_p) for mjd1_s, mjd1_p, mjd2 in zip(mjd1_s, mjd1_p, mjd2)]
                #_c stands for checked
            mjd_c=[]
            mag1_c= []
            magerr1_c= []

            mag2_c= []
            magerr2_c=[]
                
            for i in range(len(mjd1)):
                for j in range(len(mjd2)):
                    if min(mjd1_s[i], mjd1_p[i]) < mjd2[j] < max(mjd1_s[i], mjd1_p[i]):
                        mjd_c.append(mjd1[i])
                        mag1_c.append(mag1[i])
                        magerr1_c.append(magerr1[i])

                        mag2_c.append(mag2[j])
                        magerr2_c.append(magerr2[j])
                    else:
                        continue
            return mjd_c, mag1_c, magerr1_c, mag2_c, magerr2_c
        
        
        
    def color_Val_Converter(mjd, mag1, magerr1, mag2, magerr2, dist_mod):
        ''' 
        Converts mjd to days_since_detection, both mags to color, and both magerrs to a final error
        '''
        if not len(mjd) or not len(mag1) or not len(mag2):
            raise print("Not enough data on  " + snname)
        else:
            mjd_0= mjd[0]
            days_since_detection= [mjd- mjd_0 for mjd in mjd]
            color= [mag_a - mag_b for mag_a, mag_b in zip(mag1, mag2)]
            error= [np.sqrt((error1)**2 + (error2)**2) for error1, error2 in zip(magerr1, magerr2)]

            ''' 
            If days_since_detection or color are nan lists it returns not enough data
            '''
            if np.isnan(days_since_detection).all() == True or np.isnan(color).all() == True:
                raise print("Not enough data on  " + snname)
            else:
                return days_since_detection, color, error

    def filt_Value_Converter(mjd, mag, magerr, dist_mod):
        ''' 
        Converts mjd to days_since_detection, mag to absolute_mag
        '''
        if not len(mjd) or not len(mag):
            raise print("Not enough data on  " + snname)
        else:
            mjd_0= mjd[0]
            days_since_detection= [mjd- mjd_0 for mjd in mjd]
            absolute_mag= [mag - dist_mod for mag in mag]
            final_magerr= [abs(error) for error in magerr]


            ''' 
            If days_since_detection or absolute_mag are nan lists it returns not enough data
            '''
            if np.isnan(days_since_detection).all() == True or np.isnan(absolute_mag).all() == True:
                raise print("Not enough data on  " + snname)
            else:
                return days_since_detection, absolute_mag, final_magerr

    def filter_Graphs_Comb(filt_dat):

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


    #snname = 'SN2012aw'
    #Grabs data file of given supernovae
    try:
        data= [i.strip().split() for i in open("./data/"+snname+ "_uvotB15.1.dat")]
    except:
        print(str(snname) + ' is not in our data files')
        return "NULL"

    try:
        #Grabs specified data from datafile and saves as dataframe
        sndata= sn_Data(data)
        
        #Turns specified columns in dataframe to lists
        sn_filter= sndata['filter'].to_list()
        mjd= sndata['mjd'].to_list()
        mag= sndata['mag'].to_list()
        magerr= sndata['magerr'].to_list()

        #filters data based on given filter names
        filt1= ['UVW2', 'UVM2', 'UVW1', 'U', 'B', 'V']
        filt2= ['UVM2', 'UVW1', 'U', 'B', 'V']

        i=0
        color_dat=[]
        colors=[]

        while i <= len(filt1):
            if not len(filt2):
                break
            else:
                for j in range(len(filt2)):
                    color_dat.append((color_Data(sn_filter, filt1[i], filt2[j], mjd, mag, magerr)))
                    colors.append((filt1[i], filt2[j]))
                filt2.pop(0)
                i= i+1
        
    except:
        print('1c')
        return "NULL"

    try:
        filt_dat=[]
        for i in range(len(filt1)):
            filt_dat.append(filt_Data(sn_filter, filt1[i], mjd, mag, magerr))

    except:
        print('1f')
        return "NULL"
    
    #Removes NAN/Null values
    for i in range(len(color_dat)):
        try:
            color_dat[i]= null_Color_Filter(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], color_dat[i][5], i)

        except:
            color_dat[i]= "NULL"
    
        #Removes NAN/Null values
    
    for i in range(len(filt_dat)):
        try:
            filt_dat[i]= null_Filter(filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], i)

        except:
            filt_dat[i]= "NULL"
    
        #Checks if both filter mjds are in a certain range of each other
    
    for i in range(len(color_dat)):
        if color_dat[i] == "NULL":
            continue
        else:
            color_dat[i]= mjd_Check(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], color_dat[i][5])

        #Converts values for color light curve
    for i in range(len(color_dat)):
        if color_dat[i] == "NULL":
            continue
        else:
            try:
                color_dat[i]= color_Val_Converter(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], dist_mod)
            except:
                color_dat[i]= "NULL"
    
        #Coverts values for abs_mag light curve of filter1
    for i in range(len(filt_dat)):
        if filt_dat[i] == "NULL":
            continue
        else:
            try:
                filt_dat[i]= filt_Value_Converter(filt_dat[i][0], filt_dat[i][1], filt_dat[i][2], dist_mod)
            except:
                filt_dat[i]= "NULL"
    
        
    filt_dat= filter_Graphs_Comb(filt_dat)
    plots=[]
    print(filt_dat)
    print(color_dat)
        
    for i in range(len(color_dat)):
        if color_dat[i] == "NULL":
            plots.append("NULL")
        elif filt_dat[i] == "NULL":
            plots.append("NULL")
        else:
            plots.append((go.Scatter(x= filt_dat[i][0], y= filt_dat[i][1], error_y= dict(array= filt_dat[i][2], type= 'data', visible= False),mode= "lines+markers", name= snname+'_'+colors[i][0], legendgroup=snname+'_'+colors[i][0]+'-'+colors[i][1], legendgrouptitle_text=snname+'_'+colors[i][0]+'-'+colors[i][1]),
                go.Scatter(x= filt_dat[i][3], y= filt_dat[i][4], error_y= dict(array= filt_dat[i][5], type= 'data', visible= False),mode= "lines+markers", name= snname+'_'+colors[i][1], legendgroup=snname+'_'+colors[i][0]+'-'+colors[i][1]),
                go.Scatter(x= color_dat[i][0], y= color_dat[i][1], error_y= dict(array= color_dat[i][2], type= 'data', visible= False),mode= "lines+markers", name= snname+'_'+colors[i][0]+'-'+colors[i][1], legendgroup=snname+'_'+colors[i][0]+'-'+colors[i][1])))

    return plots
        

def BtnBoolList(count):
        
    count= list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in count))
        
    w2_m2= [x.replace('a', '1').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    w2_m2= list(map(bool,list(map(int, w2_m2))))

    w2_w1= [x.replace('a', '0').replace('b', '1').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    w2_w1= list(map(bool,list(map(int, w2_w1))))
        
    w2_u= [x.replace('a', '0').replace('b', '0').replace('c', '1').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    w2_u= list(map(bool,list(map(int, w2_u))))

    w2_b= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '1').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    w2_b= list(map(bool,list(map(int, w2_b))))

    w2_v= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '1').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    w2_v= list(map(bool,list(map(int, w2_v))))

    m2_w1= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '1').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    m2_w1= list(map(bool,list(map(int, m2_w1))))

    m2_u= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '1').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    m2_u= list(map(bool,list(map(int, m2_u))))

    m2_b= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '1').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    m2_b= list(map(bool,list(map(int, m2_b))))

    m2_v= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '1').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    m2_v= list(map(bool,list(map(int, m2_v))))

    w1_u= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '1').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    w1_u= list(map(bool,list(map(int, w1_u))))

    w1_b= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '1').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    w1_b= list(map(bool,list(map(int, w1_b))))

    w1_v= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '1').replace('m', '0').replace('n', '0').replace('o', '0') for x in count]
    w1_v= list(map(bool,list(map(int, w1_v))))

    u_b= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '1').replace('n', '0').replace('o', '0') for x in count]
    u_b= list(map(bool,list(map(int, u_b))))

    u_v= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '1').replace('o', '0') for x in count]
    u_v= list(map(bool,list(map(int, u_v))))

    b_v= [x.replace('a', '0').replace('b', '0').replace('c', '0').replace('d', '0').replace('e', '0').replace('f', '0').replace('g', '0').replace('h', '0').replace('i', '0').replace('j', '0').replace('k', '0').replace('l', '0').replace('m', '0').replace('n', '0').replace('o', '1') for x in count]
    b_v= list(map(bool,list(map(int, b_v))))

    return w2_m2, w2_w1, w2_u, w2_b, w2_v, m2_w1, m2_u, m2_b, m2_v, w1_u, w1_b, w1_v, u_b, u_v, b_v

def main():
    global dt

    #snname= input('Please give a Supernovae Name: ')
    snname= 'SN2005ek'
    #filter1= input('Please specify the first filter: ')
    #filter2= input('Please specify the second filter: ')
    dt= input(('Please specify a dt value or click ENTER for default dt value: '))

    #dist_mod= 34.19162132
    
    #filter1 ='UVM2'
    #filter2 = 'V'

    pd.set_option('display.max_rows', 1000)
    swift= pd.read_csv('NewSwiftSNweblist.csv')

    #Sets plotly background to dark. plotly_white gives white backgound
    pio.templates.default = "plotly_dark"

    #Sets up a 3 rowed subplot
    fig= make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0)

    #Runs data_Grapher() and returns 3 plots
    plot_data= swift.apply(lambda row: data_Grapher(row['SNname'], row['Dist_mod_cor']), axis=1)
    #data_Grapher(snname, dist_mod)
    #plot_data= data_Grapher(snname, dist_mod)

    btn_count=[]
    alpha_list=['a','b', 'c', 'd', 'e','f', 'g', 'h', 'i', 'j', 'k','l', 'm', 'n', 'o']
    for i in range(len(plot_data)):
        if plot_data[i] == "NULL":
            continue
        else:
            for j in range(len(plot_data[i])):
                if plot_data[i][j] == "NULL":
                    continue
                else:
                    fig.add_trace(plot_data[i][j][0], row=1, col=1)
                    fig.add_trace(plot_data[i][j][1], row=2, col=1)
                    fig.add_trace(plot_data[i][j][2], row=3, col=1)
                    btn_count.append(alpha_list[j])

        #w2_w1= [x.replace('g1', '0').replace('g3', '0').replace('g4', '0').replace('g5', '0').replace('g6', '0').replace('g7', '0').replace('g8', '0').replace('g9', '0').replace('g10', '0').replace('g11', '0').replace('g12', '0').replace('g13', '0').replace('g14', '0').replace('g15', '0') for x in num]
    btns= pd.DataFrame()
    btn_list= BtnBoolList(btn_count)
    for i in range(len(btn_list)):
        btns[i]= btn_list[i]
    
    #btns= pd.DataFrame(BtnBoolList(btn_count), columns=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
    '''
    #adds those plots to fig
    for i in range(len(plot_data)):
        if i <= 4:
            fig.add_trace(plot_data[i], row=1, col=1)

        elif 6< i <= 11:
            fig.add_trace(plot_data[i],row=2, col=1)
        
        elif 11< i:
            fig.add_trace(plot_data[i], row=3, col=1)
        else:
            continue

    '''
    '''
    for i in range(len(plot_data)):
        if plot_data[i][0]=="NULL":
            continue
        else:
            fig.append_trace(plot_data[i][0], row=1, col=1)
    
    for i in range(len(plot_data)):
        if plot_data[i][1]=="NULL":
            continue
        else:
            fig.append_trace(plot_data[i][1], row=2, col=1)
    
    for i in range(len(plot_data)):
        if plot_data[i][2]=="NULL":
            continue
        else:
            fig.append_trace(plot_data[i][2], row=3, col=1)
    
    '''
    
    #fig.update_layout(title_text= filter1 + '-'+ filter2 + ' Color Light Curve', legend_title= "Supernovae")

    #Edit axis orientation
    fig['layout']['yaxis']['autorange'] = "reversed"
    fig['layout']['yaxis2']['autorange'] = "reversed"

    #Edit axis labels
    fig['layout']['xaxis3']['title']='Days since first detection'

    #cycles through diffferent line types from this list
    line_styles_names = ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
    line_styles = cycle(line_styles_names)

    #pulls out specific symbol types from plotly symbol list then cycle through those symbols
    symbols = [x for x in SymbolValidator().values[::-2] if not isinstance(x, int)]
    symbols= [i for i in symbols if not i.isdigit()]
    del symbols[8:12]
    symbols_names = list(set([i.replace("-thin", "") for i in symbols]))
    markers = cycle(symbols_names)

    #once one line and marker type is used on a trace, cycle to next one
    for d in fig.data:
        d.line["dash"] = next(line_styles)
        d.marker.symbol = next(markers)
        d.marker.size = 8
    #Error Button
    #fig['layout']['yaxis']['title']= ' Absolute Mag'
    #fig['layout']['yaxis2']['title']=' Absolute Mag'
    #fig['layout']['yaxis3']['title']= 'color'
    fig.update_layout(
        updatemenus=[
            dict(
                type="dropdown", showactive=False, xanchor="left", yanchor="top", x=0.05, y=1.08,active=0, buttons=list(
                    [   

                        dict(
                            label= "All",
                            method= "update",
                            args= [{'visible': [True, True]}, {'title': "All Supernovae Color Light Curves", 'yaxis.title': 'Absolute Mag', 'yaxis2.title': 'Absolute Mag', 'yaxis3.title': 'Color'}]), 
                        dict(
                            label= "UVW2-UVM2",
                            method= "update",
                            args= [{'visible': btns[0].to_list()}, {'title': "UVW2-UVM2 Supernovae Light Curves", 'yaxis.title': 'UVW2 Absolute Mag', 'yaxis2.title': 'UVM2 Absolute Mag', 'yaxis3.title': 'UVW2-UVM2'}]),
                        dict(
                            label= "UVW2-UVW1",
                            method="update",
                            args=[{'visible': btns[1].to_list()}, {'title': "UVW2-UVW1 Supernovae Light Curves", 'yaxis.title': 'UVW2 Absolute Mag', 'yaxis2.title': 'UVW1 Absolute Mag', 'yaxis3.title': 'UVW2-UVW1'}]),
                        dict(
                            label= "UVW2-U",
                            method="update",
                            args=[{'visible': btns[2].to_list()}, {'title': "UVW2-U Supernovae Light Curves", 'yaxis.title': 'UVW2 Absolute Mag', 'yaxis2.title': 'U Absolute Mag', 'yaxis3.title': 'UVW2-U'}]),
                        dict(
                            label= "UVW2-B",
                            method= "update",
                            args= [{'visible': btns[3].to_list()}, {'title': "UVW2-B Supernovae Light Curves", 'yaxis.title': 'UVW2 Absolute Mag', 'yaxis2.title': 'B Absolute Mag', 'yaxis3.title': 'UVW2-B'}]),
                        dict(
                            label= "UVW2-V",
                            method="update",
                            args=[{'visible': btns[4].to_list()}, {'title': "UVW2-V Supernovae Light Curves", 'yaxis.title': 'UVW2 Absolute Mag', 'yaxis2.title': 'V Absolute Mag', 'yaxis3.title': 'UVW2-V'}]),
                        dict(
                            label= "UVM2-UVW1",
                            method="update",
                            args=[{'visible': btns[5].to_list()}, {'title': "UVM2-UVW1 Supernovae Light Curves", 'yaxis.title': 'UVM2 Absolute Mag', 'yaxis2.title': 'UVW1 Absolute Mag', 'yaxis3.title': 'UVM2-UVW1'}]),
                        dict(
                            label= "UVM2-U",
                            method= "update",
                            args= [{'visible': btns[6].to_list()}, {'title': "UVM2-U Supernovae Light Curves", 'yaxis.title': 'UVM2 Absolute Mag', 'yaxis2.title': 'U Absolute Mag', 'yaxis3.title': 'UVM2-U'}]),
                        dict(
                            label= "UVM2-B",
                            method="update",
                            args=[{'visible': btns[7].to_list()}, {'title': "UVM2-B Supernovae Light Curves", 'yaxis.title': 'UVM2 Absolute Mag', 'yaxis2.title': 'B Absolute Mag', 'yaxis3.title': 'UVM2-B'}]),
                        dict(
                            label= "UVM2-V",
                            method="update",
                            args=[{'visible': btns[8].to_list()}, {'title': "UVM2-V Supernovae Light Curves", 'yaxis.title': 'UVM2 Absolute Mag', 'yaxis2.title': 'V Absolute Mag', 'yaxis3.title': 'UVM2-V'}]),
                        dict(
                            label= "UVW1-U",
                            method= "update",
                            args= [{'visible': btns[9].to_list()}, {'title': "UVW1-U Supernovae Light Curves", 'yaxis.title': 'UVW1 Absolute Mag', 'yaxis2.title': 'U Absolute Mag', 'yaxis3.title': 'UVW1-U'}]),
                        dict(
                            label= "UVW1-B",
                            method="update",
                            args=[{'visible': btns[10].to_list()}, {'title': "UVW1-B Supernovae Light Curves", 'yaxis.title': 'UVW1 Absolute Mag', 'yaxis2.title': 'B Absolute Mag', 'yaxis3.title': 'UVW1-B'}]),
                        dict(
                            label= "UVW1-V",
                            method="update",
                            args=[{'visible': btns[11].to_list()}, {'title': "UVW1-V Supernovae Light Curves", 'yaxis.title': 'UVW1 Absolute Mag', 'yaxis2.title': 'V Absolute Mag', 'yaxis3.title': 'UVW1-V'}]),
                        dict(
                            label= "U-B",
                            method= "update",
                            args= [{'visible': btns[12].to_list()}, {'title': "U-B Supernovae Light Curves", 'yaxis.title': 'U Absolute Mag', 'yaxis2.title': 'B Absolute Mag', 'yaxis3.title': 'U-B'}]),
                        dict(
                            label= "U-V",
                            method="update",
                            args=[{'visible': btns[13].to_list()}, {'title': "U-V Supernovae Light Curves", 'yaxis.title': 'U Absolute Mag', 'yaxis2.title': 'V Absolute Mag', 'yaxis3.title': 'U-V'}]),
                        dict(
                            label= "B-V",
                            method="update",
                            args=[{'visible': btns[14].to_list()}, {'title': "B-V Supernovae Light Curves", 'yaxis.title': 'B Absolute Mag', 'yaxis2.title': 'V Absolute Mag', 'yaxis3.title': 'B-V'}]),
                        
                    ]
                )
            ),
            dict(
                type="buttons", showactive=False, xanchor="left", yanchor="top", x=0.18, y=1.08, buttons=list(
                    [   dict(
                            label="Error Bar",
                            method="update",
                            args=[{"error_y.visible": False}],
                            args2= [{"error_y.visible": True}],
                            
                                    )
                    ]
                ), 
            )
        ]
    )


    fig.show()
    


    
if __name__ == "__main__":  main()