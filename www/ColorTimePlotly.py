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
from itertools import cycle

def data_Grapher():
    def sn_Data(data):
        '''
        Grabs all data on the supernovae and returns it as a dataframe.
        '''
        sndata= pd.DataFrame(columns= ['filter', 'mjd', 'mag', 'magerr', 'exposure', 'elapse'])

        for sublist in data:
            if not sublist[0] == '#':
                sndata= sndata.append({'filter': sublist[0], 'mjd': sublist[1], 'mag': sublist[2], 'magerr': sublist[3],'exposure': sublist[10], 'elapse': sublist[11]}, ignore_index=True)

        sndata=sndata.sort_values(by="mjd", ascending= True, na_position='last')
        sndata.reset_index(drop=True, inplace=True)
        return sndata

    def filter_Data(sn_filter, mjd, mag, magerr):
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

    def null_Filter(mjd1_temp, mag1_temp, magerr1_temp, mjd2_temp, mag2_temp, magerr2_temp):
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
        if not len(mjd1) or not len(mag1) or not len(mjd2) or not len(mag2):
            raise print("Not enough data on " + snname)
        else:
            return mjd1, mag1, magerr1, mjd2, mag2, magerr2

    def mjd_Check(mjd1, mag1, magerr1, mjd2, mag2, magerr2):
        '''
        Checks if the mjd of bothe filters are in a certain ranage of each other, if they are
        it saves all data (mjd, mag, magerr) for them and returns them. 
        '''
        dt=.1500
        mjd1_s= [mjd-dt for mjd in mjd1]
        mjd1_p= [mjd+dt for mjd in mjd1]
        
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
        '''
        if not all(mjd_check) == True:
           ind= [i for i, x in enumerate(mjd_check) if not x]
           
           for j in ind:
               mjd1.pop(j)
               mag1.pop(j)
               magerr1.pop(j)

               mjd2.pop(j)
               mag2.pop(j)
               magerr2.pop(j)

               return mjd1, mag1, magerr1, mag2, magerr2

        else: 
            
            return mjd1, mag1, magerr1, mag2, magerr2
        '''
    def color_Val_Converter(mjd, mag1, magerr1, mag2, magerr2, dist_mod):
        ''' 
        Converts mjd to days_since_detection, both mags to color, and both magerrs to a final error
        '''
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

    def filter1_Val_Converter(mjd1, mag1, magerr1, dist_mod):
        ''' 
        Converts mjd to days_since_detection, and mag to absolute mag
        '''
        mjd_0= mjd1[0]
        days_since_detection= [mjd- mjd_0 for mjd in mjd1]
        abs_mag1= [mag - dist_mod for mag in mag1]
        magerr1= [abs(error) for error in magerr1]

        ''' 
        If days_since_detection or absolute_mag are nan lists it returns not enough data
        '''
        if np.isnan(days_since_detection).all() == True or np.isnan(abs_mag1).all() == True:
            raise print("Not enough data on  " + snname)
        else:
            return days_since_detection, abs_mag1, magerr1

    def filter2_Val_Converter(mjd2, mag2, magerr2, dist_mod):
        ''' 
        Converts mjd to days_since_detection, and mag to absolute mag
        '''
        mjd_0= mjd2[0]
        days_since_detection= [mjd- mjd_0 for mjd in mjd2]
        abs_mag2= [mag - dist_mod for mag in mag2]
        magerr2= [abs(error) for error in magerr2]

        ''' 
        If days_since_detection or absolute_mag are nan lists it returns not enough data
        '''
        if np.isnan(days_since_detection).all() == True or np.isnan(abs_mag2).all() == True:
            raise print("Not enough data on  " + snname)
        else:
            return days_since_detection, abs_mag2, magerr2


    snname = 'SN2012aw'
    #Grabs data file of given supernovae
    data= [i.strip().split() for i in open("./data/"+snname+ "_uvotB15.1.dat")]

    #Grabs specified data from datafile and saves as dataframe
    sndata= sn_Data(data)
    
    #Turns specified columns in dataframe to lists
    sn_filter= sndata['filter'].to_list()
    mjd= sndata['mjd'].to_list()
    mag= sndata['mag'].to_list()
    magerr= sndata['magerr'].to_list()

    #filters data based on given filter names
    mjd1, mag1, magerr1, mjd2, mag2, magerr2= filter_Data(sn_filter, mjd, mag, magerr)

    #Removes NAN/Null values
    mjd1, mag1, magerr1, mjd2, mag2, magerr2= null_Filter(mjd1, mag1, magerr1, mjd2, mag2, magerr2)
    

    #Makes copies of lists to more easily use later
    mjdf1= mjd1[:]
    magf1= mag1[:]
    magerrf1= magerr1[:] 
    mjdf2= mjd2[:]
    magf2= mag2[:]
    magerrf2= magerr2[:]

    #Checks if both filter mjds are in a certain range of each other
    mjd, mag1, magerr1, mag2, magerr2= mjd_Check(mjd1, mag1, magerr1, mjd2, mag2, magerr2)

    dist_mod= 30.1682355
    #Converts values for color light curve
    days, color, error= color_Val_Converter(mjd, mag1, magerr1, mag2, magerr2, dist_mod)
    #Coverts values for abs_mag light curve of filter1
    fil1_days, abs_mag1, magerr1= filter1_Val_Converter(mjdf1, magf1, magerrf1, dist_mod)
    #Coverts values for abs_mag light curve of filter2
    fil2_days, abs_mag2, magerr2= filter2_Val_Converter(mjdf2, magf2, magerrf2, dist_mod)

    #Makes filter2 light curve
    plot1= go.Scatter(x= fil1_days, y= abs_mag1, error_y= dict(array= magerr1, type= 'data', visible= False), mode= "lines+markers", name= snname+'_'+filter1)
    #Makes filter2 light curve
    plot2= go.Scatter(x= fil2_days, y= abs_mag2, error_y= dict(array= magerr2, type= 'data', visible= False), mode= "lines+markers", name= snname+'_'+filter2)
    #Makes color light curve
    plot3= go.Scatter(x= days, y= color, error_y= dict(array= error, type= 'data', visible= False), mode= "lines+markers", name= snname)

    return  plot1, plot2, plot3
        

def main():
    global filter1, filter2
    #snname= input('Please give a Supernovae Name')
    #filter1= input('Please specify the first filter')
    #filter2= input('Please specify the second filter')
    
    filter1 ='UVM2'
    filter2 = 'V'

    pd.set_option('display.max_rows', 1000)
    swift= pd.read_csv('NewSwiftSNweblist.csv')

    #Sets plotly background to dark. plotly_white gives white backgound
    pio.templates.default = "plotly_dark"

    #Sets up a 3 rowed subplot
    fig= make_subplots(rows=3, cols=1,shared_xaxes=True,vertical_spacing=0)

    #Runs data_Grapher() and returns 3 plots
    p1,p2, p3=data_Grapher()

    #adds those plots to fig
    fig.append_trace(p1, row=1, col=1)
    fig.append_trace(p2, row=2, col=1)
    fig.append_trace(p3, row=3, col=1)
    #dist_mod= 31.54400737

    fig.update_layout(legend_title= "Supernovae")

    #Edit axis orientation
    fig['layout']['yaxis']['autorange'] = "reversed"
    fig['layout']['yaxis2']['autorange'] = "reversed"

    #Edit axis labels
    fig['layout']['xaxis3']['title']='Days since first detection'

    fig['layout']['yaxis']['title']= filter1 + ' Absolute Mag'
    fig['layout']['yaxis2']['title']=filter2 + ' Absolute Mag'
    fig['layout']['yaxis3']['title']= filter1+ '-'+ filter2

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
    fig.update_layout(
        updatemenus=[
                dict(
                type="buttons", showactive=False, xanchor="left", yanchor="top", x=0.15, y=1.08, buttons=list(
                    [   dict(
                            label="Error Bar",
                            method="update",
                            args=[{"error_y.visible": False}],
                            args2= [{"error_y.visible": True}]
                                    )
                    ]
                ), 
            )
        ]
    )


    fig.show()
    '''
    def mjd_Check(mjd1, mjd2):
        dt=.1500
        mjd1_s= [mjd-dt for mjd in mjd1]
        mjd1_p= [mjd+dt for mjd in mjd1]
        
        mjd_check= [min(mjd1_s, mjd1_p) < mjd2 < max(mjd1_s, mjd1_p) for mjd1_s, mjd1_p, mjd2 in zip(mjd1_s, mjd1_p, mjd2)]
        print(mjd_check)

        if not all(mjd_check) == True:
           ind= [i for i, x in enumerate(mjd_check) if not x]
           
           for i in ind:
               mjd1.pop(i)
               mjd2.pop(i)
               
               mjd= mjd1 + mjd2
               mjd= sorted(mjd, key= float)
               return mjd

        else: 
            mjd= mjd1 + mjd2
            mjd= sorted(mjd, key= float)
            return mjd

    a= [1,2,5,6]
    b= [5, 6, 7, 8]
    c=[3, 4.8, 2,9]
    #Given 2 lists a and b we can subtract the values in the lists as such
    #We can put any operation in the function
    #print([min(a,b) < c < max(a, b) for a, b, c in zip(a, b, c)])

    '''

    
if __name__ == "__main__":  main()