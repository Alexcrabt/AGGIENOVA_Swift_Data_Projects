'''
@author: Alexander Crabtree
            acrabtree15@gmail.com
'''
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
from plotly.validators.scatter.marker import SymbolValidator
from plotly.validators.scatter.line import DashValidator
from plotly.subplots import make_subplots
from plotly import tools
import pandas as pd
import numpy as np
import itertools
import random
from itertools import cycle
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output



#Must have plotly version 5.1.0 higher installed.
#To do so use pip install plotly==5.1.0

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

def mjd_Check(mjd1, mag1, magerr1, mjd2, mag2, magerr2):
    '''
    Checks if the mjd of bothe filters are in a certain ranage of each other, if they are
    it saves all data (mjd, mag, magerr) for them and returns them. 
    '''

    if not dt=="":
        mjd1_s= [mjd-float(dt) for mjd in mjd1]
        mjd1_p= [mjd+float(dt) for mjd in mjd1]
                

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
        color= [mag_a - mag_b for mag_a, mag_b in zip(mag1, mag2)]
        error= [np.sqrt((error1)**2 + (error2)**2) for error1, error2 in zip(magerr1, magerr2)]

        ''' 
        If days_since_detection or color are nan lists it returns not enough data
        '''
        if np.isnan(mjd).all() == True or np.isnan(color).all() == True:
            raise print("Not enough data on  " + snname)
        else:
            return mjd, color, error

def data_Grapher(snname, dist_mod, type):
    global colors, filt1
    
    #Grabs data file of given supernovae
    try:
        data= [i.strip().split() for i in open("./data/"+snname+ "_uvotB15.1.dat")]
    except:
        print(str(snname) + ' is not in our data files')
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
        print(str(snname) + ' is not in our data files')
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
                    print("Not enough data on " + snname+"_"+filt1[i] +"-"+filt2[j])
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
            color_dat[i]= mjd_Check(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], color_dat[i][5])

    #Converts values for color light curve
    for i in range(len(color_dat)):
        #if list has Null in it skips it
        if color_dat[i] == "NULL":
            continue
        else:
            try:
                color_dat[i]= color_Val_Converter(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], dist_mod)
            except:
                color_dat[i]= "NULL"
    
    
    plots_dat=[]
    print(color_dat)
    for i in range(len(color_dat)):
        #Chooses random color from plotly default colors to use a legengroup color
        
        #If color list has a null value skip graphing that as well as corresponding filter data
        if color_dat[i] == "NULL":
            plots_dat.append("NULL")
        #If filter list has a null value skip graphing that as well as corresponding color data
        else:
            if type.find(SNtype) == 0:
                plots_dat.append((snname, color_dat[i][0], color_dat[i][1], color_dat[i][2]))
            else:
            #If no SNe type match append Null
                plots_dat.append("NULL")
    
    return plots_dat, colors
    


def main():
    
    global SNtype, dt

    #Asks for SNe type and dt input, defaults to dt=.1500 if no dt specified
    #SNtype= input("What type supernovae do you want to plot? (Press ENTER if you want to graph all):")
    #dt= input('Please specify a dt value or click ENTER for default dt value: ')

    SNtype=""
    dt= ""
    pd.set_option('display.max_rows',1000)
    #reads in swift data to use
    swift= pd.read_csv('NewSwiftSNweblist.csv')

    #Sets plotly background to dark. plotly_white gives white backgound
    pio.templates.default = "plotly_dark"

    #Sets up a 3 rowed subplot
    fig=go.Figure()
    
    snname= "SN2008bo"
    dist_mod= float(31.57050702)
    type=""

    plots_dat, colors= data_Grapher(snname, dist_mod, type)

    color_dat=pd.DataFrame(columns=['UVW2-UVM2', 'UVW2-UVW1', 'UVW2-U','UVW2-B','UVW2-V','UVM2-UVW1','UVM2-U','UVM2-B','UVM2-V','UVW1-U','UVW1-B','UVW1-V','U-B','U-V','B-V'])

    btns= pd.DataFrame()
    for i in range(len(plots_dat)):
        if plots_dat[i] == "NULL":
            color_dat[colors[i][0]+ '-'+ colors[i][1]]= [[[], [], [], []]]
            continue
        else:
            color_dat[colors[i][0]+ '-'+ colors[i][1]]= [[plots_dat[i][0], plots_dat[i][1], plots_dat[i][2], plots_dat[i][3]]]

    trace_color= px.colors.qualitative.Plotly
    print(color_dat)
    line_styles_names = ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']

    symbols = [x for x in SymbolValidator().values[::-2] if not isinstance(x, int)]
    symbols= [i for i in symbols if not i.isdigit()]
    del symbols[8:12]
    symbols_names = list(set([i.replace("-thin", "") for i in symbols]))
    
    
    app = dash.Dash(__name__)

    app.layout= html.Div([

        html.Div([html.Label(['Y-Axis Color:'], style={'font-weight': 'bold'}),
            dcc.Dropdown(id= 'y-dropdown', options=[
                {'label': i, 'value': i} for i in color_dat
                ], multi=False, value= 'UVW2-UVM2'

        )], style={'width': '25%', 'display': 'inline-block'}),

        html.Div([html.Label(['X-Axis Color:'], style={'font-weight': 'bold'}),
            dcc.Dropdown(id= 'x-dropdown', options=[
                {'label': i, 'value': i} for i in color_dat
                ], multi=False, value= 'UVW2-UVW1'

        )], style={'width': '25%', 'display': 'inline-block'}), 
        
        dcc.Graph(id="SNe-scatter", figure=fig), 
    ])

    @app.callback(
        Output("SNe-scatter", "figure"),
        [Input("y-dropdown", "value"), Input('x-dropdown', 'value')])
    def update_graph(y, x):

        color_traces= random.choice(trace_color)

        line_style= random.choice(line_styles_names)

        marker= random.choice(symbols_names)

        color_data= color_dat
        
        traces=[]

        for i in range(len(color_data[x])):
            if not len(color_data[x][i][1]) or not len(color_data[x][i][2]) or not len(color_data[y][i][1]) or not len(color_data[y][i][2]):
                continue
            else:
                color_c= mjd_Check(color_data[x][i][1], color_data[x][i][2], color_data[x][i][3], color_data[y][i][1], color_data[y][i][2], color_data[y][i][3])
                traces.append(go.Scatter(x= color_c[1], y= color_c[3], error_y= dict(array= [], type= 'data', visible= False), marker= {'symbol': marker, 'color': color_traces}, line={'dash': line_style, 'color': color_traces}, mode= "lines+markers", name= color_data[y][i][0]+'_'+SNtype))

        layout= go.Layout(title= y+ '-' + x +' Diagram', height=650,template= "plotly_dark", legend_title= "Supernovae", showlegend=True)

        
        return {'data': traces, 'layout':layout}
    app.run_server(debug=True)
    
  
    
    '''
    #Saves graph as html file with this name format
    
    if SNtype == "":
        fig.write_html("All_SNe_ColorLight_Curves.html")
    else:
        fig.write_html("SN_type_"+SNtype+"_ColorLight_Curves.html")
    '''
    #fig.show()
if __name__ == "__main__":  main()
