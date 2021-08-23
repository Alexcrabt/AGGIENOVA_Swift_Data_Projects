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
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
#import dash_daq as daq

app = dash.Dash(__name__)

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

def mjd_Check(mjd1, mag1, magerr1, mjd2, mag2, magerr2, dt):
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
            mjd_temp=[]
            mag1_temp= []
            magerr1_temp= []
            mag2_temp= []
            magerr2_temp= []
            for j in range(len(mjd2)):
                if min(mjd1_s[i], mjd1_p[i]) < mjd2[j] < max(mjd1_s[i], mjd1_p[i]):
                    #mjd_c.append(mjd1[i])
                    mjd_temp.append(mjd1[i])
                    #mag1_c.append(mag1[i])
                    mag1_temp.append(mag1[i])
                    #magerr1_c.append(magerr1[i])
                    magerr1_temp.append(magerr1[i])

                    #mag2_c.append(mag2[j])
                    mag2_temp.append(mag2[j])
                    #magerr2_c.append(magerr2[j])
                    magerr2_temp.append(magerr2[j])
                else:
                    continue
            
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
                if snname == 'SN2011fe':
                    print(mjd_temp)
                    print(mag1_temp)
                    print(mag2_temp)
                '''

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
            mjd_temp=[]
            mag1_temp= []
            magerr1_temp= []
            mag2_temp= []
            magerr2_temp= []
            for j in range(len(mjd2)):
                if min(mjd1_s[i], mjd1_p[i]) < mjd2[j] < max(mjd1_s[i], mjd1_p[i]):
                    #mjd_c.append(mjd1[i])
                    mjd_temp.append(mjd1[i])
                    #mag1_c.append(mag1[i])
                    mag1_temp.append(mag1[i])
                    #magerr1_c.append(magerr1[i])
                    magerr1_temp.append(magerr1[i])

                    #mag2_c.append(mag2[j])
                    mag2_temp.append(mag2[j])
                    #magerr2_c.append(magerr2[j])
                    magerr2_temp.append(magerr2[j])
                else:
                    continue
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
                if snname == 'SN2011fe':
                    print(mjd_temp)
                    print(mag1_temp)
                    print(mag2_temp)
                '''
            
        return mjd_c, mag1_c, magerr1_c, mag2_c, magerr2_c
          
def color_Val_Converter(mjd, mag1, magerr1, mag2, magerr2, dist_mod):
    ''' 
    Converts mjd to days_since_detection, both mags to color, and both magerrs to a final error
    '''
    if not len(mjd) or not len(mag1) or not len(mag2):
        raise #print("Not enough data on  " + snname)
    else:
        color= [mag_a - mag_b for mag_a, mag_b in zip(mag1, mag2)]
        error= [np.sqrt((error1)**2 + (error2)**2) for error1, error2 in zip(magerr1, magerr2)]

        ''' 
        If days_since_detection or color are nan lists it returns not enough data
        '''
        if np.isnan(mjd).all() == True or np.isnan(color).all() == True:
            raise #print("Not enough data on  " + snname)
        else:
            return mjd, color, error

def data_Grapher(snname, dist_mod, type, SNtype, dt):
    global colors, filt1
    
    #Grabs data file of given supernovae
    try:
        data= [i.strip().split() for i in open("./data/"+snname+ "_uvotB15.1.dat")]
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
                color_dat[i]= color_Val_Converter(color_dat[i][0], color_dat[i][1], color_dat[i][2], color_dat[i][3], color_dat[i][4], dist_mod)
            except:
                color_dat[i]= "NULL"
    
    plots_dat= []
    #type_temp= pd.DataFrame(columns=['type'])
    #type_temp['type']= type.str.find(SNtype)
    #print(type_temp)

    for i in range(len(color_dat)):
        #Chooses random color from plotly default colors to use a legengroup color
        
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

fig=go.Figure()
fig.update_layout(height=625,template= "plotly_dark")

app.layout= html.Div([
        html.Div([html.Label(['SNe Type:'], style={'font-weight': 'bold'}), dcc.Input(id='SNe-type', type='text', placeholder= 'Give SNe Type'), html.Button('Submit SNe Type', id='SNe-sub')], style={'width': '25%', 'display': 'inline-block'}),
        html.Div([html.Label(['dt value:'], style={'font-weight': 'bold'}), dcc.Input(id='dt', type='text', placeholder= 'empty = 0.15 value', value= ''), html.Button('Submit dt Value', id='dt-sub')], style={'right':'auto','width': '25%', 'display': 'inline-block', 'position':'absolute'}),
        html.Br(), 
        html.Div(id= 'y-axis', style={'width': '25%', 'display': 'inline-block'}),

        html.Div(id= 'x-axis', style={'width': '25%', 'display': 'inline-block'}),

        #html.Div([
            #daq.BooleanSwitch(id='error-btn', on=False, label="Error Bars:", labelPosition='top', style={'top':'9px', 'right':'auto', 'width':'6%','position': 'absolute', 'font-weight': 'bold'})], style={'display': 'inline-block'}),    
        
        dcc.Graph(id="SNe-scatter", figure=fig) 
    ])

@app.callback(
    [Output("y-axis", "children"), Output("x-axis", "children")], 
    [Input('SNe-sub', 'n_clicks'), Input('dt-sub', 'n_clicks'), State('SNe-type', 'value'), State('dt', 'value')])
def main(n_clicks1,n_clicks2, SNtype, dt):
    
    global color_dat

    if n_clicks1 is None or n_clicks2 is None:
        raise PreventUpdate
    else:
        #Asks for SNe type and dt input, defaults to dt=.1500 if no dt specified
        #SNtype= input("What type supernovae do you want to plot? (Press ENTER if you want to graph all):")
        #dt= input('Please specify a dt value or click ENTER for default dt value: ')
        #SNtype="Ib"
        #dt=""

        pd.set_option('display.max_rows',1000)
        #reads in swift data to use
        swift= pd.read_csv('NewSwiftSNweblist.csv')
        swift= swift.replace(np.nan, '', regex=True)

        #Sets plotly background to dark. plotly_white gives white backgound

        #Sets up a 3 rowed subplot
        
        
        #snname= "SN2008bo"
        #dist_mod= float(31.57050702)
        #type=""

        
        plots_data= swift.apply(lambda row: data_Grapher(row['SNname'], row['Dist_mod_cor'], row['type'], SNtype, dt), axis=1)

        color_dat=pd.DataFrame(columns=['UVW2-UVM2', 'UVW2-UVW1', 'UVW2-U','UVW2-B','UVW2-V','UVM2-UVW1','UVM2-U','UVM2-B','UVM2-V','UVW1-U','UVW1-B','UVW1-V','U-B','U-V','B-V'])

        for i in range(len(plots_data)):
            if plots_data[i] == "NULL":
                continue
            elif all(i != i for i in plots_data[i]):
                continue
            else:
                color_dat.loc[len(color_dat)]= plots_data[i]

        label1= html.Label(['Y-Axis Color:'], style={'font-weight': 'bold'}),dcc.Dropdown(id= 'y-dropdown', options=[{'label': i, 'value': i} for i in color_dat], value='UVW2-UVM2')
        label2= html.Label(['X-Axis Color:'], style={'font-weight': 'bold'}),dcc.Dropdown(id= 'x-dropdown', options=[{'label': i, 'value': i} for i in color_dat], value='UVW2-UVW1')
        return label1, label2
    
    #trace_color= px.colors.qualitative.Plotly

line_styles_names = ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']

symbols = [x for x in SymbolValidator().values[::-2] if not isinstance(x, int)]
symbols= [i for i in symbols if not i.isdigit()]
del symbols[8:12]
symbols_names = list(set([i.replace("-thin", "") for i in symbols]))
    
    

@app.callback(
    Output("SNe-scatter", "figure"),
    [Input("y-dropdown", "value"), Input('x-dropdown', 'value'), State('dt', 'value'), State('SNe-type', 'value')])
def update_graph(y, x, dt, SNtype):

    #color_traces= cycle(trace_color)

    line_style= cycle(line_styles_names)

    marker= cycle(symbols_names)

    color_data= color_dat
        
    traces=[]

    for i in range(len(color_data[x])):
        if pd.isna(color_data[x][i]) or pd.isna(color_data[y][i]):
            continue
        else:
            color_c= mjd_Check(color_data[x][i][1], color_data[x][i][2], color_data[x][i][3], color_data[y][i][1], color_data[y][i][2], color_data[y][i][3], dt)

            traces.append(go.Scatter(x= color_c[1], y= color_c[3],  error_x= dict(array= color_c[2], type= 'data', visible= False), error_y= dict(array= color_c[4], type= 'data', visible= False), marker= {'symbol': next(marker) }, line={'dash': next(line_style)}, mode= "lines+markers", name= color_data[y][i][0]+'_'+SNtype))
            
    layout= go.Layout(title= '('+y+')'+ '-' + '('+x+')' +' Diagram', yaxis=dict(title='('+y+')'), xaxis=dict(title='('+x+')'), height=625,template= "plotly_dark", legend_title= "Supernovae", showlegend=True, updatemenus=[
        dict(
            type="dropdown", showactive=False, xanchor="left", yanchor="top", x=0, y=1.08, buttons=list(
                [dict(
                        label="Hide Error Bars",
                        method="update",
                        args=[{"error_x.visible": False,"error_y.visible": False}]
                            ),
                dict(
                        label="Y & X-Axis Error Bars",
                        method="update",
                        args=[{"error_x.visible": True,"error_y.visible": True}]
                            ),
                dict(
                        label="Y-Axis Error Bars",
                        method="update",
                        args=[{"error_x.visible": False,"error_y.visible": True}]
                            ),
                dict(
                        label="X-Axis Error Bars",
                        method="update",
                        args=[{"error_x.visible": True,"error_y.visible": False}]
                            )
                ])),])
        
        
    return {'data': traces, 'layout':layout}
    
    
if __name__ == "__main__": app.run_server(port=8000, host='127.0.0.1',debug=True)

