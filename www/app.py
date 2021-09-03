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
          
def color_Val_Converter(mjd, mag1, magerr1, mag2, magerr2, dist_mod):
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

def graph_Data(snname, dist_mod, type, SNtype, dt):
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
app = dash.Dash(__name__, external_stylesheets= external_stylesheets)
fig=go.Figure()
fig.update_layout(height=625,template= "plotly_dark")
fig.add_layout_image(
    dict(
        source="https://pbrown801.github.io/SOUSA/www/sousa_galaxy.png",
        xref="paper", yref="paper",
        x=-.058, y=1.04,
        sizex=0.2, sizey=0.18,
        xanchor="left", yanchor="bottom"
    ))

app.layout= html.Div([
        html.P('Type in preferred SNe Type and dt value, then click both submit buttons. If dt value is left empty default value will be used, you must still click both buttons.', style={'color':'white'}),
        html.P('After graph appears, if you want to change SNe or dt values only click the submit button of value you changed!', style={'color':'white'}),

        #The text field where you type in SNe name
        html.Div([html.Label(['SNe Type:'], style={'font-weight': 'bold', 'color':'white', 'background':'black'}), dcc.Input(id='SNe-type', type='text', placeholder= 'Give SNe Type'), html.Button('Submit SNe Type', id='SNe-sub')], style={'width': '25%', 'display': 'inline-block'}),
        #The text field where you give a dt range
        html.Div([html.Label(['dt value:'], style={'font-weight': 'bold', 'color':'white'}), dcc.Input(id='dt', type='text', placeholder= 'empty = 0.15 value', value= ''), html.Button('Submit dt Value', id='dt-sub')], style={'right':'auto','width': '25%', 'display': 'inline-block', 'position':'absolute'}),

        html.Br(),
        html.Br(),
        #dcc loading gives the loading effect seen. Includes the color dropdown lists inside
        dcc.Loading(id='load_dat', children=[
        html.Div(id= 'y-axis', style={'width': '25%', 'display': 'inline-block'}),
        html.Div(id= 'x-axis', style={'width': '25%', 'display': 'inline-block'})], type="default"),
        
        #Makes the graph
        dcc.Graph(id="SNe-scatter", figure=fig) 
    ])

''' 
@app.callback is very important and is what makes dash run! They take inputs that are put into the function directly below it and
returns the outputs of that function to the specific dash html.Div referenced.
'''

@app.callback(
    [Output("y-axis", "children"), Output("x-axis", "children")], 
    [Input('SNe-sub', 'n_clicks'), Input('dt-sub', 'n_clicks'), State('SNe-type', 'value'), State('dt', 'value')])
def plotly_Graph(n_clicks1,n_clicks2, SNtype, dt):
    
    global color_dat

    if n_clicks1 is None or n_clicks2 is None:
        raise PreventUpdate
    else:
        #pd.set_option('display.max_rows',1000)

        #reads in swift data to use and replaces Nans with white space ''.
        swift= pd.read_csv('NewSwiftSNweblist.csv')
        swift= swift.replace(np.nan, '', regex=True)
        
        #Itterates SNname, Dist_mod_cor, and type from NewSwiftSNweblist into graph_Data function and returns the data needed to graph.
        plots_data= swift.apply(lambda row: graph_Data(row['SNname'], row['Dist_mod_cor'], row['type'], SNtype, dt), axis=1)

        #Sets up a pandas dataframe to put the data into their respective color columns.
        color_dat=pd.DataFrame(columns=['UVW2-UVM2', 'UVW2-UVW1', 'UVW2-U','UVW2-B','UVW2-V','UVM2-UVW1','UVM2-U','UVM2-B','UVM2-V','UVW1-U','UVW1-B','UVW1-V','U-B','U-V','B-V'])

        #This is the function that matches the data too the right color column in the above dataframe.
        for i in range(len(plots_data)):
            if plots_data[i] == "NULL":
                continue
            elif all(i != i for i in plots_data[i]):
                continue
            else:
                color_dat.loc[len(color_dat)]= plots_data[i]
        
        ''' 
        Here I set up the label buttons for the color dropdowns in the dash app and us a for loop to put all the color columns with data into it in the options section.
        value at the ends just specifies which color the program will start on.
        '''
        label1= html.Label(['Y-Axis Color:'], style={'font-weight': 'bold', 'color':'white'}),dcc.Dropdown(id= 'y-dropdown', options=[{'label': i, 'value': i} for i in color_dat], value='UVW2-UVM2')
        label2= html.Label(['X-Axis Color:'], style={'font-weight': 'bold', 'color':'white'}),dcc.Dropdown(id= 'x-dropdown', options=[{'label': i, 'value': i} for i in color_dat], value='UVW2-UVW1')
        return label1, label2
    
#List of line styles to cycle through in update_Graph.
line_styles_names = ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']

#List of marker styles to cycle through in update_Graph.
symbols = [x for x in SymbolValidator().values[::-2] if not isinstance(x, int)]
symbols= [i for i in symbols if not i.isdigit()]
del symbols[8:12]
symbols_names = list(set([i.replace("-thin", "") for i in symbols]))
    
    

@app.callback(
    Output("SNe-scatter", "figure"),
    [Input("y-dropdown", "value"), Input('x-dropdown', 'value'), State('dt', 'value'), State('SNe-type', 'value')])
def update_Graph(y, x, dt, SNtype):
    #Cycles through line styles
    line_style= cycle(line_styles_names)

    #Cycles through marker styles
    marker= cycle(symbols_names)

    #Makes copy of color_dat dataframe so the original is not tampered with and can be reused.
    color_data= color_dat
    
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
                ]))], 
        images=[
            dict(
        source="https://pbrown801.github.io/SOUSA/www/sousa_galaxy.png",
        xref="paper", yref="paper",
        x=-.06, y=1.012,
        sizex=0.2, sizey=0.18,
        xanchor="left", yanchor="bottom"
    )])
    
    #returns traces and layout to update graph.
    return {'data': traces, 'layout':layout}

#To turn off debug feature just remove debug variable below.
if __name__ == "__main__": app.run_server(port=7000, host='127.0.0.1',debug=True)

