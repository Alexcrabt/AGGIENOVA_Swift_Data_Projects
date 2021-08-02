'''
@author: Alexander Crabtree
            acrabtree15@gmail.com
'''

import plotly.graph_objects as go
import plotly.io as pio
from plotly.validators.scatter.marker import SymbolValidator
from plotly.validators.scatter.line import DashValidator
import pandas as pd
import numpy as np
from itertools import cycle



def Data_Grapher(snname, dist_mod, type):

    def File_Data(data, search):
        ''' 
        Uses filter name to search for mjd, mag, and magerr data for specific filter
        from the data found.
        '''
        mjd_temp= []
        mag_temp= []
        magerr_temp= []

        for sublist in data:
            if sublist[0] == search:
                mjd_temp.append(sublist[1])
                mag_temp.append(sublist[2])
                magerr_temp.append(sublist[3])

        ''' 
        If the mjd or mag lists are empty returns not enough data.
        '''
        if not len(mjd_temp) or not len(mag_temp):
            raise print("Not enough " +search+" data on " + snname)
        else:
            return mjd_temp, mag_temp, magerr_temp
    
    def Null_Filter(mjd_temp, mag_temp, magerr_temp):
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
                mjd.append(mjd_temp[i])
                mag.append(mag_temp[i])
                magerr.append(magerr_temp[i])
        
        ''' 
        If mjd or mag lists are empty after filtering or all NULL values return not enough data.
        '''
        if not len(mjd) or not len(mag):
            raise print("Not enough data on " + snname)
        else:
            return mjd, mag, magerr

    def Value_Converter(mjd, mag, magerr, dist_mod):
        ''' 
        Converts mjd to days_since_detection, mag to absolute_mag
        '''
        days_since_detection= []
        absolute_mag= []
        final_magerr= []
        
        for j in range(len(mjd)):
            days_since_detection.append(float(mjd[j]) - float(mjd[0]))
            absolute_mag.append(float(mag[j]) - dist_mod)
            final_magerr.append(abs(float(magerr[j])))

        ''' 
        If days_since_detection or absolute_mag are nan lists it returns not enough data
        '''
        if np.isnan(days_since_detection).all() == True or np.isnan(absolute_mag).all() == True:
            raise print("Not enough data on  " + snname)
        else:
            return days_since_detection, absolute_mag, final_magerr

    def dat_w2(dataContent, search1):
        ''' 
        Collects UVW2 filter data then creates plotly trace for it and returns it.
        '''
        try:
            mjd_w2, mag_w2, magerr_w2= File_Data(dataContent, search1) 

            mjd_w2, mag_w2, magerr_w2= Null_Filter(mjd_w2, mag_w2, magerr_w2)

            mjd_w2, mag_w2, magerr_w2= Value_Converter(mjd_w2, mag_w2, magerr_w2, dist_mod)
        
            if type.find(SNtype) == 0:
                uvw2= go.Scatter(x= mjd_w2, y= mag_w2, error_y= dict(array= magerr_w2, type= 'data', visible= False), mode= "lines+markers", name= snname +'_W2')
                return uvw2
            else:
                return "NULL"
        
        except:
            
            return "NULL"

    def dat_m2(dataContent, search2):
        ''' 
        Collects UVM2 filter data then creates plotly trace for it and returns it.
        '''
        try:
            mjd_m2, mag_m2, magerr_m2= File_Data(dataContent, search2)

            mjd_m2, mag_m2, magerr_m2= Null_Filter(mjd_m2, mag_m2, magerr_m2)

            mjd_m2, mag_m2, magerr_m2= Value_Converter(mjd_m2, mag_m2, magerr_m2, dist_mod)
        
            if type.find(SNtype) == 0:
                uvm2= go.Scatter(x= mjd_m2, y= mag_m2, error_y= dict(array= magerr_m2, type= 'data', visible= False), mode= "lines+markers", name= snname +'_M2')
                return uvm2
            else:
                return "NULL"
            
        except:
            
            return "NULL"
    
    def dat_w1(dataContent, search3):
        ''' 
        Collects UVW1 filter data then creates plotly trace for it and returns it.
        '''
        try:
            mjd_w1, mag_w1, magerr_w1= File_Data(dataContent, search3)

            mjd_w1, mag_w1, magerr_w1= Null_Filter(mjd_w1, mag_w1, magerr_w1)

            mjd_w1, mag_w1, magerr_w1= Value_Converter(mjd_w1, mag_w1, magerr_w1, dist_mod)
            
            if type.find(SNtype) == 0:
                uvw1= go.Scatter(x= mjd_w1, y= mag_w1, error_y= dict(array= magerr_w1, type= 'data', visible= False), mode= "lines+markers", name= snname +'_W1')
                return uvw1
            else:
                return "NULL"
            
        except:
            
            return "NULL"

    def dat_u(dataContent, search4):
        ''' 
        Collects U filter data then creates plotly trace for it and returns it.
        '''
        try:
            mjd_u, mag_u, magerr_u= File_Data(dataContent, search4)

            mjd_u, mag_u, magerr_u= Null_Filter(mjd_u, mag_u, magerr_u)

            mjd_u, mag_u, magerr_u= Value_Converter(mjd_u, mag_u, magerr_u, dist_mod)
            
            if type.find(SNtype) == 0:
                u= go.Scatter(x= mjd_u, y= mag_u, error_y= dict(array= magerr_u, type= 'data', visible= False), mode= "lines+markers", name= snname +'_U')
                return u
            else:
                return "NULL"
            
        except:
            
            return "NULL"

    def dat_b(dataContent, search5):
        ''' 
        Collects B filter data then creates plotly trace for it and returns it.
        '''
        try:
            mjd_b, mag_b, magerr_b= File_Data(dataContent, search5)

            mjd_b, mag_b, magerr_b= Null_Filter(mjd_b, mag_b, magerr_b)

            mjd_b, mag_b, magerr_b= Value_Converter(mjd_b, mag_b, magerr_b, dist_mod)
            
            if type.find(SNtype) == 0:
                b= go.Scatter(x= mjd_b, y= mag_b, error_y= dict(array= magerr_b, type= 'data', visible= False), mode= "lines+markers", name= snname +'_B')
                return b
            else:
                return "NULL"
            
        except:
            
            return "NULL"

    def dat_v(dataContent, search6):
        ''' 
        Collects V filter data then creates plotly trace for it and returns it.
        '''
        try:
            mjd_v, mag_v, magerr_v= File_Data(dataContent, search6)

            mjd_v, mag_v, magerr_v= Null_Filter(mjd_v, mag_v, magerr_v)

            mjd_v, mag_v, magerr_v= Value_Converter(mjd_v, mag_v, magerr_v, dist_mod)
            
            if type.find(SNtype) == 0:
                v= go.Scatter(x= mjd_v, y= mag_v, error_y= dict(array= magerr_v, type= 'data', visible= False), mode= "lines+markers", name= snname +'_V')
                return v
            else:
                return "NULL"
            
        except:
            
            return "NULL"

    '''
    Trys to pull data on supernovae given from data file, if it doesn't exist
    returns that the supernovae isn't in our files.
    '''
    try:
        dataContent= [i.strip().split() for i in open("./data/"+snname+ "_uvotB15.1.dat")]
    except:
        print(str(snname) + ' is not in our data files')
        return "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"

    search1= 'UVW2'
    search2= 'UVM2'
    search3= 'UVW1'
    search4= 'U'
    search5= 'B'
    search6= 'V'
    '''
    Runs the data on the given supernoava and filter name through the appropriate funcions.
    '''
    w2= dat_w2(dataContent, search1)
    m2= dat_m2(dataContent, search2)
    w1= dat_w1(dataContent, search3)
    u= dat_u(dataContent, search4)
    b= dat_b(dataContent, search5)
    v= dat_v(dataContent, search6)

    ''' 
    Returns trace data on all 6 filters
    '''
    return w2, m2, w1, u, b, v

def BtnBoolList(num):
    '''
    Depeding on the filter it converts the numbers from numList to 0 or 1 then
    converts that to a bool value and returns the filter bool lists.
    '''
    w2= []
    m2= []
    w1= []
    u= []
    b= []
    v= []

    w2= [i.replace("2", "0") for i in num]
    w2= [i.replace("3", "0") for i in w2]
    w2= [i.replace("4", "0") for i in w2]
    w2= [i.replace("5", "0") for i in w2]
    w2= [i.replace("6", "0") for i in w2]
    w2= list(map(int, w2))
    w2= list(map(bool, w2))

    m2= [i.replace("1", "0") for i in num]
    m2= [i.replace("3", "0") for i in m2]
    m2= [i.replace("4", "0") for i in m2]
    m2= [i.replace("5", "0") for i in m2]
    m2= [i.replace("6", "0") for i in m2]
    m2= list(map(int, m2))
    m2= list(map(bool, m2))

    w1= [i.replace("1", "0") for i in num]
    w1= [i.replace("2", "0") for i in w1]
    w1= [i.replace("4", "0") for i in w1]
    w1= [i.replace("5", "0") for i in w1]
    w1= [i.replace("6", "0") for i in w1]
    w1= list(map(int, w1))
    w1= list(map(bool, w1))

    u= [i.replace("1", "0") for i in num]
    u= [i.replace("2", "0") for i in u]
    u= [i.replace("3", "0") for i in u]
    u= [i.replace("5", "0") for i in u]
    u= [i.replace("6", "0") for i in u]
    u = list(map(int, u))
    u= list(map(bool, u))

    b= [i.replace("1", "0") for i in num]
    b= [i.replace("2", "0") for i in b]
    b= [i.replace("3", "0") for i in b]
    b= [i.replace("4", "0") for i in b]
    b= [i.replace("6", "0") for i in b]
    b= list(map(int, b))
    b= list(map(bool, b))

    v= [i.replace("1", "0") for i in num]
    v= [i.replace("2", "0") for i in v]
    v= [i.replace("3", "0") for i in v]
    v= [i.replace("4", "0") for i in v]
    v= [i.replace("5", "0") for i in v]
    v= list(map(int, v))
    v= list(map(bool, v))

    return w2, m2, w1, u, b, v

def main():
    #Makes the Supernovae type global for easy usage
    global SNtype

    #Reads in SwiftSNweblist
    swift= pd.read_csv("NewSwiftSNweblist.csv")

    #Takes given Supernovae type input
    SNtype= input("What type supernovae do you want to plot? (Press ENTER if you want to graph all):")

    #Sets plotly background to dark. plotly_white gives white backgound
    pio.templates.default = "plotly_dark"
    
    fig= go.Figure()

    #Uses pandas and lamda function to parse necessary swift data columns through Data_Grapher function
    data= swift.apply(lambda row: Data_Grapher(row['SNname'], row['Dist_mod_cor'], row['type']), axis=1).to_list()
    #print(data)

    numList=[]

    '''
    These for-loops assign numbers to different supernovae traces depending on filter type.
    This is for the purpose of converting to bool values to use in plotly buttons that manage
    which traces are visable when you click on them.
    
    '''
    for i in range(len(data)):
        if data[i][0]=="NULL":
            continue
        else:
            fig.add_trace(data[i][0])
            numList.append("1")
    
    for i in range(len(data)):
        if data[i][1]=="NULL":
            continue
        else:
            fig.add_trace(data[i][1])
            numList.append("2")
    
    for i in range(len(data)):
        if data[i][2]=="NULL":
            continue
        else:
            fig.add_trace(data[i][2])
            numList.append("3")

    for i in range(len(data)):
        if data[i][3]=="NULL":
            continue
        else:
            fig.add_trace(data[i][3])
            numList.append("4")

    for i in range(len(data)):
        if data[i][4]=="NULL":
            continue
        else:
            fig.add_trace(data[i][4])
            numList.append("5")

    for i in range(len(data)):
        if data[i][5]=="NULL":
            continue
        else:
            fig.add_trace(data[i][5])
            numList.append("6")

    '''
    Converts the nuumbers in the numList to bool values and retruns
    list of bool values dependent on filter.
    '''
    w2, m2, w1, u, b, v= BtnBoolList(numList)
    
    #If no SNtype given set title to a specific title
    if SNtype == "":
        fig.update_layout(title= "All Supernovae Types Band Light Curves", 
                    xaxis_title="Days since first detection", yaxis_title="Absolute Magnitude", legend_title= "Supernovae")
    else:
        fig.update_layout(title= "All Supernovae " + str(SNtype)+ " Band Light Curves", 
                    xaxis_title="Days since first detection", yaxis_title="Absolute Magnitude", legend_title= "Supernovae" + " Type "+ str(SNtype)+"'s")

    #reverses y axis so that absolute magnitude is displayed properly
    fig['layout']['yaxis']['autorange'] = "reversed"

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

    '''
    Sets the drop down buttons for the different filters. 
    The filter bool lists are used here to determine which traces are visible or not depending on button clicked.

    Error buton is also defined here, it makes error bars visible or not depending on if it is clicked.
    '''
    fig.update_layout(
        updatemenus=[
            dict(
                type="dropdown", showactive=False, xanchor="left", yanchor="top", x=0.08, y=1.08,active=0, buttons=list(
                    [   
                        dict(
                            label= "All",
                            method= "update",
                            args= [{'visible': [True, True]}, {'title': "All Supernovae " + str(SNtype)+ " Band Light Curves"}]),
                        dict(
                            label= "UVW2",
                            method= "update",
                            args= [{'visible': w2}, {'title': "Supernovae " + str(SNtype)+ " UVW2 Band Light Curves"}]),
                        dict(
                            label= "UVM2",
                            method="update",
                            args=[{'visible': m2}, {'title': "Supernovae " + str(SNtype)+ " UVM2 Band Light Curves"}]),
                        dict(
                            label= "UVW1",
                            method="update",
                            args=[{'visible': w1}, {'title': "Supernovae " + str(SNtype)+ " UVW1 Band Light Curves"}]),
                        dict(
                            label= "U",
                            method="update",
                            args= [{'visible': u}, {'title': "Supernovae " + str(SNtype)+ " U Band Light Curves"}]),
                        dict(
                            label= "B",
                            method="update",
                            args=[{'visible': b}, {'title': "Supernovae " + str(SNtype)+ " B Band Light Curves"}]),
                        dict(
                            label= "V",
                            method="update",
                            args=[{'visible': v}, {'title': "Supernovae " + str(SNtype)+ " V Band Light Curves"}])
                    ]
                )
            ),
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
    #Adds the dropdown button text
    fig.update_layout(
    annotations=[
        dict(text="Filters:", font_size=15, height=32, showarrow=False,
        x=0, y=1.085, yref="paper", align="left")
    ]
)

    #Saves graph as html file with this name format
    fig.write_html("SN_type_"+SNtype+"_Light_Curve.html")

    fig.show()
    
    
if __name__ == "__main__":  main()