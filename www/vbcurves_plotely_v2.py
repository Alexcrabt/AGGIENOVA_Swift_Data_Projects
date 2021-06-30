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

    def File_Data(data):
        mjd_temp= []
        mag_temp= []
        magerr_temp= []

        search= 'UVM2'

        for sublist in data:
            if sublist[0] == search:
                mjd_temp.append(sublist[1])
                mag_temp.append(sublist[2])
                magerr_temp.append(sublist[3])

        if not len(mjd_temp) or not len(mag_temp):
            raise print("Not enough data on " + snname)
        else:
            return mjd_temp, mag_temp, magerr_temp
    
    def Null_Filter(mjd_temp, mag_temp, magerr_temp):
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
        
        if not len(mjd) or not len(mag):
            raise print("Not enough data on " + snname)
        else:
            return mjd, mag, magerr

    def Value_Converter(mjd, mag, magerr, dist_mod):
        days_since_detection= []
        absolute_mag= []
        final_magerr= []
        
        for j in range(len(mjd)):
            days_since_detection.append(float(mjd[j]) - float(mjd[0]))
            absolute_mag.append(float(mag[j]) - dist_mod)
            final_magerr.append(abs(float(magerr[j])))

        if np.isnan(days_since_detection).all() == True or np.isnan(absolute_mag).all() == True:
            raise print("Not enough data on " + snname)
        else:
            return days_since_detection, absolute_mag, final_magerr

    try:
        dataContent= [i.strip().split() for i in open("./data/"+snname+ "_uvotB15.1.dat")]
    except:
        print(str(snname) + ' is not in our data files')
        return "NULL"

    try:
        mjd, mag, magerr= File_Data(dataContent)

        mjd, mag, magerr= Null_Filter(mjd, mag, magerr)

        mjd, mag, magerr= Value_Converter(mjd, mag, magerr, dist_mod)

        if type.find(SNtype) == 0:
            return go.Scatter(x= mjd, y= mag, error_y= dict(array= magerr, type= 'data', visible= False), mode= "lines+markers", name= snname)
        else:
            return "NULL"
    except:
        
        return "NULL"

def main():
    global SNtype

    swift= pd.read_csv("NewSwiftSNweblist.csv")

    SNtype= input("What type supernovae do you want to plot? (Press ENTER if you want to graph all):")

    pio.templates.default = "plotly_dark"
    
    fig= go.Figure()
    
    
    data= swift.apply(lambda row: Data_Grapher(row['SNname'], row['Dist_mod_cor'], row['type']), axis=1).to_list()

    for i in range(len(data)):
        if data[i]=="NULL":
            continue
        else:
            fig.add_trace(data[i])
   

    if SNtype == "":
        fig.update_layout(title= "UVW2 Band Light Curves", 
                    xaxis_title="Days since first detection", yaxis_title="Absolute Magnitude", legend_title= "Supernovae")
    else:
        fig.update_layout(title= "UVM2 Band Light Curves", 
                    xaxis_title="Days since first detection", yaxis_title="Absolute Magnitude", legend_title= "Supernovae" + " Type "+ str(SNtype)+"'s")

    fig['layout']['yaxis']['autorange'] = "reversed"

    line_styles_names = ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
    line_styles = cycle(line_styles_names)

    symbols = [x for x in SymbolValidator().values[::-2] if not isinstance(x, int)]
    symbols= [i for i in symbols if not i.isdigit()]
    del symbols[8:12]

    symbols_names = list(set([i.replace("-thin", "") for i in symbols]))
    markers = cycle(symbols_names)

    for d in fig.data:
        d.line["dash"] = next(line_styles)
        d.marker.symbol = next(markers)
        d.marker.size = 8

    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons", direction= "right", showactive=False, buttons=list(
                    [
                        dict(
                            label="Error Bar",
                            method="update",
                            args=[{"error_y.visible": False}],
                            args2= [{"error_y.visible": True}]
                           )

                    ]
                )
            )
        ]
    )

    fig.write_html("SN_type"+SNtype+"_UVM2_LC.html")

    fig.show()
    
if __name__ == "__main__":  main()