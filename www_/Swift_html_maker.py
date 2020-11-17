# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 17:17:46 2020

@author: Alexander Crabtree
"""
import numpy as np
import pandas as pd
import os
from os import path
import webbrowser



swift= pd.read_csv('NewSwiftSNweblist.csv')
swift= swift.drop(swift.index[0]) #delete header row
swift.dropna(subset=["SNname"], inplace= True) #drops all nan in SNname column
swift= swift.reset_index(drop=True) #reset the index from 0
swift = swift.replace(np.nan, '', regex=True) #turn all NANs into '' so it looks nice on webpage

#define all of the columns needed from Swift CSV
snname= swift.SNname 
host= swift.HostName 
othersn= swift.OtherName
red= swift.Redshift
#rounds the redshift numbers to the 6th decimal place
red= red.map(lambda x: round(x, 6) if isinstance(x, (int, float)) else x)
sntype= swift.type
psn= swift.PSN 

def printsnline(SNnamen, PSNn, SNname2n, galaxyn, redshiftn, typen, n): #creates new HTML code for swift_sn_middle.txt
#-----------------Paths for the differnt files, you might have to change it a bit for your computer depending on where the files are--------------------------
    datafilename= os.path.abspath("./AlexandersRepository/data/"+SNnamen+"_uvotB15.1.dat")
    
    imagefilename= os.path.abspath("./AlexandersRepository/images/"+SNnamen+"_uvot.png")

    lcfilename= os.path.abspath("./AlexandersRepository/lcplots/"+SNnamen+"_lightcurve.jpg")
    
#-----------------------------------------------------HTML text---------------------------------------------------------------
    texthead= """<tr> \n
    <td align="center" valign="top"> {n} <br></td> \n
             <td align="center" valign="top"><b> {SNnamen} </b><br> &nbsp; \n
             {PSNn} {SNname2n} \n
             <br> </td> \n
               <td align="center" valign="top"><a \n
                  href="http://ned.ipac.caltech.edu/forms/byname.html"> {galaxyn} </a><br> </td> \n
                     <td align="center" valign="top"> {redshiftn} <br> </td> \n
                        <td align="center" valign="top"><a \n
                        href="http://www.cbat.eps.harvard.edu/cbet/RecentCBETs.html"> {typen} </a><br> </td> \n
                         <td align="center" valign="top"><a \n""".format(**locals())

    #imagefile html section:
    if os.path.isfile(imagefilename)== True:
        textimage= """href="https://github.com/pbrown801/SOUSA/blob/master/images/{SNnamen}_uvot.png"><img \n
                         src="pic.jpg" alt="Swift image" style="border: 0px solid; border: 0px solid; width: 52px; height: 55px;" \n
                          border="0" height="55" width="52"></a><a \n""".format(**locals())
        texthead+= textimage

    #lcfilename html section:
    if os.path.isfile(lcfilename) == True:
        textlc= """href="https://github.com/pbrown801/SOUSA/blob/master/lcplots/{SNnamen}_lightcurve.jpg"><img \n
                    src="light.jpg" alt="Swift UVOT lightcurve" \n
                     style="border: 0px solid; border: 0px solid; width: \n
                      52px; height: 55px;" border="0" height="55" width="52"></a><a \n""".format(**locals())
        texthead+= textlc
                    
    #datafilenmesousa html section:
    if os.path.isfile(datafilename) == True:
        textdata= """href="https://github.com/pbrown801/SOUSA/blob/master/data/{SNnamen}_uvotB15.1.dat"><img \n
                      src="sousa_galaxy_small.png" alt="Swift UVOT data" \n
                       style="border: 0px solid; border: 0px solid; width: \n
                        52px; height: 89px;" border="0" height="89" width="52"></a></td> \n""".format(**locals())
        texthead+= textdata

                   
    textend= """<td align="center" valign="top"><br> \n
                 </td>\n
                  </tr>"""
    
    text= texthead + textend
      
    return(text)


swiftmid= open('new_swift_sn_middle.txt', 'w')
for i in range(0,len(snname)):
    if snname[i].startswith(('in', '?', '#'),0,4) != True:
        swiftmid.write(printsnline(str(snname[i]),str(psn[i]),str(othersn[i]),str(host[i]),str(red[i]),str(sntype[i]), i+1))
swiftmid.close

#Opens and reads the HTML code in the .txt files
import codecs
top= codecs.open("swift_sn_top.txt",'r')
mid= codecs.open("new_swift_sn_middle.txt", 'r')
bott= codecs.open("swift_sn_bottom.txt",'r')



html_t= top.read()
html_m= mid.read()
html_b= bott.read()
#defines a path for the website HTML combined code
path = os.path.abspath('SwiftSNWebsite.html')
url = 'file://' + path

#Runs all opened HTML .txt files in a browser
with open(path, 'w') as website:
    website.write(html_t)
    website.write(html_m)
    website.write(html_b)
webbrowser.open(url)


