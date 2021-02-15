# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:14:39 2020

@author: Alexander Crabtree
""" 
import requests
import json
import pandas as pd #for dataframes
import numpy as np
import math
import sys
import itertools
import threading
import time
import astropy.units as u
from astroquery.irsa_dust import IrsaDust
from astropy.coordinates import Angle,SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from bs4 import BeautifulSoup


'''A fun function Tate had in his code'''
def animate():
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done:
            break
        sys.stdout.write('\rloading ' + c)
        sys.stdout.flush()
        time.sleep(0.1)
    print(chr(27) + "[2J")
    sys.stdout.write('\rDone!')


def getLink(name):
    
    try:
        link = "http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?"
        inputs = {'objname': name,
				'extend': 'no',
				'hconst': '73',
				'omegam': '0.27',
				'omegav': '0.73',
				'corr_z': '1',
				'out_csys': 'Equatorial',
				'out_equinox': 'J2000.0',
				'obj_sort': "RA or Longitude",
				'of': 'pre_text',
				'zv_breaker': '30000.0',
				'list_limit':'10',
				'img_stamp': 'YES'}
        page = requests.get(link, params = inputs)
        soup = BeautifulSoup(page.content, 'html.parser')
        return soup
    except requests.ConnectionError:
        soup = ''
        return soup

def scrapeValues(soup):
	#----------Get Velocities----------#
    try:
        soup = soup.find("a", attrs={'name':'DerivedValues_0'})
        velocities = soup.next_sibling.next_sibling.next_sibling.find("pre")
        Helio = list(velocities.children)[2]
        VGS = list(velocities.children)[16]
        Helio = Helio.lstrip('\n')
        VGS = VGS.lstrip('\n')
        Hvals = [int(s) for s in Helio.split() if s.isdigit()]
        VGSVals = [int(s) for s in VGS.split() if s.isdigit()]
        vels = [Hvals[0],Hvals[1],VGSVals[0],VGSVals[1]]
        return vels
    except: 
        vels = [np.nan,np.nan,np.nan,np.nan]
        return vels

def get_coords(gals):
    """
        Takes list of galaxies and looks up their coordinates by name.
        If no name found: warn, skip, remove galaxy from list

        Returns:
            gals: list of galaxies minus those that weren't found
            start_coord: list of coordinates corresponding to center of galaxies in 'gals'
    """
    start_coord = []
    #bar = FillingCirclesBar('Loading galaxies', max = len(gals))
    try:
        tempCoord = SkyCoord.from_name(str(gals), frame = 'icrs')
        start_coord= tempCoord
            #bar.next()
    except NameResolveError:
        #if gals != str(''):
            #print('\nSkipping',gals,'because it couldn\'t be found.')
            #start_coord= ''
        #else:
            start_coord= np.nan
    #bar.finish()
    return(start_coord)

def coord_breakup(coord):
    """
        Breaks up a coordinate into its components
    """
    if pd.isna(coord) == True:
        ra = np.nan
        dec = np.nan
        return ra, dec
    else:
        ra = Angle(coord.ra.hour,unit = u.hour)
        dec = coord.dec
        return ra, dec


def Redshift(soup):

    try:
        soup = soup.find("a", attrs={'name':'BasicData_0'})
        soup = soup.next_sibling.next_sibling.next_sibling.find("pre")
        redshift = list(soup.children)[0]
        # print(redshift)
        redshift = str(redshift).split("Redshift",1)[1]
        # print(redshift)
        redshift = redshift.split(":",1)[1]
        redshift = redshift.split(" +",1)[0].strip(' ')
        return(redshift)
    except:
        redshift = np.nan
        return redshift
   

def Morphology(soup):
    
    try:
        soup = soup.find("a", attrs={'name':'BasicData_0'})
        morphology = soup.next_sibling.next_sibling.next_sibling.find("pre")
        morphology = list(morphology.children)[2]
        morphology = (morphology.split(": ",4)[4]).rstrip()
        return morphology
    except:
        morphology = np.nan
        return morphology


def LatLong(soup):
    
    try:
        soup = soup.find("a", attrs={'name':'Positions_0'})
        coords = soup.next_sibling.next_sibling.find("pre")
        coords = list(coords.children)[4]
        coords = (coords.split("Galactic ",1)[1]).lstrip()
        long, lat = str(coords.split()[0]), str(coords.split()[1])
 
        return long, lat
    except:
        long, lat= np.nan,np.nan
        return long, lat  

        
def getAVbest(inputcoordinates):
    '''
    Coordinates are input as a single string. Output is the recommended Av value for MW reddening, error, and reference
    '''
    
    #inputcoordinates = sys.argv[1]
    testCoords = SkyCoord(inputcoordinates,frame='fk5')

    #print('\n-----\nReading input files...')
    inFile = 'Brown_Walker_table_1.dat'
    inTable = pd.read_csv(inFile,header=None,delimiter=' ')
    ra = Angle(inTable.iloc[:,1])
    dec = Angle(inTable.iloc[:,2])
    sourceCoords = SkyCoord(ra,dec,frame='fk5')
    
    #print('Calculating separation from table coordinates')
    separations = testCoords.separation(sourceCoords).arcminute
    # compare to the distances in the table
    within = np.less(separations,inTable.iloc[:,3])
    
    # Are any of the input coordinates within the tabulated distance 
    # of the coordinates in the table?
    correctedAV = np.where(within,inTable.iloc[:,4],None) #get calculated value
    fix=any(within)
    #print('fix?',fix)
    
    if fix:
        AV = next((item for item in correctedAV if item is not None),None)
        correctedAVerr = np.where(within,inTable.iloc[:,5],None) #get calculated val
        newAVerr = next((item for item in correctedAVerr if item is not None),None)
        AVerr = math.sqrt((newAVerr)**2+(AV*0.1)**2)
        sources=np.where(within,inTable.iloc[:,6],None)
        source = next((item for item in sources if item is not None),None)+",S_F_2011"
    if not fix:
        AVtable = IrsaDust.get_extinction_table(testCoords,show_progress = False)
        AV=AVtable['A_SandF'][2]
        AVerr = AV*0.1
        source = 'S_F_2011'

    #print(AV, AVerr, source)
    return (AV, AVerr, source);

def SNdates(SNname):
    '''
        Takes Supernovae names from CSV and looks them up in the URL below
        Then parses the page for relevant information to Supernova Discovery,
        Last non-detection, and recieved date.
    '''
    url = "https://www.wis-tns.org/object/"+SNname+"/discovery-cert"

    page = requests.get(url)
    soup = BeautifulSoup(page.content, 'html.parser')
    
    sn_data= soup.find_all('span', class_='underline')
    
    rec_date= soup.find_all('em', class_='placeholder')[3].text
    
    disc_date= sn_data[0].find_next_sibling('span').text
    disc_string= 'Discovery date: '

    lnd_date= sn_data[1].find_next_sibling('span').text
    lnd_string= 'Last non-detection date: '

    if disc_string in disc_date:
        disc_date= disc_date.replace(disc_string, '')
    
    if lnd_date in ['Archival info: DSS', 'Archival info: Other', 'Archival info: SDSS']:
        lnd_date= ''
        rec_date= "'" + rec_date
        disc_date= "'" + disc_date
        
    else:
        if lnd_string in lnd_date:
            lnd_date= lnd_date.replace(lnd_string, '')
        rec_date= "'" + rec_date
        disc_date= "'" + disc_date
        lnd_date= "'" + lnd_date

    return(rec_date, disc_date, lnd_date)

# create a formatted string of the Python JSON object-------------------------
def jprint(obj):
    '''
        Used to interpret the data from astrocats catalog
    '''
    
    text = json.dumps(obj, sort_keys=True, indent=4)
    print(text)
 

'''
    Next 4 lines are important to animate function
'''
done = False
print(chr(27) + "[2J")
threader = threading.Thread(target=animate)
threader.start()


'''
    Next 2 lines read in the supernovae csv and replaces anon, Anon, AnonHost, and empty values with nan values
'''
swift= pd.read_csv('NewSwiftSNweblist.csv')
swift= swift.replace({r'anon',r'Anon',r'AnonHost', r'^\s*$'},np.nan, regex=True)


'''
    This if-else statment asks if you want to update or add new supernovae values to CSV file
'''
choice_1= input("Would you like to add new supernovae names or just update current CSV? (SN or U):")
if choice_1 == str('U'):
    print("Updating CSV now!")
else:
     choice_2=""
     while choice_2 != str('Done'):
        SN=""
        SN= input("Type supernova name:")
        new_SN= pd.DataFrame([str(SN)], columns= ['SNname'])
        swift= pd.concat([new_SN, swift]).reset_index(drop = True) 
        choice_2= input("Would you like to add another or are you done? (Y or Done):")
     print("Updating CSV now!")
      

'''
    Requests data ra, dec, and host data of all supernovae in astrocats catalog
'''
sn_data= requests.get("https://api.astrocats.space/catalog/ra+dec+host?first")
sn_catalog= sn_data.json()


def SNHost(SNname):
    ''' 
        Checks if any CSV Supernova names are in astrocats catalog
        If yes returns Host Name data, If no returns nan value
    '''
    
    if SNname in sn_catalog:
        try:
            host_names= sn_catalog[SNname]['host']['value']
        except Exception:
            host_names= np.nan
    else:
        host_names= np.nan
    return host_names

def SNRa(SNname):
    ''' 
        Checks if any CSV Supernova names are in astrocats catalog
        If yes returns Ra data, If no returns nan value
    '''
    
    if SNname in sn_catalog:
        try:
            snra= "'%s" % sn_catalog[SNname]['ra']['value']
        except Exception: 
            snra= np.nan
    else:
        snra= np.nan
    return snra

def SNDec(SNname):
    ''' 
        Checks if any CSV Supernova names are in astrocats catalog
        If yes returns Dec data, If no returns nan value
    '''
    
    if SNname in sn_catalog:
        try:
            sndec= sn_catalog[SNname]['dec']['value']
        except Exception: #if no dec value input what is given for dec 
            sndec= np.nan
    else:
        sndec= np.nan
    return sndec

'''
    The next 3 lines uses pandas apply function and lambda function to quickly
    parse supernovae names into the 3 above functions and saves data to
    relavent columns in CSV
'''
swift["HostName"]= swift.apply(lambda row: SNHost(SNname=row['SNname']), axis=1)
swift["SNra"]= swift.apply(lambda row: SNRa(SNname=row['SNname']), axis=1)
swift["SNdec"]= swift.apply(lambda row: SNDec(SNname=row['SNname']), axis=1)



def AV_empty(Ra, Dec):
    ''' 
        Returns empty nan values if condition is met
    '''
    AV, AVerr, AVsour= np.nan, np.nan, np.nan
    return pd.Series({'AV': AV, 'AVerr': AVerr, 'AVsour': AVsour})
def GrabAVbest(Ra,Dec):
    '''
        Takes Supernovae Ra,Dec data and puts in in the right format before
        runing it through the getAVbest function, if there is an Exception
        returns nan values
    '''
    ra_fix, dec_fix = Ra, Dec
    ra_fix= ra_fix.replace("'", "")
    combined= SkyCoord(ra=ra_fix, dec=dec_fix, unit=(u.hour, u.deg)).to_string('hmsdms')
    try:
        AV, AVerr, AVsour= getAVbest(combined)
    except Exception:
        AV, AVerr, AVsour= np.nan, np.nan, np.nan
    
    return pd.Series({'AV': AV, 'AVerr': AVerr, 'AVsour': AVsour})

'''
    Pandas apply function and lambda function with an embedded if else stament,
    if AV cell is empty and corresponding SNra cell is not runs SNra, SNdec data into GrabAVbest
    else runs it through AV_empty
'''
swiftAV= pd.DataFrame(swift.apply(lambda row: GrabAVbest(row['SNra'], row['SNdec']) if pd.isna(row['AV']) and pd.notna(row['SNra']) else AV_empty(row['SNra'], row['SNdec']), axis= 1))
swift= swift.fillna(swiftAV) #Fills new data into corresponding nan space in CSV


#gathers data about host galaxy ---------------------------------------------

def AllHostData():
    def host_coords(HostName):
        ra,dec= coord_breakup(get_coords(HostName))
        return ra,dec
    def host_coords_empty(HostName):
        ra,dec= np.nan, np.nan
        return ra, dec
    
    def host_lat_long(HostName):
        gal_link=getLink(HostName)
        long, lat= LatLong(gal_link)
        return(long, lat)
    def host_lat_long_empty(HostName):
        long, lat=np.nan, np.nan
        return long, lat
    
    def redshift(HostName):
        gal_link=getLink(HostName)
        red= Redshift(gal_link)
        return red
    def redshift_empty(HostName):
        red= np.nan
        return red
    
    def morphology(HostName):
        gal_link=getLink(HostName)
        morph= Morphology(gal_link)
        return morph
    def morphology_empty(HostName):
        morph= np.nan
        return morph
    
    def velocities(HostName):
        gal_link=getLink(HostName)
        vels= scrapeValues(gal_link)
        return(vels)
    def velocities_empty(HostName):
        vels= np.nan, np.nan, np.nan, np.nan
        return(vels)
    
    coord= swift.apply(lambda row: host_coords(row['HostName']) if pd.notna(row['HostName']) and pd.isna(row['HostRa']) else host_coords_empty(row['HostName']), axis=1)
    
    lon_lat= swift.apply(lambda row: host_lat_long(row['HostName']) if pd.notna(row['HostName']) and pd.isna(row['Long']) else host_lat_long_empty(row['HostName']), axis=1)
    
    reds= swift.apply(lambda row: redshift(row['HostName']) if pd.notna(row['HostName']) and pd.isna(row['Redshift']) else redshift_empty(row['HostName']), axis=1)
    
    morp= swift.apply(lambda row: morphology(row['HostName']) if pd.notna(row['HostName']) and pd.isna(row['Morphology']) else morphology_empty(row['HostName']), axis=1)
    
    hvels= swift.apply(lambda row: velocities(row['HostName']) if pd.notna(row['HostName']) and pd.isna(row['host_velocity']) else velocities_empty(row['HostName']), axis=1)
    
    
    
    return pd.Series({'HostRa':coord[0], 'HostDec':coord[1], 'Long':lon_lat[0], 'Lat':lon_lat[1], \
                      'Redshift':reds, 'Morphology':morp, 'host_velocity':hvels[0],\
                      'host_vel_err':hvels[1], 'host_vel_corr':hvels[2], 'host_vel_corr_err':hvels[3]})

'''
    Turns returned series from AllHostData into DataFrame then fills in data
    to corresponding nan spaces in CSV
'''
swiftHostData= pd.DataFrame(AllHostData())
swift= swift.fillna(swiftHostData)
    
def Dist_mod_empty(hv, hverr):
    ''' 
        Returns empty nan series if condition is met
    '''
    distance_mod_cor= np.nan
    distance_mod_cor_err= np.nan
    return pd.Series({'Dist_mod_cor': distance_mod_cor, 'Dist_mod_cor_err': distance_mod_cor_err})
def Distance_mod_cor(hv, hv_err):
    '''
        Takes host_velocity and host_vel_err data and returns distance modulous
        and distance modulous err in a series
    '''
    h0 = 72.0
    h0err = 5.0

    try:
        distance_mod_cor = 5*math.log(float(hv)/h0,10)+25 #hubble flow
        distance_mod_cor_err = math.sqrt(((5*float(hv_err))/(float(hv)*math.log(10,10)))**2+((5*200)/(float(hv)*math.log(10,10)))**2 + ((5*5.0)/(h0*math.log(10,10)))**2)
    except Exception:
        distance_mod_cor= np.nan
        distance_mod_cor_err= np.nan
    
    return pd.Series({'Dist_mod_cor': distance_mod_cor, 'Dist_mod_cor_err': distance_mod_cor_err})

'''
    Pandas apply function and lambda function with embeded if else statment. 
    If host_velocity cell is not nan then runs host_velocity, host_vel_err data
    through Distance_mod_cor, else runs it through Dist_mod_empty
'''
swiftDist_mod= pd.DataFrame(swift.apply(lambda row: Distance_mod_cor(row['host_velocity'], row['host_vel_err']) if pd.notna(row['host_velocity']) else Dist_mod_empty(row['host_velocity'], row['host_vel_err']), axis=1))
swift= swift.fillna(swiftDist_mod)

def GrabSNdates_empty(SNname):
    ''' 
        Returns empty nan series if condition is met
    '''
    rec,disc,lnd= np.nan,np.nan,np.nan
    return(pd.Series({'Date_Recieved': rec, 'Discover_date': disc, 'Last_non-detection_date': lnd}))
def GrabSNdates(SNname):
    '''
        Takes supernovae names and runs them through SNdates and returns a
        series with dates data
    '''
    try:
        rec,disc,lnd= SNdates(SNname)
    except Exception:
        rec,disc,lnd= np.nan,np.nan,np.nan
        
    return(pd.Series({'Date_Recieved': rec, 'Discover_date': disc, 'Last_non-detection_date': lnd}))

'''
    Next 3 lines creates a supernovae list from our CSV except the SN in front of 
    most of the names are removed
'''
swift_temp= pd.DataFrame({'SNname':swift['SNname'], 'Discover':swift['Discover_date'], 'Last':swift['Last_non-detection_date']})
rep= '|'.join(['SN'])
swift_temp['SNname']= swift_temp['SNname'].str.replace(rep, '')

'''
    Pandas apply function and lambda function with if else statment, If 
    Discover and Last cell is nan runs SNname through GrabSNdates, else runs it
    through GrabSNdates_empty. Then fills new data into corresponding nan spots
    in CSV
'''
swiftDates= pd.DataFrame(swift_temp.apply(lambda row: GrabSNdates(row['SNname']) if pd.notna(row['SNname']) and (pd.isna(row['Discover']) or pd.isna(row['Last'])) else GrabSNdates_empty(row['SNname']), axis=1))
swift= swift.fillna(swiftDates)

'''
    Fills nan values back with blanks then saves csv to either same file or 
    your pick of file name
'''
swift=swift.fillna('')
swift.to_csv('TestSwiftSNweblist.csv', index= False)

'''Last segment of animate function that makes the animation run'''
done= True
