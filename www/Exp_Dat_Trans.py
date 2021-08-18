'''
@author: Alexander Crabtree
            acrabtree15@gmail.com
'''
import pandas as pd
import numpy as np


def Grab_AED_In(SNname):
    idex= AEDswift.loc[AEDswift.isin([SNname]).any(axis=1)].index.tolist()
    if not idex:
        return np.nan
    else:
        return idex[0]
def empty(SNname):
    return np.nan

def GrabData(ind):
    sn= AEDswift['SNname'][ind]
    ed= AEDswift['ExplosionDate'][ind]
    r= AEDswift['R'][ind]
    ref= AEDswift['Reference'][ind]
    meth= AEDswift['Method'][ind]

    return sn, ed, r,ref, meth

def emptyData(ind):
    return np.nan, np.nan, np.nan, np.nan, np.nan

def PlaceData(sn,ed,r,ref,meth):
    idex= swift.loc[swift.isin([str(sn)]).any(axis=1)].index.tolist()
    if not idex:
        return 
    else: 
        try:
            swift['Explosion_Date'][idex[0]]= AEDdat['ExplosionDate'][idex[0]]
            swift['Explosion_Date'][idex[1]]= AEDdat['ExplosionDate'][idex[1]]

            swift['Reference'][idex[0]]= AEDdat['R'][idex[0]]
            swift['Reference'][idex[1]]= AEDdat['R'][idex[1]]

            swift['ReferenceLinks'][idex[0]]= AEDdat['Reference'][idex[0]]
            swift['ReferenceLinks'][idex[1]]= AEDdat['Reference'][idex[1]]

            swift['Method'][idex[0]]= AEDdat['Method'][idex[0]]
            swift['Method'][idex[1]]= AEDdat['Method'][idex[1]]
            return 
        except Exception:
            swift['Explosion_Date'][idex[0]]= ed
            swift['Reference'][idex[0]]= r
            swift['ReferenceLinks'][idex[0]]= ref
            swift['Method'][idex[0]]= meth
            return 

def emptyGrave(sn, ed, r,ref, meth):
    return 

def Grab_Swift_In(SNname):
    idex= swift.loc[swift.isin([SNname]).any(axis=1)].index.tolist()
    if not idex:
        return np.nan
    else:
        return idex[0]
def empty(SNname):
    return np.nan

def main():
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    global swift, AEDswift, AEDdat
    swift= pd.read_csv("NewSwiftSNweblist.csv")
    AEDswift= pd.read_csv("AED_SwiftSNweblist.csv")

    ind1= pd.Series(swift.apply(lambda row: Grab_AED_In(row['SNname']) if pd.notna(row['SNname']) else empty(row['SNname']), axis=1))
    AEDin= pd.DataFrame(ind1.tolist(), columns=['ind'], index= ind1.index)

    ind2= pd.Series(AEDin.apply(lambda row: GrabData(row['ind']) if pd.notna(row['ind']) else emptyData(row['ind']), axis=1))
    AEDdat= pd.DataFrame(ind2.tolist(), columns=['SNname', 'ExplosionDate', 'R', 'Reference', 'Method'], index= ind2.index)

    AEDdat.apply(lambda row: PlaceData(row['SNname'], row['ExplosionDate'],row['R'], row['Reference'], row['Method']) if pd.notna(row['SNname']) else emptyGrave(row['SNname'],  row['ExplosionDate'], row['R'], row['Reference'],  row['Method']), axis=1)


    ind3= pd.Series(AEDswift.apply(lambda row: Grab_Swift_In(row['SNname']) if pd.notna(row['SNname']) else empty(row['SNname']), axis=1))
    swiftin= pd.DataFrame(ind3.tolist(), columns=['ind'], index= ind3.index)

    swift=swift.fillna('')

    swift.to_csv('NewSwiftSNweblist.csv', index= False)


if __name__ == "__main__":  main()