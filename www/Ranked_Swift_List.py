'''
@author: Alexander Crabtree
            acrabtree15@gmail.com
'''
import pandas as pd
import numpy as np

''' 
    Since the list is saved as it's own csv I am commenting out this function.
    If you wanted a list in a single column of a csv this would be the way to go.

def Ranked_List(ranList):
    def Combine_list(obs, SN, type, Disc_date, Expl_date):
        comb= str(obs) + ", " + str(SN) + ", " + str(type) + ", " + str(Disc_date) + ", " + str(Expl_date)

        return comb

    sortList=ranList.sort_values(by="DatPoints", ascending= False, na_position='last')
    #print(sortList.head(20))
    combList= pd.Series(sortList.apply(lambda row: Combine_list(row['DatPoints'], row['SNname'], row['SNtype'], row['Disc'], row['Expl']), axis=1))
    combList= pd.DataFrame(combList.tolist(), columns=['Ranked_List'], index= combList.index)

    return combList
'''

def main():
    
    swift= pd.read_csv("NewSwiftSNweblist.csv")#reads in swift csv

    '''
    Saves relavent swift columns to a new dataframe
    '''
    ranList= pd.DataFrame({"DatPoints": swift['UVM2_DatPoints'], "SNname": swift['SNname'],"SNtype": swift['type'], "Disc_Date": swift['Discover_date'], "Expl_Date":swift['Explosion_Date']})
    ranList.drop(ranList.head(1).index, inplace=True)#Drops first row since it didn't have any info
    ranList.reset_index(drop=True, inplace=True)#Resets the index after dropping first row


    sortList=ranList.sort_values(by="DatPoints", ascending= False, na_position='last')#Sorts the ranList by the DatPoints in decending order
    sortList.reset_index(drop=True, inplace=True)#Resets index after sorting
    #print(sortList.head(20))

    text=['In descending order of number of UVM2_DatPoints']#Info on how list is sorted
    text= pd.DataFrame({"DatPoints":text})#Making it a dataframe allows for easy merging with sortList
    #sortList= pd.DataFrame(np.insert(sortList.values, 1, 'Observations, SNname, type, Disc_date, Expl_date', axis=0))
    #print(sortList.head(20))

    sortList= pd.concat([text, sortList])

    sortList.to_csv('SwiftRankedList.csv', index=False, na_rep='nan')#Saves sortList and has nan values print as nan

if __name__ == "__main__":  main()