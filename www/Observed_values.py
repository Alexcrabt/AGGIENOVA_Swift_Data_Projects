'''
@author: Alexander Crabtree
            acrabtree15@gmail.com
'''
import pandas as pd
import numpy as np
def W1_Counter(snname):

    def count(mag):
        magcount=[]
        for i in range(len(mag)):
            if mag[i] == "NULL":
                continue
            else:
                magcount.append(mag[i])
        return magcount

    try:
        mag= []
        dataContent= [i.strip().split() for i in open("./data/"+snname+ "_uvotB15.1.dat")]
        search= 'UVW1'

        for sublist in dataContent:
            if sublist[0] == search:
                mag.append(sublist[2])
        
        final_count= count(mag)

        return len(final_count)

    except:
        return np.nan

    return

def dataPoint_Check(sortList):

    for i in range(len(sortList)):
        if sortList['W1_DatPoints'][i] == int(3) and sortList['W1_DatPoints'][i+1] == int(2):
                sortList= sortList.iloc[:i+1]
                return sortList

        else:
            continue

def main():

    swift= pd.read_csv("NewSwiftSNweblist.csv")

    #print(swift.head(40))

    w1_num= pd.Series(swift.apply(lambda row: W1_Counter(row['SNname']), axis=1))

    ranList= pd.DataFrame({"W1_DatPoints": w1_num, "SNname": swift['SNname'],"SNtype": swift['type'], "Expl_Date":swift['Explosion_Date']})
    #print(ranList.head(40))
    ranList= ranList.iloc[2:]

    ranList.reset_index(drop=True, inplace=True)#Resets the index after dropping first row
    #print(ranList.head(40))

    pd.set_option('display.max_rows', 1000)
    sortList=ranList.sort_values(by="W1_DatPoints", ascending= False, na_position='last')#Sorts the ranList by the DatPoints in decending order
    sortList.reset_index(drop=True, inplace=True)#Resets index after sorting
    
    sortList= dataPoint_Check(sortList)
    print(sortList)



if __name__ == "__main__":  main()