import pandas as pd
import numpy as np

def Mag_Counter(snname):

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
        search= 'UVM2'

        for sublist in dataContent:
            if sublist[0] == search:
                mag.append(sublist[2])
        
        final_count= count(mag)

        return len(final_count)

    except:
        return np.nan

    return


def main():
    #file = "data/"+snname + '_uvotB15.1.dat'

    swift=pd.read_csv("NewSwiftSNweblist.csv")

    num= pd.Series(swift.apply(lambda row: Mag_Counter(row['SNname']), axis=1))
    
    num= pd.DataFrame(num.to_list(), columns=['UVM2_DatPoints'], index=num.index)

    swift= pd.concat([swift, num], axis=1)

    swift.to_csv('NewSwiftSNweblist.csv', index= False)



if __name__ == "__main__":  main()