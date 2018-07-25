
#@Author: Rohan Hundia 

import pandas as pd
import glob
import numpy as np
import math
import random

all_files  = glob.glob("*.csv")
df_list = []

for filepathname in sorted(all_files):
    print filepathname
    df_list.append(pd.read_csv(filepathname))

full_df = pd.concat(df_list)

#print(full_df)

sorteddata = full_df 

#sorteddata = sorted(full_df, key=lambda row: row[8], reverse=True) #specify row value according to which you want to sort data

subject = sorteddata['subjectId']
blockId = sorteddata['blockId']
trial = sorteddata['trialId']
mu1 = sorteddata['mu1']
mu2 = sorteddata['mu2']
choice = list(sorteddata['responseKey.keys'])
reward = (sorteddata['reward'])
RTs = list(sorteddata['responseKey.rt'])
cond = list(sorteddata['condition'])


runId = sorteddata['runId']
reward = sorteddata['reward']

blockId = list(blockId[~np.isnan(blockId)]) #done tbd changes
subject = subject[~pd.isnull(subject)] #done
trial = list(trial[~pd.isnull(trial)]) #done
mu1 =  list(mu1[~pd.isnull(mu1)]) #done
mu2 =  list(mu2[~pd.isnull(mu2)]) #done
reward =  list(reward[~pd.isnull(reward)]) #done
runId = runId[~np.isnan(runId)]
block = list(((runId - 1) * 4) + blockId)
runId = list(runId)
subject = list(subject)



print(block[0])

# RT same as choice 

# left block, RT, choice 


new_cond, new_choice, new_RTs  = [], [], []


for i in xrange(len(cond)):
    cond[i] = str(cond[i])
    if('nan' not in cond[i]):
        new_cond.append(cond[i])
        try:

            if(choice[i] == 'None'):
                #new_choice.append(-1)
                #print new_choice
                #new_RTs.append(-1)
                new_choice.append(1)
                new_RTs.append(0.56)
            else:
                new_choice.append(choice[i])
                new_RTs.append(RTs[i])
                
                
        except(KeyError):
            pass
            new_choice.append(1)
            new_RTs.append(0.56)
        

for i in xrange(len(new_cond)):
    
    if('RS' in new_cond[i]):
        new_cond[i] = 1
    elif('SR' in new_cond[i]):
        new_cond[i] = 2
    elif('RR' in new_cond[i]):
        new_cond[i] = 3
    elif('SS' in new_cond[i]):
        new_cond[i] = 4


new_RTs = [i * 1000 for i in new_RTs]
#print(new_choice)
#np.savetxt('final_data.csv', ('subject', 'block','trial','mu1','mu2','choice','reward', 'RT', 'cond'), delimiter = ',')



datafile = open('final_data.csv','a')
headers = 'subject, block, trial, mu1, mu2, choice, reward, RT, cond \n'
datafile.write(headers)

for i in xrange(len(subject)):
    #print(subject[i])
    numeric_subject = int(filter(str.isdigit, subject[i]))
    vals = str(numeric_subject)+ ","+ str(block[i]) + "," + str(trial[i]) + "," + str(mu1[i]) + "," + str(mu2[i]) + "," + str(new_choice[i])+ ","+ str(reward[i])+ "," + str(new_RTs[i]) + "," + str(new_cond[i]) + "\n"
    datafile.write(vals)
datafile.close()
   


