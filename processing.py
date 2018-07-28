
#@Author: Rohan Hundia 

import pandas as pd
import glob
import numpy as np
import math
import random

all_files  = glob.glob("data/*.csv")
df_list, file_paths = [],[]

for filepathname in sorted(all_files):
    file_paths.append(filepathname)

file_paths.sort()
print(file_paths)

for i in xrange(len(file_paths)):
    df_list.append(pd.read_csv(file_paths[i]))
    

#print(df_list)
full_df = pd.concat(df_list)

#print(full_df)

sorteddata = full_df 

#sorteddata = sorted(full_df, key=lambda row: row[8], reverse=True) #specify row value according to which you want to sort data

subject = sorteddata['subjectId']
blockId = sorteddata['blockId']
trial = sorteddata['trialId']
mu1 = sorteddata['mu1']
mu2 = sorteddata['mu2']
r1 = sorteddata['leftAnswer']
r2 = sorteddata['rightAnswer']
choice = sorteddata['responseKey.keys']
reward = sorteddata['reward']
RTs = sorteddata['responseKey.rt']
cond = sorteddata['condition']
trial_onset = sorteddata['actualChoiceOnset']
choice_onset = sorteddata['actualChoiceOffset']
feedback_onset = sorteddata['actualFeedbackOnset']


runId = sorteddata['runId']
reward = sorteddata['reward']

which_rows = ~np.isnan(blockId)

subject = list(subject[which_rows])
runId = runId[which_rows]
run = list(runId)
blockId = list(blockId[which_rows])
block = list(((runId - 1) * 4) + blockId)
trial = list(trial[which_rows])
mu1 =  list(mu1[which_rows])
mu2 =  list(mu2[which_rows])
r1 =  list(r1[which_rows])
r2 =  list(r2[which_rows])
choice =  list(choice[which_rows])
reward =  list(reward[which_rows])
RTs = list(RTs[which_rows])
cond = list(cond[which_rows])
trial_onset = list(trial_onset[which_rows]) 
choice_onset = list(choice_onset[which_rows]) 
feedback_onset = list(feedback_onset[which_rows]) 



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
                new_choice.append(1) # MOMCHIL TODO: what's this?
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



datafile = open('data.csv','w') # we process all subjects each time
headers = 'subject, run, block, trial, mu1, mu2, r1, r2, choice, reward, RT, cond, trial_onset, choice_onset, feedback_onset \n'
datafile.write(headers)

for i in xrange(len(subject)):
    #print(subject[i])
    numeric_subject = int(filter(str.isdigit, subject[i]))
    vals = str(numeric_subject)+ "," + str(run[i]) + "," + str(block[i]) + "," + str(trial[i]) + "," + str(mu1[i]) + "," + str(mu2[i]) + "," + str(r1[i]) + "," + str(r2[i]) + ","  + str(new_choice[i])+ ","+ str(reward[i])+ "," + str(new_RTs[i]) + "," + str(new_cond[i]) + "," + str(trial_onset[i]) + "," + str(choice_onset[i]) + "," + str(feedback_onset[i]) + "\n"
    datafile.write(vals)
datafile.close()
   


