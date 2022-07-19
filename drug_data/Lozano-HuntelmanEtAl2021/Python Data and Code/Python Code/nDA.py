# [1] --------- SETUP ----------
# Import libraries
import pandas as pd #use "pd" to call for pandas
import itertools # for combinations
import os

# [2] --- 3 combo ---
storage = {}
fileList = os.listdir(r"/u/.../HD/C3") # for csvs dealing with potential hidden suppression

x = 0
for file in fileList:
	fl_store = r"/u/.../HD/C3/%s"%(file)
	storage[x] = pd.read_csv(fl_store) # csvs now in storage dict
	x = x+1

COMB3 = "nC3_" # naming purposes
number3 = 0

for num in range (0,1428): ##### CHANGE VALUE DEPENDING ON AMOUNT OF COMBOS IN QUESTION #####
	look = storage[num] # looping through dictionary
	lowFit = look.iloc[0]["F2"] # temporary value
	
	for index, rows in look.iterrows(): # looking for smallest fitness at lower orders
		nf2 = rows["F2"]
		if nf2 < lowFit:
			lowFit = nf2

	EDIT3 = pd.DataFrame(columns = ['C3','F3','C2','F2','C1','F1','newDA','INTR'])   
	for indx, rws in look.iterrows(): # adding adjusted DA row
		c3 = rws["C3"]
		f3 = rws["F3"] #
		c2 = rws["C2"]
		f2 = rws["F2"]
		c1 = rws["C1"]
		f1 = rws["F1"] # don't need to check bc already gone over via DA3 >= 1.3 check
		i3 = rws["INTR"]
		if lowFit == 0:
			lowFit = 0.000000000001 # so no division by 0
		nd = f3/lowFit

		EDIT3 = EDIT3.append({'C3':c3,'F3':f3,'C2':c2,'F2':f2,'C1':c1,'F1':f1,'newDA':nd,'INTR':i3}, ignore_index=True)

	f_name = "/u/.../HD/nC3/%s%s.csv"%(COMB3,str(number3))
	EDIT3.to_csv(f_name, index = False)
	
	number3 = number3 + 1

# [3] --- 4 combo ---
# structure same as above
storage4 = {}
fileList4 = os.listdir(r"/u/.../HD/C4")

x4 = 0
for file in fileList4:
	fl_store4 = r"/u/.../HD/C4/%s"%(file)
	storage4[x4] = pd.read_csv(fl_store4)
	x4 = x4+1
	
COMB4 = "nC4_"
number4 = 0

for num in range (0,4849): ##### CHANGE VALUE DEPENDING ON AMOUNT OF COMBOS IN QUESTION #####
	look4 = storage4[num] # looping through dictionary
	lowFit4 = look4.iloc[0]["F2"] # temporary value
	
	for index, rows in look4.iterrows(): # looking for smallest fitness at lower orders
		nf3 = rows["F3"]
		nf2 = rows["F2"]
		if nf3 < lowFit4:
			lowFit4 = nf3
		if nf2 < lowFit4:
			lowFit4 = nf2
	 
	EDIT4 = pd.DataFrame(columns = ['C4','F4','C3','F3','C2','F2','C1','F1','newDA','INTR'])   
	for indx, rws in look4.iterrows(): # adding adjusted DA row
		c4 = rws["C4"]
		f4 = rws["F4"] #
		c3 = rws["C3"]
		f3 = rws["F3"] 
		c2 = rws["C2"]
		f2 = rws["F2"]
		c1 = rws["C1"]
		f1 = rws["F1"] # don't need to check bc already gone over via DA3 >= 1.3 check
		i4 = rws["INTR"]
		if lowFit4 == 0:
			lowFit4 = 0.000000000001
		nd = f4/lowFit4

		EDIT4 = EDIT4.append({'C4':c4,'F4':f4,'C3':c3,'F3':f3,'C2':c2,'F2':f2,'C1':c1,'F1':f1,'newDA':nd,'INTR':i4}, ignore_index=True)

	f_name = "/u/.../HD/nC4/%s%s.csv"%(COMB4,str(number4))
	EDIT4.to_csv(f_name, index = False)
	
	number4 = number4 + 1

# [4] --- 5 combo ---
storage5 = {}
fileList5 = os.listdir(r"/u/.../HD/C5")

x5 = 0
for file in fileList5:
	fl_store5 = r"/u/.../HD/C5/%s"%(file)
	storage5[x5] = pd.read_csv(fl_store5)
	x5 = x5+1
	
COMB5 = "nC5_"
number5 = 0

for num in range (0,10972): ##### CHANGE VALUE DEPENDING ON AMOUNT OF COMBOS IN QUESTION #####
	look5 = storage5[num] # looping through dictionary
	lowFit5 = look5.iloc[0]["F2"] # temporary value
	
	for index, rows in look5.iterrows(): # looking for smallest fitness at lower orders
		nf4 = rows["F4"]
		nf3 = rows["F3"]
		nf2 = rows["F2"]
		if nf4 < lowFit5:
			lowFit5 = nf4
		if nf3 < lowFit5:
			lowFit5 = nf3
		if nf2 < lowFit5:
			lowFit5 = nf2

	EDIT5 = pd.DataFrame(columns = ['C5','F5','C4','F4','C3','F3','C2','F2','C1','F1','newDA','INTR'])   
	for indx, rws in look5.iterrows(): # adding adjusted DA row
		c5 = rws["C5"]
		f5 = rws["F5"] #
		c4 = rws["C4"]
		f4 = rws["F4"] 
		c3 = rws["C3"]
		f3 = rws["F3"] 
		c2 = rws["C2"]
		f2 = rws["F2"]
		c1 = rws["C1"]
		f1 = rws["F1"] # don't need to check bc already gone over via DA3 >= 1.3 check
		i5 = rws["INTR"]
		if lowFit5 == 0:
			lowFit5 = 0.000000000001
		nd = f5/lowFit5

		EDIT5 = EDIT5.append({'C5':c5,'F5':f5,'C4':c4,'F4':f4,'C3':c3,'F3':f3,'C2':c2,'F2':f2,'C1':c1,'F1':f1,'newDA':nd,'INTR':i5}, ignore_index=True)

	f_name = "/u/.../HD/nC5/%s%s.csv"%(COMB5,str(number5))
	EDIT5.to_csv(f_name,index = False)
	
	number5 = number5 + 1

