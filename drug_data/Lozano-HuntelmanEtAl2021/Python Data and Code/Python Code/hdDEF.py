# [1] --------- SETUP ----------
# Import libraries
import pandas as pd #use "pd" to call for pandas
import itertools # for combinations
import os

## NOTE: very similar to spDEF, classification is slightly different
# [2] --- 3 combo ---
YN3 = {}
fileList3 = os.listdir(r"/u/.../HD/nC3") ##### FOR CSVS THAT HAD DA<1.3, LOOKING FOR HIDDEN SUPPRESSION #####

x3 = 0
for file in fileList3:
	store3 = r"/u/.../HD/nC3/%s"%(file)
	YN3[x3] = pd.read_csv(store3)
	x3 = x3+1
rng3 = x3 + 1
#----------------- # category counts
types3 = pd.DataFrame(columns = ['FSH','PH','NA'])
FSH3 = 0 # fully hidden suppression
PH3 = 0 # partially hidden
NA3 = 0 # no "case"
#----------------- # combo lists
fsh3 = pd.DataFrame(columns = ['combos','INTR']) # intr for "original" interaction
ph3 = pd.DataFrame(columns = ['combos','2 only','INTR'])
na3 = pd.DataFrame(columns = ['combos','INTR'])
#----------------- # nDA under 1.3
nDAF3 = pd.DataFrame(columns = ['number']) # "new" DA doesn't pass check
n3fail = pd.DataFrame(columns = ['combos','INTR'])
ndaf3 = 0


for num in range (0,1428): ##### CHANGE VALUE DEPENDING ON CSVS APPLICABLE #####
	look3 = YN3[num] # looping through dictionary - 1 combo per dictionary key
	name = look3.iloc[0]["C3"]
	i3 = look3.iloc[0]["INTR"]
	#----------------- overall counts
	nested_count = 0
	LRsml_count = 0
	f32_count = 0
	#----------------- PH counts
	h2 = 0
	
	newDA = look3.iloc[0]["newDA"]
	if newDA >= 1.3:
		for index, rows in look3.iterrows(): # counts
			C3 = rows["C3"]
			F3 = rows["F3"]
			C2 = rows["C2"]
			F2 = rows["F2"]
			C1 = rows["C1"]
			F1 = rows["F1"]

			# -------------- overall counts
			if F3 > F2:    
				nested_count = nested_count + 1
			if F3 < F1:
				LRsml_count = LRsml_count + 1
			if F3 > F2:
				f32_count = f32_count + 1

			# -------------- PH counts
			if F3 > F2:
				h2 = h2 + 1

		# classification
		if LRsml_count == 6 and f32_count == 6: # fully suppressed hidden
			FSH3 = FSH3 + 1
			fsh3 = fsh3.append({'combos':name,'INTR':i3}, ignore_index=True)
		elif LRsml_count == 6 and f32_count > 0: # partially hidden
			PH3 = PH3 + 1
			ph3 = ph3.append({'combos':name,'2 only':h2,'INTR':i3}, ignore_index=True)
		elif LRsml_count != 6:
			NA3 = NA3 + 1
			na3 = na3.append({'combos':name,'INTR':i3}, ignore_index=True)
	else:
		ndaf3 = ndaf3 + 1
		n3fail = n3fail.append({'combos':name,'INTR':i3}, ignore_index=True)

types3 = types3.append({'FSH':FSH3,'PH':PH3,'NA':NA3}, ignore_index=True)
types3.to_csv(r"/u/.../HD/hdSORT/H3_Counts.csv",index=False)
fsh3.to_csv(r"/u/.../HD/hdSORT/fsh3_Combos.csv",index=False)
ph3.to_csv(r"/u/.../HD/hdSORT/ph3_Combos.csv",index=False)
nDAF3 = nDAF3.append({'number':ndaf3}, ignore_index=True)
nDAF3.to_csv(r"/u/.../HD/hdSORT/nDA3_fail.csv",index=False)
n3fail.to_csv(r"/u/.../HD/hdSORT/nDA3_fail_Combos.csv",index=False)
na3.to_csv(r"/u/.../HD/hdSORT/NA3_Combos.csv",index=False)


# [3] --- 4 combo ---
# similar structure as above
YN4 = {}
fileList4 = os.listdir(r"/u/.../HD/nC4")

x4 = 0
for file in fileList4:
	store4 = r"/u/.../HD/nC4/%s"%(file)
	YN4[x4] = pd.read_csv(store4)
	x4 = x4+1
rng4 = x4 + 1

#----------------- # category counts
types4 = pd.DataFrame(columns = ['FNH','NH','FSH','PH','NA'])
FNH4 = 0 # fully nested hidden suppression
NH4 = 0 # nested hidden suppression
FSH4 = 0
PH4 = 0
NA4 = 0
#----------------- # combo lists
fnh4 = pd.DataFrame(columns = ['combos','INTR'])
nh4 = pd.DataFrame(columns = ['combos','INTR'])
fsh4 = pd.DataFrame(columns = ['combos','INTR'])
ph4 = pd.DataFrame(columns = ['combos','2 only','3 only','2 and 3','INTR'])
na4 = pd.DataFrame(columns = ['combos','INTR'])
#----------------- # nDA under 1.3
nDAF4 = pd.DataFrame(columns = ['number'])
n4fail = pd.DataFrame(columns = ['combos','INTR'])
ndaf4 = 0

for num in range (0,4849): ##### CHANGE VALUE DEPENDING ON CSVS APPLICABLE #####
	look4 = YN4[num] # looping through dictionary - 1 combo per dictionary key
	name = look4.iloc[0]["C4"]
	i4 = look4.iloc[0]["INTR"]

	#----------------- overall counts
	nested_count = 0
	LRsml_count = 0
	f432_count = 0
	#----------------- PH counts
	h2 = 0
	h3 = 0
	h23 = 0
	
	newDA = look4.iloc[0]["newDA"]
	if newDA >= 1.3:    
		for index, rows in look4.iterrows(): # counts
			C4 = rows["C4"]
			F4 = rows["F4"]
			C3 = rows["C3"]
			F3 = rows["F3"]
			C2 = rows["C2"]
			F2 = rows["F2"]
			C1 = rows["C1"]
			F1 = rows["F1"]

			# -------------- overall counts
			if F4>F3 and F4>F2 and F3 > F2:    
				nested_count = nested_count + 1
			if F4 < F1:
				LRsml_count = LRsml_count + 1
			if F4 > F3 and F4 > F2:
				f432_count = f432_count + 1

			# -------------- PH counts
			if F4 > F2 and F4 < F3:
				h2 = h2 + 1
			if F4 < F2 and F4 > F3:
				h3 = h3 + 1
			if F4 > F2 and F4 > F3: 
				h23 = h23 + 1

		# classification
		if nested_count == 24: # fully nested suppression
			FNH4 = FNH4 + 1
			fnh4 = fnh4.append({'combos':name,'INTR':i4}, ignore_index=True)
		elif nested_count > 0: # nested suppression
			NH4 = NH4 + 1
			nh4 = nh4.append({'combos':name,'INTR':i4}, ignore_index=True)
		elif nested_count == 0:
			if LRsml_count == 24 and f432_count == 24: # fully suppressed
				FSH4 = FSH4 + 1
				fsh4 = fsh4.append({'combos':name,'INTR':i4}, ignore_index=True)
			elif LRsml_count == 24 and f432_count > 0: # partially hidden
				PH4 = PH4 + 1
				ph4 = ph4.append({'combos':name,'2 only':h2,'3 only':h3,'2 and 3':h23,'INTR':i4}, ignore_index=True)
			elif LRsml_count != 24:
				NA4 = NA4 + 1
				na4 = na4.append({'combos':name,'INTR':i4}, ignore_index=True)

	else:
		ndaf4 = ndaf4 + 1
		n4fail = n4fail.append({'combos':name,'INTR':i4}, ignore_index=True)
			
types4 = types4.append({'FNH':FNH4,'NH':NH4,'FSH':FSH4,'PH':PH4,'NA':NA4}, ignore_index=True)
types4.to_csv(r"/u/.../HD/hdSORT/H4_Counts.csv",index=False)
fnh4.to_csv(r"/u/.../HD/hdSORT/fnh4_Combos.csv",index=False)
nh4.to_csv(r"/u/.../HD/hdSORT/nh4_Combos.csv",index=False)
fsh4.to_csv(r"/u/.../HD/hdSORT/fsh4_Combos.csv",index=False)
ph4.to_csv(r"/u/.../HD/hdSORT/ph4_Combos.csv",index=False)
nDAF4 = nDAF4.append({'number':ndaf4}, ignore_index=True)
nDAF4.to_csv(r"/u/.../HD/hdSORT/nDA4_fail.csv",index=False)
n4fail.to_csv(r"/u/.../HD/hdSORT/nDA4_fail_Combos.csv",index=False)
na4.to_csv(r"/u/.../HD/hdSORT/NA4_Combos.csv",index=False)


# [4] --- 5 combo ---
YN5 = {}
fileList5 = os.listdir(r"/u/.../HD/nC5")

x5 = 0
for file in fileList5:
	store5 = r"/u/.../HD/nC5/%s"%(file)
	YN5[x5] = pd.read_csv(store5)
	x5 = x5+1
rng5 = x5 + 1

#----------------- # category counts
types5 = pd.DataFrame(columns = ['FNH','NH','FSH','PH','NA'])
FNH5 = 0
NH5 = 0
FSH5 = 0
PH5 = 0
NA5 = 0
#----------------- # combo lists
fnh5 = pd.DataFrame(columns = ['combos','INTR'])
nh5 = pd.DataFrame(columns = ['combos','INTR'])
fsh5 = pd.DataFrame(columns = ['combos','INTR'])
ph5 = pd.DataFrame(columns = ['combos','2 only','3 only','4 only','2 and 3','2 and 4','3 and 4','234','INTR'])
na5 = pd.DataFrame(columns = ['combos','INTR'])
#----------------- # nDA under 1.3
nDAF5 = pd.DataFrame(columns = ['number'])
n5fail = pd.DataFrame(columns = ['combos','INTR'])
ndaf5 = 0

for num in range (0,10972): ##### CHANGE VALUE DEPENDING ON CSVS APPLICABLE #####
	look5 = YN5[num] # looping through dictionary - 1 combo per dictionary key
	name = look5.iloc[0]["C5"]
	i5 = look5.iloc[0]["INTR"]

	#----------------- overall counts
	nested_count = 0
	LRsml_count = 0
	f5432_count = 0
	#----------------- PH counts
	h2 = 0
	h3 = 0
	h4 = 0
	h23 = 0
	h24 = 0
	h34 = 0
	h234 = 0
	
	newDA = look5.iloc[0]["newDA"]
	if newDA >= 1.3:    
		for index, rows in look5.iterrows(): # counts
			C5 = rows["C5"]
			F5 = rows["F5"]
			C4 = rows["C4"]
			F4 = rows["F4"]
			C3 = rows["C3"]
			F3 = rows["F3"]
			C2 = rows["C2"]
			F2 = rows["F2"]
			C1 = rows["C1"]
			F1 = rows["F1"]

			# -------------- overall counts
			if F5>F4 and F5>F3 and F5>F2 and F4>F3 and F4>F2 and F3 > F2:    
				nested_count = nested_count + 1
			if F5 < F1:
				LRsml_count = LRsml_count + 1
			if F5 > F4 and F5 > F3 and F4 > F2:
				f5432_count = f5432_count + 1

			# -------------- PH counts
			if F5 > F2 and F5 < F3 and F5 < F4:
				h2 = h2 + 1
			if F5 < F2 and F5 > F3 and F5 < F4:
				h3 = h3 + 1
			if F5 < F2 and F5 < F3 and F5 > F4:
				h4 = h4 + 1
			if F5 > F2 and F5 > F3 and F5 < F4: 
				h23 = h23 + 1
			if F5 > F2 and F5 < F3 and F5 > F4:
				h24 = h24 + 1
			if F5 < F2 and F5 > F3 and F5 > F4:
				h34 = h34 + 1
			if F5 > F2 and F5 > F3 and F5 > F4:
				h234 = h234 + 1

		# classification
		if nested_count == 120: # fully nested suppression
			FNH5 = FNH5 + 1
			fnh5 = fnh5.append({'combos':name,'INTR':i5}, ignore_index=True)
		elif nested_count > 0: # nested suppression
			NH5 = NH5 + 1
			nh5 = nh5.append({'combos':name,'INTR':i5}, ignore_index=True)
		elif nested_count == 0:
			if LRsml_count == 120 and f5432_count == 120: # fully suppressed
				FSH5 = FSH5 + 1
				fsh5 = fsh5.append({'combos':name,'INTR':i5}, ignore_index=True)
			elif LRsml_count == 120 and f5432_count > 0: # partially hidden
				PH5 = PH5 + 1
				ph5 = ph5.append({'combos':name,'2 only':h2,'3 only':h3,'4 only':h4,'2 and 3':h23,'2 and 4':h24,'3 and 4':h34,'234':h234,'INTR':i5}, ignore_index=True)
			elif LRsml_count != 120:
				NA5 = NA5 + 1
				na5 = na5.append({'combos':name,'INTR':i5}, ignore_index=True)

	else:
		ndaf5 = ndaf5 + 1
		n5fail = n5fail.append({'combos':name,'INTR':i5}, ignore_index=True)
				
types5 = types5.append({'FNH':FNH5,'NH':NH5,'FSH':FSH5,'PH':PH5,'NA':NA5}, ignore_index=True)
types5.to_csv(r"/u/.../HD/hdSORT/H5_Counts.csv",index=False)
fnh5.to_csv(r"/u/.../HD/hdSORT/fnh5_Combos.csv",index=False)
nh5.to_csv(r"/u/.../HD/hdSORT/nh5_Combos.csv",index=False)
fsh5.to_csv(r"/u/.../HD/hdSORT/fsh5_Combos.csv",index=False)
ph5.to_csv(r"/u/.../HD/hdSORT/ph5_Combos.csv",index=False)
nDAF5 = nDAF5.append({'number':ndaf5}, ignore_index=True)
nDAF5.to_csv(r"/u/.../HD/hdSORT/nDA5_fail.csv",index=False)
n5fail.to_csv(r"/u/.../HD/hdSORT/nDA5_fail_Combos.csv",index=False)
na5.to_csv(r"/u/.../HD/hdSORT/NA5_Combos.csv",index=False)


