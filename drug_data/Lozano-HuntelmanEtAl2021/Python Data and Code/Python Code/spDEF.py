# [1] --------- SETUP ----------
# Import libraries
import pandas as pd #use "pd" to call for pandas
import itertools # for combinations
import os

# [2] --- 3 combo ---
YN3 = {}
fileList3 = os.listdir(r"/u/.../SUPR/yC3") ##### FOR CSVS THAT PASSED DA>=1.3 CHECK #####

x3 = 0
for file in fileList3:
	store3 = r"/u/.../SUPR/yC3/%s"%(file)
	YN3[x3] = pd.read_csv(store3)
	x3 = x3+1
rng3 = x3 + 1
#----------------- # category counts
types3 = pd.DataFrame(columns = ['FNS','NS','FS','PH','NHS','NA'])
FNS3 = 0 # fully nested suppression
NS3 = 0 # nested suppression
FS3 = 0 # fully suppressed
PH3 = 0 # partially hidden suppression
NHS3 = 0 # no hidden suppression
NA3 = 0 # no "case"
#----------------- # combo lists
fns3 = pd.DataFrame(columns = ['combos','INTR']) # intr is "original" interaction
ns3 = pd.DataFrame(columns = ['combos','INTR'])
fs3 = pd.DataFrame(columns = ['combos','INTR'])
ph3 = pd.DataFrame(columns = ['combos','2 only','INTR'])
nhs3 = pd.DataFrame(columns = ['combos','INTR'])
na3 = pd.DataFrame(columns = ['combos','INTR'])

for num in range (0,84): ##### CHANGE VALUE DEPENDING ON CSVS APPLICABLE #####
	look3 = YN3[num] # looping through dictionary - 1 combo per dictionary key
	name = look3.iloc[0]["C3"]
	i3 = look3.iloc[0]["INTR"]
	#----------------- overall counts
	nested_count = 0
	LRsml_count = 0
	f32_count = 0
	#----------------- PH counts
	h2 = 0
	
	for index, rows in look3.iterrows(): # counts
		C3 = rows["C3"]
		F3 = rows["F3"]
		C2 = rows["C2"]
		F2 = rows["F2"]
		C1 = rows["C1"]
		F1 = rows["F1"]

		# -------------- overall counts
		if F3 > F2 and F3 > F1 and F2 > F1:    
			nested_count = nested_count + 1
		if F3 > F1:
			LRsml_count = LRsml_count + 1
		if F3 > F2:
			f32_count = f32_count + 1
		
		# -------------- PH counts
		if F3 > F2:
			h2 = h2 + 1
		
	# classification
	if nested_count == 6: # fully nested suppression
		FNS3 = FNS3 + 1
		fns3 = fns3.append({'combos':name,'INTR':i3}, ignore_index=True)
	elif nested_count > 0: # nested suppression
		NS3 = NS3 + 1
		ns3 = ns3.append({'combos':name,'INTR':i3}, ignore_index=True)
	elif nested_count == 0:
		if LRsml_count > 0 and f32_count == 6: # fully suppressed
			FS3 = FS3 + 1
			fs3 = fs3.append({'combos':name,'INTR':i3}, ignore_index=True)
		elif LRsml_count > 0 and f32_count > 0: # partially hidden
			PH3 = PH3 + 1
			ph3 = ph3.append({'combos':name,'2 only':h2,'INTR':i3}, ignore_index=True)
		elif LRsml_count > 0 and f32_count == 0: # no hidden
			NHS3 = NHS3 + 1
			nhs3 = nhs3.append({'combos':name,'INTR':i3}, ignore_index=True)
		else:
			NA3 = NA3 + 1
			na3 = na3.append({'combos':name,'INTR':i3}, ignore_index=True)

			
types3 = types3.append({'FNS':FNS3,'NS':NS3,'FS':FS3,'PH':PH3,'NHS':NHS3,'NA':NA3}, ignore_index=True)
types3.to_csv(r"/u/.../SUPR/suprSORT/Supr3_Counts.csv",index=False)
fns3.to_csv(r"/u/.../SUPR/suprSORT/fns3_Combos.csv",index=False)
ns3.to_csv(r"/u/.../SUPR/suprSORT/ns3_Combos.csv",index=False)
fs3.to_csv(r"/u/.../SUPR/suprSORT/fs3_Combos.csv",index=False)
ph3.to_csv(r"/u/.../SUPR/suprSORT/ph3_Combos.csv",index=False)
nhs3.to_csv(r"/u/.../SUPR/suprSORT/nhs3_Combos.csv",index=False)
na3.to_csv(r"/u/.../SUPR/suprSORT/NA3_Combos.csv",index=False)

# [3] --- 4 combo ---
# follows same structure
YN4 = {}
fileList4 = os.listdir(r"/u/.../SUPR/yC4")

x4 = 0
for file in fileList4:
	store4 = r"/u/.../SUPR/yC4/%s"%(file)
	YN4[x4] = pd.read_csv(store4)
	x4 = x4+1
rng4 = x4 + 1

#----------------- # category counts
types4 = pd.DataFrame(columns = ['FNS','NS','FS','PH','NHS','NA'])
FNS4 = 0
NS4 = 0
FS4 = 0
PH4 = 0
NHS4 = 0
NA4 = 0
#----------------- # combo lists
fns4 = pd.DataFrame(columns = ['combos','INTR'])
ns4 = pd.DataFrame(columns = ['combos','INTR'])
fs4 = pd.DataFrame(columns = ['combos','INTR'])
ph4 = pd.DataFrame(columns = ['combos','2 only','3 only','2 and 3','INTR'])
nhs4 = pd.DataFrame(columns = ['combos','INTR'])
na4 = pd.DataFrame(columns = ['combos','INTR'])


for num in range (0,821): ##### CHANGE VALUE DEPENDING ON CSVS APPLICABLE #####
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
		if F4>F3 and F4>F2 and F4>F1 and F3 > F2 and F3 > F1 and F2 > F1:    
			nested_count = nested_count + 1
		if F4 > F1:
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
		FNS4 = FNS4 + 1
		fns4 = fns4.append({'combos':name,'INTR':i4}, ignore_index=True)
	elif nested_count > 0: # nested suppression
		NS4 = NS4 + 1
		ns4 = ns4.append({'combos':name,'INTR':i4}, ignore_index=True)
	elif nested_count == 0:
		if LRsml_count > 0 and f432_count == 24: # fully suppressed
			FS4 = FS4 + 1
			fs4 = fs4.append({'combos':name,'INTR':i4}, ignore_index=True)
		elif LRsml_count > 0 and f432_count > 0: # partially hidden
			PH4 = PH4 + 1
			ph4 = ph4.append({'combos':name,'2 only':h2,'3 only':h3,'2 and 3':h23,'INTR':i4}, ignore_index=True)
		elif LRsml_count > 0 and f432_count == 0: # no hidden
			NHS4 = NHS4 + 1
			nhs4 = nhs4.append({'combos':name,'INTR':i4}, ignore_index=True)
		else:
			NA4 = NA4 + 1
			na4 = na4.append({'combos':name,'INTR':i4}, ignore_index=True)
			
types4 = types4.append({'FNS':FNS4,'NS':NS4,'FS':FS4,'PH':PH4,'NHS':NHS4,'NA':NA4}, ignore_index=True)
types4.to_csv(r"/u/.../SUPR/suprSORT/Supr4_Counts.csv",index=False)
fns4.to_csv(r"/u/.../SUPR/suprSORT/fns4_Combos.csv",index=False)
ns4.to_csv(r"/u/.../SUPR/suprSORT/ns4_Combos.csv",index=False)
fs4.to_csv(r"/u/.../SUPR/suprSORT/fs4_Combos.csv",index=False)
ph4.to_csv(r"/u/.../SUPR/suprSORT/ph4_Combos.csv",index=False)
nhs4.to_csv(r"/u/.../SUPR/suprSORT/nhs4_Combos.csv",index=False)
na4.to_csv(r"/u/.../SUPR/suprSORT/NA4_Combos.csv",index=False)

# [4] --- 5 combo ---
YN5 = {}
fileList5 = os.listdir(r"/u/.../SUPR/yC5")

x5 = 0
for file in fileList5:
	store5 = r"/u/.../SUPR/yC5/%s"%(file)
	YN5[x5] = pd.read_csv(store5)
	x5 = x5+1
rng5 = x5 + 1

#----------------- # category counts
types5 = pd.DataFrame(columns = ['FNS','NS','FS','PH','NHS','NA'])
FNS5 = 0
NS5 = 0
FS5 = 0
PH5 = 0
NHS5 = 0
NA5 = 0
#----------------- # combo lists
fns5 = pd.DataFrame(columns = ['combos','INTR'])
ns5 = pd.DataFrame(columns = ['combos','INTR'])
fs5 = pd.DataFrame(columns = ['combos','INTR'])
ph5 = pd.DataFrame(columns = ['combos','2 only','3 only','4 only','2 and 3','2 and 4','3 and 4','234','INTR'])
nhs5 = pd.DataFrame(columns = ['combos','INTR'])
na5 = pd.DataFrame(columns = ['combos','INTR'])


for num in range (0,2630): ##### CHANGE VALUE DEPENDING ON CSVS APPLICABLE #####
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
		if F5>F4 and F5>F3 and F5>F2 and F5>F1 and F4>F3 and F4>F2 and F4>F1 and F3 > F2 and F3 > F1 and F2 > F1:    
			nested_count = nested_count + 1
		if F5 > F1:
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
		FNS5 = FNS5 + 1
		fns5 = fns5.append({'combos':name,'INTR':i5}, ignore_index=True)
	elif nested_count > 0: # nested suppression
		NS5 = NS5 + 1
		ns5 = ns5.append({'combos':name,'INTR':i5}, ignore_index=True)
	elif nested_count == 0:
		if LRsml_count > 0 and f5432_count == 120: # fully suppressed
			FS5 = FS5 + 1
			fs5 = fs5.append({'combos':name,'INTR':i5}, ignore_index=True)
		elif LRsml_count > 0 and f5432_count > 0: # partially hidden
			PH5 = PH5 + 1
			ph5 = ph5.append({'combos':name,'2 only':h2,'3 only':h3,'4 only':h4,'2 and 3':h23,'2 and 4':h24,'3 and 4':h34,'234':h234,'INTR':i5}, ignore_index=True)
		elif LRsml_count > 0 and f5432_count == 0: # no hidden
			NHS5 = NHS5 + 1
			nhs5 = nhs5.append({'combos':name,'INTR':i5}, ignore_index=True)
		else:
			NA5 = NA5 + 1
			na5 = na5.append({'combos':name,'INTR':i5}, ignore_index=True)
			
types5 = types5.append({'FNS':FNS5,'NS':NS5,'FS':FS5,'PH':PH5,'NHS':NHS5,'NA':NA5}, ignore_index=True)
types5.to_csv(r"/u/.../SUPR/suprSORT/Supr5_Counts.csv",index=False)
fns5.to_csv(r"/u/.../SUPR/suprSORT/fns5_Combos.csv",index=False)
ns5.to_csv(r"/u/.../SUPR/suprSORT/ns5_Combos.csv",index=False)
fs5.to_csv(r"/u/.../SUPR/suprSORT/fs5_Combos.csv",index=False)
ph5.to_csv(r"/u/.../SUPR/suprSORT/ph5_Combos.csv",index=False)
nhs5.to_csv(r"/u/.../SUPR/suprSORT/nhs5_Combos.csv",index=False)
na5.to_csv(r"/u/.../SUPR/suprSORT/NA5_Combos.csv",index=False)

