# [1] --------- SETUP ----------
# Import libraries
import pandas as pd # use "pd" to call for pandas
import itertools # for combinations
import os

# [2] --- IMPORT ---
DAF2 = pd.read_csv(r"/u/.../2Drug.csv") # reading in csvs with fitness info
DAF3 = pd.read_csv(r"/u/.../3Drug.csv")
DAF4 = pd.read_csv(r"/u/.../4Drug.csv")
DAF5 = pd.read_csv(r"/u/.../5Drug.csv")
E1 = pd.read_csv(r"/u/.../F1.csv")

# [3] --- 3 combo ---
COMB3 = "C3_" # for naming purposes
number3 = 0 # for naming purposes
for index, row in DAF3.iterrows():
	SORTING3 = pd.DataFrame(columns = ['C3','F3','C2','F2','C1','F1','INTR'])

	if row["DA"] < 1.30: ##### SEARCH FOR HIDDEN SUPPRESSION #### >= 1.30 USED FOR SUPPRESSIVE COMBINATIONS #####
		d1 = row["Drug1"] # NOTE: uses title of columns
		d2 = row["Drug2"]
		d3 = row["Drug3"]
		f3 = row["fitness"] # goes into PD
		c3 = frozenset((d1,d2,d3))
		print3 = d1+d2+d3 # goes into PD
		intr3 = row["interaction"]
		
		f2 = 0
		print2 = 0

		combo_gen2 = [d1,d2,d3]
		combo2 = list(itertools.combinations(combo_gen2,2)) # combos within higher level combo
		for i in combo2:
			dr1, dr2 = i
			c2 = frozenset((dr1,dr2)) # for comparison
			print2 = dr1+dr2 # rewrites above; goes into PD
			
			for a,b in DAF2.iterrows():
				r1 = b["Drug1"]
				r2 = b["Drug2"]
				comp2 = frozenset((r1,r2)) # for comparison (current combo w/ what's generated from above)
				if comp2 == c2:
					f2 = b["fitness"] # rewrites above; goes into PD
					break
		
			for ind, rw in E1.iterrows():
				drug1 = rw["DRUG"]
				if drug1 == dr1:
					f1 = rw["MEDIAN"]
					SORTING3 = SORTING3.append({'C3':print3,'F3':f3,'C2':print2,'F2':f2,'C1':drug1,'F1':f1,'INTR':intr3}, ignore_index=True) # a single path

				if drug1 == dr2:
					f1 = rw["MEDIAN"]
					SORTING3 = SORTING3.append({'C3':print3,'F3':f3,'C2':print2,'F2':f2,'C1':drug1,'F1':f1,'INTR':intr3}, ignore_index=True)

		f_name = "/u/.../HD/C3/%s%s.csv"%(COMB3,str(number3))
		SORTING3.to_csv(f_name) # export csv
	
		number3 = number3 + 1

# [4] --- 4 COMBO ---
# the format is largely the same
COMB4 = "C4_"
number4 = 0
for index, row in DAF4.iterrows():
	SORTING4 = pd.DataFrame(columns = ['C4','F4','C3','F3','C2','F2','C1','F1','INTR'])

	if row["DA"] < 1.30:
		d1 = row["Drug1"]
		d2 = row["Drug2"]
		d3 = row["Drug3"]
		d4 = row["Drug4"]
		f4 = row["fitness"] #goes into PD
		c4 = frozenset((d1,d2,d3,d4))
		print4 = d1+d2+d3+d4 # goes into PD
		intr4 = row["interaction"]
		
		f3 = 0
		print3 = 0
		f2 = 0
		print2 = 0

		combo_gen3 = [d1,d2,d3,d4]
		combo3 = list(itertools.combinations(combo_gen3,3))
		for i in combo3:
			dr1, dr2, dr3 = i
			c3 = frozenset((dr1,dr2,dr3)) # for comparison
			print3 = dr1+dr2+dr3 # rewrites above; goes into PD
			
			for a,b in DAF3.iterrows():
				r1 = b["Drug1"]
				r2 = b["Drug2"]
				r3 = b["Drug3"]
				comp3 = frozenset((r1,r2,r3)) # for comparison
				if comp3 == c3:
					f3 = b["fitness"] # rewrites above; goes into PD
					break
		
			combo_gen2 = [dr1,dr2,dr3]
			combo2 = list(itertools.combinations(combo_gen2,2))
			for j in combo2:
				du1, du2 = j
				c2 = frozenset((du1,du2)) # for comparison
				print2 = du1+du2 # rewrites above; goes into PD
				
				for w, x in DAF2.iterrows():
					u1 = x["Drug1"]
					u2 = x["Drug2"]
					comp2 = frozenset((u1,u2)) # for comparison
					if comp2 == c2:
						f2 = x["fitness"] # rewrites above; goes into PD
						break

				for ind, rw in E1.iterrows():
					drug1 = rw["DRUG"]
					if drug1 == du1:
						f1 = rw["MEDIAN"]
						SORTING4 = SORTING4.append({'C4':print4,'F4':f4,'C3':print3,'F3':f3,'C2':print2,'F2':f2,'C1':drug1,'F1':f1,'INTR':intr4}, ignore_index=True) # single path
						
					if drug1 == du2:
						f1 = rw["MEDIAN"]
						SORTING4 = SORTING4.append({'C4':print4,'F4':f4,'C3':print3,'F3':f3,'C2':print2,'F2':f2,'C1':drug1,'F1':f1,'INTR':intr4}, ignore_index=True)

		f_name = "/u/.../HD/C4/%s%s.csv"%(COMB4,str(number4))
		SORTING4.to_csv(f_name)
	
		number4 = number4 + 1

# [5] --- 5 COMBO ---
COMB5 = "C5_"
number5 = 0
for index, row in DAF5.iterrows():
	SORTING = pd.DataFrame(columns = ['C5','F5','C4','F4','C3','F3','C2','F2','C1','F1','INTR'])
	
	if row["DA"] < 1.30:
		d1 = row["Drug1"]
		d2 = row["Drug2"]
		d3 = row["Drug3"]
		d4 = row["Drug4"]
		d5 = row["Drug5"]
		f5 = row["fitness"] #goes into PD
		c5 = frozenset((d1,d2,d3,d4,d5))
		print5 = d1+d2+d3+d4+d5 # goes into PD
		intr5 = row["interaction"]
		
		f4 = 0
		print4 = 0
		f3 = 0
		print3 = 0
		f2 = 0
		print2 = 0

		combo_gen4 = [d1,d2,d3,d4,d5]
		combo4 = list(itertools.combinations(combo_gen4,4))
		for i in combo4:
			dr1, dr2, dr3, dr4 = i
			c4 = frozenset((dr1,dr2,dr3,dr4)) # for comparison
			print4 = dr1+dr2+dr3+dr4 # rewrites above; goes into PD
			
			for a,b in DAF4.iterrows():
				r1 = b["Drug1"]
				r2 = b["Drug2"]
				r3 = b["Drug3"]
				r4 = b["Drug4"]
				#f4 = b["Fitness"]
				comp4 = frozenset((r1,r2,r3,r4)) # for comparison
				if comp4 == c4:
					f4 = b["fitness"] # rewrites above; goes into PD
					break
		
			combo_gen3 = [dr1,dr2,dr3,dr4]
			combo3 = list(itertools.combinations(combo_gen3,3))
			for j in combo3:
				du1, du2, du3 = j
				c3 = frozenset((du1,du2,du3)) # for comparison
				print3 = du1+du2+du3 # rewrites above; goes into PD
				
				for w, x in DAF3.iterrows():
					u1 = x["Drug1"]
					u2 = x["Drug2"]
					u3 = x["Drug3"]
					comp3 = frozenset((u1,u2,u3)) # for comparison
					if comp3 == c3:
						f3 = x["fitness"] # rewrites above; goes into PD
						break
						
				combo_gen2 = [du1,du2,du3]
				combo2 = list(itertools.combinations(combo_gen2,2))
				for k in combo2:
					dg1, dg2 = k
					c2 = frozenset((dg1,dg2)) # for comparison
					print2 = dg1+dg2 # rewrites above; goes into PD
					
					for y, z in DAF2.iterrows():
						g1 = z["Drug1"]
						g2 = z["Drug2"]
						comp2 = frozenset((g1,g2)) # for comparison
						if comp2 == c2:
							f2 = z["fitness"] # rewrites above; goes into PD
							break
							
					for ind, rw in E1.iterrows():
						drug1 = rw["DRUG"]
						if drug1 == dg1:
							f1 = rw["MEDIAN"]
							SORTING = SORTING.append({'C5':print5,'F5':f5,'C4':print4,'F4':f4,'C3':print3,'F3':f3,'C2':print2,'F2':f2,'C1':drug1,'F1':f1,'INTR':intr5}, ignore_index=True) # single path

						if drug1 == dg2:
							f1 = rw["MEDIAN"]
							SORTING = SORTING.append({'C5':print5,'F5':f5,'C4':print4,'F4':f4,'C3':print3,'F3':f3,'C2':print2,'F2':f2,'C1':drug1,'F1':f1,'INTR':intr5}, ignore_index=True)
		
		f_name = "/u/.../HD/C5/%s%s.csv"%(COMB5,str(number5))
		SORTING.to_csv(f_name)
	
		number5 = number5 + 1



