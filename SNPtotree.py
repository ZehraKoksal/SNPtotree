print("STARTING SNPTOTREE")
import argparse
import time
import shutil
import glob
import os, sys
import pandas as pd
import numpy as np
import pprint
import re
import itertools
from itertools import product, groupby
from itertools import count
from itertools import combinations
from itertools import chain

start = time.time()

#User input paths
parser = argparse.ArgumentParser()
parser.add_argument("input", help="Tab-delimited csv file of a dataframe with individuals as columns and variants as rows. "
  "Insert variant data as 'A' for ancestral, 'D' for derived allele, and 'X' for missing data. Add variant names as first column and sequence name as header row. "
  "Include only polymorphic variants! Refrain from using comma separted marker or sequence names.")
parser.add_argument("output", help="Path to output folder")
# parser.add_argument("newick_output", help="Path to newick output folder")
parser.add_argument("-contradictory_variants", action="store_true", help="Optional: print contradictory variants based on pairwise comparisons in case contradictory variants exist for dataset. ")
parser.add_argument("-ambiguous_variants", action="store_true", help="Optional: print removed variants due to their ambiguous positions in the tree. More informative sequencing information could help find the correct tree position of these variants. ")
parser.add_argument("-metadata_individuals", action="store_true", help="Optional: print which individuals showed derived alleles for the variants in each tree row. For a set of multiple equal variants, the sequences of this row were derived for at least one of the variants. ")


args = parser.parse_args()


df = pd.read_csv(args.input, sep="\t")#, header=None)  #

#Check if there are elements other than A,D and X
df_flat = set(df.iloc[:,1:].to_numpy().flatten().tolist())
unexpect_character = df_flat - {"A","D","X"}
if len(unexpect_character) > 0:
    print("\n"+"Unexpected character(s) were found in the input file: "+str(unexpect_character)+" . Only ancestral (A), derived (D), and missing (X) alleles are accepted in the input file!")
    sys.exit()


#Search input file for variants that have no ancestral alleles
query_A_allele = df.apply(lambda x: x != "A" , axis=1)
non_poly = df[query_A_allele.agg(min, axis=1)]

if len(non_poly) > 0:
    non_poly_list = non_poly.iloc[:,0].values.tolist()
    print()
    print("Input file contains one or more variants that are not polymorphic. We will remove variant(s): "+str(non_poly_list)+" !")
    df=df[~df.iloc[:,0].isin(non_poly_list)]
    print()


#Insert root marker at top
root = ["D"] * len(df.axes[1])
root[0] = "Root"
df.loc[-1]=root
df.index = df.index + 1
df = df.sort_index()
print(df)

df_copy = df.copy()
marker_nr = list(range(len(df)))
all_samples = df.columns.values[1:].tolist()
marker_name= df.iloc[:,0].values.tolist() #

#make a marker dictionary based on the initial input overview instead of a separate dictionary to prevent naming the markers wrong because the input file has a different order e.g.
df_dic = pd.DataFrame()
marker_nr = [str(x) for x in marker_nr]
df_dic["marker_nr"]=marker_nr
df_dic["marker_name"]=marker_name


marker_count = len(df)
header_int = np.arange(marker_count).tolist()
header_int =[str(x) for x in header_int]

df_T =df.T 

df["Marker_nr"] = list(marker_nr)
df_T = df_T.replace(r"\(.*\)","X",regex=True) #r"\(.*\)" takes anything within parentheses --> bad quality results are also just "X"
# all_samples = df_T.iloc[0,1:].values.tolist()
# print(all_samples)


df_T = df_T.iloc[1:] 
df_T.columns = header_int


Output= []
Parallel_output = []
Remove_markers_list = []

#loop through the columns of markers
ds_comb = tuple(itertools.combinations(header_int,2))
print("Pairwise comparison is starting ...")
for pair in ds_comb:
    print(pair)
    mn_D_df=df_T[df_T[pair[0]].str.contains("D")]
    # observations of the paired column/marker
    mc_D=mn_D_df[mn_D_df[pair[1]].str.contains("D")]
    mc_A=mn_D_df[mn_D_df[pair[1]].str.contains("A")]
    mc_X=mn_D_df[mn_D_df[pair[1]].str.contains("X")]
    if len(mc_D) == 0 and len(mc_A) >= 1 and len(mc_X) == 0: #there is only A in alternative column
        Result_list=["parallel","upstream"]
    elif len(mc_D) >= 1 and len(mc_A) == 0 and len(mc_X) == 0: #there is only D in alternative column
        Result_list=["equal","downstream","upstream"]
    elif len(mc_D) >= 1 and len(mc_A) >= 1 and len(mc_X) == 0: #there is D and A in alternative column
        Result_list=["upstream"]
    elif len(mc_D) == 0 and len(mc_A) == 0 and len(mc_X) >= 1: #there is only X in alternative column
        Result_list=["equal","downstream","upstream","parallel"]
    elif len(mc_D) >= 1 and len(mc_A) == 0 and len(mc_X) >= 1: #there is D and X in alternative column
        Result_list=["equal","downstream","upstream"]
    elif len(mc_D) == 0 and len(mc_A) >= 1 and len(mc_X) >= 1: #there is A and X in alternative column
        Result_list=["parallel","upstream"]
    elif len(mc_D) >= 1 and len(mc_A) >= 1 and len(mc_X) >= 1: #there is D and A and X in alternative column
        Result_list=["upstream"]
    else:
        Result_list=["ERROR"]
    #comparison in opposite direction
    mn_D_df=df_T[df_T[pair[1]].str.contains("D")]
    mn_D=mn_D_df[mn_D_df[pair[0]].str.contains("D")]
    mn_A=mn_D_df[mn_D_df[pair[0]].str.contains("A")]
    mn_X=mn_D_df[mn_D_df[pair[0]].str.contains("X")]
    if len(mn_D) == 0 and len(mn_A) >= 1 and len(mn_X) >= 1: #there is A and X in alternative column
        Result2_list=["parallel","downstream"]
    elif len(mn_D) == 0 and len(mn_A) >= 1 and len(mn_X) == 0: #there is only A in alternative column
        Result2_list=["parallel","downstream"]
    elif len(mn_D) >= 1 and len(mn_A) == 0 and len(mn_X) == 0: #there is only D in alternative column
        Result2_list=["equal","downstream","upstream"]
    elif len(mn_D) >= 1 and len(mn_A) >= 1 and len(mn_X) == 0: #there is D and A in alternative column
        Result2_list=["downstream"]
    elif len(mn_D) == 0 and len(mn_A) == 0 and len(mn_X) >= 1: #there is only X in alternative column
        Result2_list=["equal","downstream","upstream","parallel"]
    elif len(mn_D) >= 1 and len(mn_A) == 0 and len(mn_X) >= 1: #there is D and X in alternative column
        Result2_list=["equal","downstream","upstream"]
    elif len(mn_D) >= 1 and len(mn_A) >= 1 and len(mn_X) >= 1: #there is D and A and X in alternative column
        Result2_list=["downstream"]
    else:
        Result_list=["ERROR"] #Error means there is no derived allele of this marker
    
    #Combine both comparisons
    same_value = set(Result_list) & set(Result2_list)
    #only used to take into account uncertainty
    if same_value == {"parallel"}: 
        parallel = (str(pair[0])+"/"+ str(pair[1]))
        Parallel_output = Parallel_output +[parallel]    
    elif "parallel" in same_value:
        speed_up = "speed_up" #to speed up in case of informative 
    elif same_value == {"upstream"}:
        upstr = (str(pair[0])+"-"+ str(pair[1]))
        Output = Output +[upstr]
    elif same_value == {"downstream"}:
        downstr = (str(pair[1])+"-"+ str(pair[0]))
        Output = Output +[downstr]
    elif "equal" in same_value:
        equal = (str(pair[0])+"&"+ str(pair[1]))
        Output = Output +[equal]
    ###
    elif len(same_value) == 0:
        rmv = (str(pair[0])+"_"+ str(pair[1]))
        Remove_markers_list = Remove_markers_list + [rmv]


Parallel_output = Parallel_output + Output #for statistical power determination --> is continued later on after tree generation
Pairwise_complete  = pd.DataFrame ()
Pairwise_complete['relationships'] = Parallel_output
Pairwise_complete[['M1','M2']] = Pairwise_complete.relationships.str.split("&|-|/", expand=True) 

Pairwise_complete_filter=Pairwise_complete[~Pairwise_complete["relationships"].str.contains("&")]

Pairwise_complete = Pairwise_complete.iloc[:,1:3] #filter to show the last 2 columns
Pairwise_complete_filter = Pairwise_complete_filter.iloc[:,1:3]


print("Pairwise relationships determined!")
Remove_markers_df = pd.DataFrame()
to_remove = []
#Generating contradicting marker list that will be removed from future analysis and can be saved 
if len(Remove_markers_list) != 0:
    Remove_markers_df["contr"]=Remove_markers_list
    Remove_markers_df[['parted1','parted2']] = Remove_markers_df["contr"].str.split("_", expand=True) #update if changing _part2.py script
    #Join all entries in one list
    first_col = Remove_markers_df["parted1"].values.tolist()
    second_col = Remove_markers_df["parted2"].values.tolist()
    # print(second_col)
    combined = first_col + second_col
    combined_unique = list(set(combined))
    #Count how many times each element is observed in both columns
    combined_copy=combined[:]
    combined.sort(key=lambda x:combined_copy.count(x)) #gives from least common to most common entries
    combined = list(dict.fromkeys(combined))
    combined.reverse() #needs to be reversed in the order, so the most common markers are in the beginning 
    #remove until dataframe gets empty
    n = 0
    while len(Remove_markers_df) != 0:
        Remove_markers_df=Remove_markers_df[~Remove_markers_df["parted1"].isin([combined[n]])]
        Remove_markers_df=Remove_markers_df[~Remove_markers_df["parted2"].isin([combined[n]])]
        to_remove.append(combined[n])
        first_col = Remove_markers_df["parted1"].values.tolist()
        second_col = Remove_markers_df["parted2"].values.tolist()
        combined = first_col + second_col
        combined_copy=combined[:]
        combined.sort(key=lambda x:combined_copy.count(x)) #gives from least common to most common entries
        combined = list(dict.fromkeys(combined))
        combined.reverse() #needs to be reversed in the order, so the most common markers are in the beginning 

   
<<<<<<< HEAD

if args.contradictory_variants:
    path = args.output + "contradictory_variants.csv"
=======
if type(args.contradictory_variants) ==str: # != False:
    path = args.contradictory_variants + "contradictory_variants.csv"
>>>>>>> b36a4c3f45f06f30125b0c39c479333dc2f3cab7
    to_remove_df = pd.DataFrame()
    to_remove_df[0] = to_remove
    to_remove_df[0] = to_remove_df[0].replace(marker_nr,marker_name)
    np.savetxt(path, to_remove_df, delimiter=",", fmt="%s", comments="")
print("Contradictory markers removed!")
#Preparation: Combining equal markers



#Tree with upstream, downstream markers
Tree_raw = pd.DataFrame()
Tree_raw['output'] = Output
Tree_raw[['parted1','parted2']] = Tree_raw.output.str.split("&|-", expand=True) 

#Filtering contradictory markers out
if len(to_remove) != 0:
    markers = set(header_int) - set(to_remove)
    markers = list(markers)
    Tree_raw_parted1=Tree_raw[~Tree_raw["parted2"].isin(to_remove)] #filter in table for entries from the first marker set in parted1 and in parted2 column
    Tree_raw=Tree_raw_parted1[~Tree_raw_parted1["parted1"].isin(to_remove)]
else:
    markers = header_int

print("Combining equal markers...")
Tree_raw_hyp=Tree_raw[Tree_raw["output"].str.contains("-")]
Tree_raw_hyp = Tree_raw_hyp.sort_values(by=["parted1"])
Tree_raw_hyp = Tree_raw_hyp.reset_index(drop=True)
Tree_raw_eq=Tree_raw[Tree_raw["output"].str.contains("&")]

#Combining equal markers

to_remove_2 = []
true_equals= []
eq_eq=[]
dwn_list = []


Equal_ds=Tree_raw_hyp[Tree_raw_hyp["parted1"].isin(markers)]
Equal_ds=Equal_ds[Equal_ds["parted2"].isin(markers)]
Equal_ds=Equal_ds.drop(columns=["output"])
Equal_ds.columns = ["equal_marker","downstream"]
Equal_ds = Equal_ds.reset_index(drop=True)

#Identifying downstream markers of each marker
eq_m=itertools.groupby(Equal_ds["equal_marker"])
for ds,key in eq_m:
    grouped = Equal_ds.groupby(['equal_marker'])
    df_eq = grouped.get_group(ds)
    dwn = df_eq.iloc[:,1].values.tolist()
    dwn.sort()
    dwn_list.append(dwn)
    eq_eq.append(ds) 

Equal_ds_2 = pd.DataFrame()
Equal_ds_2["equal_marker"] = eq_eq
Equal_ds_2["downstream"] = dwn_list


#turn list to strings for comparison
Equal_ds_strings =[]
for e in Equal_ds_2["downstream"]:
    m = ','.join(e)
    Equal_ds_strings.append(m)   
Equal_ds_2["downstream_str"] = Equal_ds_strings


#groupby downstream markers

Equal_ds_2 = Equal_ds_2.sort_values(by=["downstream_str"])
Equal_ds_2 = Equal_ds_2.reset_index(drop=True)

equal_ds_group=itertools.groupby(Equal_ds_2["downstream_str"])

equal_marker_names_list = []
equal_marker_ds_list = []
true_equals_raw=[]
for ds,key in equal_ds_group:
    grouped = Equal_ds_2.groupby(['downstream_str'])
    df_eq = grouped.get_group(ds)
    equal_marker_names = df_eq.iloc[:,0].values.tolist()
    equal_marker_names_list.append(equal_marker_names)
    equal_marker_ds = df_eq.iloc[0,1] #could save strings; save just first list entry, because they are all the same
    equal_marker_ds_list.append(equal_marker_ds)
true_equals_raw.append(equal_marker_names_list)
true_equals_raw = true_equals_raw[0]

# #make pairwise sets of the downstream markers to see if there are overlapping values in the lists
ds_comb = tuple(itertools.combinations(equal_marker_ds_list,2))
mk_comb = tuple(itertools.combinations(equal_marker_names_list,2))
cnt = 0
for c in ds_comb:
    overlapping = set(c[0])&set(c[1])
    current_marker_set = mk_comb[cnt]       #which marker pair combination goes with ds_comb
    if len(overlapping) > 0:
        ups_list = []
        eqs_list = []
        m1 = list(set(c[0])-set(c[1]))
        m2 = list(set(c[1])-set(c[0]))
        #if one is longer than the other:
        if len(m1)==0 or len(m2)==0:
            if len(m1)==0: #if m1 is zero, c[0] does not contain markers unique compared to c[1] --> the second marker, c[1] compared to c[0] is the upstream one, sicne it has additonal downstream markers
                #is the marker list with no unique ds markers among the downstream markers of the other marker?
                hyphen = set(current_marker_set[0])&set(m2)
                if len(hyphen)==0: #current_marker_set[0] is in m2: 
                    if current_marker_set[0] in true_equals_raw:
                        true_equals_raw.remove(current_marker_set[0])
                    to_remove_2.append(current_marker_set[0])  
            if len(m2)==0: #if m1 is zero, c[0] does not contain markers unique compared to c[1] --> the second marker, c[1] compared to c[0] is the upstream one, sicne it has additonal downstream markers
                hyphen = set(current_marker_set[1])&set(m1)
                if len(hyphen)==0: #current_marker_set[0] in m2: #c[1]:
                    if current_marker_set[1] in true_equals_raw:
                        true_equals_raw.remove(current_marker_set[1])
                    to_remove_2.append(current_marker_set[1])   

        elif len(m1)!=0 and len(m2)!=0:
            #remove that marker and their unique downstream markers, which has the smaller downstream marker list
            hyphen = set(current_marker_set[0])&set(m2)
            hyphen_1 = set(current_marker_set[1])&set(m1)
            if len(hyphen)==0 and len(hyphen_1)==0:
                if len(m1) < len(m2):
                    rmv = current_marker_set[0]                 
                    to_remove_2.append(rmv)
                    if current_marker_set[0] in true_equals_raw:
                        true_equals_raw.remove(current_marker_set[0])
                elif len(m2) < len(m1):
                    rmv = current_marker_set[1]                 
                    to_remove_2.append(rmv)
                    if current_marker_set[1] in true_equals_raw:
                        true_equals_raw.remove(current_marker_set[1])
                else:
                    rmv = current_marker_set[1]                 
                    to_remove_2.append(rmv)
                    if current_marker_set[1] in true_equals_raw:
                        true_equals_raw.remove(current_marker_set[1])
    cnt = cnt + 1

for l in true_equals_raw:
    if len(l)>1:
        l = '&'.join(map(str,l))
        true_equals.append(l)
        # print(true_equals)

#remove the removed markers from Equal_ds dataframe
prelim_to_remove_list = list(chain.from_iterable(to_remove_2))
prelim_to_remove_unique = list(set(prelim_to_remove_list))

Equal_ds=Equal_ds[~Equal_ds["equal_marker"].isin(prelim_to_remove_unique)]
Equal_ds=Equal_ds[~Equal_ds["downstream"].isin(prelim_to_remove_unique)]
Equal_ds = Equal_ds.reset_index(drop=True)


eq_m=itertools.groupby(Equal_ds["equal_marker"])

eq_eq=[]
dwn_list = []
for ds,key in eq_m:
    grouped = Equal_ds.groupby(['equal_marker'])
    df_eq = grouped.get_group(ds)
    dwn = df_eq.iloc[:,1].values.tolist()
    dwn_list.append(dwn)
    eq_eq.append(ds) 

markers = list(set(markers)-set(prelim_to_remove_unique))
identical = list(set(markers)-set(eq_eq))
while len(identical) != 0:
    initial_len = len(Equal_ds)
    #join equal markers without downstream markers with &
    #group identical markers, if they were found in one sample_list
    sample_list = []
    sample_list_str = []
    upstream_counter = []
    for e in identical: #range(sel_marker_counter): 
        df_m=df[df["Marker_nr"].isin([e])]
        df_m = df_m.set_index(list(df_m.columns[0:1]))    #"Sample")
        df_m = df_m.T
        df_m = df_m.iloc[:-1,:]
        mn_D_df=df_m[df_m.iloc[:,0].str.contains("D")]
        samples = mn_D_df.index.values.tolist() #samples that have this variant
        sample_list.append(samples) #str(samples))
        sample_list_str.append(str(samples))
        #count how many upstream markers
        Upstream_m =Tree_raw_hyp[Tree_raw_hyp["parted2"].isin([e])]
        Upstream_marker_list = len(Upstream_m["parted1"].values.tolist())
        upstream_counter.append(Upstream_marker_list)
    Ident_Output = pd.DataFrame()
    Ident_Output['marker']= identical
    Ident_Output['samples']=sample_list_str
    Ident_Output['upstream_nr']=upstream_counter

    #Get list of samples that have the markers
    sample_id = []
    for sample in sample_list:
        for element in sample:
            sample_id.append(element)
    sample_id = list(set(sample_id))

    Ident_Output = Ident_Output.sort_values(by=["samples"])
    Ident_Output = Ident_Output.reset_index(drop=True)

    equal_marker_names_list = []
    equal_marker_ds_list = []
    #filter for rows where each sample is included and remove all entries except the ones with the max number of upstream markers

    for s in sample_id:
        ds_per_sample_df=Ident_Output[Ident_Output["samples"].str.contains(s)]
        if len(ds_per_sample_df) > 1:
            upmax = ds_per_sample_df["upstream_nr"].max()
            upmax_filter = ds_per_sample_df[ds_per_sample_df["upstream_nr"] ==upmax]
            #remove the markers without downstream markers of one sample, if they have fewer upstream markers
            fail_upmax_filter = ds_per_sample_df[ds_per_sample_df["upstream_nr"] < upmax]
            fail_upmax = fail_upmax_filter.iloc[:,0].values.tolist()
            if len(upmax_filter) > 1:
                ident_markers = upmax_filter.iloc[:,0].values.tolist()
                ident_markers = '&'.join(ident_markers)
                true_equals.append(ident_markers)
            if len(fail_upmax) > 0:
                to_remove_2.append(fail_upmax)
                Equal_ds=Equal_ds[~Equal_ds["downstream"].isin(fail_upmax)] #remove the just removed branch end markers from the initial dataframe

    Equal_ds = Equal_ds.reset_index(drop=True)
    eq_m=itertools.groupby(Equal_ds["equal_marker"])
    eq_eq_2=[]
    dwn_list_2 = []
    for ds,key in eq_m:
        grouped = Equal_ds.groupby(['equal_marker'])
        df_eq = grouped.get_group(ds)
        dwn = df_eq.iloc[:,1].values.tolist()
        dwn_list_2.append(dwn)
        eq_eq_2.append(ds) 
    #which markers are now only downstream and never upstream?
    left_markers = Equal_ds["equal_marker"].values.tolist()
    right_markers = Equal_ds["downstream"].values.tolist()
    identical = list(set(right_markers)-set(left_markers))
    current_len = len(Equal_ds)
    if initial_len == current_len:
        break

#there are duplicate entries in true_equals:
true_equals_sorted = []
for trues in true_equals:
    true_1 = trues.split("&")
    true_1.sort()
    true_1 = '&'.join(true_1)
    true_equals_sorted.append(true_1)
true_equals_sorted = list(set(true_equals_sorted))

b = tuple(itertools.combinations(true_equals_sorted,2))

for elements in b:
    c1 = set(elements[0].split("&"))
    c2 = set(elements[1].split("&"))
    c = c1&c2
    if len(c) > 0:
        if len(c1) > len(c2):
            c2 = list(c2)
            c2.sort()
            c2 = '&'.join(c2)
            if c2 in true_equals_sorted:
                true_equals_sorted.remove(c2)
        else:
            c1 = list(c1)
            c1.sort()
            c1 = '&'.join(c1)
            if c1 in true_equals_sorted:
                true_equals_sorted.remove(c1)

to_remove_list = list(chain.from_iterable(to_remove_2))
to_remove_unique = list(set(to_remove_list))

<<<<<<< HEAD
if args.ambiguous_variants:
    path = args.output + "ambiguous_variants.csv"
=======
if type(args.ambiguous_variants) == str: #!= False:
    path = args.ambiguous_variants + "ambiguous_variants.csv"
>>>>>>> b36a4c3f45f06f30125b0c39c479333dc2f3cab7
    to_remove_df = pd.DataFrame()
    to_remove_df[0] = to_remove_unique
    to_remove_df[0] = to_remove_df[0].replace(marker_nr,marker_name)
    np.savetxt(path, to_remove_df, delimiter=",", fmt="%s", comments="")

#SCRIPT:Combining equal markers in list of markers
#1.remove the ambiguous markers from list

if len(to_remove_unique) != 0:
    Tree_raw_parted1=Tree_raw[~Tree_raw["parted2"].isin(to_remove_unique)] #(non_equal_markers[:])] #filter in hyphen list for entries from the first marker set in parted1 and afterwards in parted2 column
    Tree_raw=Tree_raw_parted1[~Tree_raw_parted1["parted1"].isin(to_remove_unique)]
#updated hyphenated tree list
Tree_raw_hyp=Tree_raw[Tree_raw["output"].str.contains("-")]
#turn hyphen separated columns into lists
Tree_raw_hyp_p1 = Tree_raw_hyp.iloc[:,1].values.tolist()
Tree_raw_hyp_p2 = Tree_raw_hyp.iloc[:,2].values.tolist()
#replace entries with true equal markers in these two lists

for equal in true_equals_sorted:
    single_markers_eq = equal.split("&")
    #turn list entries from strings to int
    single_markers_eq = [str(j) for j in single_markers_eq]
    for i in single_markers_eq:
        Tree_raw_hyp_p1 = [equal if item ==i else item for item in Tree_raw_hyp_p1]
        Tree_raw_hyp_p2 = [equal if item ==i else item for item in Tree_raw_hyp_p2]

Result = pd.DataFrame(list(zip(Tree_raw_hyp_p1,Tree_raw_hyp_p2)))
output=Result[0]+"-"+Result[1]
Result.insert(0,"output",output)
Result.rename(columns={0:"parted1",1:"parted2"},inplace=True)

#remove rows that are duplicate (because we have replaced the equal markers with same name, there should be some rows that have the identical information
Result = Result.drop_duplicates( keep='last')
Result = Result.reset_index(drop=True)
Result_ed = Result.copy() #new
print("Equal markers are combined!")

#edit Result table to make markers separable
Result_ed['parted1'] = '_' + Result_ed['parted1'].astype(str) + '_'
Result_ed['parted2'] = '_' + Result_ed['parted2'].astype(str) + '_'
Result_ed['output'] = '_' + Result_ed['output'].astype(str) + '_'
Result_ed = Result_ed.replace(r"&","_&_",regex=True)
Result_ed = Result_ed.replace(r"-","_-_",regex=True)

to_remove_3 = []

for sa in all_samples:
    df_T = df.T
    allele = df_T.loc[sa].values.tolist() 
    marker_allele = df_T.loc["Marker_nr"].values.tolist()
    marker_allele = [str(x) for x in marker_allele]
    current_df = pd.DataFrame()
    current_df["Marker_nr"]=marker_allele
    current_df["Allele"]=allele    
    mn_D_df=current_df[current_df["Allele"].str.contains("D")]    
    if len(mn_D_df) > 2:
        m_D = mn_D_df["Marker_nr"].values.tolist()

        #add characters around every marker name to make them easier to filter
        m_D = ["_" + suit + "_" for suit in m_D]

        Tree_hyp_rmv_2=Result_ed.copy()
        Tree_hyp_rmv_2=Tree_hyp_rmv_2[Tree_hyp_rmv_2.parted1.apply(lambda x: True if any(i in x for i in m_D) else False)]
        Tree_hyp_rmv_2=Tree_hyp_rmv_2[Tree_hyp_rmv_2.parted2.apply(lambda x: True if any(i in x for i in m_D) else False)]
  
        if len(Tree_hyp_rmv_2) > 1:
            Tree_hyp_rmv_parted2 = Tree_hyp_rmv_2["parted2"].values.tolist()
            Tree_hyp_rmv_parted1 = Tree_hyp_rmv_2["parted1"].values.tolist()
            Tree_hyp_rmv_parted2_1 = set(Tree_hyp_rmv_parted2) - set(Tree_hyp_rmv_parted1)
            Tree_hyp_rmv_parted2_1 = list(Tree_hyp_rmv_parted2_1)

            #filter for these only downstream markers as parted1 for all samples
            if len(Tree_hyp_rmv_parted2_1) > 1:
                Result_ed_fltr = Result_ed[Result_ed["parted1"].isin(Tree_hyp_rmv_parted2_1)]

                if len(Result_ed_fltr) > 0: #one of the markers has actually downstream markers
                    Result_ed_nodwn = Result_ed_fltr["parted1"].values.tolist()
                    removable = set(Tree_hyp_rmv_parted2_1)-set(Result_ed_nodwn)
                    removable = list(removable)
                    to_remove_3.append(removable)
                else: #if none of them have downstream markers, rmv randomly; this should not be an option, because dwnstream markers within one sample were already filtered out
                    to_remove_3.append(Tree_hyp_rmv_parted2_1[1:])


to_remove_3 = list(chain.from_iterable(to_remove_3))
to_remove_3 = list(set(to_remove_3))
to_remove_3 = [sub.replace('_', '') for sub in to_remove_3]

#filter markers out
Result=Result[~Result["parted1"].isin(to_remove_3)]
Result=Result[~Result["parted2"].isin(to_remove_3)]

print("Generate final pairwise marker relationships...")
#SCRIPT:generate final pairwise markers and phyloxml format of file
Result_up = [] #pd.DataFrame() #maybe change to DF afterwards
Result_down = []

Result = Result.sort_values(by=["parted1"]) #needs to be sorted, so the groupby works
Result = Result.reset_index(drop=True)
#groupby the upstream marker
Result["parted1"]=Result["parted1"].astype(str)

Tree_updated = Result.copy()
p1_group=itertools.groupby(Result["parted1"])
markers_with_depth = []

for key,group in p1_group:
    markers_with_depth.append(key)

index_list =["0"]
for up in index_list:
    Tree_0=Result[Result["parted1"].isin([up])] 
    Tree_0_ds = Tree_0['parted2'].values.tolist()   #list to compare: are downstream of 0, now: find out which are downstream of each other
    ##the updated Tree which will be used continuously
    Tree_updated=Tree_updated[~Tree_updated["parted1"].isin([up])] 
    Tree_upd_ds = Tree_updated['parted2'].values.tolist()     #second list to compare
    next_layer = (list(set(Tree_0_ds)-set(Tree_upd_ds)))    #get what is unique to the previous list, Tree_0_ds compared to Tree_upd_ds
    ##the most immediate downstream markers belonging to the up marker
    Result_up.append(up)
    Result_down.append(next_layer)
    for down in next_layer:
        if down in markers_with_depth:
            index_list.append(down)

Result = pd.DataFrame(list(zip(Result_up, Result_down)))
Final = "("+str(Result.loc[0,0])+",("+str(Result.loc[0,1])+")"
flatten_list = Result.iloc[0,1]
for i in range(len(Result)): #we start with one because we skip the base of the tree (see "Final")
    Result_0=Result[Result[0].isin(flatten_list)]
    Result_0 = Result_0.reset_index(drop=True)
    if len(Result_0) == 0:
        break
    for one_layer in range(len(Result_0)):
        upstream = Result_0[0]
        downstream = Result_0[1]
        to_replace = "'"+upstream[one_layer]+"'"
        replace_with = "'"+upstream[one_layer]+"',("+str(downstream[one_layer])+")"
        Final = re.sub(to_replace,replace_with,Final)
    k = downstream.values.tolist()
    flatten_list = [j for sub in k for j in sub]

Final=Final.replace(',([])', '') 
Final=Final.replace(',(', ',(\n') 
Final=Final.replace('),', '),\n')
Final=Final.replace(' ', '') 
Final=Final.replace(',(', '(') 
Final=Final.replace('),', ')') 
Final=Final.replace(',', '\n') 
Final=Final.replace('[', '') 
Final=Final.replace(']', '') 
Final=Final.replace("'", "") 

 
path = args.output + "temporary_Output_tree.txt"
text_file = open(path,"w")
text_file.write(Final)
text_file.close()

df_final = pd.read_csv(path, sep=";", header=None)
os.remove(path)



df_raw = df_final.iloc[:,0].values.tolist()
location = [0,1]
#ASSUMPTIONS: there cannot be '(' and ')' in one entry
xaxis= 1
for elements in df_raw[1:]:
    count_dwn =elements.find("(")
    count_up =elements.find(")")
    closed_bracket_counter = elements.count(")")
    if count_dwn > 0:
        xaxis = xaxis + 1
        location.append(xaxis)
    elif count_up > 0:
        xaxis = xaxis - closed_bracket_counter
        location.append(xaxis)
    else: #if there is no parantheses, the position on x axis stays the same
        location.append(xaxis)

#df['delta_xaxis'] = parentheses_cnt[:-1] #Don't save the last entry
df_final['xaxis'] = location[:-1] 
tabs_xml = location[1:-1] #Save this for generating phylxml file


#move "delta_xaxis" column as first column
df_final = df_final[['xaxis',0]]
#give highest x axis value, i.e. most resolved SNP and add as many new columns
max_xaxis = df_final['xaxis'].max()

column_nr =int(int(max_xaxis)+1)
for e in range(max_xaxis):
    #add empty column 
    column = "layer"+str(e)
    df_final[column] = df_final[0]

#replace brackets in dataframe
df_final = df_final.replace(r"\(","",regex=True)
df_final = df_final.replace(r"\)","",regex=True)


#move deeper located SNPs to the other columns
cnt=1
for depth in range(max_xaxis+1):
    for index,row in df_final.iterrows():
        c_name = df_final.columns[depth+1]
        if row['xaxis'] != depth:
            df_final.at[index,c_name]="" # #np.nan #instead of 1 put name of columns
    cnt = cnt +1
            

index_counter = 0
cnt=0
for index,row in df_final.iterrows():
  #  print(row) #gives for this row, for each column the column name and element
    column_counter = 0 #because we skip the xaxis column
    for element in row[1:]: #from 1 on, so we skip the xaxis column
        index = index_counter
        column_counter = column_counter + 1
        if element.find(",") != -1: #meaning found a comma/parallel SNPs at the element of this index and column
            #transform element in row[1:] into a list separated by commas, and insert them element by element in the current cell
            parallel = element.split(",") #turn string with parallel SNPs into list of parallel SNPs
            #add empty row to dataframe depending on len of parallel SNPs list
            for parallele_elements in range(len(parallel)): #-1): #CHANGE
                index = index_counter
                added_row = [df_final.iloc[index,0]]+[""]*(column_nr)
                #remove element in list at certain position and insert parallel snp in list added_row at that same position
                added_row.pop(column_counter)
                added_row.insert(column_counter,parallel[parallele_elements])
                df_final.loc[index+0.5]=added_row #need to add index+0.5, so it is in the row after the current index, and we don't overwrite current rows 
                df_final = df_final.sort_index().reset_index(drop=True)
                index_counter = index_counter + 1 #important to keep up with the inserted columns and shift a
    cnt=cnt+1
    index_counter = index_counter + 1
    column_counter = column_counter + 1

#Remove rows that contain commas
m = 1 #not zero, we skip the xaxis column

df_final = df_final.sort_index().reset_index(drop=True)
#remove x axis column:
df_final = df_final.iloc[:,1:]


#split equal markers in all columns:
for e in range(column_counter-1):
    df_final.iloc[:,e] = df_final.iloc[:,e].str.split("&")

Final = pd.DataFrame()


#looping through the columns:
for e in range(column_counter-1):
    Rows_final = [] #reset for new column
    for m in range(len(df_final)):
        liste=df_final.iloc[m,e] #
        if liste[0] != "":#df.iloc[m,e]!="[]":
            update_df_2 = pd.Series(df_final.iloc[m,e])
            update_df_2 = update_df_2.replace(marker_nr,marker_name)
            update_df_2 = update_df_2.tolist()
            Rows_final.append(update_df_2)
        else:
            Rows_final.append('['']')
    Final[e] = Rows_final


Final = Final.replace(r"\[","",regex=True)
Final = Final.replace(r"\]","",regex=True)

#combine all entries in each row into a single entry of the dataframe
label_variants = Final.apply(lambda row: ''.join(map(str,row)), axis=1)
label_variants_list = label_variants.tolist()
coordinate_xml = pd.DataFrame({"1":tabs_xml, "2":label_variants_list[1:]})
#remove ' and [, ] from variant names
last_column = coordinate_xml.columns[-1]
coordinate_xml[last_column] = coordinate_xml[last_column].str.replace("'", "").str.replace("[", "").str.replace("]", "")


#create the phylxml file for the tree
path = args.output + "Output_tree.xml"
open_clade_counter = 1
close_clade_counter = 0

with open(path, 'w') as file:
    file.write('<phyloxml xmlns="http://www.phyloxml.org" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd">\n<phylogeny rooted="true">\n<clade>\n<color>\n<red>200</red>\n<green>20</green>\n<blue>51</blue>\n</color>')
    tree_layers = len(coordinate_xml)
    end_clade = '\n</clade>'

    for i in range(tree_layers):
        current_variant = coordinate_xml.iloc[i,1]
        insert_text = '\n<clade>\n<name>'+current_variant+'</name>'
        current_position = coordinate_xml.iloc[i,0]
        previous_position = coordinate_xml.iloc[i-1,0]
        if i == 0:
            file.write(insert_text)
            open_clade_counter += 1
        else:
            if current_position > previous_position:
                file.write(insert_text)
                open_clade_counter += 1
            if current_position == previous_position:
                file.write(end_clade)
                close_clade_counter += 1
                file.write(insert_text)
                open_clade_counter += 1
            if current_position < previous_position:
                delta = previous_position - current_position + 1
                for m in range(delta):
                    file.write(end_clade)
                    close_clade_counter += 1
                file.write(insert_text)
                open_clade_counter += 1
    end_close_clade = open_clade_counter - close_clade_counter
    for n in range(end_close_clade):
        file.write(end_clade)
    file.write('\n</phylogeny>\n</phyloxml>')

print()
print("The Phyloxml tree file was successfully generated in "+str(path)+" !")

#CSV file
path = args.output + "Output_tree.csv"
np.savetxt(path, Final, delimiter="\t", fmt="%s", comments="")
print()
print("The csv tree file was successfully generated in "+str(path)+" !")

##To account for uncertainty Part 2:

#1: how many variants are contained in the final tree:
if len(to_remove) != 0:
    Pairwise_complete = Pairwise_complete[~Pairwise_complete['M1'].isin(to_remove) & ~Pairwise_complete['M2'].isin(to_remove)]
if len(to_remove_unique) != 0:
    Pairwise_complete = Pairwise_complete[~Pairwise_complete['M1'].isin(to_remove_unique) & ~Pairwise_complete['M2'].isin(to_remove_unique)]
Pairwise_complete = Pairwise_complete.replace(marker_nr,marker_name) #replace marker number by marker names in dataframe
# print(Pairwise_complete)

#put all markers of both columns in the dataframe into a list, combine both lists and keep unique entries
strings_from_column1 = Pairwise_complete["M1"].tolist()
strings_from_column2 = Pairwise_complete["M2"].tolist()
all_strings = set(strings_from_column1 + strings_from_column2)
all_strings = list(all_strings)
all_strings.remove("Root") #Remove root entry
all_strings_counter = len(all_strings)-1 #so many pairwise relationships to OTHER variants are possible for each variant
# print(all_strings_counter)
#loop through dataframe of all pairwise variants to count how many relationships were reported to other variants
# print(all_strings)


#2: get informative (=downstream,upstream, parallel) pairwise relationships 

if len(to_remove) != 0:
    Pairwise_complete_filter = Pairwise_complete_filter[~Pairwise_complete_filter['M1'].isin(to_remove) & ~Pairwise_complete_filter['M2'].isin(to_remove)]
if len(to_remove_unique) != 0:
    Pairwise_complete_filter = Pairwise_complete_filter[~Pairwise_complete_filter['M1'].isin(to_remove_unique) & ~Pairwise_complete_filter['M2'].isin(to_remove_unique)]
Pairwise_complete_filter = Pairwise_complete_filter.replace(marker_nr,marker_name) #replace marker number by marker names in dataframe
Pairwise_complete_filter = Pairwise_complete_filter[~Pairwise_complete_filter['M1'].isin(["Root"])]


#3: get uncertainty
certainty = []


for var in all_strings:
    #when filtering using .str.contains, var needs to be escaped, for cases where it contains +,;; or similar in the string name. Escaping however introduces an error when filtering the columns, when there are other special characters, e.g. for Q1b1a1a1i2~_Z35616
    # var= re.escape(var)
    # filtered_df = Pairwise_complete_filter[(Pairwise_complete_filter['M1'].str.contains(var)) | (Pairwise_complete_filter['M2'].str.contains(var))]
    filtered_df = Pairwise_complete_filter[(Pairwise_complete_filter['M1'] == var) | (Pairwise_complete_filter['M2'] == var)]
    current_string_len = len(filtered_df)
    #get current certainty value
    current_certainty = current_string_len/all_strings_counter
    current_certainty_exp = "{:e}".format(current_certainty)
    certainty = certainty + [current_certainty_exp]
    

Result = pd.DataFrame()
Result["variables"] = all_strings
Result["certainty values"]=certainty 
Result.columns = ["variables","certainty values"]
path_stat = args.output + "certainty_values.csv"
Result.to_csv(path_stat, index = False, sep="\t")
# np.savetxt(path, Result, delimiter="\t", fmt="%s", comments="")
print("\n"+"The certainty values for variants in the phylogenetic tree were successfully generated in "+str(path_stat)+" !")

##

#METADATA
<<<<<<< HEAD
if args.metadata_individuals:
=======
if type(args.metadata_individuals) == str: # != False:
>>>>>>> b36a4c3f45f06f30125b0c39c479333dc2f3cab7
    with open(path, "r") as file:
    # with open("/mnt/ngs/scratch_areas/nxd426/1.5_SNPtotree/output5_newTree_Testdata1.txt", "r") as file:
        newText=file.read()
        newText=newText.replace("['", '') 
        newText=newText.replace("']", '') 
        newText=newText.replace("', '", ' , ')
        newText=newText.replace(']', '') 
    with open(path, "w") as file:
        file.write(newText)

    df_tree = pd.read_csv(path, sep="\t")
    # df_tree = df_tree.replace(r"\[\'","",regex=True) #takes anything within parentheses --> bad quality results are also just "X"
    # df_tree = df_tree.replace(r"\'\]","",regex=True) #takes anything within parentheses --> bad quality results are also just "X"
    # df_tree = df_tree.replace(r"\', '","&",regex=True) #takes anything within parentheses --> bad quality results are also just "X"
    # df_tree = df_tree.replace(r"\]","",regex=True) #takes anything within parentheses --> bad quality results are also just "X"

    marker_selection = []
    marker_grouping = []
    df_tree_T = df_tree.T
    col = len(df_tree_T)
    

    #Get marker list in the tree from top to bottom
    for index,row in df_tree.iterrows():
        for element in df_tree.iloc[index,1:col]:
            if type(element) == str:
            # if element != "":
                marker_grouping.append(element)
                element = element.split(" , ")
                # for e in element:
                    # marker_selection.append(e)
    # sel_marker_counter = len(marker_selection)

    df = df.iloc[:,:-1]

    marker_grouping_samples = []
    for m in marker_grouping:
        m = m.split(" , ")
        entries = len(m)
        df_m=df[df.iloc[:,0].isin(m)]
        df_m = df_m.set_index(list(df_m.columns[0:1]))
        df_m = df_m.T
        m = 0
        samples_per_layer = []
        for e in range(entries):
            mn_D_df=df_m[df_m.iloc[:,m].str.contains("D")]
            samples = mn_D_df.index.values.tolist()
            samples_per_layer.append(samples)
            m = m+1
        samples_per_layer = list(chain.from_iterable(samples_per_layer))
        samples_per_layer = list(set(samples_per_layer))
        marker_grouping_samples.append(samples_per_layer)

    Result = pd.DataFrame()
    Result["marker"] = marker_grouping
    Result["samples"]=marker_grouping_samples  
    path = args.output + "Variant_individuals_tree.csv"
    np.savetxt(path, Result, delimiter="\t", fmt="%s", comments="")
    print("\n"+"The metadata file was successfully generated in "+str(path)+" !")

if len(non_poly) > 0:
    print("\n"+"USER MESSAGE: Your input file contained one or more variants that were not polymorphic. For the analysis we removed variant(s): "+str(non_poly_list)+" !")
#Give execution time
end = time.time()
total_time = end - start
print("\n"+"This run took a total of "+ str(total_time)+" sec!")
