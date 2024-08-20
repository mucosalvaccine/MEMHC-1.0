#Creator; Mohammad Arabpour, mohammad.arabpour.sinior@gmail.com
## License

This software is licensed under the [CC BY-NC 4.0 License](LICENSE) for non-commercial use.

For commercial use, please contact Mohammad Arabpour at mohammad.arabpour.sinior@gmail.com to obtain a commercial license.


# minimal-epitope-for-maximum-MHC-coverage
The code  provided is a Python script that processes a CSV file containing MHC (Major Histocompatibility Complex) ranking data for a set of peptides. It performs various operations such as filtering, binary conversion of Kd scores, calculating coverage, and generating output files.
Using this code you would be able to use a csv file including columnes with 'id' = prptide ID that could be peptide index and or any attribution , 'peptide' = peptide sequence in AA,'peptide_length'= lenght of peptide, 'start'= starting amino acid position in protein,  the rest f columns are the MHC alleles headings for correspondent peptides and kd(nM) score for each peptide.

The output is a ranking file that give the minimal epitopes(peptides) needed to cover maximum HLA alleles provided.
the applications encompases but not limited to  any immunoassay and vaccine developments need the minimal epitope with maximal HLA coverage.

Here is an example of input csv file :

A part of  protein inflrunza A/Sendai/TU66/2008(H1N1)):   mkvkllvllctft
peptide sequences 8-11 mers
HLA type I supertype binding affinity Kd(nM): HLA-A01:01,HLA-B15:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01

id	peptide	peptide_length	start	HL+D1:BM18A-A01:01	HLA-A02:01	HLA-A03:01	HLA-A24:02	HLA-A26:01	HLA-B07:02	HLA-B08:01	HLA-B15:01	HLA-B27:05	HLA-B39:01	HLA-B40:01	HLA-B58:01
Sequence1	MKVKLLVLLC	10	0	37290.2266	24836.6055	29509.4492	38739.0938	40020.6367	36766.1914	17256.5215	24450.7969	12728.3242	15750.3096	31257.5508	25245.1699
.
.
Sequence14	LVLLCTFT	8	5	36421.332	19980.8965	33236.9688	44220.8711	45617.0742	38217.8594	26305.9082	33364.5117	38968.2148	42882.7188	44307.5664	34264.1367

the output csv  file for above would  like :
peptide sequence	# HLA allel hits	accumulativcoverage(%)	Peptide lenght
MKVKLLVLL	1	33.33333333	9
MKVKLLVL	1	66.66666667	8
KLLVLLCTF	1	100	9


import matplotlib
import pandas as pd
import numpy as np

df1 = pd.read_csv("MHC kd input .csv")

# Group rows by the same value in the 'start' column
grouped = df1.groupby(['start'])

# Create an empty dataframe to store the filtered values
filtered_df = pd.DataFrame(columns=df1.columns)

# Loop through each group
for group_name, group_data in grouped:
    # Get the minimum value of each column (excluding peptide_length, start, and id)
    min_values = group_data.drop(['peptide_length', 'start', 'id', 'peptide'], axis=1).min()
    # Create a new dataframe with the same columns as group_data
    filtered_group = pd.DataFrame(columns=df1.columns)
    # Copy the peptide_length, start, and id columns to the filtered_group dataframe
    filtered_group['peptide_length'] = group_data['peptide_length']
    filtered_group['start'] = group_data['start']
    filtered_group['id'] = group_data['id']
    filtered_group['peptide'] = group_data['peptide']
    # Set values greater than the minimum to 1000
    for col in min_values.index:
        filtered_group[col] = group_data[col].apply(lambda x: x if x == min_values[col] else 1000)
    # Append the filtered_group dataframe to the filtered_df dataframe
    filtered_df = pd.concat([filtered_df, filtered_group], ignore_index=True)
df1 = filtered_df
""""
# format I for going through individual Kd scores for binary conversion using a for loop(time consuming)
#defining the Kd_treshold for peptides in nM
Kd_treshold_nM = 500

for i in range(0, (df1.shape[0])):
    for j in range(0, (df1.shape[1])):
        p = df1.iloc[i][j]
        if type(p) !=  str:
            if int(p) <= Kd_treshold_nM:
                p = 1
            elif int(p) > Kd_treshold_nM:
                p = 0
        elif type(p) == str:
            p = p
        df1.iloc[i,j] = p
"""
Kd_treshold_nM = 500
# format II for going through individual Kd scores for binary conversion using np method(time efficient)
exclude_cols = ['peptide_length', 'start', 'id', 'peptide']  # list of column names to exclude
numeric_cols = df1.columns.drop(exclude_cols)  # list of numeric column names

# convert numeric columns to numeric dtype
df1[numeric_cols] = df1[numeric_cols].apply(pd.to_numeric, errors="coerce")
# Convert Kd_treshold_nM to 1 and 0, excluding NaN values

# apply the condition to the binary numeric cells where it is  a number
#df1[numeric_cols] = np.where(df1[numeric_cols] <= Kd_treshold_nM, 1, 0)
#df1[numeric_cols] = np.where(np.isan(df1[numeric_cols]) & np.where(df1[numeric_cols] <= Kd_treshold_nM, 1, 0))
df1[numeric_cols] = np.where(~np.isnan(df1[numeric_cols]), np.where(df1[numeric_cols] <= Kd_treshold_nM, 1, 0), np.nan)
print(df1)  # print the updated dataframe
# Initialize the total coverage to 0
relative_HLA_coverage = 0
# Loop through each column (excluding peptide_length, start, id, and peptide)
for col in df1.columns[4:]:
    # Sum the values of the column and set the sum value to 1 if the sum is greater than 0, and 0 otherwise
    sum_value = df1[col].sum()
    if sum_value > 0:
        sum_value = 1
    else:
        sum_value = 0

    # Add the sum_value to the total_coverage
    relative_HLA_coverage += sum_value
print(relative_HLA_coverage)

df1.to_csv('filtered MHC output file.csv', index=False)
print(df1)


#make the sume coverage value column for each peptide
df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)
#print(df1['sum_weight'])
df1.to_csv('result_table I.csv', index=False)
total1 = (df1.shape[1])-1
#print(total1)
df1 = df1.sort_values(by='sum_weight',axis=0, ascending =False, inplace=False, kind='quicksort', na_position='last', ignore_index=False, key=None)
maximum1 = max(df1['sum_weight'])
lead_peptide = (df1.iloc[0,1])
Lead_p = len(df1.iloc[0,1])

leadcoverage = (maximum1 / relative_HLA_coverage) * 100
row_number = df1[df1['sum_weight'] == maximum1].index
df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)
# print ('totalscore',total1)
# print('lead pep score',df1.iloc[0]['sum_weight'])
# print('lead pep',lead_peptide)
# print('leadcoverage',leadcoverage)
# print('rowindex',row_number)
# print(df1)
totalcoverage = leadcoverage
total = total1
df1.to_csv('result_table II.csv', index=False)
max_num = maximum1
#print(maximum1)
#####going to second line
#####------------------------------------------------------------------------------------------------------------------
totalscore =[total1]
#lead_pep_score = [df1.iloc[0]['sum_weight']]
lead_pep_score = [max_num]
lead_pep = [lead_peptide]
Lead_peptide_lenght = [int(Lead_p)]
rowindex = [row_number]
accuulativcoverage = [leadcoverage]
while total > 0:
    df= df1.iloc[0]
    s = ((df.index))
    g = []
    for i in range (0, len(s)):
        h = df.iloc[i]
        if type(h) != str:
            if h == 0:
             g = [s[i]] + g
    ali = ['peptide'] + g
    df1 = df1[ali]
    df1.shape[0]

    df1 = df1[1:(df1.shape[0])]
    df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)
    total = df1['sum_weight'].sum()
    #print(total)
    df1 = df1.sort_values(by='sum_weight',axis=0, ascending =False, inplace=False, kind='quicksort', na_position='last', ignore_index=False, key=None)
    maximum = max(df1['sum_weight'])
    max_num += maximum
    lead_peptide = (df1.iloc[0,0])
    Lead_peptide_lenght += [len(df1.iloc[0, 0])]
    leadcoverage = [(maximum / total1) * 100]
    #print (max_num)
    row_number = df1[df1['sum_weight'] == maximum].index
    #calculate accumulative coverage from all the HLA allotypes (including those that are not recognized by any peptide in peptide pool)
    absolut_acumulativelcoverage = (max_num/total1) * 100
    #calculate acumulative coverage from all the recognized hits of HLA allotypes (excluding those that are not recognized by any peptide in peptide pool)
    acumulativelcoverage = (max_num/relative_HLA_coverage) * 100
    #print(max_num)
    # print('rowindex', row_number[0])
    # print('lead pep score', df1.iloc[0]['sum_weight'])
    # print('lead pep', lead_peptide)
    # print('leadcoverage', leadcoverage)
    # print ('accuulativcoverage', acumulativelcoverage)
    # print(df1)
    accuulativcoverage += [acumulativelcoverage]
    lead_pep_score += [df1.iloc[0]['sum_weight']]

    #print(lead_pep_score)
    lead_pep += [lead_peptide]
    rowindex += [(row_number[0])+1]


print(Lead_peptide_lenght)
data ={'peptide sequence':lead_pep,'# HLA allel hits': lead_pep_score, 'accumulativcoverage(%)':accuulativcoverage ,"Peptide lenght":Lead_peptide_lenght}
df3 = pd.DataFrame(data)

df3.to_csv('minimal-epitop-max_coverage_MHC.csv', index=False)
