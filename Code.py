import requests
from rdkit import Chem
import pandas as pd
from molmass import Formula
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

########################################################

# Specify all file paths involved in the program. Please change the following to your actual paths.
# Metabolite namelist path
METABOLITE_NAMELIST_PATH="C:/Users/celis/Downloads/29387495/metabolites.csv"
# MS data file path
MS_DATA_PATH="C:/Users/celis/Downloads/29387495/MS_data.csv"
# Output file folder path
OUT_PATH="C:/Users/celis/Downloads/29387495/"
# Standard background (non-related) peak file path
NR_PATH="C:/Users/celis/Downloads/29387495/NR.csv"

# Specify parameters used in the program
# Set the minimum intensity threshold for peak identification
int_min_thred=0
# Set the maximum tolerated m/z difference for peak annotation
mzdiff=0.05

########################################################

def name_to_smiles_pubchem(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        smiles = response.json()['PropertyTable']['Properties'][0]['ConnectivitySMILES']
        return smiles
    else:
        return None

# Read names of metabolites
data=pd.read_csv(METABOLITE_NAMELIST_PATH,header=None, encoding='utf-8-sig')
data.columns=["Name"]
for i in range(len(data)):
    data.loc[i,'Name']=(data.loc[i,'Name']).lower()

print("Original molecule names loaded")   

# Write SMILES and molecular weight
data.insert(data.shape[1],'SMILES',None)
data.insert(data.shape[1],'MW',None)

for i in range(len(data)):
    compound_name=data.loc[i,'Name']
    smiles = name_to_smiles_pubchem(compound_name)
    if smiles:
        data.loc[i,'SMILES']=smiles
    else:
        data.loc[i,'SMILES']="#"

# Find properties
data.insert(data.shape[1],'Formula',None)
data.insert(data.shape[1],'Properties',None)
amine_pattern = Chem.MolFromSmarts("[NH2]")
COOH_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
CHO_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
ammoniumsalt_pattern = Chem.MolFromSmarts("[N+]")
OH_pattern = Chem.MolFromSmarts("[#6][OX2H]")
alkene_pattern = Chem.MolFromSmarts("[$([CX3]=[CX3])]")
phenol_pattern = Chem.MolFromSmarts("[OX2H][cX3]:[c]")

amine_list = pd.DataFrame(columns=['Name','MW'])
amine_index = 0
ketone_list = pd.DataFrame(columns=['Name','MW'])
ketone_index = 0
COOH_list = pd.DataFrame(columns=['Name','MW'])
COOH_index = 0
ammoniumsalt_list = pd.DataFrame(columns=['Name','MW'])
ammoniumsalt_index = 0
OH_list = pd.DataFrame(columns=['Name','MW'])
OH_index = 0
alkene_list = pd.DataFrame(columns=['Name','MW'])
alkene_index = 0
phenol_list = pd.DataFrame(columns=['Name','MW'])
phenol_index = 0

for i in range(len(data)):
    data.loc[i,'Properties']=""

BasicProperties="H-Na-K-Z-N-"

print("Masses and SMILES codes of original molecules generated")     

for i in range(len(data)):
    smiles=data.loc[i,'SMILES']
    if smiles!="#":
        mol = Chem.MolFromSmiles(smiles)
#        data.loc[i,'MW']=round(MolWt(mol),3)
        ChemF = Formula(CalcMolFormula(mol))
        data.loc[i,'Formula'] = CalcMolFormula(mol)
        data.loc[i,'MW'] = round(ChemF.monoisotopic_mass,3)
        data.loc[i,'Properties']=data.loc[i,'Properties']+BasicProperties
        if mol and mol.HasSubstructMatch(amine_pattern):
            data.loc[i,'Properties']=data.loc[i,'Properties']+"A-"
            amine_list.loc[amine_index,'Name']=data.loc[i,'Name']
            amine_list.loc[amine_index,'MW']=data.loc[i,'MW']
            amine_index += 1
        if mol and mol.HasSubstructMatch(COOH_pattern):
            data.loc[i,'Properties']=data.loc[i,'Properties']+"C-"
            COOH_list.loc[COOH_index,'Name']=data.loc[i,'Name']
            COOH_list.loc[COOH_index,'MW']=data.loc[i,'MW']
            COOH_index += 1
        if mol and mol.HasSubstructMatch(ketone_pattern):
            data.loc[i,'Properties']=data.loc[i,'Properties']+"Ke-"
            ketone_list.loc[ketone_index,'Name']=data.loc[i,'Name']
            ketone_list.loc[ketone_index,'MW']=data.loc[i,'MW']
            ketone_index += 1
        if mol and mol.HasSubstructMatch(ammoniumsalt_pattern):
            data.loc[i,'Properties']=data.loc[i,'Properties']+"N+-"
            ammoniumsalt_list.loc[ammoniumsalt_index,'Name']=data.loc[i,'Name']
            ammoniumsalt_list.loc[ammoniumsalt_index,'MW']=data.loc[i,'MW']
            ammoniumsalt_index += 1
        if mol and mol.HasSubstructMatch(OH_pattern):
            data.loc[i,'Properties']=data.loc[i,'Properties']+"OH-"
            OH_list.loc[OH_index,'Name']=data.loc[i,'Name']
            OH_list.loc[OH_index,'MW']=data.loc[i,'MW']
            OH_index += 1
        if mol and mol.HasSubstructMatch(alkene_pattern):
            data.loc[i,'Properties']=data.loc[i,'Properties']+"C=C-"
            alkene_list.loc[alkene_index,'Name']=data.loc[i,'Name']
            alkene_list.loc[alkene_index,'MW']=data.loc[i,'MW']
            alkene_index += 1
        if mol and mol.HasSubstructMatch(phenol_pattern):
            data.loc[i,'Properties']=data.loc[i,'Properties']+"ArOH-"
            alkene_list.loc[alkene_index,'Name']=data.loc[i,'Name']
            alkene_list.loc[alkene_index,'MW']=data.loc[i,'MW']
            alkene_index += 1

print("Original molecule information written to file")
data_index=len(data)

# Write radical products
for i in range(len(data)):
    if data.loc[i,'MW'] != None:
        data.loc[data_index,'Name']="adduct of " + data.loc[i,'Name'] + " with OH radical"
        data.loc[data_index,'MW']=round(data.loc[i,'MW'] + 15.995 + 1.008, 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1

        data.loc[data_index,'Name']="adduct of " + data.loc[i,'Name'] + " with NO radical"
        data.loc[data_index,'MW']=round(data.loc[i,'MW'] + 14.003 + 15.995, 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1  

        data.loc[data_index,'Name']="adduct of " + data.loc[i,'Name'] + " with Cl radical"
        data.loc[data_index,'MW']=round(data.loc[i,'MW'] + 34.969, 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1  
        
        data.loc[data_index,'Name']="adduct of " + data.loc[i,'Name'] + " with CH3 radical"
        data.loc[data_index,'MW']=round(data.loc[i,'MW'] + 12.000 + 1.008 * 3, 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1 

        data.loc[data_index,'Name']="adduct of " + data.loc[i,'Name'] + " with CH3O radical"
        data.loc[data_index,'MW']=round(data.loc[i,'MW'] + 12.000 + 1.008 * 3 + 15.995, 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1  
print("Adduct with radicals generated")

# Write reductive hydrogenation products
for j in range(len(ketone_list)):
    data.loc[data_index,'Name']=ketone_list.loc[j,'Name']+" reductive hydrogenation product"
    data.loc[data_index,'MW']=round(ketone_list.loc[j,'MW'] + 1.008 * 2, 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1
print("Reductive hydrogenation products generated")

# Write reductive amination products
for j in range(len(ketone_list)):
    for i in range(len(amine_list)):
        data.loc[data_index,'Name']="reductive amination product of " + amine_list.loc[i,'Name'] + " and " + ketone_list.loc[j,'Name']
        data.loc[data_index,'MW']=round(amine_list.loc[i,'MW'] + ketone_list.loc[j,'MW'] - 15.995, 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1
print("Reductive amination products generated")

# Write -H+2Na products
for i in range(len(COOH_list)):
    data.loc[data_index,'Name']=COOH_list.loc[i,'Name'] + " -H +2Na product"
    data.loc[data_index,'MW']=round(COOH_list.loc[i,'MW'] + 21.982 , 3)
    data.loc[data_index,'Properties']="Z-"
    data_index += 1
print("-H+2Na products generated")

# Write decarboxylation products
for i in range(len(COOH_list)):
    data.loc[data_index,'Name']="decarboxylation product of " + COOH_list.loc[i,'Name']
    data.loc[data_index,'MW']=round(COOH_list.loc[i,'MW'] - 43.990 , 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1
print("Decarboxylation products generated")

# Write decarboxylation coupling (R-R') products
for i in range(len(COOH_list)):
    for j in range(i,len(COOH_list)):
        data.loc[data_index,'Name']="decarboxylation coupling (R-R') product of " + COOH_list.loc[i,'Name'] + " and " + COOH_list.loc[j,'Name']
        data.loc[data_index,'MW']=round(COOH_list.loc[i,'MW'] + COOH_list.loc[j,'MW'] - 43.990 * 2 - 2 * 1.008, 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1
print("Decarboxylation coupling (R-R') products generated")

# Write decarboxylation coupling (RCOOR') products
for i in range(len(COOH_list)):
    for j in range(i,len(COOH_list)):
        data.loc[data_index,'Name']="decarboxylation coupling (RCOOR') product of " + COOH_list.loc[i,'Name'] + " and " + COOH_list.loc[j,'Name']
        data.loc[data_index,'MW']=round(COOH_list.loc[i,'MW'] + COOH_list.loc[j,'MW'] - 43.990 - 2 * 1.008, 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1
print("Decarboxylation coupling (RCOOR') products generated")

# Write decarboxylation C-N coupling products            
for i in range(len(COOH_list)):
    for j in range(len(amine_list)):
        data.loc[data_index,'Name']="decarboxylation C-N coupling product of " + COOH_list.loc[i,'Name'] + " and " + amine_list.loc[j,'Name']
        data.loc[data_index,'MW']=round(COOH_list.loc[i,'MW'] + amine_list.loc[j,'MW'] - 43.990 - 2 * 1.008 , 3)
        data.loc[data_index,'Properties']=BasicProperties
        data_index += 1
print("Decarboxylation C-N coupling products generated")

# Write -TMA products
for i in range(len(ammoniumsalt_list)):
    data.loc[data_index,'Name']=ammoniumsalt_list.loc[i,'Name'] + " de-trimethylamine product"
    data.loc[data_index,'MW']=round(ammoniumsalt_list.loc[i,'MW'] - 59.073, 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1   
print("De-trimethylamine product generated")

# Write de-NH2 products
for i in range(len(amine_list)):
    data.loc[data_index,'Name']=amine_list.loc[i,'Name'] + " de-NH2 product"
    data.loc[data_index,'MW']=round(amine_list.loc[i,'MW'] - 16.019, 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1   
print("De-NH2 product generated")

# Write dehydroxyl products
for i in range(len(OH_list)):
    data.loc[data_index,'Name']=OH_list.loc[i,'Name'] + " de-OH product"
    data.loc[data_index,'MW']=round(OH_list.loc[i,'MW'] - 17.003, 3)
    data.loc[data_index,'Properties']="Z-"
    data_index += 1   
print("De-OH product generated")

# Write phenol dehydrogen products
for i in range(len(phenol_list)):
    data.loc[data_index,'Name']=phenol_list.loc[i,'Name'] + " phenol dehydrogenation product"
    data.loc[data_index,'MW']=round(phenol_list.loc[i,'MW'] - 2.016, 3)
    data.loc[data_index,'Properties']="Z-"
    data_index += 1   
print("Phenol oxidation(dehydrogenation) product generated")

# Write alkene oxidation and add H2O products
for i in range(len(alkene_list)):
    data.loc[data_index,'Name']=alkene_list.loc[i,'Name'] + " add 1 oxygen"
    data.loc[data_index,'MW']=round(alkene_list.loc[i,'MW'] + 15.995, 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1   
    
    data.loc[data_index,'Name']=alkene_list.loc[i,'Name'] + " add 2 oxygens"
    data.loc[data_index,'MW']=round(alkene_list.loc[i,'MW'] + 15.995 * 2, 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1  

    data.loc[data_index,'Name']=alkene_list.loc[i,'Name'] + " add H2O"
    data.loc[data_index,'MW']=round(alkene_list.loc[i,'MW'] + 18.011, 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1   
print("Addition of O/2O/H2O products of alkene generated")

# Write CH3CN products (with CH3CO+)
for i in range(len(OH_list)):
    data.loc[data_index,'Name']=OH_list.loc[i,'Name'] + " add CH3CO+ (from CH3CN) product"
    data.loc[data_index,'MW']=round(OH_list.loc[i,'MW'] + 43.0446, 3)
    data.loc[data_index,'Properties']="Z-"
    data_index += 1 
for i in range(len(amine_list)):
    data.loc[data_index,'Name']=amine_list.loc[i,'Name'] + " add CH3CO+ (from CH3CN) product"
    data.loc[data_index,'MW']=round(amine_list.loc[i,'MW'] + 43.0446, 3)
    data.loc[data_index,'Properties']="Z-"
    data_index += 1
print("Reactions with CH3CO+ (from CH3CN) products generated")

# Write CH3CN products (with CH2CN-)
for i in range(len(ketone_list)):
    data.loc[data_index,'Name']=ketone_list.loc[i,'Name'] + " add CH2CN- (from CH3CN) product"
    data.loc[data_index,'MW']=round(ketone_list.loc[i,'MW'] + 41.0519, 3)
    data.loc[data_index,'Properties']="Na-K-Z-N-"
    data_index += 1 

    data.loc[data_index,'Name']=ketone_list.loc[i,'Name'] + " add -CH2COOH (from CH3CN) product"
    data.loc[data_index,'MW']=round(ketone_list.loc[i,'MW'] + 60.0520, 3)
    data.loc[data_index,'Properties']="Na-K-Z-N-"
    data_index += 1 

    data.loc[data_index,'Name']=ketone_list.loc[i,'Name'] + " add CH2CN- (from CH3CN) product (H2O eliminated)"
    data.loc[data_index,'MW']=round(ketone_list.loc[i,'MW'] + 41.0519 - 18.011, 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1 

    data.loc[data_index,'Name']=ketone_list.loc[i,'Name'] + " add CH2COOH (from CH3CN) product (H2O eliminated)"
    data.loc[data_index,'MW']=round(ketone_list.loc[i,'MW'] + 60.0520 - 18.011, 3)
    data.loc[data_index,'Properties']=BasicProperties
    data_index += 1 

print("Reactions with CH2CN- (from CH3CN) products generated")

data.to_csv(OUT_PATH+'Generated species list.csv', index=False, encoding='utf-8-sig')
print("Species list generating from metabolite names completed. Saved to\""+OUT_PATH+'Generated species list.csv'+"\"")

# Adduct ion modes
modes=pd.DataFrame({
    'Mode': ['H', 'H', 'Na', 'Z', 'K', 'N'],
    'Description':['(add H+)', '(add H3O+)', '(add Na+)', '(self)', '(add K+)', '(add NH4+)'],
    'Multiplier': ['1', '1', '1', '1', '1', '1'],
    'Addition': ['1.008', '19.023', '22.99', '0', '39.098', '18.035']
})
PeakList=pd.DataFrame(columns=['Mass','Description'])

# Generating standard peaks ("the dictionary")
for i in range(len(data)):
    if data.loc[i,'SMILES']!='#':
        MW=data.loc[i,'MW']
        Properties=data.loc[i,'Properties'].split('-')
        for m in range(len(Properties)):
            for j in range(len(modes)):
                if (Properties[m] == modes.loc[j,'Mode']) and Properties[m]:
                    PeakMass=MW * float(modes.loc[j,'Multiplier']) + float(modes.loc[j,'Addition'])
                    Description = str(data.loc[i,'Name']) + str(modes.loc[j,'Description'])
                    PeakList.loc[len(PeakList.index)] = [round(PeakMass,3), Description]
    if i % 3000==0:
        print("Peak generated:"+str(i)+"/"+str(len(data)))

print("Peak generation completed")

# Sort peaks
PeakList=PeakList.sort_values(by='Mass')
PeakList=PeakList.reset_index(drop=True)

print("Peak values sorted")

# Combine peaks with the same m/z
Is_OK=False
while Is_OK == False:
    Is_OK = True
    LEN=len(PeakList)
    i=1
    while i<LEN:
        if PeakList.loc[i-1,'Mass'] == PeakList.loc[i,'Mass']:
            PeakList.loc[i-1,'Description']=PeakList.loc[i-1,'Description'] + "; " + PeakList.loc[i,'Description']
            PeakList=PeakList.drop([i])
            PeakList=PeakList.reset_index(drop=True)
            Is_OK=False
            LEN-=1
        i+=1

print("Peak values with same m/z values combined")
to_file_path=OUT_PATH+'Generated standard peaks list.csv'
PeakList.to_csv(to_file_path, index=False, encoding='utf-8-sig')
print("Standard peak list generating completed. Saved to\""+to_file_path+"\"")

# Data extraction 
MSdata=pd.read_csv(MS_DATA_PATH,header=None, encoding='utf-8-sig')
MSdata.columns=['Mass','Intensity']
#MSdata=MSdata.drop([0,1,2,3,4,5,6,7]) # For removing headings, add if necessary
MSdata=MSdata.reset_index(drop=True)
In=MSdata['Intensity']

# Find maximum intensity
Max=0
for i in range(0,len(MSdata)):
    if Max<float(In[i]):
        Max=float(In[i])

# Find peaks
Peaks=pd.DataFrame(columns=['Mass','Intensity'])
for i in range(1,len(MSdata)):
    if float(In[i])>float(In[i-1]) and float(In[i])>float(In[i+1]):
        if float(In[i])>int_min_thred:
            Peaks=pd.concat([Peaks,MSdata[i:i+1]])

Peaks=Peaks.reset_index(drop=True)
print('Peak value identification completed. Detected '+str(len(Peaks))+' MS peaks.')

# Eliminate non-related peaks according to file input
NR=pd.read_csv(NR_PATH,header=None, encoding='utf-8-sig')
if NR.empty == False:
    NR.columns = ['Mass','Descriptions']
    DropNRList=[]
    NRNum=0
    for i in range(1,len(Peaks)):
        for j in range(len(NR)):
            ERROR=abs(float(Peaks.loc[i,'Mass'])-float(NR.loc[j,'Mass']))
            if ERROR<mzdiff:
                DropNRList.append(i)
                NRNum+=1
    if NRNum!=0:
        Peaks=Peaks.drop(DropNRList)
        Peaks=Peaks.reset_index(drop=True)
    
    print("Eliminated "+str(NRNum)+" background peaks")

# Calculate relative intensity
Peaks.insert(Peaks.shape[1],'Relative%',None)
for i in range(len(Peaks)):
    Peaks.loc[i,'Relative%']=str(round(float(Peaks.loc[i,'Intensity'])/Max*100,2))+"%"

Peaks=Peaks.reset_index(drop=True)
print("Relative intensity calculation completed")

# Comparison
PeaksSTD=PeakList
PeaksSTD.columns=['Mass','Description']
PeaksRes=Peaks
PeaksRes=PeaksRes.drop(columns='Intensity')
PeaksRes.insert(PeaksRes.shape[1],'Description',None)

INDEX=mzdiff/0.001+1
MATCHED=0

for k in range(len(PeaksRes)):
    # Initialization, all descriptions are set to 'No matching standard peak' first.
    PeaksRes.loc[k,'Description']="No matching standard peak" 

    # Show matching progress
    if k%1000==0:
        print("Peak matching progress: "+str(k)+"/"+str(len(PeaksRes)))

    start_index = max(0, int(INDEX - mzdiff/0.001 - 1))

    # Matching
    for i in range(start_index, len(PeaksSTD)): 
        if float(PeaksSTD.loc[i,'Mass'])-float(PeaksRes.loc[k,'Mass'])>mzdiff:
            break
        # because the "resolution" is 0.001, if STD is smaller than INDEX-mzdiff/0.001-1, then no comparison is needed. 
        # for example, in extreme conditions, STD=[50.000, 50.001, ..., 50.008, 50.009], Res=[50.007, 50.008], then we can start from 50.007-6*0.001=50.001 when matching 50.007.

        # Comparison
        ERROR=abs(float(PeaksRes.loc[k,'Mass'])-float(PeaksSTD.loc[i,'Mass']))
        if ERROR<mzdiff:
            if PeaksRes.loc[k,'Description'] == "No matching standard peak":
                PeaksRes.loc[k,'Description']=PeaksSTD.loc[i,'Description'] + ", with the error of " + str(round(ERROR,4)) + ". "
            else:
                PeaksRes.loc[k,'Description']=PeaksRes.loc[k,'Description'] + PeaksSTD.loc[i,'Description'] + ", with the error of " + str(round(ERROR,4))+ "."

            INDEX=i+1

    if PeaksRes.loc[k,'Description'] != "No matching standard peak":
        MATCHED+=1

print("Peak matching completed, matching rate {:.2%}".format(MATCHED/len(PeaksRes)))

PeaksRes.to_csv(OUT_PATH+"Results.csv", index=False, encoding='utf-8-sig')

print("Output to file completed. Saved to \""+OUT_PATH+"Results.csv\"")
print("Process completed")