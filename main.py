#%%
import json 
import function
# %%
fname = input("file name: ")
with open('./'+fname, 'r') as f:
    Q = json.load(f)

print(json.dumps(Q) )

name = Q['name']
fnum = Q['fnum']
m1 = Q['m1']
m2 = Q['m2']
a = Q['a']
e = Q['e']
Nperiod = Q['Nperiod']
timestepmax = Q['timestepmax']
PNon = Q['PNon']
#%%
function.TwoBody(fnum, name, m1, m2, e, a, Nperiod, timestepmax, PNon)

# %%
