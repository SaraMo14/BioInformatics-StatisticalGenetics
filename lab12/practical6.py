import numpy as np
import pandas as pd

#%%
# Replace 'dtype' with the actual data type and 'filename.raw' with your file name
data = pd.read_csv("yri.csv")

np_dat = data.iloc[:, 5:].values
#%%
n_ind = np_dat.shape[0]
n_snps = np_dat.shape[1]

#%%


mean_dist = np.zeros((n_ind, n_ind))
var_dist = np.zeros((n_ind, n_ind))

from tqdm import tqdm

for i in tqdm(range(n_ind)):
    for j in range(i+1, n_ind):
        shared = 2 - abs(np_dat[i,:] - np_dat[j,:])
        mean_dist[i, j] = shared.mean()
        mean_dist[j, i] = shared.mean()
        var_dist[i, j] = shared.std()
        var_dist[j, i] = shared.std()


#if you want the rownames, we would have to map them from data
print(mean_dist[:5])
print(var_dist[:5])

#%%
print(pd.DataFrame(mean_dist[:5, :5]).to_latex(index=False, caption="Mean", label="tab:mean"))
print(pd.DataFrame(var_dist[:5, :5]).to_latex(index=False, caption="Mean", label="tab:mean"))


#%%
fract_of_zeros = np.zeros((n_ind, n_ind))
fract_of_sharing_both = np.zeros((n_ind, n_ind))

for i in tqdm(range(n_ind)):
    for j in range(i+1, n_ind):
        fract_of_sharing_both[i,j] = np.sum((np_dat[i,:] - np_dat[j,:]) == 0) / n_snps
        fract_of_sharing_both[j, i] = fract_of_sharing_both[i, j]
        fract_of_zeros[i, j] = np.sum(abs((np_dat[i, :] - np_dat[j, :])) == 2) / n_snps
        fract_of_zeros[j, i] = fract_of_zeros[i, j]

print(pd.DataFrame(fract_of_zeros[:5, :5]).to_latex(index=False, caption="Mean", label="tab:mean"))
print(pd.DataFrame(fract_of_sharing_both[:5, :5]).to_latex(index=False, caption="Mean", label="tab:mean"))

#%%
m = mean_dist
s = var_dist
p_0 = fract_of_zeros
p_2 = fract_of_sharing_both
check_m = m - (1 - p_0 + p_2)
np.fill_diagonal(check_m, 0)

#%%



#%%
parent_offspring = np.zeros((n_ind, n_ind))

data_indexed = data.copy()
data_indexed["indexer"] = data_indexed.index
data_indexed.set_index("IID", inplace=True)

for idx, row in data.iterrows():
    if row["PAT"] != '0':
        pat = row["PAT"]
        try:
            idx_pat = data_indexed.loc[pat, "indexer"]
            parent_offspring[idx, idx_pat] = 1
            parent_offspring[idx_pat, idx] = 1
        except:
            pass
    if row["MAT"] != '0':
        pat = row["MAT"]
        try:
            idx_pat = data_indexed.loc[pat, "indexer"]
            parent_offspring[idx, idx_pat] = 1
            parent_offspring[idx_pat, idx] = 1
        except:
            pass

#%%

import matplotlib.pyplot as plt
plt.scatter(m, s, c=parent_offspring, cmap='bwr')
plt.xlabel("mean of shared alleles")
plt.ylabel("standard deviation of shared alleles")
plt.xlim([1.2,1.55])
plt.ylim([0.45,0.7])

plt.show()
plt.scatter(p_0, p_2, c=parent_offspring, cmap='bwr')
plt.xlabel("fract of 0 shared alleles")
plt.ylabel("fract of 2 shared alleles")
plt.xlim([-0.01,0.14])
plt.ylim([0.3,0.55])
plt.legend()

plt.show()