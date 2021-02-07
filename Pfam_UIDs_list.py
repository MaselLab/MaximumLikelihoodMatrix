import pandas as pd

pfam_df = pd.read_csv("PfamUIDs.csv", sep="\t")

pfamID_list = list(pfam_df["PfamUID"].values)

check_pfamID_list = ['PF00457', 'PF07806', 'PF08854', 'PF10056', 'PF11397', 'PF12720', 'PF14527', 'PF15897', 'PF17182', 'PF17343']

result = []

i = 0
a = 0
while i < len(pfamID_list) and a < len(check_pfamID_list):
    if pfamID_list[i] == check_pfamID_list[a]:
        result.append(check_pfamID_list[a])
        a += 1
    elif pfamID_list[i] > check_pfamID_list[a]:
        a += 1
    i += 1

print(result)
'''
result_1 = []
for id in check_pfamID_list:
    if id in pfamID_list:
        result_1.append(id)
print(result_1)
'''