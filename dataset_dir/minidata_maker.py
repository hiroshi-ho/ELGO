L_num=1000

with open('minidataAll_Data.gene_info','w') as m, open('All_Data.gene_info', 'r') as go:
    counter = 0
    for line in go:
        m.write(line)
        counter += 1
        if counter == L_num:
            break

print('done')