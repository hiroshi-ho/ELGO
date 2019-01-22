L_num=1000

with open('minidata','w') as m, open('go.obo', 'r') as go:
    counter = 0
    for line in go:
        m.write(line)
        counter += 1
        if counter == L_num:
            break

print('done')