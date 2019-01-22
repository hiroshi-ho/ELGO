import json, re

def GOprocessor(go_obo_filepath):
    GO_json = {}

    another_tag = []

    with open(go_obo_filepath) as gof:
        one_gene_data = {'synonym': [],
                                     'is_a': [],
                                     'relationship': [],
                                     'xref': [],
                                     'intersection_of': []}

        for line in gof:
            line = line.strip()
            if line == '':
                continue

            if re.match(r'\[Term\]',line):
                if 'id' not in one_gene_data:
                    continue

                GO_json.update({one_gene_data['id'] : one_gene_data})

                one_gene_data = {'synonym': [],
                                 'is_a': [],
                                 'relationship': [],
                                 'xref': [],
                                 'intersection_of': []}


            if line[0:3] == 'id:':
                one_gene_data['id'] = line[7:].strip()

            elif line[0:5] == 'name:':
                one_gene_data['name'] = line[6:].strip()

            elif line[0:10] == 'namespace:':
                one_gene_data['namespace'] = line[11:].strip()

            elif line[0:4] == 'def:':
                one_gene_data['def'] = re.search(pattern=r'def: "(.+)"', string=line.strip()).group(1)

            elif line[0:8] == 'synonym:':
                one_gene_data['synonym'].append(re.search(pattern=r'synonym: "(.+)"', string=line.strip()).group(1))


            elif line[0:5] == 'is_a:':
                one_gene_data['is_a'].append(re.search(pattern=r'is_a:\s(GO:.+)\s!', string=line.strip()).group(1))

            elif line[0:13] == 'relationship:':
                one_gene_data['relationship'].append(re.search(pattern=r'relationship: (.+) !', string=line.strip()).group(1))

            elif line[0:5] == 'xref:':
                one_gene_data['xref'].append(line[6:].strip())

            elif line[0:15] == 'intersection_of':
                one_gene_data['intersection_of'].append(re.search(pattern=r'intersection_of: (.+) !', string=line.strip()).group(1).split())

            else:
                another_tag.append(line[0:15])


        if len(one_gene_data) != 0:
            GO_json.update({one_gene_data['id'] : one_gene_data})

    return GO_json, another_tag








if __name__ == '__main__':
    dataset_dir = './dataset_dir/'

    go_obo_file = 'go.obo'
    go_obo_json = 'go.obo.json'

    words_file = 'types.txt'
    vectors_file = 'vectors.txt'

    GO_json, another_tag = GOprocessor(go_obo_filepath=dataset_dir+go_obo_file)

    with open(dataset_dir+go_obo_json,'w') as dj:
        json.dump(GO_json,dj,ensure_ascii=False, indent=4)

    print('exception tags',another_tag)