import json, re, nltk

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
                one_gene_data['is_a'].append(re.search(pattern=r'is_a: (.+) !', string=line.strip()).group(1))

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


def all_definition_getter_from_ontology(go_jsonpath):
    definition_sentence_list = list()
    with open(go_jsonpath,'r') as gj:
        go_json = json.load(gj)

        for go_id, its_data in go_json.items():
            definition_sentence_list.append(its_data['def'])

    return definition_sentence_list

def make_vocabulary_dict_from_all_sentences_in_ontology_definition(definition_set_list):
    vocab_list = list()

    for one_definition_paragraph in definition_set_list:
        sentence_list = nltk.sent_tokenize(one_definition_paragraph)

        for one_sentence in sentence_list:
            vocab_list += nltk.word_tokenize(one_sentence)

    return list(set(vocab_list))

def dictmaker(all_of_vocab):

    word2iddict = {}
    character2iddict = {}

    id2word_dict = {}
    id2character_dict = {}


    word2iddict.update({'<pad>':len(word2iddict)})
    character2iddict.update({'<pad>':len(character2iddict)})

    id2word_dict.update({len(id2word_dict):'<pad>'})
    id2character_dict.update({len(id2character_dict):'<pad>'})


    for one_vocab in all_of_vocab:
        one_word_revised = str((one_vocab)).lower()
        if not one_word_revised in word2iddict:
            word2iddict.update({one_word_revised: len(word2iddict)})
            id2word_dict.update({len(id2word_dict): one_word_revised})

        char_list = list(one_word_revised)
        for one_char in char_list:
            if not one_char in character2iddict:
                character2iddict.update({str(one_char): len(character2iddict)})
                id2character_dict.update({len(id2character_dict): str(one_char)})

    if not '<unknown_word>' in word2iddict.keys():
        word2iddict.update({'<unknown_word>':len(word2iddict)})
        id2word_dict.update({len(id2word_dict):'<unknown_word>'})

    if not '<unknown_char>' in character2iddict.keys():
        character2iddict.update({'<unknown_char>':len(character2iddict)})
        id2character_dict.update({len(id2character_dict):'<unknown_char>'})

    return word2iddict,character2iddict,id2word_dict,id2character_dict

def when_you_input_one_wordsplited_sentence_then_return_one_intnized_sentence_and_one_intnized_char_sentence(one_raw_wordsplitted_sentence,word2iddic,char2iddic):
    intnized_sentence = []
    intnized_char_sentence = []
    for one_word in one_raw_wordsplitted_sentence:
        if one_word in word2iddic:
            intnized_sentence.append(word2iddic[one_word])
        else:
            intnized_sentence.append(word2iddic['<unknown_word>'])

        characters_in_one_word = list(one_word)
        intnized_characters_in_one_word = []
        for one_raw_character in characters_in_one_word:

            if one_raw_character in char2iddic:
                intnized_characters_in_one_word.append(char2iddic[one_raw_character])
            else:
                intnized_characters_in_one_word.append(char2iddic['<unknown_char>'])

    return intnized_sentence, intnized_char_sentence


if __name__ == '__main__':
    dataset_dir = './dataset_dir/'

    go_obo_file = 'go.obo'
    go_obo_json = 'go.obo.json'

    words_file = 'types.txt'
    vectors_file = 'vectors.txt'

    GO_json, another_tag = GOprocessor(go_obo_filepath=dataset_dir+go_obo_file)

    with open(dataset_dir+go_obo_json,'w') as dj:
        json.dump(GO_json,dj,ensure_ascii=False, indent=4)

    # print('exception tags',another_tag)