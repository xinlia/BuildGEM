# env python 3.9
# this is a module for KEGG mapping to BIGG

import json
with open('KEGG\\KEGG_BIGG_dir_selected.json', 'r') as f:
    KEGG_BIGG_dir_selected = json.load(f)
f.close()
with open('KEGG\\KEGG_pathways_name.json', 'r') as f:
    pathways_name_dir = json.load(f)
f.close()

import pandas as pd
import re
KEGG_pathways_entrys = pd.read_csv('KEGG\\KEGG_pathways_entrys.txt', sep='\t',dtype={'pw_ID':str,'pw_name':str,'href':str,'title':str}).drop_duplicates()
KEGG_pathways_reactions = pd.read_csv('KEGG\\KEGG_pathways_reactions.txt', sep='\t',dtype={'pw_ID':str,'pw_name':str,'coords':str,'href':str,'title':str,'R_id':str})
KEGG_pathways_kos = pd.read_csv('KEGG\\KEGG_pathways_kos.txt', sep='\t',dtype={'pw_ID':str,'pw_name':str,'coords':str,'href':str,'title':str,'KO_id':str})


def printmyrxn(model, rxn):
    if type(rxn) == str:
        rxn = model.reactions.get_by_id(rxn)
    print(rxn.id,',', rxn.name)
    #print('\t https://modelseed.org/biochem/reactions/'+rxn.id.split('_')[0])
    print('\t',rxn.reaction)
    for meta, sto in rxn.metabolites.items():
        metablic = model.metabolites.get_by_id(meta.id)
        print('\t\t',metablic.id,'\t', metablic.name,'\t', sto)


def findentry(entry,KEGG_pathways_entrys=KEGG_pathways_entrys,globalsearch=False):
    if globalsearch:
        search = KEGG_pathways_entrys.loc[KEGG_pathways_entrys["title"].str.contains(entry)]
    else:
        search = KEGG_pathways_entrys.iloc[17994:].loc[KEGG_pathways_entrys["title"].str.contains(entry)]     
    if len(search) == 0:
        raise TypeError('Cannot find %s!'%(entry))
    else:
        return search.reset_index()

def getPathways(entry,globalsearch=False):
    df =findentry(entry,globalsearch=globalsearch)
    pw_list = df['pw_ID'].tolist()
    return_list = list(set(pw_list))
    return_list.sort(key=pw_list.index)
    return return_list
          
def printurl(entry,globalsearch=False):
    pws = getPathways(entry,globalsearch=globalsearch)
    for pw in pws:
        print('https://www.kegg.jp/kegg-bin/show_pathway?map=%s&multi_query=%s+yellow'%(pw,entry))

def getRfromKO(entry,globalsearch=False):
    df =findentry(entry,globalsearch=globalsearch)
    R_list=[]
    for i in range(len(df)):
        string = df.loc[i,'title']
        R = re.findall(r"R\d{5}", string)
        R_list+=R
    if len(R_list) > 0:
        return_list = list(set(R_list))
        return_list.sort(key=R_list.index)
        if len(return_list) == 1:
            return return_list[0]
        else:
            return return_list
    else:
        print('Cannot find R for %s!'%(entry))
        return []

def getECfromKO(entry,globalsearch=False):
    df =findentry(entry,globalsearch=globalsearch)
    EC_list=[]
    for i in range(len(df)):
        string = df.loc[i,'title']
        elements = string.split(',')
        for element in elements:
            if len(element.split('.'))>3 and (not'(' in element):
                EC_list.append(element.strip())
    if len(EC_list) > 0:
        return_list = list(set(EC_list))
        return_list.sort(key=EC_list.index)
        if len(return_list) == 1:
            return return_list[0]
        else:
            return return_list
    else:
        print('Cannot find EC for %s!'%(entry))
        return []


def getRfromPathway(pw_ID,KEGG_pathways_reactions=KEGG_pathways_reactions):
    df = KEGG_pathways_reactions[KEGG_pathways_reactions['pw_ID']==pw_ID]
    R_list = df['R_ID'].tolist()
    return_list = list(set(R_list))
    return return_list





def comparemodel(rxnlist1,rxnlist2):
    both = [r for r in rxnlist1 if r in rxnlist2]
    only1 = [r for r in rxnlist1 if r not in rxnlist2]
    only2 = [r for r in rxnlist2 if r not in rxnlist1]
    print('both: %d, only1: %d, only2: %d'%(len(both),len(only1),len(only2)))
    return [both,only1,only2]


def pathway_count(rxnlist,*args):

    if args:
        refdict={}
        for arg in args:
            refdict.update(arg.to_dict())

    pathways =[]
    for r in rxnlist:
        try:
            pathways+=getPathways(r)
        except:
            pass
    import pandas as pd
    result = pd.value_counts(pathways)
    for index, value in result.items():
        name = pathways_name_dir.get(index)
        if args:
            if refdict.get(index):
                total_value = refdict.get(index)
            else:
                total_value = 'Inf'
            print(index, value,'/',total_value, name)
        else:
            print(index, value, name)

    return result


def printbypathway(rxnlist,countresult,KEGG_pathways_reactions=KEGG_pathways_reactions,pathways_name_dir=pathways_name_dir):
    refdict=countresult.to_dict()
    for map, value in refdict.items():
        name = pathways_name_dir.get(map)
        print(map, value, name)
        Rxn_in_pathway = KEGG_pathways_reactions[KEGG_pathways_reactions['pw_ID']==map]['R_ID'].to_list()
        for r in Rxn_in_pathway:
            if r in rxnlist:
                print('     ',r)
        print('\n')


def color_compare(rxnlist1,rxnlist2):
    both = [r for r in rxnlist1 if r in rxnlist2]
    only1 = [r for r in rxnlist1 if r not in rxnlist2]
    only2 = [r for r in rxnlist2 if r not in rxnlist1]
    for r in both:
        print(r,'yellow','pink')
    for r in only1:
        print(r,'yellow','white')
    for r in only2:
        print(r,'white','pink')

def get_key(val,KEGG_BIGG_dir_selected=KEGG_BIGG_dir_selected):
    for key, value in KEGG_BIGG_dir_selected.items():
         if val == value:
             return key
    return "key doesn't exist"

def get_pathway_from_bigg(bigg_id):
    R_id = get_key(bigg_id)
    print(R_id)
    try:
        print(getPathways(R_id))
        return getPathways(R_id)
    except:
        print('Cannot find pathway for %s!'%(bigg_id))
        return []