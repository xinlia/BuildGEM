# env python 3.9
# this is a module for KEGG download

import requests
from lxml import etree
import os
import xml.dom.minidom as xo
import pandas as pd


def getpwlist(org_name):
    url ='https://www.kegg.jp/kegg-bin/show_organism?menu_type=pathway_maps&org='+org_name
    headers = {'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/108.0.0.0 Safari/537.36 Edg/108.0.1462.54'}
    html = requests.get(url,headers = headers)
    if not html.ok:
        raise TypeError('wrong KEGG name')
    e = etree.HTML(html.text)
    strain_name = e.xpath('//font[@class="title1"]/text()')[0]
    print('Organism: '+strain_name+'\n')
    DATA = e.xpath('//ul/a/@href')
    pathway_ID_list = [i.replace('/pathway/','') for i in DATA]
    print(str(len(pathway_ID_list))+' pathway maps are found')
    return pathway_ID_list

def downloadxml(org_name,pathway_ID_list):
    # create folder
    path=".\\"+org_name+"\\xml"
    if not os.path.exists(path):
        os.makedirs(path)

    #download xml   
    i = 1
    for ID in pathway_ID_list:
        downloadurl = 'https://www.kegg.jp/kegg-bin/download?entry='+ID+'&format=kgml'
        headers={'Referer':'https://www.kegg.jp/pathway/'+ID,'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/108.0.0.0 Safari/537.36 Edg/108.0.1462.54'}
        xml=requests.get(downloadurl,headers = headers)
        if not xml.ok:
            print(ID+'.xml can not be found.   '+str(i)+'/'+ str(len(pathway_ID_list)))
            i+=1
            continue
        with open(path+'\\'+ID+".xml","wb") as code:
            code.write(xml.content)
        print(ID+'.xml write secceed.   '+str(i)+'/'+ str(len(pathway_ID_list)))
        i+=1

def getentriesfromxml(org_name):
    path=".\\"+org_name+"\\xml"
    xmls = [path+'\\'+i for i in os.listdir(path) if i.endswith('.xml')]
    if len(xmls) == 0:
        raise TypeError('xml files not found')

    df_xmlentry = pd.DataFrame(columns=["pathway_ID", "pathway_title", "entry_id", "entry_name","entry_type","entry_reaction","entry_link"])
    df_xmlreaction = pd.DataFrame(columns=["reaction_ID", "reaction_name", "reaction_reversible", "reaction_substrates","reaction_products"])
    i=1
    for xml in xmls[0:1]:
        domtree = xo.parse(xml)
        pathway = domtree.documentElement
        pathway_ID = pathway.getAttribute("name")
        pathway_org = pathway.getAttribute("org")
        pathway_title = pathway.getAttribute("title")
        entrys = pathway.getElementsByTagName("entry")
        reactions = pathway.getElementsByTagName("reaction")
        for entry in entrys:
            entry_id= pathway_ID+'-'+entry.getAttribute("id")
            entry_name = entry.getAttribute("name")
            entry_type = entry.getAttribute("type")
            entry_reaction = entry.getAttribute("reaction")
            entry_link = entry.getAttribute("link")
            data = [pathway_ID,pathway_title,entry_id,entry_name,entry_type,entry_reaction,entry_link]
            df_xmlentry.loc[len(df_xmlentry.index)]=data
            #graphics = entry.getElementsByTagName('graphics')[0]
            #graphics_name = graphics.getAttribute("name")
        for reaction in reactions:
            reaction_id = pathway_ID+'-'+reaction.getAttribute("id")
            reaction_name = reaction.getAttribute("name")
            reaction_reversible = reaction.getAttribute("type")
            substrates = reaction.getElementsByTagName("substrate")
            reaction_substrates = [substrate.getAttribute("name") for substrate in substrates]
            products = reaction.getElementsByTagName("product")
            reaction_products = [product.getAttribute("name") for product in products]
            data = [reaction_id, reaction_name, reaction_reversible, str(reaction_substrates),str(reaction_products)]
            df_xmlreaction.loc[len(df_xmlreaction.index)]=data
        print('File:',xml,'done.  ',i,'/',len(xmls))
        i+=1

    df_r = df_xmlreaction.rename(columns={'reaction_name':'entry_reaction'}).drop(columns='reaction_ID')
    df_org_reaction = pd.merge(df_xmlentry.query('entry_type == "gene" and entry_reaction != ""'),df_r,how='left', on='entry_reaction').drop_duplicates(subset=['entry_reaction'])
    df_org_reaction.to_excel(org_name+'\\'+org_name+'_KEGG_pathway_reaction.xlsx', sheet_name='Sheet1', header=True)

    titlename= ['entry_name','pathway_ID','pathway_title','reaction_reversible','reaction_substrates','reaction_products']
    df_gene = df_org_reaction[['entry_reaction']+titlename]
    df_reaction_list = pd.concat([df_gene['entry_reaction'].str.split('rn:',expand=True).iloc[:,1:],df_gene[titlename]],axis=1)
    df_reaction_melt = df_reaction_list.melt(id_vars = titlename, value_vars = [1,2,3,4],var_name = 'r',value_name = 'Kegg_ID').drop(['r'],axis = 1).dropna(subset=['Kegg_ID'])
    df_reaction_melt['Kegg_ID'] = df_reaction_melt['Kegg_ID'].str.strip()
    df_reaction = df_reaction_melt.drop_duplicates(subset=['Kegg_ID'])[['Kegg_ID']+titlename]
    df_reaction.to_excel(org_name+'\\'+org_name+'_clean_kegg_reaction.xlsx', sheet_name='Sheet1', header=True)

    return df_reaction