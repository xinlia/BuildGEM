# env python 3.9
# this is a module for building model


import cobra
biggunimodel = cobra.io.load_json_model('BIGG\\universal_model_modified.json')

def get_metabolite_from_unibigg(metabolite_id,biggunimodel=biggunimodel):
    metabolite = biggunimodel.metabolites.get_by_id(metabolite_id)
    cobra_metabolite = cobra.Metabolite(metabolite_id)
    cobra_metabolite.name = metabolite.name
    compartmentid = metabolite.id.split('_')[-1]
    if compartmentid in ['c','e','p']:
        cobra_metabolite.compartment = compartmentid
    else:
        raise ValueError(metabolite.id+' is not in c,e,p')
    cobra_metabolite.annotation = metabolite.annotation
    return cobra_metabolite

def add_rxn_from_unibigg(new_model,rxn_id_list,biggunimodel=biggunimodel):
    if not type(rxn_id_list) is list:
        rxn_id_list = [rxn_id_list]

    for rxn_id in rxn_id_list:
        rxn = biggunimodel.reactions.get_by_id(rxn_id)
        cobra_reaction = cobra.Reaction(rxn_id)
        cobra_reaction.name = rxn.name
        cobra_reaction.lower_bound = rxn.lower_bound
        cobra_reaction.upper_bound = rxn.upper_bound
        cobra_reaction.annotation = rxn.annotation

        for metabolite, stoich in rxn.metabolites.items():
            cobra_metabolite = get_metabolite_from_unibigg(metabolite.id)
            cobra_reaction.add_metabolites({cobra_metabolite: stoich})
        new_model.add_reactions([cobra_reaction])

               
def add_rxn_from_othermodel(new_model,rxn_id,drift_model):
    rxn = drift_model.reactions.get_by_id(rxn_id)
    cobra_reaction = cobra.Reaction(rxn_id)
    cobra_reaction.name = rxn.name
    cobra_reaction.lower_bound = rxn.lower_bound
    cobra_reaction.upper_bound = rxn.upper_bound
    cobra_reaction.annotation = rxn.annotation

    for metabolite, stoich in rxn.metabolites.items():
        try:
            new_model.metabolites.get_by_id(metabolite.id)
            cobra_reaction.add_metabolites({new_model.metabolites.get_by_id(metabolite.id): stoich})
        except:
            cobra_metabolite = cobra.Metabolite(metabolite.id)
            cobra_metabolite.name = metabolite.name
            compartmentid = metabolite.id.split('_')[-1]
            if compartmentid in ['c','e','p']:
                cobra_metabolite.compartment = compartmentid
            else:
                raise ValueError(metabolite.id+' is not in c,e,p')
            cobra_metabolite.annotation = metabolite.annotation
            cobra_reaction.add_metabolites({cobra_metabolite: stoich})
    new_model.add_reactions([cobra_reaction])


import re
def check_reversible(rxn_string):
    if '<' in rxn_string and '>' in rxn_string:
        return -1000,1000
    elif '>' in rxn_string:
        return 0,1000
    elif '<' in rxn_string:
        return -1000,0
    else:
        raise ValueError('rxn_string error')

def get_rxn_from_string(rxn_string):

    rxn_string_list = re.split('<=>|<-->|=>|-->|<=|<--|<-|<->|->',' '+rxn_string)
    if not len(rxn_string_list) == 2:
        raise ValueError('rxn_string error')
    
    stoichiometry = {}
    pattern = r'(\d*\.*\d*\s*)(\w+\_\w+)'
    matches1 = re.findall(pattern, rxn_string_list[0])
    matches2 = re.findall(pattern, rxn_string_list[1])
    for match in matches1:
        coefficient = match[0].strip()
        compound = match[1]
        if coefficient == '':
            coefficient = 1
        coefficient = -1*float(coefficient)
        stoichiometry[compound] = coefficient
    for match in matches2:
        coefficient = match[0].strip()
        compound = match[1]
        if coefficient == '':
            coefficient = 1
        coefficient = float(coefficient)
        stoichiometry[compound] = coefficient
    
    return stoichiometry

def add_rxn_by_string(model,rxn_id,rxn_name,rxn_string,biggunimodel=biggunimodel):
    try:
        if model.reactions.get_by_id(rxn_id):
            model.remove_reactions([model.reactions.get_by_id(rxn_id)])
    except:
        pass
    rxn = cobra.Reaction(rxn_id)
    rxn.name = rxn_name
    rxn.lower_bound,rxn.upper_bound = check_reversible(rxn_string)
    stoichiometry = get_rxn_from_string(rxn_string)
    for compound_id,coefficient in stoichiometry.items():
        try:
            meta = model.metabolites.get_by_id(compound_id)
            if meta:
                rxn.add_metabolites({meta:coefficient})
        except:
            try:
                meta = get_metabolite_from_unibigg(compound_id,biggunimodel=biggunimodel)
                if meta:
                    rxn.add_metabolites({meta:coefficient})
            except:
                compound = cobra.Metabolite(compound_id)
                compound.name = compound_id
                print(compound_id + ' is not in model and unibiggmodel,it need to add name and annotation manually')
                compound.compartment = compound_id.split('_')[-1]
                rxn.add_metabolites({compound:coefficient})
    model.add_reactions([rxn])

def set_rxn_direction(model,rxn_id,lb,ub):
    model.reactions.get_by_id(rxn_id).lower_bound = lb
    model.reactions.get_by_id(rxn_id).upper_bound = ub

def print_rxn(model,*rxn_id):
    if not rxn_id:
        rxn_ids = [rxn.id for rxn in model.reactions]
    else:
        rxn_ids = [i for i in rxn_id]
    for rxn_id in rxn_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        print(rxn.id+'\t'+rxn.name+'\t'+rxn.reaction)

def print_metabolite(model,*metabolite_id):
    if not metabolite_id:
        metabolite_ids = [metabolite.id for metabolite in model.metabolites]
    else:
        metabolite_ids = [i for i in metabolite_id]
    for metabolite_id in metabolite_ids:
        metabolite = model.metabolites.get_by_id(metabolite_id)
        print(metabolite.id+'\t'+metabolite.name)