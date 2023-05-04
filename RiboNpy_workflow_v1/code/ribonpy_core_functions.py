# -*- coding: utf-8 -*-
# This code is used to introduce catalytic constraints in GEMs
# by COBRApy and to calculate the parameters that need to be entered
# during the construction of the catalysis-constrained model.
#from warnings import warn

import pandas as pd
import numpy as np
import json
import cobra
import math
import re
import random
import statistics
import os
import shutil
import random
from cobra.core import Reaction
from cobra.io.dict import model_to_dict
from cobra.util.solver import set_objective
from xml.dom import minidom
from optlang.symbolics import Zero, add

def convert_to_irreversible(model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.

    """
    #warn("deprecated, not applicable for optlang solvers", DeprecationWarning)
    reactions_to_add = []
    coefficients = {}
    for reaction in model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[
                reverse_reaction] = reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                             for k, v in reaction._metabolites.items()}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)


def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value
                    model.add_reaction(r_add)
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model


def get_genes_and_gpr(model,gene_outfile,gpr_outfile):
    """Retrieving genes and gene_reaction_rule from GEM.

    Arguments
    ----------
    * model: cobra.Model ~ A genome scale metabolic network model for
        constructing the enzyme-constrained model.

    :return: all genes and gpr in model.
    """
    model_dict = model_to_dict(model, sort=False)
    genes = pd.DataFrame(model_dict['genes']).set_index(['id'])
    genes.to_csv(gene_outfile)
    all_gpr = pd.DataFrame(model_dict['reactions']).set_index(['id'])
    all_gpr.to_csv(gpr_outfile)
    return [genes, all_gpr]


def get_reaction_gene_subunit_MW(reaction_gene_subunit_file,gene_molecular_weight_file,save_file):
    """Retrieving genes,subunits and MW, and split 'or' type of reaction

    Arguments
    ----------
    * reaction_gene_subunit_file: gene-molecular_weight file eg. b3500,48771.94
    * gene_molecular_weight_file: manually get the subunit of each protein from EcoCy  

    :return: all reaction with gene, subunit, and MW.
    """
    reaction_gene_subunit_MW_new = pd.DataFrame()
    reaction_gene_subunit = pd.read_csv(reaction_gene_subunit_file, index_col=0)
    protein_mw=pd.read_csv(gene_molecular_weight_file, index_col=0)
    for reaction, data in reaction_gene_subunit.iterrows():
        if re.search(" or ", data['gene_reaction_rule']):
            gene = enumerate(data['gene_reaction_rule'].split(" or "))
            subunit_num = data['subunit_num'].split(" or ")
            for index, value in gene:
                if index == 0:
                    reaction_new = reaction + "_num1"
                    reaction_gene_subunit_MW_new.loc[reaction_new,'name'] = data['name']
                    reaction_gene_subunit_MW_new.loc[reaction_new, 'gene_reaction_rule'] = value
                    reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_num'] = subunit_num[index]
                    if re.search(" and ", value):
                        reaction_gene_subunit_MW = []
                        gene2 = enumerate(value.replace('(', '').replace(")", '').replace(" ", '').split('and'))
                        for index2, value2 in gene2:
                            reaction_gene_subunit_MW.append(str(protein_mw.loc[value2,'mw']))
                        reaction_gene_subunit_MW=' and '.join(reaction_gene_subunit_MW)
                        if re.search('\(',value):
                            reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = '( '+reaction_gene_subunit_MW+ ' )'
                        else:
                            reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = reaction_gene_subunit_MW
                    else:
                        reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = protein_mw.loc[value,'mw']
                else:
                    reaction_new = reaction + "_num" + str(index+1)
                    reaction_gene_subunit_MW_new.loc[reaction_new, 'name'] = data['name']
                    reaction_gene_subunit_MW_new.loc[reaction_new, 'gene_reaction_rule'] = value
                    reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_num'] = subunit_num[index]
                    if re.search(" and ", value):
                        reaction_gene_subunit_MW = []
                        gene3 = enumerate(value.replace('(', '').replace(")", '').replace(" ", '').split('and'))
                        for index3, value3 in gene3:
                            reaction_gene_subunit_MW.append(str(protein_mw.loc[value3,'mw']))
                        reaction_gene_subunit_MW=' and '.join(reaction_gene_subunit_MW)
                        if re.search('\(',value3):
                            reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = '( '+reaction_gene_subunit_MW+ ' )'
                        else:
                            reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = reaction_gene_subunit_MW
                    else:
                        reaction_gene_subunit_MW_new.loc[reaction_new,'subunit_mw'] = protein_mw.loc[value,'mw']                                     
        elif re.search(" and ", data['gene_reaction_rule']):
            reaction_gene_subunit_MW = []
            gene4 = enumerate(data['gene_reaction_rule'].replace('(', '').replace(")", '').replace(" ", '').split('and'))
            reaction_gene_subunit_MW_new.loc[reaction, 'name'] = data['name']
            reaction_gene_subunit_MW_new.loc[reaction, 'gene_reaction_rule'] = data['gene_reaction_rule']
            reaction_gene_subunit_MW_new.loc[reaction,'subunit_num'] = data['subunit_num']
            for index4, value4 in gene4:
                reaction_gene_subunit_MW.append(str(protein_mw.loc[value4,'mw']))
            reaction_gene_subunit_MW=' and '.join(reaction_gene_subunit_MW)
            if re.search('\(',value4):
                reaction_gene_subunit_MW_new.loc[reaction,'subunit_mw'] = '( '+reaction_gene_subunit_MW+ ' )'
            else:
                reaction_gene_subunit_MW_new.loc[reaction,'subunit_mw'] = reaction_gene_subunit_MW
        else:
            reaction_gene_subunit_MW_new.loc[reaction, 'name'] = data['name']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'gene_reaction_rule'] = data['gene_reaction_rule']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'subunit_mw'] = protein_mw.loc[data['gene_reaction_rule'],'mw']
            reaction_gene_subunit_MW_new.loc[reaction,
                                             'subunit_num'] = data['subunit_num']
    reaction_gene_subunit_MW_new.to_csv(save_file)
    return reaction_gene_subunit_MW_new


def calculate_ribozyme_MW(nucleotide_MW_file, ribozyme_data_file, sequence_length, nucleotide_proportion=[0.25,0.25,0.25,0.25], random_proportion=False, sequence_motifs_file=None, new_ribozyme=True):
  """Calculate the molecular weight of an RNA based on an inputted or created RNA sequence

  Arguments
  ----------
  * nucleotide_MW_file: csv file with the MW of each RNA nucleotide (ATP, CTP, GTP, UTP) and the MW of H2O.
  * sequence_length: float, defines the length of the RNA sequence.
  * nucleotide_proportion: the proportions of AGCU (in order) to be used in the sequence. Default: [0.25,0.25,0.25,0.25]
  * random_proportion: boolean, if true then the proportions are randomised according to the given proportions, if false the exact given proportions are used.
  * sequence_motifs_file: optional, csv file with RNA sequence motifs. If inputted, motifs will be added to specifics part of sequence. Default: None
  * new_ribozyme: boolean condition:
    - True: runs the function and generates a ribozyme sequence and calculates its weight.
    - False (default): runs the function and generates a ribozyme sequence unless the outfile already contains a MW value.

  :return: ribozyme_data_file, a csv with sequence, length, MW and actual proportion of nucleotides in the sequence.
  """

  nucleotides = ["A", "G", "C", "U"]
  ribozyme_MW = pd.DataFrame()
  
  if not new_ribozyme and 'MW' in ribozyme_MW.columns and not ribozyme_MW['MW'].isnull().values.any():
    ribozyme_MW = pd.read_csv(ribozyme_data_file, index_col=0)
    print("MW column already exists and is not empty. Skipping...")
    return ribozyme_MW
  
  elif new_ribozyme: 
    # Read nucleotide molecular weight data from file
    nucleotide_MW = pd.read_csv(nucleotide_MW_file, index_col=0) # could we retrieve this from the GEM?

    # Create a new RNA sequence
    if random_proportion:
      sequence = ''.join(np.random.choice(nucleotides, sequence_length, nucleotide_proportion))

    elif not random_proportion:
      # Create a list of nucleotides with the specified proportions
      nucleotide_list = []
      for i in range(len(nucleotides)):
          nucleotide_list.extend([nucleotides[i]] * int(sequence_length * nucleotide_proportion[i]))
      
      # Shuffle the nucleotide list
      random.shuffle(nucleotide_list)
      
      # Create the final sequence by joining the shuffled nucleotide list
      sequence = ''.join(nucleotide_list)


    # Insert motifs if provided
    if sequence_motifs_file:
        motif_df = pd.read_csv(sequence_motifs_file)
        for motif, index in motif_df.iterrows():
            sequence = sequence[:index] + motif + sequence[index:]
        
    # Calculate molecular weight
    sequence_MW = sum(nucleotide_MW.loc[n, 'MW'] for n in sequence) - (sequence_length - 1) * nucleotide_MW.loc['H2O', 'MW']
                    
    # Calculate percentage of each nucleotide
    nucleotide_count = {'A': sequence.count('A'), 'G': sequence.count('G'), 'C': sequence.count('C'), 'U': sequence.count('U')}
    nucleotide_percentage = {k: v / sequence_length * 100 for k, v in nucleotide_count.items()}

    # Add results to dataframe
    ribozyme_MW = pd.DataFrame({'MW': [sequence_MW], 'length': [sequence_length], 'sequence': [sequence]})
    for nuc, percentage in nucleotide_percentage.items():
        ribozyme_MW[nuc] = round(percentage, 2)

    # write results to CSV
    ribozyme_MW.to_csv(ribozyme_data_file)
    return ribozyme_MW
  
  else:
    raise ValueError("new_ribozyme must be 'True' or 'False'")


def calculate_reaction_mw(reaction_gene_subunit_MW,reaction_mw_outfile):
    """Calculate the molecular weight of the enzyme that catalyzes each
    reaction in GEM based on the number of subunits and
    molecular weight of each gene.

    Arguments
    ----------
    * reaction_gene_subunit_MW: A CSV file contains the GPR relationship
     for each reaction in the GEM model,the number of subunit components 
     of each gene expressed protein, and the molecular weight of each 
     gene expressed protein.

    :return: The molecular weight of the enzyme that catalyzes each reaction
     in the GEM model.
    """
    reaction_gene_subunit_MW = pd.read_csv(
        reaction_gene_subunit_MW, index_col=0)
    reaction_mw = pd.DataFrame()
    for reaction_id in reaction_gene_subunit_MW.index:
        subunit_mw_list = reaction_gene_subunit_MW.loc[reaction_id, 'subunit_mw'].\
            replace('(', '').replace(")", '').replace(" ", '').split('or')
        subunit_num_list = reaction_gene_subunit_MW.loc[reaction_id, 'subunit_num'].\
            replace('(', '').replace(")", '').replace(" ", '').split('or')

        mw_s = ''
        for mw_i in range(0, len(subunit_mw_list)):
            mw_list = np.array(subunit_mw_list[mw_i].split('and'))
            num_list = np.array(subunit_num_list[mw_i].split('and'))
            mw_list = list(map(float, mw_list))
            num_list = list(map(float, num_list))
            mw_s = mw_s + \
                str(round(np.sum(np.multiply(mw_list, num_list)), 4)) + ' or '

        reaction_mw.loc[reaction_id, 'MW'] = mw_s.rstrip(' or ')
    reaction_mw.to_csv(reaction_mw_outfile)
    return reaction_mw


def modify_reaction_MW(reaction_MW_file, modified_reaction_MW_file, reaction_fraction, ribozyme_data_file, selection_type='random', reactions_subset_file=None, specific_reactions=None, new_reaction_selection=True):
    """Selects a set of reactions and changes their MW value with the ribozyme MW value to convert them into ribozymes.

    Arguments
    ----------
    * reaction_MW_file: str, the path to a CSV file that contains the MW values for each reaction in the model
    * reaction_fraction: float, a numerical percentage that determines how many values will be selected from the file (e.g. 0.5)
    * selection_type: str, the type of selection to use (default: 'random')
      - 'random': select values randomly from the entire CSV
      - 'subset': select values randomly from a specified subset of the CSV.
      - 'reactions': only modify the efficiency of specified reactions
    * reactions_subset_file: str, the path to a CSV file containing a subset of the reactions (only used if selection_type is 'subset')
    * specific_reactions: list of str, the names of reactions to modify (only used if selection_type is 'reactions')
    * new_reaction_selection: boolean condition:
      - False: reads the csv file with modified MW data.
      - True (default): runs the code and modifies the reaction MW for the selected reactions. Change to 'False' after running the function, otherwise
      the next time the code is run, the modified values will be randomised again.
    
    :return: 
    * modified_reaction_MW_file: A CSV file containing the reaction MW values for the reactions in the model, 
    with an additional column to indicate which values were modified (proteozyme - original or ribozyme - modified).
    """

    if not new_reaction_selection:
        # Read CSV file with modified MW data
        reaction_MW_df = pd.read_csv(modified_reaction_MW_file, index_col=0)
        return reaction_MW_df
    
    elif new_reaction_selection: 
        # Read the CSV file into a Pandas dataframe
        reaction_MW_df = pd.read_csv(reaction_MW_file, index_col=0)

        # Select reactions based on the specified selection type
        if selection_type == 'random':
            selected_reactions = reaction_MW_df.sample(frac=reaction_fraction)
        elif selection_type == 'subset':
            subset_df = pd.read_csv(reactions_subset_file)
            selected_reactions = subset_df.sample(frac=reaction_fraction)
        elif selection_type == 'reactions':
            selected_reactions = reaction_MW_df[reaction_MW_df[0].isin(specific_reactions)]
        else:
            raise ValueError("Invalid selection type. Must be 'random', 'subset', or 'reactions'.")

        # Check for additional reactions matching the selected reactions (applies to subunits and reversible reactions)
        additional_reactions = pd.DataFrame() # (index=reaction_MW_df.index, columns=reaction_MW_df.columns)
        for reaction in selected_reactions.index:
            selector = "^" + re.escape(reaction.split("_")[0]) # + "_"
            additional_matches = reaction_MW_df[reaction_MW_df.index.str.contains(selector) & ~reaction_MW_df.index.isin(selected_reactions.index)]
            additional_reactions = pd.concat([additional_reactions, additional_matches])

        # Combine selected_reactions with additional_reactions
        selected_reactions = pd.concat([selected_reactions, additional_reactions])

        # Select ribozyme MW value
        ribozyme_MW_value = pd.read_csv(ribozyme_data_file)['MW'].iloc[0]

        # Modify the MW of the selected reactions
        reaction_MW_df.loc[selected_reactions.index, 'MW'] = ribozyme_MW_value 
        reaction_MW_df.loc[selected_reactions.index, 'type'] = 'ribozyme'
        
        # Label the unmodified reactions as 'proteozyme'
        reaction_MW_df.loc[~reaction_MW_df.index.isin(selected_reactions.index), 'type'] = 'proteozyme'

        # Count the number of reactions labeled as ribozyme and the total number of reactions
        num_ribozyme_reactions = len(reaction_MW_df[reaction_MW_df['type'] == 'ribozyme'])
        total_num_reactions = len(reaction_MW_df)
        fraction_ribozyme_reactions = num_ribozyme_reactions/total_num_reactions*100
        print(f"{num_ribozyme_reactions} reactions ({round(fraction_ribozyme_reactions, 2)}%) have been modified out of {total_num_reactions} total reactions in the model.")
        
        # Save the modified reaction kcat values to a new CSV file
        reaction_MW_df.to_csv(modified_reaction_MW_file, index=True)
        
        return reaction_MW_df, fraction_ribozyme_reactions
    
    else:
        raise ValueError("new_reaction_selection must be 'True' or 'False'")


def modify_reaction_biosynthesis_MW(modified_reaction_MW_file, modified_biosynthesis_reaction_MW_file, biosynth_costs_penalty, biosynth_costs=False):
    """Calculates the biosynthetic costs associated with the synthesis of an individual proteozyme or ribozyme and adds them to the MW of the reaction.

    Arguments
    ----------
    * modified_reaction_MW_file: str, the path to a csv file containing the MW values for the reactions in the model.
    * modified_biosynthesis_reaction_MW_file: str, the CSV outfile path
    * biosynth_costs penalty: num, penalty imposed
    * biosynth_costs: bool, turns on and off the function

    :return: modified_biosynthesis_reaction_MW_file, str, the path to a csv file containing the modified MW values for the proteozyme reactions in the model
    """

    if not biosynth_costs:
        # Read CSV file with modified MW data and save to biosynthesis CSV file
        modified_biosynthesis_reaction_MW_df = pd.read_csv(modified_reaction_MW_file, index_col=0)
        modified_biosynthesis_reaction_MW_df.to_csv(modified_biosynthesis_reaction_MW_file)
        return modified_biosynthesis_reaction_MW_df
    
    elif biosynth_costs: 
        # Read the CSV file to a Pandas dataframe
        modified_reaction_MW_df = pd.read_csv(modified_reaction_MW_file, index_col=0)

        # Create an empty DataFrame for storing the modified data
        modified_biosynthesis_reaction_MW_df = pd.DataFrame()

        # Iterate through the rows of the original DataFrame
        for index, row in modified_reaction_MW_df.iterrows():
            if row['type'] == 'proteozyme':
                # Multiply the MW value by biosynth_costs_penalty for proteozyme reactions
                new_MW = row['MW'] * biosynth_costs_penalty
            else:
                # Leave the MW unmodified for ribozyme reactions
                new_MW = row['MW']

            # Save the modified data to the new DataFrame
            modified_biosynthesis_reaction_MW_df.loc[index, 'MW'] = new_MW
            modified_biosynthesis_reaction_MW_df.loc[index, 'type'] = row['type']

        # Save the new DataFrame to a CSV file
        modified_biosynthesis_reaction_MW_df.to_csv(modified_biosynthesis_reaction_MW_file)
        return modified_biosynthesis_reaction_MW_df
    
    else:
        raise ValueError("biosynth_costs must be 'True' or 'False'")
    

def calculate_modified_reaction_kcat_mw(reaction_kcat_file, modified_biosynthesis_reaction_MW_file, modified_reaction_kcat_mw_file, select_key, ef_factor):
    """Calculating kcat/MW

    Arguments
    ----------
    * reaction_kcat_file: A json file contains the kcat values for each
     reaction in the model.
    * modified_biosynthesis_reaction_MW_file: The molecular weight of the proteozyme that catalyzes
     each reaction in the GEM model, some modified to be transformed into ribozymes, biosynthetic costs accounted as well.
    * select_key:  ####3
    * ef_factor: float, efficiency factor applied to the proteozymes to be transformed in ribozymes, eg 0.01 (100 times less efficient).

    :return: modified_reaction_kcat_mw_file, a CSV file containing the kcat/MW value of the 
    proteozyme/ribozyme catalyzing each reaction in the GEM model.
    """
    reaction_kcat = json_load(reaction_kcat_file)
    modified_reaction_kcat_mw = pd.DataFrame()
    modified_reaction_MW = pd.read_csv(modified_biosynthesis_reaction_MW_file, index_col=0)
    for modified_reaction_idmw in modified_reaction_MW.index:
        if re.search('_num', modified_reaction_idmw):
            modified_reaction_id = modified_reaction_idmw.split('_num')[0]
        else:
            modified_reaction_id = modified_reaction_idmw
        for key,value in reaction_kcat.items():
            if re.search('_b',key):
                reaction_kcat_id = key.split('_b')[0]+'_reverse'
                if modified_reaction_id == reaction_kcat_id and str(value[select_key]) != 'nan':
                    modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'MW'] = modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                    if modified_reaction_MW.loc[modified_reaction_idmw, 'type'] == 'ribozyme':
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw,'kcat'] = value[select_key]*3600*ef_factor
                        kcat_MW = value[select_key]*3600*ef_factor/modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'type'] = 'ribozyme'
                    else:
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw,'kcat'] = value[select_key]*3600
                        kcat_MW = value[select_key]*3600/modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'type'] = 'proteozyme'
                    modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'kcat_MW'] = kcat_MW
            elif re.search('_f',key):
                reaction_kcat_id = key.split('_f')[0]
                if modified_reaction_id == reaction_kcat_id and str(value[select_key]) != 'nan':
                    modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'MW'] = modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                    if modified_reaction_MW.loc[modified_reaction_idmw, 'type'] == 'ribozyme':
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw,'kcat'] = value[select_key]*3600*ef_factor
                        kcat_MW = value[select_key]*3600*ef_factor/modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'type'] = 'ribozyme'
                    else:
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw,'kcat'] = value[select_key]*3600
                        kcat_MW = value[select_key]*3600/modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'type'] = 'proteozyme'
                    modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'kcat_MW'] = kcat_MW
            elif modified_reaction_id == key and str(value[select_key]) != 'nan':
                modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'MW'] = modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                if modified_reaction_MW.loc[modified_reaction_idmw, 'type'] == 'ribozyme':
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw,'kcat'] = value[select_key]*3600*ef_factor
                        kcat_MW = value[select_key]*3600*ef_factor/modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                        modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'type'] = 'ribozyme'
                else:
                    modified_reaction_kcat_mw.loc[modified_reaction_idmw,'kcat'] = value[select_key]*3600
                    kcat_MW = value[select_key]*3600/modified_reaction_MW.loc[modified_reaction_idmw, 'MW']
                    modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'type'] = 'proteozyme'
                modified_reaction_kcat_mw.loc[modified_reaction_idmw, 'kcat_MW'] = kcat_MW                
    modified_reaction_kcat_mw.to_csv(modified_reaction_kcat_mw_file)
    return modified_reaction_kcat_mw


def trans_model2enz_json_model_split_isoenzyme(model_file, modified_reaction_kcat_mw_file, f, ptot, sigma, lowerbound, upperbound, json_output_file):
    """Tansform cobra model to json model with  
    enzyme concentration constraint.

    Arguments
    ----------
    * model_file:   The path of sbml model
    * reaction_kcat_mw_file: The path of storing kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint. 

    """
    model = cobra.io.read_sbml_model(model_file)
    convert_to_irreversible(model)
    model = isoenzyme_split(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = "./model/%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    dictionary_model['enzyme_constraint'] = {'enzyme_mass_fraction': f, 'total_protein_fraction': ptot,
                                             'average_saturation': sigma, 'lowerbound': lowerbound, 'upperbound': upperbound}
    # Reaction-kcat_mw file.
    # eg. AADDGT,49389.2889,40.6396,1215.299582180927
    reaction_kcat_mw = pd.read_csv(modified_reaction_kcat_mw_file, index_col=0)
    reaction_kcay_mw_dict = {}
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
            reaction_kcay_mw_dict[reaction_id] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    dictionary_model['enzyme_constraint']['kcat_MW'] = reaction_kcay_mw_dict
    json_write(json_output_file, dictionary_model)


def get_enzyme_constraint_model(json_model_file):
    """using enzyme concentration constraint
    json model to create a COBRApy model.

    Arguments
    ----------
    * json_model_file: json Model file.

    :return: Construct an enzyme-constrained model.
    """

    dictionary_model = json_load(json_model_file)
    model = cobra.io.json.load_json_model(json_model_file)

    coefficients = dict()
    for rxn in model.reactions:
        if rxn.id in dictionary_model['enzyme_constraint']['kcat_MW'].keys():
            coefficients[rxn.forward_variable] = 1 / \
                float(dictionary_model['enzyme_constraint']['kcat_MW'][rxn.id])

    lowerbound = dictionary_model['enzyme_constraint']['lowerbound']
    upperbound = dictionary_model['enzyme_constraint']['upperbound']
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model


def get_fluxes_detail_in_model(model, fluxes_outfile, modified_reaction_kcat_mw_file):
    """Get the detailed information of each reaction

    Arguments
    ----------
    * model: cobra.Model.
    * fluxes_outfile: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.

    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """
    model_pfba_solution = cobra.flux_analysis.pfba(model)
    model_pfba_solution = model_pfba_solution.to_frame()
    reaction_kcat_mw = pd.read_csv(modified_reaction_kcat_mw_file, index_col=0)
    model_pfba_solution_detail = pd.DataFrame()
    for index, row in model_pfba_solution.iterrows():
        reaction_detail = model.reactions.get_by_id(index)
        model_pfba_solution_detail.loc[index, 'fluxes'] = row['fluxes']
        if index in reaction_kcat_mw.index:
            model_pfba_solution_detail.loc[index,
                                           'kcat'] = reaction_kcat_mw.loc[index, 'kcat']
            model_pfba_solution_detail.loc[index,
                                           'MW'] = reaction_kcat_mw.loc[index, 'MW']
            model_pfba_solution_detail.loc[index,
                                           'kcat_MW'] = reaction_kcat_mw.loc[index, 'kcat_MW']
            if 'source' in reaction_kcat_mw.columns:
                model_pfba_solution_detail.loc[index,
                                               'source'] = reaction_kcat_mw.loc[index, 'source']
        model_pfba_solution_detail.loc[index, 'equ'] = reaction_detail.reaction
    model_pfba_solution_detail.to_csv(fluxes_outfile)
    return model_pfba_solution


def get_enzyme_usage(enz_total, reaction_flux_file, reaction_kcat_mw_file, reaction_enz_usage_file):
    """Get the enzyme usage of each reaction

    Arguments
    ----------
    * enz_total: total enzyme amount(e.g. 0.227).
    * reaction_flux_file: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_enz_usage_file: the output file contain reaction and enzyme usage of each reaction.

    :return: reaction_enz_usage_file.
    """
    reaction_fluxes = pd.read_csv(reaction_flux_file, index_col=0)
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    reaction_enz_usage_df = pd.DataFrame()
    for index, row in reaction_kcat_mw.iterrows():
        if index in reaction_fluxes.index:
            reaction_enz_usage_df.loc[index, 'type'] = row['type']
            reaction_enz_usage_df.loc[index, 'kcat_mw'] = row['kcat_MW']
            reaction_enz_usage_df.loc[index,
                                    'flux'] = reaction_fluxes.loc[index, 'fluxes']
            reaction_enz_usage_df.loc[index,
                                    'enz useage'] = reaction_fluxes.loc[index, 'fluxes']/row['kcat_MW']
            reaction_enz_usage_df.loc[index, 'enz ratio'] = reaction_fluxes.loc[index,
                                                                                'fluxes']/row['kcat_MW']/enz_total
    reaction_enz_usage_df = reaction_enz_usage_df.sort_values(
        by="enz ratio", axis=0, ascending=False)
    reaction_enz_usage_df.to_csv(reaction_enz_usage_file)

    return reaction_enz_usage_df


def change_reaction_kcat_by_database(json_model_path,select_reactionlist, kcat_database_combined_file, reaction_kcat_mw_file, reaction_kapp_change_file):
    """Use the kcat in database to change reaction kcat in model

    Arguments
    ----------
    * json_model_path: The file storing json model.
    * select_reactionlist: reaction list need to change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_kapp_change_file: changed file stored reaction kcat/mw.

    :return: a dataframe stored new reaction kcat/mw .
    """
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0) #reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0).query("type == 'proteozyme'")
    Brenda_sabio_combined_select = json_load(kcat_database_combined_file)
    json_model=cobra.io.load_json_model(json_model_path)
    reaction_change_accord_database = []
    for eachreaction in select_reactionlist:
        if reaction_kcat_mw.loc[eachreaction, 'type'] == 'proteozyme':
            select_reaction = json_model.reactions.get_by_id(eachreaction)
            if "ec-code" in select_reaction.annotation.keys():
                ec_number = select_reaction.annotation["ec-code"]
                kcat_max_list = []
                if isinstance(ec_number, str):
                    if ec_number in Brenda_sabio_combined_select.keys():
                        reaction_kcat_max = Brenda_sabio_combined_select[ec_number]['kcat']
                        if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max * 3600:
                            reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max * 3600
                            reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600/reaction_kcat_mw.loc[eachreaction, 'MW']
                            reaction_change_accord_database.append(eachreaction) 
                else:
                    for eachec in ec_number:
                        if eachec in Brenda_sabio_combined_select.keys():
                            kcat_max_list.append(Brenda_sabio_combined_select[eachec]['kcat'])
                    reaction_kcat_max = np.max(kcat_max_list)     
                    if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max * 3600:
                        reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max * 3600
                        reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600/reaction_kcat_mw.loc[eachreaction, 'MW']
                        reaction_change_accord_database.append(eachreaction)
    reaction_kcat_mw.to_csv(reaction_kapp_change_file)
    return(reaction_change_accord_database)


def get_enz_model_use_enz_usage_by_eckcat(enz_ratio,json_model_path, reaction_flux_file, reaction_kcat_mw_file, reaction_enz_usage_file, kcat_database_combined_file, model_file, f, ptot, sigma, lowerbound, upperbound, json_output_file, reaction_kcat_mw_outfile):

    """Get new enzyme model using enzyme usage to calibration

    Arguments
    ----------
    * enz_ratio: enzyme ratio which needed change.
    * json_model_path: The file storing json model.
    * reaction_flux_file: reaction-flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_enz_usage_file： enzyme usage of each reaction.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * model_file: cobra model.
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint.  
    * json_output_file: json file store json model
    * reaction_kcat_mw_outfile: changed file stored reaction kcat/mw.

    :return: new enzyme model
    """
    reaction_enz_usage_df = get_enzyme_usage(
        upperbound, reaction_flux_file, reaction_kcat_mw_file, reaction_enz_usage_file)

    select_reaction = list(
        reaction_enz_usage_df[(reaction_enz_usage_df['type'] == 'proteozyme') & (reaction_enz_usage_df['enz ratio'] > enz_ratio)].index) # more than 1%

    print('need changing proteozyme reaction: ')
    print(select_reaction)
    change_reaction_list_round1 = change_reaction_kcat_by_database(json_model_path,select_reaction, kcat_database_combined_file, reaction_kcat_mw_file, reaction_kcat_mw_outfile)
    print('changed proteozyme reaction: ')
    print(change_reaction_list_round1)

    trans_model2enz_json_model_split_isoenzyme(
        model_file, reaction_kcat_mw_outfile, f, ptot, sigma, lowerbound, upperbound, json_output_file)

    enz_model = get_enzyme_constraint_model(json_output_file)
    return enz_model


def select_calibration_reaction_by_c13(reaction_kcat_mw_file, c13reaction_file, enzyme_amount, percentage, sigma):
    """Get reaction list need change kcat using c13 data

    Arguments
    ----------
    * reaction_kcat_mw_file: original file stored kcat/mw of reaction.
    * c13reaction_file: The file contained reaction list whcih had c13 flux.
    * enzyme_amount:  total enzyme amount(e.g. 0.227).
    * percentage:  percentage which needed change.
    * sigma: The approximated average saturation of enzyme. 
    
    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """    
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    c13reaction = pd.read_csv(c13reaction_file, index_col=0)
    c13reaction_select = []
    for index, row in c13reaction.iterrows():
        if index in reaction_kcat_mw.index and reaction_kcat_mw.loc[index, 'type'] == 'proteozyme':
            RiboNpy_c13_reaction_flux = reaction_kcat_mw.loc[index,
                                                           'kcat_MW']*enzyme_amount*percentage*sigma
            if RiboNpy_c13_reaction_flux < row['Flux norm']:
                c13reaction_select.append(index)
    return(c13reaction_select)


def get_enz_model_use_c13(reaction_kcat_mw_file,json_model_path, c13reaction_file, c13_percentage, kcat_database_combined_file, model_file, f, ptot, sigma, lowerbound, upperbound, json_output_file, reaction_mw_outfile):
    """Get new enzyme model using C13 reaction to calibration

    Arguments
    ----------
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * json_model_path: The file storing json model.
    * c13reaction_file: The file contained reaction list whcih had c13 flux.
    * c13_percentage:  percentage which needed change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.  
    * model_file: cobra model.
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint.  
    * json_output_file: json file store json model
    * reaction_mw_outfile: changed file stored reaction kcat/mw.

    :return: new enzyme model.
    """
    c13reaction_select = select_calibration_reaction_by_c13(
        reaction_kcat_mw_file, c13reaction_file, upperbound, c13_percentage, sigma)
    print('need changing proteozyme reaction: ')
    print(c13reaction_select)

    #if isinstance(df_reaction_select, pd.DataFrame):
    #    reaction_kcat_mw_file = "./analysis/toy/reaction_change_by_biomass.csv"

    reaction_kapp_change_file = reaction_mw_outfile
    #c13reaction_selecet=['CS','ACONTa','ACONTb','ICDHyr','MALS', 'MDH', 'ICL', 'SUCOAS_reverse', 'SUCDi', 'AKGDH']
    change_reaction_list_round1 = change_reaction_kcat_by_database(json_model_path,
        c13reaction_select, kcat_database_combined_file, reaction_kcat_mw_file, reaction_kapp_change_file)
    print('changed proteozyme reaction: ')
    print(change_reaction_list_round1)

    reaction_kcat_mw_file = reaction_mw_outfile
    trans_model2enz_json_model_split_isoenzyme(
        model_file, reaction_kcat_mw_file, f, ptot, sigma, lowerbound, upperbound, json_output_file)
    enz_model = get_enzyme_constraint_model(json_output_file)
    return enz_model


def get_diff_reaction_use_c13(c13reaction_file, model_fluxes):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * c13reaction_file: The file contained reaction list whcih had c13 flux.
    * model_fluxes: calulated fluxes
    
    :return: reaction list.
    """   
    c13reaction = pd.read_csv(c13reaction_file, index_col=0)
    c13reaction = list(c13reaction.index)
    enz_model_pfba_solution_select = model_fluxes[model_fluxes['fluxes'] > 0]
    enz_model_pfba_solution_select_id = []
    for eachreaction in enz_model_pfba_solution_select.index:
        if re.search('_num', eachreaction):
            enz_model_pfba_solution_select_id.append(
                eachreaction.split('_num')[0])
        else:
            enz_model_pfba_solution_select_id.append(eachreaction)
    c13reaction_2_enz_model_diff = list(
        set(c13reaction).difference(set(enz_model_pfba_solution_select_id)))
    return(c13reaction_2_enz_model_diff)


def get_ribozyme_constraint_model_pfba_obj_solution(model_name, obj_reaction, sequence_length, nucleotide_proportion, random_proportion, new_ribozyme, reaction_fraction, selection_type, new_reaction_selection, biosynth_costs, biosynth_costs_penalty, ef_factor):
    """Add ribozyme constraint to model and get the pFBA solution of the objective function.
    For a more detailed explanation on how the function works, see the jupyter notebook.

    :return: float, pFBA solution flux of the objective function of the model
    """

    # Step 1: construct the raw ribozyme-constrained model

    ### 2.1: convert to irreversible and split isoenzyme

    model = cobra.io.read_sbml_model(model_name)

    convert_to_irreversible(model)

    model = isoenzyme_split(model)

    #### 2.2: Retrieving enzyme and ribozyme kinetics and proteomics data

    get_genes_and_gpr(model,gene_outfile,gpr_outfile)

    get_reaction_gene_subunit_MW(reaction_gene_subunit_file,gene_molecular_weight_file,reaction_gene_subunit_MW_file)

    calculate_ribozyme_MW(nucleotide_MW_file, ribozyme_data_file, sequence_length, nucleotide_proportion, random_proportion, sequence_motifs_file, new_ribozyme)

    modify_reaction_MW(reaction_MW_file, modified_reaction_MW_file, reaction_fraction, ribozyme_data_file, selection_type, reactions_subset_file, specific_reactions, new_reaction_selection)

    modify_reaction_biosynthesis_MW(modified_reaction_MW_file, modified_biosynthesis_reaction_MW_file, biosynth_costs_penalty, biosynth_costs)

    calculate_modified_reaction_kcat_mw(reaction_kcat_file, modified_biosynthesis_reaction_MW_file, modified_reaction_kcat_MW_file, select_key, ef_factor)

    ### 2.3: Save ribozyme concentration constraint model as json file
    #The enzyme mass fraction 
    f = 0.406
    # The total protein fraction in cell.
    ptot = 0.56 
    # The approximated average saturation of enzyme.
    sigma = 1 
    lowerbound = 0   
    upperbound = round(ptot * f * sigma, 3)

    #create enzyme concentration constraint model
    trans_model2enz_json_model_split_isoenzyme(model_name, modified_reaction_kcat_MW_file, f, ptot, sigma, lowerbound, upperbound, json_output_file)

    enz_model=get_enzyme_constraint_model(json_output_file)

    enz_model_pfba_solution = get_fluxes_detail_in_model(enz_model,RiboNpy_fluxes_outfile,reaction_kcat_MW_file) # no needed for the main function?
    
    print(enz_model_pfba_solution.fluxes[obj_reaction])

    # Step 2: construct final ribozyme-constrained model

    ### 3.1: Calibration enzyme kcat according enzyme usage 
    enz_model=get_enz_model_use_enz_usage_by_eckcat(enz_ratio,json_model_path,fluxes_infile_ori,reaction_kcat_MW_file,\
                                      reaction_enz_usage_file,kcat_database_combined_file, model_name, \
                                      f, ptot, sigma, lowerbound, upperbound, json_round1_output_file, \
                                      reaction_kcat_MW_round1_outfile)

    enz_model_pfba_solution = get_fluxes_detail_in_model(enz_model,round1_fluxes_outfile,reaction_kcat_MW_round1_outfile)
    print(enz_model_pfba_solution.fluxes[obj_reaction])

    c13reaction_2_enz_model_diff = get_diff_reaction_use_c13(c13reaction_file,enz_model_pfba_solution)
    print (c13reaction_2_enz_model_diff)

    ### 3.2: Calibration enzyme kcat according c13 reaction list
    enz_model=get_enz_model_use_c13(reaction_kcat_MW_round1_outfile, json_model_path, c13reaction_file, c13_percentage, \
                                kcat_database_combined_file,model_name, f, ptot, sigma, lowerbound, \
                                upperbound, json_round2_output_file,reaction_kcat_MW_round2_outfile)

    enz_model_pfba_solution = get_fluxes_detail_in_model(enz_model,round2_fluxes_outfile,reaction_kcat_MW_round2_outfile)
    print(enz_model_pfba_solution.fluxes[obj_reaction])

    c13reaction_2_enz_model_diff = get_diff_reaction_use_c13(c13reaction_file,enz_model_pfba_solution)
    print (c13reaction_2_enz_model_diff)

    ### 3.3: Solving ribozyme concentration constraint by COBRApy
    enz_model=get_enzyme_constraint_model(json_round2_output_file)
    pfba_solution = cobra.flux_analysis.pfba(enz_model)
    pfba_solution_df = pfba_solution.to_frame()
    pfba_solution_df.to_csv(RiboNpy_solution_df_pfba_file)
    print(pfba_solution.fluxes[obj_reaction])

    return pfba_solution.fluxes[obj_reaction], pfba_solution_df
