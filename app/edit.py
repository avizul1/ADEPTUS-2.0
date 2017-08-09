# -*- coding: utf-8 -*-

import subprocess
from subprocess import Popen
import random, string, re, os
import numpy as np
import time, math
from shutil import copyfile

#Define pathes for the R scrits working dir
#home_dir_1 = '/specific/scratch/davidama/django/app/RScripts/'
#home_dir_1 = '/specific/scratch/amirvizel/django/app/RScripts/'
home_dir_1 = '/specific/netapp5/gaga/Spike/adeptus/app/RScripts/'

#home_dir_1 = '/specific/scratch/davidama/django/app/'
#home_dir_2 = '/specific/scratch/amirvizel/django/app/'
home_dir_2 = '/specific/netapp5/gaga/Spike/adeptus/app/'

'''
input: gene list as string line, seperate by white-spaces(" ") and disease name (from close select list at the broswer)
output: string that represant the nodes variables section in the cytoscape js code on the browser and the nodes name as a list.
'''
def create_nodes(query, disease, tissue, qVal, rep_sig = 0.5, spec_sig = 0.7, bcg_sig = 0.7):
    # Each node line at the html page will contain those str + relevant data
    nodes_str = 'nodes: [\n'
    start_node = '                        { data: { id: \''
    node_dr = 'drblack: 0, drwhite: 2.5,'
    nodes_end = '                    ],\n\n'

    if disease != None:
        disease = str(disease)
    # check if the query is on specific genes or on the whole disease
    if(query == 'disease'): # 'go' == magic word for whole disease query
        query_lst = get_all_signific_disease_genes(disease, tissue, qVal, rep_sig, spec_sig, bcg_sig)   # get the lst of the disease most significant genes
    else: # or create lst from the user query
        #reunion the string query:        
        temp_str = '' #str for combining the query from the splited list
        for q in query: #add each part from the list to the string
            temp_str += q
        query = temp_str
        query_lst = re.split("[ \t\n,]+",query.strip())        
        
    #remove genes that not exist in the data, based on BCG and NEG data
    query_lst = remove_unknown_genes(query_lst)
    
    # for each gene in the query - build is html line
    for i in range(len(query_lst)):     
        query_lst[i] = query_lst[i].strip().upper() # to uppercase and remove whitespaces
        neg_val,bcg_val = get_exp_values(disease, tissue, query_lst[i]) # check their score
        # define their color
        node_ur = define_up_colors(bcg_val,"right", bcg_sig)
        node_ul = define_up_colors(neg_val,"left", bcg_sig)
        node_dl = define_dl_colors(query_lst[i])
        node_dr = define_dr_colors(query_lst[i],disease)
        # build the string line        
        nodes_str += start_node + query_lst[i] + node_ur + node_dr + node_dl + node_ul
    
    nodes_str += nodes_end
    
    return nodes_str,query_lst  # also return the query list for use in the edges build


'''
input: genes list and edge type (from close select list at the broswer, default - ppi intact)
output: string that represant the edges variables section in the cytoscape js code on the browser.
'''
def create_edges(query_lst,edge_type):
    # Each edge line at the html page will contain those str + relevant data    
    edges_str = 'edges: [\n'
    start_edge = '                        { data: { id: \''
    edge0 = '\', weight: '
    edge1 = ', source: \''
    edge2 = '\', target: \''
    edge3 = '\' } },\n'
    edges_end = '\' } }\n                    ]\n'
    
    # create edges list (tuples) by the edge type and the genes list
    lst = create_edges_lst(query_lst,edge_type)
    if (len(lst) == 0): 
        return (edges_str + ']\n', ('No', 'edges'))

    for i in range(len(lst) - 1): #for each edge tuple, build the line
        one_edge_str = start_edge + lst[i][0] + lst[i][1] + edge0 + str(i) + edge1 + lst[i][0] + edge2 + lst[i][1] + edge3
        edges_str += one_edge_str

    edges_str += start_edge + lst[len(lst) - 1][0] + lst[len(lst) - 1][1] + edge0 + str(len(lst) - 1) + edge1 + lst[len(lst) - 1][0] + edge2 + lst[len(lst) - 1][1] + edges_end
    
    return edges_str, lst


'''
the function return string that will placed in the html page.
the string create select form on the browser, 
and will contain all the disease that supported by the data.
the disease list based on the BCG data file.
'''
def create_disease_select_form():
    select_str = ''
    
    for name in disease_lst: # for each disease - build the select line
        select_str += "<option value=\"" + name + "\">" + name[0].capitalize() + "</option>\n"
    
    return select_str

'''
Input: user selection of tissue and sel_type(control/general/default = with tissue)
Output: string of cy elements (nodes and edges)
'''
def create_disease_onto_graph(sel_type,tissue):
    # node line strings
    nodes_str = 'nodes: [\n'
    node_start = '{ data: { id: \''
    node_end = '\' } },\n'
    nodes_str_end = '     ],\n'
    
    if sel_type != 'control' and sel_type != 'default': #sel_type is general, no specific tissue
        tmp = list(disease_ontology_dic[tissue])
        dis_lst = list(set([t[0]for t in tmp] + [t[1]for t in tmp]))
        key_lst = list(disease_ontology_dic[tissue])

    elif sel_type == 'default':
#########        
        keys = list(disease_ontology_dic.keys())
        dis_lst = []
        key_lst = []
        for key in keys:
        #per key=tissue, get all his disease, if our input is in this list, return the key
            dis_lst += list(set([str(tup[0])+'('+str(key)+')' for tup in disease_ontology_dic[key]] + [str(tup[1])+'('+str(key)+')' for tup in disease_ontology_dic[key]]))
            key_lst += [(str(tup[0])+'('+str(key)+')',str(tup[1])+'('+str(key)+')') for tup in disease_ontology_dic[key]]
#        tmp = list(disease_ontology_dic.values())
#        key_lst = set([])        
#        for lst in tmp:
#            key_lst.update(lst)
#        key_lst = list(key_lst)
#        dis_lst = list(set([t[0] for t in key_lst] + [t[1] for t in key_lst]))
    
#    if 'control' in dis_lst:
#        dis_lst.remove('control')
    for dis in dis_lst:
        if dis.find('control') != -1:
            dis_lst.remove(dis)

    for i in range(len(dis_lst)):
        if dis_lst[i].find('None') != -1:
            dis_lst[i] = dis_lst[i].split('(')[0]
 
    for i in range(len(key_lst)):
        if key_lst[i][0].find('None') != -1 or key_lst[i][1].find('None') != -1 :
            key_lst[i] = (key_lst[i][0].split('(')[0], key_lst[i][1].split('(')[0])
###########
    for disease in dis_lst: #create each node string - each disease got node
        nodes_str += node_start + disease + node_end
    
    nodes_str += nodes_str_end  # close the nodes string

    #edge line strings
    edges_str = 'edges: [\n'
    edge_start = '{ data: { source: \''
    edge_mid = '\', target: \''
    edge_end = '\' } },\n'
    edges_str_end = '    ]\n'
    
    for key in key_lst: # create each edge string
        if key[1] != key[0]:
            edges_str += edge_start + key[1] + edge_mid + key[0] + edge_end
    
    edges_str += edges_str_end  # close edges string
    
    return nodes_str + edges_str    # return the comon string - the cy elements

    
#---- HELP FUNCTIONS ----#

# -- PATHES --#
ppi_intact_file = r'app/data/ppi_intact.txt'
genemania_gi_file = r'app/data/genemania_gi.txt'
genemania_ppi_file = r'app/data/genemania_ppi.txt'
background_roc_file = r'app/data/background_roc.txt'
negative_roc_file = r'app/data/negative_roc.txt'
disease_ontology_file = r'app/data/selected_terms_network_formatted.txt'
rep_score_file = r'app/data/replicability_scores.txt'
spec_score_file = r'app/data/specificity_scores.txt'
smqs_score_file = r'app/data/smqs_scores.txt'
drugs_file = r'app/data/drug2genes_table.txt'
genes_id_file = r'app/genes_id_table.txt'
DOID_to_name_file = r'app/data/DOID_to_name.txt'
Lawrance_mut_file = r'app/data/mutated_genes.txt'
selected_terms = r'app/data/selected_terms_formatted.txt'
label2metaData_file = r'app/data/label2metadata.txt'

# -- replace with those for internal tests: -- #
#ppi_intact_file = r'data/ppi_intact.txt'
#genemania_gi_file = r'data/genemania_gi.txt'
#genemania_ppi_file = r'data/genemania_ppi.txt'
#background_roc_file = r'data/background_roc.txt'
#negative_roc_file = r'data/negative_roc.txt'
#disease_ontology_file = r'data/disease_ontology_subnetwork.txt'
#rep_score_file = r'data/gene_smqs.txt'
#drugs_file = r'data/drug2genes_table.txt'
#genes_id_file = r'genes_id_table.txt'

# -- GLOBAL VARIABLES -- #
ppi_intact_dic = {}
genemania_ppi_dic = {}
genemania_gi_dic = {}
disease_name_dic = {}
background_genes_name_dic = {}
negative_genes_name_dic = {}
disease_ontology_dic = {}
have_drugs_dic = {}
DOID_to_name_dic = {}
name_to_DOID_dic = {}
Lawrance_mut_dic = {}
dis_reports_dic = {}
tissue_DO_set = set([])
general_DO_set = set([])

##################### Numpy functions ##################
# Added by David Amar, Jan 2017
def read_matrix_obj_using_numpy(path):
    try:
        f = open(path,'r')
        header = f.readline().rstrip().split("\t")
        f.close()
        mat = np.loadtxt(path,skiprows=1,usecols=range(1,1+len(header)))
        genes = np.loadtxt(path,skiprows=1,usecols=[0],dtype=str)
        for i in range(len(genes)):
            #note: the first replace is because we made manipulate on the disease names in the data file
            #the second replcae is for js problem with - in vars
            #the third replace is for names with ' - make confused with strings
            genes[i] = genes[i].upper().replace("_"," ").replace("-","_").replace("'","")
        if len(header)-mat.shape[1] == 1:
            header = header[1:]
        gene2ind = {}
        for i in range(len(genes)):
            gene2ind[genes[i]] = i
        col2ind = {}
        for i in range(len(header)):
            col2ind[header[i]] = i
        obj = [col2ind,gene2ind,mat,header,genes]  
        return obj
    except Exception as e:
        return None
        

def get_matrix_column(term,matrix_obj):
    if not matrix_obj[0].has_key(term):
        return None
    return matrix_obj[2][:,matrix_obj[0][term]]

def get_matrix_value(term,gene,matrix_obj):
    gene = gene.replace("-","_")    
    if not matrix_obj[0].has_key(term):
        return None
    if not matrix_obj[1].has_key(gene):
        return None
    gene_ind = matrix_obj[1][gene]
    col_ind = matrix_obj[0][term]
    return matrix_obj[2][gene_ind,col_ind]

########################################################
'''
input: string (usually from text file)
output: the string with spaces(' ') instead of any other regular seperator (' ','\t','\n',',')
'''
def parse_gene_file (file_text):
    file_text = file_text.decode()
    lst = re.split("[ \t\n,]+",file_text)    
    return ' '.join(lst) 

'''
this function called by the run.py function when the server up.
its initialized all the the global variables - 
actually, in that way we processed the data files once, when the server upload. 
'''
def init_data_vars():
    # more GLOBAL variables
    global background_mat
    global negative_mat
    global rep_score_mat
    global spec_score_mat
    global smqs_score_mat
    global rep_score_lines_lst
    global background_lines_lst
    global negative_lines_lst         
    global disease_lst
    global gene_to_id_lst
    global DOID_to_name
    global Lawrance_mut
    global general_DO
    global tissue_DO
    global control_DO
    global selected_terms_lst
    
    general_DO = ''
    tissue_DO = '' 
    control_DO = ''
    #need forcreate_disease_name_dic():
    background_lines_lst = create_line_lst(background_roc_file)
    selected_terms_lst = create_line_lst(selected_terms)
    
    DOID_to_name = create_line_lst(DOID_to_name_file)
    for line in DOID_to_name:
        tmp = line.split("\t")        
        name_to_DOID_dic[tmp[1].lower().strip()] = tmp[0].strip()
        DOID_to_name_dic[tmp[0].strip()] = tmp[1].lower().strip()
    
    #init vars for disease   
    create_disease_name_dic() # -> init disease_name_dic
    disease_lst = disease_name_dic.keys()

    #init data vars    
    negative_lines_lst = create_line_lst(negative_roc_file)
    rep_score_lines_lst = create_line_lst(rep_score_file)
    spec_score_lines_lst = create_line_lst(spec_score_file)
    smqs_score_lines_lst = create_line_lst(smqs_score_file)
    
    create_genes_name_dic(background_lines_lst,background_genes_name_dic)
    create_genes_name_dic(negative_lines_lst,negative_genes_name_dic)

    background_mat = read_matrix_obj_using_numpy(background_roc_file)
    negative_mat = read_matrix_obj_using_numpy(negative_roc_file)
    rep_score_mat = read_matrix_obj_using_numpy(rep_score_file)
    spec_score_mat = read_matrix_obj_using_numpy(spec_score_file)
    smqs_score_mat = read_matrix_obj_using_numpy(smqs_score_file)
    meta_data_mat = read_matrix_obj_using_numpy(label2metaData_file)
    
    create_dis_reports_dic(meta_data_mat)
   
    parse_disease_ontology()   # create the disease main page graph 
    
    Lawrance_mut = create_line_lst(Lawrance_mut_file)
    for line in Lawrance_mut:
        tmp = line.split("\t")
        for i in range(len(tmp)):
            tmp[i] = tmp[i].strip()
        Lawrance_mut_dic[tmp[0].strip()] = tmp[1:]

    #init edges data
    parse_ppi_intact()
    parse_genemania_gi()
    parse_genemania_ppi()
    #init the gene-drug info
    parse_drugs()
    #init table for the go enrichment
    #create_gene_id_table() #need only for first use in new location, this source update seldom
    gene_to_id_lst = create_line_lst(genes_id_file)
    #init DO terms for disease query
    general_DO = create_general_do()
    tissue_DO = create_tissue_do()
    control_DO = create_control_do()
    

'''
the function parse the ppi intact file
each line parsed into two keys in the dic - for both directions
each key is edge
'''
def parse_ppi_intact():
    lst = create_line_lst(ppi_intact_file)
    for line in lst:
        edge = line.split("\t")
        ppi_intact_dic[(edge[0].strip(),edge[1].strip())] = 1
        ppi_intact_dic[(edge[1].strip(),edge[0].strip())] = 1


'''
the function parse the genemania gi file
each line parsed into two keys in the dic - for both directions
each key is edge
'''
def parse_genemania_gi():
    lst = create_line_lst(genemania_gi_file)
    for line in lst:
        edge = line.split("\t")
        genemania_gi_dic[(edge[0].strip(),edge[1].strip())] = 1
        genemania_gi_dic[(edge[1].strip(),edge[0].strip())] = 1


'''
the function parse the genemania ppi file
each line parsed into two keys in the dic - for both directions
each key is edge
'''
def parse_genemania_ppi():
    lst = create_line_lst(genemania_ppi_file)
    for line in lst:
        edge = line.split("\t")
        genemania_ppi_dic[(edge[0].strip(),edge[1].strip())] = 1
        genemania_ppi_dic[(edge[1].strip(),edge[0].strip())] = 1


'''
the function parse the disease ontology file
each key is tissue, if there is no tissue the key is None, 
the values are list of tuples, first var in tuple is parent, second is son
'''
def parse_disease_ontology():

    disease_is_a_lst = create_line_lst(disease_ontology_file)
    for line in disease_is_a_lst:
        edge = line.split("\t")
        part_one = edge[0].strip().split(";")
        part_two = edge[1].strip().split(";")
        if len(part_one) > 1:
            tissue = part_one[-1]
        elif len(part_two) > 1 :
            tissue = part_two[-1]
        else:
            tissue = None

        name_one = part_one[0].strip().replace("'","")
        name_two = part_two[0].strip().replace("'","")
        
        if tissue in disease_ontology_dic:
            disease_ontology_dic[tissue].append((name_one,name_two))
        else:
            disease_ontology_dic[tissue] = [(name_one,name_two)]



'''
the function parse the drugs file
in this file each line contain target gene and his drugs, if exist
could be more than one drug for gene and deug could be for more then one gene
each key is gene name, and the value is the list of drugs that work on him
'''
def parse_drugs():
    drugs_lines_lst = create_line_lst(drugs_file)
    #for each line in the file
    #we look on columns 1/2. column 0 is line number
    for i in range(1,len(drugs_lines_lst)): #"jump" on the first line - titles
        line = drugs_lines_lst[i].split("\t")
        gene = line[1].strip().upper()
        drug = line[2].strip()
        
        #add/create drug list as the dic value
        if gene in have_drugs_dic:
            have_drugs_dic[gene].append(drug)
        else:
            have_drugs_dic[gene] = [drug]


'''
the function get the query genes list
and return potential edgelist:
for each gene in the list we couple all the other genes.
return all the existance couples.
'''
def create_potential_lst(query_lst):
    lst = []
    #for each gene
    for i in range(len(query_lst)):
        #remove from the query gene spaces and change to upper
        first = query_lst[i].strip().upper()
        #take all the remains genes, as the first gene couple
        for j in range(i+1,len(query_lst)):
            second = query_lst[j].strip().upper()
            lst.append((first,second))
    
    return lst
 

'''
the function get the query potential edge list and the edge type.
return actual edge list - edge that could be in the graph and exist in the data
'''
def create_edges_lst(query_lst,edge_type):

    potential_lst = create_potential_lst(query_lst)
    edges_lst = []
    dic = {}
    #copy the dic of the specific edge type
    if edge_type == 'ppi_intact':
        dic = ppi_intact_dic
    elif edge_type == 'genemania_ppi':
        dic = genemania_ppi_dic
    elif edge_type == 'genemania_gi':
        dic = genemania_gi_dic
    
    #take only the edges that exist in this type
    for edge in potential_lst:
        if edge in dic:
            edges_lst.append(edge)

    return edges_lst


'''
the function get file path and return list of his data
each cell in the list is a line in the file.
'''
def create_line_lst(file_path):
    with open (file_path, 'r') as data:
        lst = data.readlines()
    
    return lst


'''
each disease is the key and the column is the value (begin at 1)
based on background roc file - 
assuming that in the all files there are the same diseases.
'''
def create_disease_name_dic():
    names = background_lines_lst[0].split('\t') #take first line and split to list by tabs
    index = 1 #column number - the 0 is the genes name
    tmp_dic = {}
    
    #create dic of all terms as keys
    for line in names:
        data = line.split(";")
        name = data[0].strip().replace("'","")
        
        if (len(data) > 1):
            tissue = data[1].strip()
        else:
            tissue = None

        tmp_dic[(name,tissue)] = index
        index += 1
    
    #remain the selected terms only
    for line in selected_terms_lst:
        data = line.split(";")
        name = data[0].strip().replace("'","")
        if (len(data) > 1):
            tissue = data[1].strip()
        else:
            tissue = None
        term = (name,tissue)
        
        if term in tmp_dic:
            disease_name_dic[term] = tmp_dic[term]


'''
create genes names dic
each gene is the key and the column is the value
the dic and genes name list are the input, and dic build base on the list.
'''
def create_genes_name_dic(lst,dic):
    index = 0
    
    for line in lst:
        dic[line.split('\t')[0].strip().upper().replace("-","_")] = index
        #note: the replacing of "-" with "_" is because js have problem with elemnts with "-" in their names
        index += 1


'''
the function get an diseae and qVal.
return list of the disease most significant genes - by two conditions:
maximum 200 genes, their q-value is under the qval.
if there last than 10, take 10 highest PN-ROC score
'''
def get_all_signific_disease_genes(disease, tissue, qVal, rep_sig = 0.5, spec_sig = 0.7, bcg_sig = 0.7):
    global background_mat
    global negative_mat
    global rep_score_mat
    global spec_score_mat
    global smqs_score_mat
#    rep_genes_dic = {} #variables of old implemntation
#    spec_genes_dic = {}
#    smqs_genes_dic = {}
#    bcg_genes_dic = {}
#    disease_genes_lst = [] 
    high_limit = 200
    low_limit = 10
    smqs_sig = 0.01  
    #convert the term
    if tissue != None:
        term = disease.strip().replace("'","") + ";" + tissue.strip()
    else:
        term = disease.strip().replace("'","")
    #get scores
    pnroc = get_matrix_column(term,negative_mat)
    pnroc2 = np.maximum(pnroc,np.ones(pnroc.shape)-pnroc)
    pbroc = get_matrix_column(term,background_mat)
    pbroc2 = np.maximum(pbroc,np.ones(pbroc.shape)-pbroc)
    specroc = get_matrix_column(term,spec_score_mat)    
    if specroc is None and re.search("control",term):
        specroc = np.ones(pnroc.shape)
    if specroc is None and not re.search("control",term):
        specroc = np.zeroes(pnroc.shape)+0.5
    specroc2 = np.maximum(specroc,np.ones(specroc.shape)-specroc)
    reps = get_matrix_column(term,rep_score_mat)
    smqs = get_matrix_column(term,smqs_score_mat)

    # first - cut using our thresholds
    selection = np.logical_and(pnroc2 >= bcg_sig , pbroc2 >= bcg_sig)
    selection = np.logical_and(selection,smqs<=smqs_sig)
    selection = np.logical_and(selection,reps>=rep_sig)
    selection = np.logical_and(selection,specroc2 >= spec_sig)

    # add the min and max number filters
    if selection.sum() < low_limit:
        min_thr = sorted(pnroc2)[-low_limit]
        selection = np.logical_or(selection , pnroc2 >= min_thr)
    if selection.sum() > high_limit:
        scores = pnroc2[selection]
        max_thr = sorted(scores)[-high_limit]
        selection = np.logical_and(selection , pnroc2 >= max_thr)

    inds = np.where(selection)
    lst = negative_mat[4][inds]
    final_lst = lst.tolist()
    return final_lst


'''
input: gene list
output: the most significant PNROC score (high or low) genes, up to 200
'''
def get_best_NP_roc(disease_genes_lst):
    high_limit = 100
    disease_genes_lst = sorted(disease_genes_lst, key=lambda x: x[1])
    min_lst = disease_genes_lst[:high_limit]
    length = len(disease_genes_lst)    
    max_lst = [(gene[0],1-gene[1]) for gene in disease_genes_lst[length-high_limit:]]
    com_lst = sorted(min_lst + max_lst)
    return com_lst


#unused function
#'''
#the function get list of lines of data
#split each line by tab 
#and insert the line data as list to list
#return data matrix
#'''
#def create_data_mat(lines_lst):
#    mat = []    
#    for line in lines_lst:
#        mat.append(line.split('\t'))
#    return mat


'''
the function get a genes query list
remove all the genes that not exist in at least one of BCG or NEG datafiles.
convert the exist genes to upper and remove spaces from it
return the edited list.
'''
def remove_unknown_genes(query_lst):
    edit_lst = []

    for gene in query_lst:
        gene = str(gene).strip().upper().replace("-","_")
        #note: the replacing of "-" with "_" is because js have problem with elemnts with "-" in their names
        if ((gene in background_genes_name_dic) and (gene in negative_genes_name_dic)):
            edit_lst.append(gene)
    
    return edit_lst


'''
the function return the expression values of the genes in the BCG and NEG data
'''
def get_exp_values(disease, tissue, gene):
    #if there is no selection disease, 
    #the color will be not gene expression specific   
    if disease == None:
        return 0.5, 0.5
    #otherwise, color by the specific gene expression in the disease
    if tissue != None:
        tissue = ";" + str(tissue)
    if tissue == None:
        tissue = ""
    bcg_val = get_matrix_value(str(disease)+tissue,gene,background_mat)
    neg_val = get_matrix_value(str(disease)+tissue,gene,negative_mat)

    return neg_val,bcg_val


'''
the function return the string part of the for the js node up color definition
based on the expression value and the user ROCscore definition
'''
def define_up_colors(val, side, score):
    low_sig = 1 - score
    high_sig = score
    
    if side == 'right':
        if val >= high_sig:
            return '\', urred: 2.5, urcyan: 0, urgreen: 0, '
        elif val <= low_sig:
            return '\', urred: 0, urcyan: 0, urgreen: 2.5, '
        else:
            return '\', urred: 0, urcyan: 2.5, urgreen: 0, '

    else:
        if val >= high_sig:
            return ', ulred: 2.5, ulcyan: 0, ulgreen: 0} },\n'
        elif val <= low_sig:
            return ', ulred: 0, ulcyan: 0, ulgreen: 2.5} },\n'
        else:
            return ', ulred: 0, ulcyan: 2.5, ulgreen: 0} },\n'


'''
return the down-left node color string part for the js
based on the gene-drug data
'''
def define_dl_colors(gene):
    if gene in have_drugs_dic:
        return ' dlblue: 2.5 , dlwhite: 0'
    
    return ' dlblue: 0 , dlwhite: 2.5'


'''
return the down-right node color string part for the js
based on the mutation data
'''
def define_dr_colors(gene,disease):
    if disease == None:
         return 'drblack: 0, drwhite: 2.5,'
         
    if disease.lower() in name_to_DOID_dic:

        DO = name_to_DOID_dic[disease.lower()]
        if DO in Lawrance_mut_dic:
            if gene in Lawrance_mut_dic[DO]:
                return 'drblack: 2.5, drwhite: 0,'
    
    return 'drblack: 0, drwhite: 2.5,'
                

'''
return converted string for the html page of all the tissue terms in our DO
'''
def create_tissue_do():
    tissue_DO_lst = list(disease_ontology_dic.keys())
    tissue_DO_lst.remove(None)
    tissue_DO_lst.sort()

    result =''
    
    for do in tissue_DO_lst:
        result += '<option value="' + do + '">' + do  + '</option>'

    return result


'''
getter for the tissue DO terms
'''
def get_tissue_do():
    return tissue_DO


'''
return converted string for the html page of all the NONE tissue terms in our DO
'''
def create_general_do():
    lst = disease_ontology_dic[None]
    
    for tup in lst:
        general_DO_set.add(tup[0])
        general_DO_set.add(tup[1])

    general_DO_lst = list(general_DO_set)
    result =''
    
    for do in general_DO_lst:
        result += '<option value="' + do + '">' + do  + '</option>'
    
    return result


'''
getter for the general DO terms
'''
def get_general_do():
    return general_DO


'''
return converted string for the html page of all the control terms in our DO
'''
def create_control_do():
    lst = background_lines_lst[0].split("\t")
    result = ''
    
    for name in lst:
        name = name.split(";")
        if name[0] == "control":
            result += '<option value="' + name[1] + '">' + name[1] + '</option>'
    
    return result


'''
getter for the control DO terms
'''
def get_control_do():
    return control_DO

        
'''
input: disease name
output: the label (tissue/general) that tissue remain of
'''
def find_tissue_per_disease(disease):
    #get from the disease name his tissue
    #by spliting by '(', get the second term and left the term before ')'
    tup = disease.split('(')
    if len(tup) == 1:
        return None, ''
    #otherwise
    else:
        tissue = tup[1][:-1]
        return tissue, tissue
#    keys = list(disease_ontology_dic.keys())
#    
#    for key in keys:
#        #per key=tissue, get all his disease, if our input is in this list, return the key
#        values = [tup[0] for tup in disease_ontology_dic[key]] + [tup[1] for tup in disease_ontology_dic[key]]
#        if disease in values:
#            return key
    
'''
input: lst of genes, and term (disease and tissue)
output: table of all the gene scores for that term
'''
#Importanat notes:
# lables in data are in template <diease>;<tissue> (look for the ; !)
# so we must add it for the name in search
def get_gene_lst_reports(gene_lst, disease, tissue):
    reports = '{'
    if tissue != None:
        term = str(disease) + ";" + str(tissue)
    else:
        term = str(disease)

    if tissue == None and disease == None:
        for gene in gene_lst:
            reports += str(gene) + ': "No data to present for this kind of query",'
    else:
        for gene in gene_lst:
            PNROC = format_pvalue_to_string(get_matrix_value(term, gene, negative_mat))
            PBROC = format_pvalue_to_string(get_matrix_value(term, gene, background_mat))
            rep = get_matrix_value(term, gene, rep_score_mat)
            if rep != None:            
                rep = format_pvalue_to_string(rep)
            else:
                rep = "NULL"
                
            if disease == "control":
                spec = "NULL"
            else:
                spec = format_pvalue_to_string(get_matrix_value(term, gene, spec_score_mat))
            SMQ = get_matrix_value(term, gene, smqs_score_mat)
            if SMQ != None:
                SMQ = format_pvalue_to_string(SMQ)
            else:
                SMQ = "NULL"
                
            table_str = '\"<table>'
            table_str += '\t<tr>'
            table_str += '\t\t<td>' + "" + '</td>'
            table_str += '\t\t<td>' + "Gene Stats" + '</td>'
            table_str += '\t</tr>'
            table_str += '\t<tr>'        
            table_str += '\t\t<td>' + "PNROC" + '</td>'
            table_str += '\t\t<td>' + PNROC + '</td>'
            table_str += '\t</tr>'
            table_str += '\t<tr>'        
            table_str += '\t\t<td>' + "PBROC" + '</td>'
            table_str += '\t\t<td>' + PBROC + '</td>'
            table_str += '\t</tr>'
            table_str += '\t<tr>'        
            table_str += '\t\t<td>' + "Replicability" + '</td>'
            table_str += '\t\t<td>' + rep + '</td>'
            table_str += '\t</tr>'
            table_str += '\t<tr>'        
            table_str += '\t\t<td>' + "Specifity" + '</td>'
            table_str += '\t\t<td>' + spec + '</td>'
            table_str += '\t</tr>'
            table_str += '\t<tr>'        
            table_str += '\t\t<td>' + "Meta-analysis q-value" + '</td>'
            table_str += '\t\t<td>' + SMQ + '</td>'
            table_str += '\t</tr>'
            table_str += '</table>\"\n'
            
            reports += str(gene) + ':' + table_str + ','
    
    reports = reports[:-1] + '}'

    return reports



def get_drugs_report(gene_list):
    
    reports = '{'
    
    for gene in gene_list:
        if gene in have_drugs_dic:
            reports += str(gene) + ':' + '"Drugs targeting the gene (durgbank ids)\: ' + ', '.join(have_drugs_dic[gene]) + '",'
        else:
            reports += str(gene) + ': "",'
    
    reports = reports[:-1] + '}'
    return reports


'''
input: term (disease and tissue)
output: table of the term properties
'''
def get_dis_reports(disease, tissue):
    #build the term key    
    if disease == None or '*' in disease:
        term = '0' #key for None term in the disease dic
    elif tissue == None:
        term = str(disease).replace("-","_")
    else:
        term = str(disease).replace("-","_") + ";" + str(tissue).replace("-","_")
    
    term = term.upper()
    return dis_reports_dic[term]


'''
create for each term (disease and tissue) his propertis table
return: dic with term as key and his properties table as value
'''
def create_dis_reports_dic(mat):
    diseases = list(mat[len(mat) - 1])
    for dis in diseases:
        TotalDatasets = str(int(get_matrix_value('TotalDatasets',dis,mat)))
        RepAnalysisDatasets = str(int(get_matrix_value('RepAnalysisDatasets',str(dis),mat)))
        pos = str(int(get_matrix_value('P',str(dis),mat)))
        neg = str(int(get_matrix_value('N',str(dis),mat)))
        bgc = str(int(get_matrix_value('BGC',str(dis),mat)))

        table_str = '\"<table>'
        table_str += '\t<tr>'
        table_str += '\t\t<td>' + "" + '</td>'
        table_str += '\t\t<td>' + "Disease Stats" + '</td>'
        table_str += '\t</tr>'
        table_str += '\t<tr>'        
        table_str += '\t\t<td>' + "TotalDatasets" + '</td>'
        table_str += '\t\t<td>' + TotalDatasets + '</td>'
        table_str += '\t</tr>'
        table_str += '\t<tr>'        
        table_str += '\t\t<td>' + "RepAnalysisDatasets" + '</td>'
        table_str += '\t\t<td>' + RepAnalysisDatasets + '</td>'
        table_str += '\t</tr>'
        table_str += '\t<tr>'        
        table_str += '\t\t<td>' + "P" + '</td>'
        table_str += '\t\t<td>' + pos + '</td>'
        table_str += '\t</tr>'
        table_str += '\t<tr>'        
        table_str += '\t\t<td>' + "N" + '</td>'
        table_str += '\t\t<td>' + neg + '</td>'
        table_str += '\t</tr>'
        table_str += '\t<tr>'        
        table_str += '\t\t<td>' + "BGC" + '</td>'
        table_str += '\t\t<td>' + bgc + '</td>'
        table_str += '\t</tr>'
        table_str += '</table>\"\n'
        
        dis_reports_dic[dis] = table_str        
    #add empty "table" for None term query, his key is '0'
    dis_reports_dic['0'] = '""'

    return 


'''
input: term (disease, tissue) and gene lst (default is empty lst)
output: js dictionary when each gene is key and the value is his bar plots of his score on different datasets
'''
def get_bar_plots(disease,tissue,lst = []):
    charts_js_dic = '{'
    if tissue == None and disease == None:
        for gene in lst:
            charts_js_dic += str(gene) + ': "undefined",'###
    else:
        mat = read_file_for_bar_plot(disease, tissue)
        if mat == None:
            for gene in lst:
                charts_js_dic += str(gene) + ': " ",'
            charts_js_dic = charts_js_dic + '};'
            return charts_js_dic
            
        terms = mat[0]
        genes = list(mat[len(mat) - 1])
    
        for gene in genes:
    
            chart = '"'
            for term in terms:       
                val = format_pvalue_to_string(-math.log10(get_matrix_value(term,gene,mat)))
                px = str(int(15 * float(val)))
                chart += '<p>' + str(term) + '</p><div style=\'width:' + px + 'px;\'>' + val.replace("-","") + '</div>'
            
            charts_js_dic += str(gene) + ':' + chart +'",'
    
    charts_js_dic = charts_js_dic[:-1] + '};'
    
    return charts_js_dic 
    

'''
return data from term relevant file wich contain the term pvalues of all datasets
'''
def read_file_for_bar_plot(disease, tissue):
    path = r'app/data/pvalue_matrices/'

    disease = disease.lower().replace("*","_star")
    if tissue == None:
        path += disease + '.txt'
    else:
        path += disease + ';' + tissue + '.txt'
    
    return read_matrix_obj_using_numpy(path)


def create_gene_id_table():
    command = ['R','CMD','BATCH']#r'C:\Program Files\R\R-3.3.1\bin\x64\Rscript.exe' #Rscript
    path2script = 'app/gene_id_table.R'
    command.append(path2script)#[command, path2script]
    subprocess.call(command, universal_newlines=True)

#=======================#
#---- GO Enrichment ----#
#=======================#
'''
input: gene
output: his entrez id, if exist. "NA" otherwise.
'''
def gene_to_id(gene):
    
    for line in gene_to_id_lst:
        data = line.split("\t")
        if gene == data[1].strip().upper().replace('"',''):
            return data[2]
    
    return "NA";


'''
convert from ENTREZ id to gene
'''
def id_to_gene(gene_id):
    
    for line in gene_to_id_lst[1:]:
        data = line.split("\t")
        if gene_id == data[2].strip():
            return data[1].strip().upper().replace('"','')


'''
GO enrichment use batch file of Expander
We build input file and batch file that:run the command of the Expander batch file
'''

'''
create from the query gene list input file for the GO batch
'''
def create_go_input_file(query_lst, file_name):
    
    file_text = ''
    file_name = 'app/Expander/input/' + file_name + '.txt'
    for gene in query_lst:
        gene_id = gene_to_id(gene).strip()
        file_text += str(gene_id) + '\t\n'

    with open(file_name,'w') as file:
        file.write(file_text)


'''
create the batch for GO command
'''
def create_go_bat_file(file_name):
    
    bat_start_text = "cd 'app/Expander'\necho 'hello'>check.txt\njava -jar -Xms128m -Xmx128m executable7.0Batch.jar 5 human " #for linux
    bat_mid_text = ' organisms/human/allGenes.txt 1000 3000 0.05 '
    bat_end_text = 'sleep 3\nexit' #'sleep 7\nexit'    
    bat_start_text += 'input/' + file_name + '.txt' + bat_mid_text + 'temp/' + file_name + '_tmp_result.txt\n' + bat_end_text
    
    file_name = 'app/Expander/' + file_name + '.sh'#'.bat'
    
    with open(file_name,'w') as file:
        file.write(bat_start_text)


'''
convert the GO result file
'''
def create_res_file(file_name):
    tmp_res_file = 'app/Expander/temp/' + file_name + '_tmp_result.txt'
    with open(tmp_res_file, 'r') as file_data:
        tmp_res_lst = file_data.readlines()
    
    gene_names = ''
    text = ''
    
    for line in tmp_res_lst:
        split_line = line.split('[')
        if len(split_line) > 1: #for convert only the lines with IDs
            ids_lst = split_line[1].split(', ')
            for gene_id in ids_lst:
                gene_id = gene_id.replace(']','')
                name = id_to_gene(gene_id)
                if name != None:
                    gene_names += id_to_gene(gene_id) + ','
            text += split_line[0] + '[' + gene_names[:-1] + ']\r\n'
        else:
            text += split_line[0] 
    
    os.remove(tmp_res_file)
    
    res_file = 'app/Expander/result/' + file_name + '_tmp_result.txt'
    
    with open(res_file, 'w') as file:
        file.write(text)
    

'''
main GO func, manage all the functions above
'''
def run_go(query_lst):
    #generate file name    
    file_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    
    create_go_input_file(query_lst, file_name)
    create_go_bat_file(file_name)

    run_file = file_name + '.sh'
    path = r'app/Expander/' + run_file
    command = 'chmod 777 ' + path
    
    os.system(command)
    os.system(path)

    #create the res file for being shown in the browser
    create_res_file(file_name)
    return file_name

#============================#
#---- Pathway Enrichment ----#
#============================#
'''
formating the pvalue string, for present only 3 numbers after .
'''
def format_pvalue_to_string(pval):
    if pval>0.001:
        return str('%.3f'%round(float(pval),3))
    return("{:.2E}".format(float(pval)))


'''
input: query gene list
output: pathes to output file and user download file
'''
# Important notes:
#   The function assumes that the working dir is the RScripts dir
#   This diectory should contain all R runnables.
#   This function also assumes that the working dir has an object called bio_pathways.RData
def run_pathway_enrichment_analysis(query_test,
        working_dir = home_dir_1):
    #generate file name
    file_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    #create shell command
    in_file = working_dir + "input/" + file_name
    download = "static/downloads/" + file_name + ".txt"    
    out_file = working_dir + "output/" + file_name
    create_simple_text_file_mutation_analysis(query_test,in_file)
    command = "R --no-save --no-restore --args " + in_file +" NULL " + working_dir + " " + out_file +" kegg <" + working_dir+ "pathway_enrichment_analysis.R"
    #run R script    
    os.system(command)

    copyfile(out_file, "app/" + download)    
    
    os.remove(in_file)
    return out_file, download


'''
create string for html table from the output file
'''
def create_pathway_enrichment_results_table(file_name):
    res_file =  file_name
    with open(res_file, 'r') as file_text:
        text = file_text.readlines()
    table_str = '<table>\n'
    if len(text) <= 1:
        return None
    i = 0
    is_number = [False,False,False,True,True,False]
    for line in text:
        if line != '\n':
            columns = line.split('\t')
            if i==0:
                table_str += '\t<tr>\n'
                table_str += '\t\t<td>' + "Pathway name" + '</td>\n'
                table_str += '\t\t<td>' + "Number of genes" + '</td>\n'
                table_str += '\t\t<td>' + "Pathway size" + '</td>\n'
                table_str += '\t\t<td>' + "p-value" + '</td>\n'
                table_str += '\t\t<td>' + "q-value" + '</td>\n'
                table_str += '\t\t<td>' + "Gene IDs" + '</td>\n'
                table_str += '\t</tr>\n'
                i = i + 1
                continue
            for jj in range(1,len(columns)):
                if is_number[jj]:
                    columns[jj] = format_pvalue_to_string(float(columns[jj]))
            table_str += '\t<tr>\n'
            for col in columns:
                table_str += '\t\t<td>' + col + '</td>\n'
            table_str += '\t</tr>\n'
        i = i + 1
    table_str += '</table>\n'
    os.remove(file_name)
    return table_str


#=================================#
#---- Classification analysis ----#
#=================================#
'''
same logic as the pathway enrichment part
'''

def parse_gene_file_mutation_analysis(file_text):
    file_text = file_text.decode()
    return file_text


def create_simple_text_file_mutation_analysis(file_text, file_name):

    with open(file_name,'w') as file:
        file.write(file_text)


def run_classification_from_input_string(file_text,working_dir = home_dir_2,
                    rscript_path = 'RScripts/',
                    classifier_path = 'data/classifiers/cancer_mutations/',
                    dir_path = 'RScripts/input/', script_name= "run_somatic_mutation_ranger_classifiers_batch.R"):

    rscript_path = working_dir + rscript_path
    classifier_path = working_dir + classifier_path
    dir_path = working_dir + dir_path
    file_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    in_file = dir_path + file_name + "_input.txt"
    create_simple_text_file_mutation_analysis(file_text, in_file)

    download = "static/downloads/" + file_name + "_output.txt"
    out_file = dir_path + file_name + "_output.txt"
    command = "R --no-save --no-restore --args " + in_file + " " + out_file + " " + working_dir +" "+ rscript_path+" " +classifier_path +  " < " + rscript_path + script_name
    run_command(command)
    
    copyfile(out_file, "app/" + download)
    
    os.remove(in_file)
    return out_file, download


'''
allowing to run shell command with os.system or with popen
we build it for tests...
'''
def run_command(command,use_os_system=False):
    if(use_os_system):
        val = os.system(command)
        return (val)
    else:
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while proc.poll() is None:
            proc.stdout.readline() #give output from your execution/your own message
        commandResult = proc.wait() #catch return code
        return(commandResult)


def create_predictions_browser_table(file_name):
    res_file = file_name
    with open(res_file, 'r') as file_text:
        text = file_text.readlines()
    table_str = '<table>\n'
    if len(text) <= 1:
        return None
    gene_list = []
    i = 0
    for line in text:
        if line != '\n':
            columns = line.split('\t')
            if i==0:
                table_str += '\t<tr>\n'
                for colind in range(len(columns)):
                    if colind==0:
                        table_str += '\t\t<td>' + "Disease name" + '</td>\n'
                    else:
                        table_str += '\t\t<td>' + "Sample " + str(colind) + '</td>\n'
                table_str += '\t</tr>\n'
                
            for jj in range(1,len(columns)):
                columns[jj] = str('%.4f'%round(float(columns[jj]),4)) #for check
            table_str += '\t<tr>\n'
            for col in columns:
                table_str += '\t\t<td>' + col + '</td>\n'
            table_str += '\t</tr>\n'
        i = i + 1
    table_str += '</table>\n'
    #os.remove(file_name)
    return table_str

###########################

## Help functions for GO ##

'''
edit the presenting of long gene lists at GO table
'''
def edit_gene_column(lst):
    res = ''    
    if len(lst) > 10:
        lst = lst[:10]

    for i in range(len(lst) - 1):
        res += lst[i] + ", " 
    
    if len(lst) == 10:
        res += lst[len(lst) - 1] + ",..."
    else:
        res += lst[len(lst) - 1]
    
    return res


'''
formating the go table for presenting at the html page
and create the string of the html table
'''
def cerate_go_browser_table(file_name):
    res_file = 'app/Expander/result/' + file_name + '_tmp_result.txt'
    download = '/static/downloads/' + file_name + '_tmp_result.txt'
    copyfile(res_file, "app/" + download)
    
    with open(res_file, 'r') as file_text:
        text = file_text.readlines()
    
    table_str = '<table>\n'

    if len(text) <= 1:
        return None, download

    table_list = []

    for line in text:
        if line != '\n':
            columns = line.split('\t')
            columns.pop(4) 
            if columns[2] != 'p-value':
                columns[2] = format_pvalue_to_string(float(columns[2]))
            if columns[3] != 'Raw Pvalue':
                columns[3] = columns[3][:5] + columns[3][-3:]

            if columns[4] != 'Gene IDs':
                columns[4] = columns[4].replace('[','').replace(']','')
                #save it before cutting off the gene list, for download
                table_list .append(columns)
                
                columns[4] = edit_gene_column(columns[4].split(','))

            columns.pop(0)
            
            table_str += '\t<tr>\n'
            for col in columns:
                table_str += '\t\t<td>' + col + '</td>\n'
            table_str += '\t</tr>\n'
    
    table_str += '</table>\n'
    os.remove(res_file)
    
    return table_str, download
  
#=======================#
#---- DO Enrichment ----#
#=======================#

'''
same logic as the pathway enrichment
'''
# Important notes:
#   The function assumes that the working dir is the RScripts dir
#   This diectory should contain all R runnables
def run_DO(query_lst,
           working_dir = home_dir_1,
           genescores_file = r'app/data/negative_roc.txt',
           labels_file = r'app/data/selected_terms_formatted.txt',
           run_method = "0"):
    #create random file name
    file_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    in_file = working_dir + "input/" + file_name
    out_file = working_dir + "output/" + file_name
    download = "static/downloads/" + file_name + ".txt"
    create_simple_text_file_mutation_analysis(query_lst, in_file)

    # correct the genescores file path for the R command
    # R needs the full path and its running dir is the RScripts dir
    # We go back two directories since all global file names
    # start with 'app'
    genescores_file = working_dir + "../../" + genescores_file
    labels_file = working_dir + "../../" + labels_file

    command = "R --no-save --no-restore --args " + genescores_file + " "
    command += in_file +" " + run_method + " " + working_dir + " "
    command += "-out="+ out_file + " -labels="+labels_file + " "
    command += "<" + working_dir+ "gsea_batch.R"     
    os.system(command)
    copyfile(out_file, "app/" + download)
    table_text = create_do_browser_table(out_file)
    os.remove(in_file)
    return table_text, download
    

def create_do_browser_table(file_name):
    with open(file_name, 'r') as file_text:
        text = file_text.readlines()  
    table_str = '<table>\n'
    if len(text) <= 1:
        return None
    i = 0
    is_number = [False,True,True,False]
    for line in text:
        if line != '\n':
            columns = line.split('\t')
            if i==0:
                table_str += '\t<tr>\n'
                table_str += '\t\t<td>' + "Label name" + '</td>\n'
                table_str += '\t\t<td>' + "p-value" + '</td>\n'
                table_str += '\t\t<td>' + "q-value" + '</td>\n'
                table_str += '\t\t<td>' + "Up/Down regulation" + '</td>\n'
                table_str += '\t</tr>\n'
                i = i + 1
                continue
            for jj in range(1,len(columns)):
                if is_number[jj]:
                    columns[jj] = format_pvalue_to_string(float(columns[jj]))
            table_str += '\t<tr>\n'
            for col in columns:
                table_str += '\t\t<td>' + col + '</td>\n'
            table_str += '\t</tr>\n'
        i = i + 1
    table_str += '</table>\n'
    os.remove(file_name)
    return table_str    
