# -*- coding: utf-8 -*-

from flask import render_template, redirect, request, session
from app import app
from .forms import LoginForm
from .edit import *

import tempfile
import pickle
import shutil

app.debug = True

#default forwarding to home page
@app.route('/')
@app.route('/index')
def index():
    return redirect("/home")


#---- main pages ----#
@app.route('/home')
def home():
    global dis_query #declare on global variable, this var is used on the graph pages, but need to declered bafore
    return render_template('home.html')


@app.route('/tutorials')
def tutorials():    
    return render_template('tutorials.html')
    

@app.route('/downloads')
def downloads():    
    return render_template('downloads.html')


@app.route('/contact') #==contact us page
def contact():    
    return render_template('contact.html')


@app.route('/profile_analysis_selection')
def profile_analysis_selection():    
    return render_template('profile_analysis_selection.html')


#---- query pages ----#
j = 0   #graph countr, to avoid diffrent gephes with the same name
@app.route('/geneList', methods=['GET', 'POST'])
def geneList():
    global query
    
    form = LoginForm()
    error = ''
    
    if form.validate_on_submit():   #if user submit his query
        session['var'] = form.openid.data #get query list
        if len(session['var']) == 0: #if no user input at the table, check for file
            f = request.files['upload']
            text = f.read()
            if text == b'': #if no user upload file, return same page with error massege
                error = 'Please insert a gene list or upload file'
                return render_template('geneList.html', 
                           title='Main',
                           form=form,
                           error = error,
                           providers=app.config['OPENID_PROVIDERS'])
            #otherwise, parse file into list
            session['var'] = parse_gene_file(text)
        
        #the length of the seesion string is to big - we split him to a list of string
        session['split_var'] = split_session(session['var'])
        session['disease'] = None  #define no disease query
        session['score'] = 0.7    #define ROC score
        session['edge_type'] = 'genemania_ppi' #default edge type
        
        # create graph nodes part and get edges list
        nodes, lst = create_nodes(session['split_var'], session['disease'], None, session['score'])       
        # check for user choosen option
        if request.form['btn'] == "View Gene Network":
            #create edges part for th graph
            edges, edges_lst = create_edges(lst, session['edge_type'])  
            query = nodes + edges

            global j
            j += 1
            return redirect('/graph/' + str(j))

        if request.form['btn'] == "Gene Ontology Enrichment":
            session['var'] = lst
            return redirect('/tango')
        
        if request.form['btn'] == "Pathway Enrichment":
            session['var'] = lst
            return redirect('/pathway_enrichement_result')
        
        if request.form['btn'] == "Disease Ontology Enrichment":
            session['var'] = lst            
            return redirect('/DO')


    #otherwise return the genelist page
    return render_template('geneList.html', 
                           title='Main',
                           form=form,
					  error = error,
                           providers=app.config['OPENID_PROVIDERS'])


@app.route('/disease', methods=['GET', 'POST'])
def disease():
    form = LoginForm()
    #get all choose options for the right bar radio
    general = get_general_do()
    tissue = get_tissue_do()
    control = get_control_do()
        
    if form.validate_on_submit():
        #create disease graph - for user selection
        # check for radio user selection
        radio = request.form.get('analyze')
        if radio == "tissue":
            session["tissue"] = request.form.get('selTissue')
        elif radio == "control":
            session["tissue"] = request.form.get('selControl')
        else:
            session["tissue"] = None
        
        if request.form['btn'] == "View Gene Network":
            return redirect('/graph_dis/control/' + session["tissue"])

        selected = session["tissue"]
        if selected == None:
            selected = " for general selection (none tissue)" 

        elements = create_disease_onto_graph(radio, session["tissue"])
        session["elements"] = elements

        return render_template('disease.html', 
                           title='Main',
                           form=form,
                           general = general,
                           selected = selected,
                           tissue = tissue,
                           control = control,
                           elements = elements,
                           con = 1,
                           providers=app.config['OPENID_PROVIDERS'])
    #create default - all disease - network
    elements = create_disease_onto_graph('default',None)
    selected = ''
    return render_template('disease.html', 
                           title='Main',
                           form=form,
                           general = general,
                           selected = selected,
                           tissue = tissue,
                           control = control,
                           elements = elements,
                           con = 0,
                           providers=app.config['OPENID_PROVIDERS'])     


@app.route('/profile', methods=['GET', 'POST'])
def profile():
    form = LoginForm()
    error = ''
    
    if form.validate_on_submit():   #if user submit his query
        f = request.files['upload']
        text = f.read()        
        if text == b'':
            error = 'You must insert gene expression profile file'
            return render_template('profile.html',
                          form=form,
                          error = error,
                          providers=app.config['OPENID_PROVIDERS'])

        return redirect('profileGraph')
        
    return render_template('profile.html',
                           form=form,
                           error = error,
                           providers=app.config['OPENID_PROVIDERS'])


@app.route('/profile_cancer',methods=['GET', 'POST'])
def profile_cancer():    
    form = LoginForm()
    error = ''
    if form.validate_on_submit():   #if user submits his query
        
        session['var'] = form.openid.data #get query list
        if len(session['var']) == 0: #if there is no query on input table check for upload file
            f = request.files['upload']
            text = f.read()
            if text == b'':#if there is no upload file, return to the same page and show input error
                error = 'Please insert a gene list or upload file'
                return render_template('profile_cancer.html', 
                           title='Main',
                           form=form,
                           error = error,
                           providers=app.config['OPENID_PROVIDERS'])      
            #else, parse the file to list
            session['var'] = parse_gene_file_mutation_analysis(text) 
        #the length of the session string is to big - we split him to a list of string
        session['split_var'] = split_session(session['var'])
        
        if request.form['btn'] == "Run classification":
            return redirect('/profile_cancer_results')     
    #otherwise(==no submition == default page) return the profile cancer page
    return render_template('profile_cancer.html', 
                           title='Main',form=form,error = error,
                           providers=app.config['OPENID_PROVIDERS'])


@app.route('/profile_expression',methods=['GET', 'POST'])
def profile_expression():    
    form = LoginForm()
    error = ''
    if form.validate_on_submit():   #if user submits his query

        f = request.files['upload']        
        text = f.read()
        if text == b'':#check if no input file
            error = 'Please upload a file'
            return render_template('profile_expression.html', 
                    title='Main',form=form,error = error)      
        #if there is user upload file, save his data in tmp file
        session['tempdir'] = tempfile.mkdtemp()
        outfile = open(session['tempdir'] + '/filename', 'wb')
        pickle.dump(text, outfile)
        outfile.close()        

        if request.form['btn'] == "Run classification":
            return redirect('/profile_expression_results')     
    #otherwise return the default page
    return render_template('profile_expression.html', 
                           title='Main',form=form,error = error)
    

#---- query result pages ----# 
@app.route('/graph_dis/<disease>/<tissue>', methods=['GET','POST'])
def graph_dis(disease, tissue):
    form = LoginForm()
    try:
        os.remove(session["remove"])
    finally:
        # define default values
        qVal = '0.01'
        edge_type = 'genemania_ppi'
    
        session['sel_enrich'] = request.form.get('select_enrich')
        session['roc_score'] = request.form.get('PNroc')
        if session['roc_score'] != None:
            session['roc_score'] = float(session['roc_score'])
        else:
            session['roc_score'] = 0.7
        
        session['spec_score'] = request.form.get('specifity')
        if session['spec_score'] != None:
            session['spec_score'] = float(session['spec_score'])
        else:
            session['spec_score'] = 0.7
    
        session['rep_score'] = request.form.get('replicability')
        if session['rep_score'] != None:
            session['rep_score'] = float(session['rep_score'])
        else:
            session['rep_score'] = 0.5
    
        # get the user disease selection
        session['disease'] = disease
        # if the user select 'disease' node - return the current page
        if disease == 'disease':
            return redirect('/disease')
    
        #get the the user definitions
        session['qVal'] = request.form.get('qVal')    
        session['edge_type'] = request.form.get('edge_type') 
        # if there are changes - update
        if session['qVal'] != None and session['qVal'] != qVal:
            qVal = session['qVal'] 
                    
        if session['edge_type'] != None and session['edge_type'] != edge_type:
            edge_type = session['edge_type']
        
        # create the new graph and return
        txt_tissue = tissue
        if tissue == " for general selection (none tissue)":
            tissue = None
            txt_tissue = ''
    
        nodes, lst = create_nodes('disease', disease, tissue, qVal, session['rep_score'], session['spec_score'], session['roc_score'])    
        session['var'] = lst
        
        edges, edges_lst = create_edges(lst, edge_type) 
        dis_query = nodes + edges

        #extra data
        gene_reports = get_gene_lst_reports(session['var'], disease, tissue)
        dis_reports = get_dis_reports(disease, tissue)
        charts = get_bar_plots(disease,tissue)
        
        #GO enrichment option
        if form.validate_on_submit():
            if request.form['btn'] == "GO Enrichment":
                if session['sel_enrich'] != '':
                    session['var'] = str(session['sel_enrich']).split(",")
                return redirect('/tango')
        
        #user download option
        file_name = '/static/downloads/' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8)) + ".txt"
        
        with open('app/' + file_name,'w') as file_text:
            text = "Nodes:\n" + "\n".join(lst) + "\n\nEdges\n" + "\n".join(map(str,edges_lst))
            file_text.write(text)
            session["remove"] = 'app/' + file_name
            
        return render_template('graph.html',
                               title='Graph Result', 
                               result = dis_query, 
                               form = form, 
                               gene_reports = gene_reports,
                               dis_reports = dis_reports,
                               charts = charts,
                               disease = disease, 
                               tissue = txt_tissue,
                               download = file_name,
                               rep_val = session['rep_score'],
                               spec_val = session['spec_score'],
                               roc_val = session['roc_score'],)#, select_part = select_part)


@app.route('/graph/<i>', methods=['GET', 'POST'])
def graph(i):

    form = LoginForm()
    try:
        os.remove(session["remove"])
    finally:
        # define default values    
        qVal = '0.01'
        edge_type = 'genemania_ppi'
    
        session['sel_enrich'] = request.form.get('select_enrich')
        session['roc_score'] = 0.7        
        session['spec_score'] = 0.7    
        session['rep_score'] = 0.5    
        session['disease'] = None
        session['qVal'] = request.form.get('qVal')    

        #get the the user definitions
        session['edge_type'] = request.form.get('edge_type') 

        # if there are changes - update
        if session['qVal'] != None and session['qVal'] != qVal:
            qVal = session['qVal'] 

        if session['edge_type'] != None and session['edge_type'] != edge_type:
            edge_type = session['edge_type']

        # create the new graph and return
        disease = None
        tissue = None
        txt_tissue = ''
    
        nodes, lst = create_nodes(session['split_var'], session['disease'], tissue, qVal, session['rep_score'], session['spec_score'], session['roc_score'])    
        session['var'] = lst
        
        gene_reports = get_gene_lst_reports(session['var'], disease, tissue)
        charts = get_bar_plots(None, None, session['var'])
        edges, edges_lst = create_edges(lst, edge_type) 
        dis_query = nodes + edges
        
        #GO enrichment
        if form.validate_on_submit():#GO enrichment
            if request.form['btn'] == "GO Enrichment":
                if session['sel_enrich'] != '':
                    session['var'] = str(session['sel_enrich']).split(",")
                return redirect('/tango')
        
        #data (genes&edges) file for user download
        file_name = '/static/downloads/' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8)) + ".txt"
        
        with open('app/' + file_name,'w') as file_text:
            text = "Nodes:\n" + "\n".join(lst) + "\n\nEdges\n" + "\n".join(map(str,edges_lst))
            file_text.write(text)
            session["remove"] = 'app/' + file_name
            
        return render_template('graph.html',
                               title='Graph Result', 
                               result = dis_query, 
                               form = form, 
                               gene_reports = gene_reports,
                               dis_reports = '" ";',
                               charts = charts,
                               disease = disease, 
                               tissue = txt_tissue,
                               download = file_name,
                               rep_val = session['rep_score'],
                               spec_val = session['spec_score'],
                               roc_val = session['roc_score'],)#, select_part = select_part)


@app.route('/profile_cancer_results')
def profile_cancer_results():
    try:#clean session last file, if exist
        os.remove(session["remove"])
    finally:
        lst = session['var']
        file_name, download = run_classification_from_input_string(lst) #run the R scripts
        #note: the downloat var is path for the output copy file, for download file user option
        table = create_predictions_browser_table(file_name) #create results output table
        if table == None:
            table = "No prediction was found."
        
        os.remove(file_name)
        session["remove"] = 'app/' + download #save file name, for future dwlwtion.
        return render_template('profile_cancer_results.html', table = table, download = download)


@app.route('/profile_expression_results')
def profile_expression_results():
    try:
        os.remove(session["remove"])
    finally:
        #get data from tmp file and remove him
        infile = open(session['tempdir'] + '/filename', 'rb')
        lst = pickle.load(infile)
        infile.close()
        shutil.rmtree(session['tempdir'])
        session.pop('tempdir', None)
        #from here - same as profile_cancer_results()    
        file_name, download = run_classification_from_input_string(lst, script_name= "run_gene_expression_classifiers_batch.R", classifier_path = 'data/classifiers/gene_expression/')
        table = create_predictions_browser_table(file_name)
        if table == None:
            table = "No prediction was found."
        os.remove(file_name)        
        session["remove"] = 'app/' + download
        return render_template('profile_expression_results.html', table = table, download = download)


@app.route('/tango')
def tango():
    # same as profile_cancer_results()
    try:
        os.remove(session["remove"])
    finally:
        lst = remove_unknown_genes(session['var']) 
        file_name = run_go(lst)
        table, download = cerate_go_browser_table(file_name)
        session["remove"] = "app/" + download
        if table == None:
            table = "No significant enrichment was found."
            
        return render_template('tango.html', table = table, download = download)


@app.route('/pathway_enrichement_result')
def pathway_enrichement_result():
    # same as profile_cancer_results()
    try:
        os.remove(session["remove"])
    finally:
        file_name, download = run_pathway_enrichment_analysis(" ".join(session['var']))
        table = create_pathway_enrichment_results_table(file_name)
        session["remove"] = "app/" + download        
        if table == None:
            table = "No significant enrichment was found."
        return render_template('pathway_enrichement_result.html', table = table, download = download)


@app.route('/DO')
def do_page():
    # same as profile_cancer_results()
    try:
        os.remove(session["remove"])
    finally:
        lst = remove_unknown_genes(session['var']) 
        table, download = run_DO(" ".join(lst))
        session["remove"] = "app/" + download        
        if table == None:
            table = "No DO enrichment was found."

    return render_template('do.html', table = table, download = download)

### to remove?
@app.route('/profileGraph', methods=['GET', 'POST'])
def profileGraph():
    form = LoginForm()
    #create disease graph - for user selection
    elements = create_disease_onto_graph()

    return render_template('profileGraph.html', 
                    title='Main',
                    form=form,
                    elements = elements,
                    providers=app.config['OPENID_PROVIDERS'])


#---- View Help Function ----#
def split_session(param):
    length = len(param)
    split_var = int(length/130) + 1
    chunk_size = int(length/split_var)
    lst = [ param[i:i+chunk_size] for i in range(0, length, chunk_size) ]
   
    return lst
