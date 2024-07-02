# display pyvis graph with streamlit

import streamlit as st 
import pyvis.network as net
import streamlit.components.v1 as components
# from stvis import pv_static

import os
# from decouple import AutoConfig # Install python-decouple
import requests
import pandas as pd
import plotly as plt
import matplotlib.pyplot as matplt
import pyvis
# from pyvis import network as net
from pyvis.network import Network
import networkx as nx

from copy import deepcopy
from community import community_louvain
from collections import Counter

from itertools import tee, combinations

# config = AutoConfig(".env") # Create a file called .env file in the same directory.
#                             # This is just a text file that contains the BEARER TOKEN so that we don't 
#                             # Have to include it in the code.
#                             # It will have one line like the following (exclude the angle brackets):
#                             # BEARER_TOKEN=<MY_BEARER_TOKEN>
                
BASE = "https://api-ibis.rilm.org/200/haverford/"

# # don't share the following!
BEARER_TOKEN="47E96109-DC23-41E1-8295-377667F587DC"

# BEARER_TOKEN = st.secrets['tk']
# /Users/rfreedma/Documents/CRIM_Python/RILM/rilm_streamlit_dev

URLS = {
    "year": BASE + "rilm_index_RYs",
    "terms": BASE + "rilm_index_top_terms",
    "index": BASE + "rilm_index"
}

HEADERS = {
    "Authorization": f"Bearer {BEARER_TOKEN}"
}

# a function to deal with pairs
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return tuple(zip(a, b))

def simple_search(term):
    
    params = {
            "termName": search_term,
            "includeAuthors": True
        }

    # and get the response
    response = requests.get(
        URLS["index"], 
        headers=HEADERS, 
        params=params
    )

    # get the data
    data = response.json()
    results = pd.DataFrame(data)
    results = results.fillna('')
    # combines year and accession number to make unique id for each item
    results['full_acc'] = results.ry.apply(str) + "-"  + results.ac.apply(str)
    results.rename(columns = {'ry':'year', 
                              'ac': 'item', 
                              'ent' : 'entry', 
                              'lvl': 'level', 
                              'name': 'term', 
                              'cat': 'category', 
                              'full_acc': 'full_id'}, inplace=True)
    return results

# This function removes pairs of terms that are just the same term 2x, 
# or ones that are reverses of each other
# it's used below
def _clean_pairs(list_of_pairs): 
    pairs_no_reverse = [] 
    pairs_cleaned = []
    for item in list_of_pairs: 
        if set(item) not in pairs_no_reverse: 
            pairs_no_reverse.append(set(item))
    for pair in pairs_no_reverse:
        if len(set(pair)) > 1:
            pairs_cleaned.append(pair)
    pairs_cleaned = [tuple(s) for s in pairs_cleaned]
    return pairs_cleaned

# this function gets the top terms for the given author, along with counts and groups
# the returned values are then used to create nodes (the terms) and the edges (the pairs of terms)
def _get_author_terms_and_values(author_name, results):
    # select the author
    selected_results = results[results['author'] == author_name]
    # a dictionary of the nodes and their counts
    term_node_values = selected_results['term'].value_counts().to_dict()
    # limiting the dictionary to counts above X (5, for instance), in order to avoid a dense graph
    top_term_node_values = {key: value for (key, value) in term_node_values.items() if value > 2 }
    # and just the keys of that subset
    top_keys = top_term_node_values.keys()
    # now narrow the results so that we only see the top keys for our author
    top_keys_for_author = selected_results[selected_results['term'].isin(top_keys)]
    # and group them according to the bibliographical item and terem
    author_terms_grouped = top_keys_for_author.groupby(['full_id'])['term'].apply(list).reset_index()
    return author_terms_grouped, top_keys, top_term_node_values

def _get_pairs(author_terms_grouped): 
    pairs = author_terms_grouped['term'].apply(lambda x: list(combinations(x, 2)))
    unique_pairs = pairs.explode().dropna().unique()
    final_pairs = _clean_pairs(unique_pairs)
    return final_pairs

def add_communities(G):
    G = deepcopy(G)
    partition = community_louvain.best_partition(G)
    nx.set_node_attributes(G, partition, "group")
    return G

# Initialize size and color

def one_author_graph(author_name, results):
    author_terms_and_values = _get_author_terms_and_values(author_name, results)
    author_terms = author_terms_and_values[1]
    author_term_values = author_terms_and_values[2]
    author_terms_grouped = author_terms_and_values[0]
    final_pairs = _get_pairs(author_terms_grouped)

    pyvis_graph = Network(bgcolor="black", 
                          font_color="white")
    G = nx.Graph()

    # add nodes, and sizes, one at a time
    for node in author_terms:
        G.add_node(node, size=author_term_values[node])
    # add the edges
        G.add_edges_from(final_pairs)
    # visualize with pyvis
    G = add_communities(G)
    pyvis_graph.from_nx(G)
    pyvis_graph.save_graph('graph.html')
    HtmlFile = open("graph.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    with open('graph.html', "r") as file:
        btn = st.download_button(
                label="Download graph",
                data=file,
                file_name="graph.html")
    components.html(source_code,height=1500, width=1500)
    

def _author_pairs(results, author_impact_ratio):
#     author_impact_ratio = author_impact_ratio # <== put your ratio here .3 is about right to start.
    # a dictionary of the nodes and their counts
    author_node_values = results['author'].value_counts().to_dict()
    # find number of unique writings
    number_unique_ids = results['full_id'].nunique()
    # dictionary that tells us how many items each author wrote
    items_per_author = results.groupby(['author'])['full_id'].nunique().to_dict()
    # finding the 'top' authors, according to author_impact_ratio set above
    top_authors = []
    author_ratios = {}
    for author in items_per_author.items():
        author_ratio = (author[1]/number_unique_ids)*100
        author_ratios.update({author[0]:author_ratio})
        if author_ratio > author_impact_ratio:
            top_authors.append(author[0])
    # now narrow the results so that we only see the top keys for our author
    items_for_top_authors = results[results['author'].isin(top_authors)]
    # and group them:  for each 'term' in the results, make a list of authors connected with that term
    authors_by_term = items_for_top_authors.groupby(['term'])['author'].apply(list).reset_index()
    # create pairs for edges
    pairs = authors_by_term['author'].apply(lambda x: list(combinations(x, 2)))
    unique_pairs = pairs.explode().dropna().unique()
    final_author_pairs = _clean_pairs(unique_pairs)
    return top_authors, final_author_pairs, author_ratios
    
    
def _add_author_communities(G):
    G = deepcopy(G)
    partition = community_louvain.best_partition(G)
    nx.set_node_attributes(G, partition, "group")
    return G

# Initialize size and color
def graph_author_communities(results, author_impact_ratio, graph_name):
    top_author_pairs_and_names = _author_pairs(results, author_impact_ratio)
    top_authors = top_author_pairs_and_names[0]
    top_author_pairs = top_author_pairs_and_names[1]
    author_ratios = top_author_pairs_and_names[2]
    pyvis_graph = Network(notebook=True, width="900", 
                          height="900", 
                          bgcolor="black", 
                          font_color="white")
    G = nx.Graph()
    # add nodes, and sizes, one at a time
    for node in top_authors:
        G.add_node(node, size=(author_ratios[node]*25))
    # add the edges
        G.add_edges_from(top_author_pairs)
    # add communities
    G = _add_author_communities(G)
    pyvis_graph.from_nx(G)

    display(pyvis_graph.show(graph_name + "_graph.html"))    


def get_query_data(search_term):
    """ Returns the results of an API query for the given search term """
    # query the API
    params = {
        "termName": search_term,
        "includeAuthors": True
    }

    # and get the response
    response = requests.get(
        URLS["index"], 
        headers=HEADERS, 
        params=params
    )

    # get the data
    data = response.json()
    results = pd.DataFrame(data)
    results = results.fillna('')
    if len(results) > 0:
    # # combines year and accession number to make unique id for each item
        results['full_acc'] = results.ry.apply(str) + "-"  + results.ac.apply(str)
        results.rename(columns = {'ry':'year', 'ac': 'item', 'ent' : 'entry', 'lvl': 'level', 'name': 'term', 'cat': 'category', 'full_acc': 'full_id'}, inplace=True)
        return results
    else:
        return(print("SORRY! There were no results for the folowing term: " + search_term))


def clean_query_data(results, year_list, categories):
    """
    Cleans the query results for a given a search term, list of years, and list of categories
    results : pandas dataframe containing the results from the API query
    year : list of ints
    category : list of strings
    """
    # parse results for corresponding entries
    if year_list is not None:
        results = results[results['year'].isin(year_list)]
    if categories is not None:
        # results = results[results['category'] == category]
        results = results[results['category'].isin(categories)]
    results = results.drop_duplicates(['term', 'full_id'])
    return results


def create_concept_map(results, weight_threshold=1):
    """
    Creates a concept map given cleaned query results
    results : pandas dataframe containing the results from the API query cleaned for a given a search term, list of 
                years, and list of categories
    """
    
    # get dictionary with key=full_id, value=list of unique terms
    terms_dict = {}

    past_id = results.iloc[0]['full_id']
    terms_list = []
    for index, row in results.iterrows():
        curr_id = row['full_id']
        if curr_id == past_id:
            terms_list.append(row['term'])
        else:
            terms_dict[past_id] = terms_list
            past_id = curr_id
            terms_list = [row['term']]
    terms_dict[past_id] = terms_list

    # get list of all combinations of pairs for each entry
    pairs_list = []
    for key, value in terms_dict.items():
        pairs_list += list(combinations(value, 2))
        
    # get edge weights and unique nodes
    for i, p in enumerate(pairs_list):
        pairs_list[i] = tuple(sorted(p))
        
    if weight_threshold == 0:
        weighted_plist = [[elem, count] for elem, count, in Counter(pairs_list).items() if count >= weight_threshold]
        nodes = results['term'].unique()
    else:
        nodes = set()
        weighted_plist = []
        for ele, count in Counter(pairs_list).items():
            if count >= weight_threshold:
                weighted_plist.append([ele, count])
                if weight_threshold > 0:
                    nodes.add(ele[0])
                    nodes.add(ele[1])
                    
    # get the information about each unique node [category, list of years, number of years]
    nodes_dict = {}
    for node in nodes:
        node_info = []
        node_info.append(results[results['term'] == node]['category'].unique()[0])
        node_info.append(results[results['term'] == node]['year'].unique())
        node_info.append(len(node_info[1]))
        nodes_dict[node] = node_info

    # create network
    G = nx.Graph()
    cmap = Network(notebook=True, width=1000, height = 1000)
    
    for name, info in nodes_dict.items():
        years = f"years: {*info[1],}"
        
        G.add_node(name, value=info[2], group=info[0], title=years)
        
    for pair, weight in weighted_plist:

        G.add_edge(pair[0], pair[1], value=weight, title=str(weight))
    cmap.from_nx(G)
    return cmap

def term_hist(cleaned_df, num_terms=5):
    """ Creates, shows, and returns a histogram showing the number of times each term appears in the DataFrame
    
    @param cleaned_df: the cleaned DataFrame to count the term occurences in
    @param num_terms: the number of terms to show on the histogram
    @return: the AxesSubplot object of the histogram
    """
    counts = dict(cleaned_df['term'].value_counts())
    counts_frame = pd.DataFrame({'term' : counts.keys(), 'occurences' : counts.values()})
    my_plot = counts_frame.head(num_terms).plot(x='term', 
                                                y='occurences', 
                                                kind='bar', 
                                                ylabel='occurences', 
                                                xlabel='term', 
                                                title='Number of Occurences of Terms', 
                                                legend=False)
    matplt.show()
    return my_plot

# Create a basic graph

# st.setFrameHeight(1000)
st.set_page_config(layout="centered")
st.title("Richard's new App")

# height = st.sidebar.slider("Height", 200, 1500, 300, 100)
# width = st.sidebar.slider("Width", 200, 1500, 600, 100)

search_term = "Renaissance"
results = simple_search(search_term)
# select author
author_name = "Vendrix, Philippe"
st.subheader("Graph of " + author_name + " and " + search_term)
out = one_author_graph(author_name, results)
# with open('graph.html', "r") as file:
#     btn = st.download_button(
#             label="Download graph",
#             data=file,
#             file_name="graph.html")



# components.html(f'Frühauf, Tina graph.html', height=1000, width=600)