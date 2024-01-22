import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path
import requests
# from decouple import AutoConfig # Install python-decouple
import requests # Install requests
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as matplt
import pyvis
from pyvis import network as net
from pyvis.network import Network
import networkx as nx
# import SessionState
# from streamlit import caching

from copy import deepcopy
from community import community_louvain
from collections import Counter

from itertools import tee, combinations

from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_object_dtype,
)

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return tuple(zip(a, b))

# config = AutoConfig(".env") # Create a file called .env file in the same directory.
#                             # This is just a text file that contains the BEARER TOKEN so that we don't 
#                             # Have to include it in the code.
#                             # It will have one line like the following (exclude the angle brackets):
#                             # BEARER_TOKEN=<MY_BEARER_TOKEN>
                
BASE = "https://api-ibis.rilm.org/200/haverford/"

# SECRET_TOKEN =  st.secrets["db_password"]



URLS = {
    "year": BASE + "rilm_index_RYs",
    "terms": BASE + "rilm_index_top_terms",
    "index": BASE + "rilm_index",
    "biocards": "https://ibis2.rilm.org/api/bio_cards/musicologist_places"
}

# HEADERS = {
#     "Authorization": st.secrets["SECRET_TOKEN"]
# }
BEARER_TOKEN = st.secrets["SECRET_TOKEN"]
HEADERS = {
    "Authorization": f"Bearer {BEARER_TOKEN}"
}
# st.cache speeds things up by holding data in cache

# the functions
@st.cache_data()
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return tuple(zip(a, b))

def simple_search(term):
    """
    Returns results for a single term.
    term : string
    """
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
    simple_search_results = pd.DataFrame(data)
    simple_search_results = simple_search_results.fillna('')
    # combines year and accession number to make unique id for each item
    simple_search_results['full_acc'] = simple_search_results.ry.apply(str) + "-"  + simple_search_results.ac.apply(str)
    simple_search_results.rename(columns = {'ry':'year', 
                              'ac': 'item', 
                              'ent' : 'entry', 
                              'lvl': 'level', 
                              'name': 'term', 
                              'cat': 'category', 
                              'full_acc': 'full_id'}, inplace=True)
    return simple_search_results

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

def _get_one_author_terms_and_values(author_name, results):
    """
    Helper function for graph of ONE author for given set of results.
    selected_authors : list 
    results : dataframe
    """
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
    """
    helper function for multi-author graph
    """
    pairs = author_terms_grouped['term'].apply(lambda x: list(combinations(x, 2)))
    unique_pairs = pairs.explode().dropna().unique()
    final_pairs = _clean_pairs(unique_pairs)
    return final_pairs

def add_communities(G):
    """
    adds Louvain Community Detection to graphs
    """
    G = deepcopy(G)
    partition = community_louvain.best_partition(G)
    nx.set_node_attributes(G, partition, "group")
    return G

# Initialize size and color

def one_author_graph(selected_author, results):
    author_terms_and_values = _get_one_author_terms_and_values(selected_author, results)
    author_terms = author_terms_and_values[1]
    author_term_values = author_terms_and_values[2]
    author_terms_grouped = author_terms_and_values[0]
    final_pairs = _get_pairs(author_terms_grouped)
    # size must be in px.  Should also match streamlit widget with components below
    pyvis_graph = Network(width="1000px",
                          height="1000px",
                          bgcolor="black", 
                          font_color="white")
    # Set the physics layout of the network
    pyvis_graph.set_options("""
    {
    "physics": {
    "enabled": true,
    "forceAtlas2Based": {
        "springLength": 1
    },
    "solver": "forceAtlas2Based"
    }
    }
    """)
    # correct node label size
   
    G = nx.Graph()
    # add nodes, and sizes, one at a time
    for node in author_terms:
        G.add_node(node, 
                   size=author_term_values[node],
                   font={"size": 40})
    # add the edges
        G.add_edges_from(final_pairs)
    # visualize with pyvis
    G = add_communities(G) 
    author_graph = pyvis_graph.from_nx(G)
    author_graph = pyvis_graph.show("author_graph.html")
    return author_graph
    
def _author_pairs(results, author_impact_ratio):
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
    """ Helper function for author community """
    G = deepcopy(G)
    partition = community_louvain.best_partition(G)
    nx.set_node_attributes(G, partition, "group")
    return G

# Initialize size and color
def graph_author_communities(results, author_impact_ratio=.3):
    """ Louvain Community Detection network for a group of authors allied with a given term """
    top_author_pairs_and_names = _author_pairs(results, author_impact_ratio)
    top_authors = top_author_pairs_and_names[0]
    top_author_pairs = top_author_pairs_and_names[1]
    author_ratios = top_author_pairs_and_names[2]
    # size must be in px.  Should also match streamlit widget with components below
    pyvis_graph = Network(width="1000px", 
                          height="1000px", 
                          bgcolor="black", 
                          font_color="white") 
    # Set the physics layout of the network
    pyvis_graph.set_options("""
    {
    "physics": {
    "enabled": true,
    "forceAtlas2Based": {
        "springLength": 1
    },
    "solver": "forceAtlas2Based"
    }
    }
    """)
    G = nx.Graph()
    # add nodes, and sizes, one at a time
    for node in top_authors:
        G.add_node(node, 
                   size=(author_ratios[node]),
                   font={"size": 10})
    # add the edges
        G.add_edges_from(top_author_pairs)
    # visualize with pyvis
    G = _add_author_communities(G)
    pyvis_graph.from_nx(G)
    author_community_graph = pyvis_graph.show('author_community_graph.html')
    return author_community_graph

def get_query_data(search_term):
    """ 
    Returns the results of an API query for the given search term 
    search_term : the given single RILM term
    """
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
    
    clean_pair_list = []
    for pair in pairs_list:
        unique_set = set(pair)
        if len(unique_set) > 1:
            clean_pair_list.append(pair)
        
    # get edge weights and unique nodes
    for i, p in enumerate(clean_pair_list):
        clean_pair_list[i] = tuple(sorted(p))
        
    if weight_threshold == 0:
        weighted_plist = [[elem, count] for elem, count, in Counter(clean_pair_list).items() if count >= weight_threshold]
        nodes = results['term'].unique()
    else:
        nodes = set()
        weighted_plist = []
        for ele, count in Counter(clean_pair_list).items():
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
    # size must be in px.  Should also match streamlit widget with components below
    pyvis_graph = Network(width="1000px",
                          height="1000px",
                          bgcolor="black", 
                          font_color="white")
    
    # Set the physics layout of the network
    pyvis_graph.set_options("""
    {
    "physics": {
    "enabled": true,
    "forceAtlas2Based": {
        "springLength": 1
    },
    "solver": "forceAtlas2Based"
    }
    }
    """)
    G = nx.Graph()
    
    for name, info in nodes_dict.items():
        years = f"years: {*info[1],}"
        
        G.add_node(name, value=info[2], group=info[0], title=years)
        
    for pair, weight in weighted_plist:

        G.add_edge(pair[0], pair[1], value=weight, title=str(weight))
    pyvis_graph.from_nx(G)
    concept_map = pyvis_graph.show('concept_graph.html')
    # return the network
    return concept_map


def term_hist_px(cleaned_df, num_terms=20):
    """ Creates, shows, and returns a histogram showing the number of times each term appears in the DataFrame
    
    @param cleaned_df: the cleaned DataFrame to count the term occurences in
    @param num_terms: the number of terms to show on the histogram
    @return: the Figure object of the histogram
    """
    counts = dict(cleaned_df['term'].value_counts())
    counts_frame = pd.DataFrame({'term' : counts.keys(), 'occurences' : counts.values()})
    term_hist = px.histogram(counts_frame.head(num_terms), 
                           x='term', y='occurences', 
                           title='Distribution of Terms',
                          labels = {'term': "Search Term", 'occurences': "Number of Occurences"})

    return term_hist

def scatter_plot(final_results, term_threshold=10, legend=False):
    
    """Creates and shows a scatterplot of terms, years and authors from the given dataframe.
    The terms are the Y axis, and the years are the X axis
    The size of the marker shows the number of times the term appears in that year
    The color indicates the author (as shown in legend)
    The hover data reveals author 
    
    @ param final_results:  the dataframe of RILM search results
    @ term_threshold:  the minimum number of occurences of the given term
    """
    filtered_df = final_results.groupby('term').filter(lambda x: len(x) > term_threshold)
    filtered_df = filtered_df.fillna('-')
    value_counts = filtered_df.groupby(['term', 'year']).size()

    filtered_df['marker_size'] = filtered_df.apply(lambda row: value_counts.get((row['term'], row['year']), 0), axis=1)

    fig_scatter = px.scatter(filtered_df,
                             x = 'year', y = 'term',
                            hover_data = ['author'],
                            color='author',
                            labels = {'term': "Term", "author": "Author", "year": "Publication Year"},
                            size='marker_size',
                            height=750)
    # fig_scatter.update_layout(height=800) 
    fig_scatter.update_layout(showlegend=legend) 
    fig_scatter.update_yaxes(categoryorder='category descending')
    # fig_scatter.show()
    return fig_scatter

# User Enters Search Term
st.header("RILM Data Visualizations")
st.write("Enter a search term in the box below.  Then use the options at left to see histograms and scatterplots of terms.  You can also view networks of terms and authors.  Use the Filter option to limit your results in various ways.")
st.header("Subject Search")
search_term = st.text_input("Enter Search Term")
if len(search_term) == 0:
    st.write("Please enter a search term")
    if 'simple_search_results' not in st.session_state:
        st.session_state.simple_search_results = pd.DataFrame()
    # st.session_state.simple_search_results = simple_search_results
else:
    simple_search_results = simple_search(search_term)
    # convert year to integer
    simple_search_results['year'] = simple_search_results['year'].astype(int)
    simple_search_results = simple_search_results.drop(['langTransFrom'], axis=1)
    # display the term
    st.subheader("Search Results for " + "'" + search_term + "'")
    length = simple_search_results['full_id'].nunique()
    st.write("Your simple search returned " + str(length) + " unique RILM items")
    if 'simple_search_results' not in st.session_state:
        st.session_state.simple_search_results = pd.DataFrame()
    st.session_state.simple_search_results = simple_search_results

# st.subheader("Search Results")
if st.session_state.simple_search_results is not None and not st.session_state.simple_search_results.empty:
    st.subheader("Filter Your Search")
    filter_option = st.checkbox("Filter Your Search")
    if filter_option:
        with st.form('my_form'):
            min_year, max_year = int(st.session_state.simple_search_results['year'].min()), int(st.session_state.simple_search_results['year'].max())
            # selected_year = st.slider('Select Year Range', min_year, max_year)
            selected_years = st.slider('Select Year Range', min_year, max_year, (min_year, max_year))

            min_level, max_level = int(st.session_state.simple_search_results['level'].min()), int(st.session_state.simple_search_results['level'].max())
            # selected_year = st.slider('Select Year Range', min_year, max_year)
            selected_levels = st.slider('Select RILM Index Level Range', min_level, max_level, (min_level, max_level))

            unique_categories = st.session_state.simple_search_results['category'].unique().tolist()
            selected_categories = st.multiselect('Select Categories', unique_categories, default=unique_categories)
            
            unique_authors = st.session_state.simple_search_results['author'].unique().tolist()
            selected_authors = st.multiselect('Select Authors', unique_authors, default=unique_authors)

            unique_terms = st.session_state.simple_search_results['term'].unique().tolist()
            # selected_terms = st.multiselect('Select terms', unique_terms, default=unique_terms)

            submitted = st.form_submit_button('Apply Filters')

            if submitted:
                filtered_results = st.session_state.simple_search_results[(st.session_state.simple_search_results['year'] >= selected_years[0]) 
                                                        & (st.session_state.simple_search_results['year'] <= selected_years[1]) & 
                                (st.session_state.simple_search_results['level'].isin(selected_levels)) & 
                                (st.session_state.simple_search_results['author'].isin(selected_authors)) &
                                # (simple_search_results['term'].isin(selected_terms)) &
                                (st.session_state.simple_search_results['category'].isin(selected_categories))
                                ]
                if 'filtered_results' not in st.session_state:
                    st.session_state.filtered_results = pd.DataFrame()
                st.session_state.filtered_results = filtered_results

                # st.dataframe(st.session_state.filtered_results)

                st.write("Your filtered search has " + str(st.session_state.filtered_results['full_id'].nunique()) + " unique RILM items")


if st.sidebar.checkbox("Show Histogram of Results by Term"):
    num_terms = st.sidebar.slider('Adjust Number of Terms in Histogram', 
                                    min_value=5, 
                                    max_value=25, 
                                    value=10, 
                                    step=1)
    st.plotly_chart(term_hist_px(st.session_state.filtered_results, num_terms=num_terms))
# 
if st.sidebar.checkbox("Show Scatterplot of Results by Term"):
    # check max count of any term, which is used in slider for scatter plot
    counts = st.session_state.filtered_results['term'].value_counts()
    most_common_value = counts.idxmax()
    max_count = counts.max()
    # create slider for plot threshold of term counts
    count_threshold = st.sidebar.slider('Adjust Threshold for Minimum Number of Results for Each Term in Chart', 
                                        min_value=5, 
                                        max_value=int(max_count), 
                                        value=round(max_count/5), 
                                        step=10)
    st.plotly_chart(scatter_plot(st.session_state.filtered_results, term_threshold=count_threshold))

# concept map
if st.sidebar.checkbox("Show Concept Map for Given Search Results"):

    weight_threshold = st.sidebar.slider('Adjust Weight Treshold for Concept Map', 
                                    min_value=0, 
                                    max_value=10, 
                                    value=5, 
                                    step=1)
    # call the function
    concept_map  = create_concept_map(st.session_state.filtered_results, weight_threshold=weight_threshold)
        # Display HTML in Streamlit size of components should match network function above
    with open("concept_graph.html", "r") as f:
        concept_html_string = f.read()
        components.html(concept_html_string, height=1000, width=1000)

# now single author graph
if st.sidebar.checkbox("Show Single Author Graph"):
    author_list = st.session_state.filtered_results['author'].unique()
    selected_author = st.sidebar.selectbox('Select an author', author_list)
    if len(selected_author) == 0:
        st.sidebar.write("Please select an author from the list of results")
    else:
    # Call the function
        author_graph = one_author_graph(selected_author, st.session_state.filtered_results)
    # Display HTML in Streamlit
    with open("author_graph.html", "r") as f:
        author_html_string = f.read()
    # "st." not needed with components
    # here the height must be integer!
    # Display HTML in Streamlit size of components should match network function above
    if len(selected_author) > 1:
        components.html(author_html_string, height=1000, width=1000)


if st.sidebar.checkbox("Show Multi-Author Graph"):
    author_list = st.session_state.filtered_results['author'].unique()
    # selected_authors = st.sidebar.multiselect('Select authors', author_list)
    author_impact_ratio = st.sidebar.slider('Adjust Author Impact Ratio', 
                                    min_value=.1, 
                                    max_value=float(10), 
                                    value=.3, 
                                    step=.1)
    # Call the function
    author_community_graph = graph_author_communities(st.session_state.filtered_results, author_impact_ratio)
    # Display HTML in Streamlit

    with open("author_community_graph.html", "r") as f:
        html_string = f.read()
    # "st." not needed with components
    # here the height must be integer!
    # Display HTML in Streamlit size of components should match network function above
        components.html(html_string, height=1000, width=1000)

    
