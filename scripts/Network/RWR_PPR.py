#!/usr/bin/env python3
#-*- coding:utf-8 -*-
import sys,os
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
from RWR_PPR_function import run_rwr,run_modified_rwr,run_personalized_pagerank


trait_name = sys.argv[1]
OUTPUT = sys.argv[2]
PPI_file = sys.argv[3]

# Step 1: Load Data ---------------
# Load your gene scores
GAMMA_V2G_scores = pd.read_csv(OUTPUT + '/V2G/score/' + trait_name + '_GAMMA_V2G.summary', sep='\t', header=0, low_memory=False)
GAMMA_L2G_scores = pd.read_csv(OUTPUT + '/L2G/score/' + trait_name + '_GAMMA_xQTL.summary', sep='\t', header=0, low_memory=False)

gene_scores_df = pd.merge(GAMMA_V2G_scores[['gene_id', 'gene_name', 'GAMMA_V2G']],
                          GAMMA_L2G_scores[['gene_id', 'gene_name', 'GAMMA_xQTL']],
                          on=['gene_id', 'gene_name'],
                          how='outer')  

gene_scores_df['GAMMA'] = gene_scores_df['GAMMA_V2G'].fillna(0) + gene_scores_df['GAMMA_xQTL'].fillna(0)
print(gene_scores_df.head())


# Step 2: Prepare the Network ---------------
# Load PPI data ---------------
ppi_df = pd.read_csv(PPI_file, sep='\t', header = 0)

G = nx.Graph()
for _, row in ppi_df.iterrows():
    G.add_edge(row['gene_1'], row['gene_2'], weight=row['evidence_score'])
# G = nx.from_pandas_edgelist(ppi_df, 'gene_1', 'gene_2', edge_attr='max_evidence_score')



#### Step 3: Run RWR and PPR analysis
# Step 3_1: Run RWR analysis ---------------
# convergence = False
# restart_prob = 0.8
# tolerance = 1e-6

GAMMA_threshold = 0
restart_prob = 0.8
tolerance = 1e-10

seed_genes = gene_scores_df[gene_scores_df['GAMMA'] > GAMMA_threshold]['gene_id'].tolist()
probabilities = run_rwr(G, seed_genes, restart_prob, tolerance)


# Add the RWR scores to the gene_scores_df DataFrame
gene_scores_df['RWR_Score'] = gene_scores_df['gene_id'].map(probabilities).fillna(0)
# gene_scores_df.sort_values(by='RWR_Score', ascending=False, inplace=True)

# Step 4_2: Run modified RWR analysis ---------------
# Set the restart probability and tolerance
# restart_prob = 0.8
# tolerance = 1e-6
# Run the modified RWR

node_importance = {row['gene_id']: row['GAMMA'] for _, row in gene_scores_df.iterrows() if row['gene_id'] in G}
rwr_scores = run_modified_rwr(G, node_importance, restart_prob, tolerance)

# Add the PageRank scores to the gene_scores_df DataFrame
gene_scores_df['modified_RWR_Score'] = gene_scores_df['gene_id'].map(rwr_scores).fillna(0)
# gene_scores_df.sort_values(by='modified_RWR_Score', ascending=False, inplace=True)


# Step 4_3: Run modified RWR analysis ---------------
# Create a personalization vector

personalization = {row['gene_id']: row['GAMMA'] for _, row in gene_scores_df.iterrows() if row['gene_id'] in G}
# Set the alpha and tolerance values
# alpha = 0.8
# tolerance = 1e-6

# Run Personalized PageRank
pagerank_scores = run_personalized_pagerank(G, personalization, restart_prob, tolerance)

# Add the PageRank scores to the gene_scores_df DataFrame
gene_scores_df['PPR_Score'] = gene_scores_df['gene_id'].map(pagerank_scores).fillna(0)
gene_scores_df.sort_values(by='PPR_Score', ascending=False, inplace=True)

# Output the results
# print(gene_scores_df)

gene_scores_df.to_csv(OUTPUT+"/RWR_PPR/summary/"+trait_name+"_"+str(GAMMA_threshold)+"_RWR_PPR.txt", sep="\t", index=False)




