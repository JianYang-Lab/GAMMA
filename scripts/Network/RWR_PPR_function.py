#!/usr/bin/env python3
#-*- coding:utf-8 -*-

import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random

#### RWR (equal important, probabilities sum to 1) --------

def run_rwr(G, seed_genes, restart_prob, tolerance):
    """
    Run Random Walk with Restart on a given graph.

    Parameters:
    G (networkx.Graph): The graph on which to run RWR.
    seed_genes (list): List of seed nodes.
    restart_prob (float): The probability of restarting at a seed node.
    tolerance (float): The convergence tolerance.

    Returns:
    dict: Node importance scores from RWR.
    """
    # Initialize probabilities
    probabilities = {node: 1.0 / len(seed_genes) if node in seed_genes else 0 for node in G.nodes()}
    # probabilities = {node: 1.0 if node in seed_genes else 0 for node in G.nodes()}  # Same the above, but more iterations.
    
    
    # RWR algorithm
    iteration_count = 0  # Initialize the counter
    convergence = False
    while not convergence:
        iteration_count += 1  # Increment the counter at the start of each iteration
        new_probabilities = probabilities.copy()
        
        for node in G.nodes():
            neighbor_sum = sum(probabilities[neighbor] * G[node][neighbor]['weight'] for neighbor in G.neighbors(node))
            weighted_degree = sum(G[node][neighbor]['weight'] for neighbor in G.neighbors(node))
            new_probabilities[node] = restart_prob * (1.0 / len(seed_genes) if node in seed_genes else 0) + (1 - restart_prob) * neighbor_sum / weighted_degree

        convergence = np.all(np.abs(np.array(list(new_probabilities.values())) - np.array(list(probabilities.values()))) < tolerance)
        probabilities = new_probabilities
        
        # Logging for each iteration
        print(f"Iteration {iteration_count}: Convergence not reached")
   
    # After the loop completes, print the number of iterations
    print(f"Total number of iterations: {iteration_count}")
    
    return probabilities




# Example usage
# G = your graph
# seed_genes = your list of seed genes
# restart_prob = 0.8
# tolerance = 1e-6
# rwr_scores = run_rwr(G, seed_genes, restart_prob, tolerance)


#### Modified RWR (With node importance considered) --------

def run_modified_rwr(G, node_importance, restart_prob, tolerance):
    """
    Run Modified Random Walk with Restart on a given graph, incorporating node importance.

    Parameters:
    G (networkx.Graph): The graph on which to run RWR.
    node_importance (dict): Dictionary of node importance scores.
    restart_prob (float): The probability of restarting at a seed node.
    tolerance (float): The convergence tolerance.

    Returns:
    dict: Node importance scores from RWR.
    """
    # Normalize node importance scores to sum to 1
    total_importance = sum(node_importance.values())
    normalized_importance = {node: score / total_importance for node, score in node_importance.items()}

    # Initialize probabilities
    probabilities = {node: normalized_importance.get(node, 0) for node in G.nodes()}
    
    # Precompute restart component for each node
    restart_component = {node: restart_prob * normalized_importance.get(node, 0) for node in G.nodes()}

    # RWR algorithm
    iteration_count = 0  # Initialize the counter
    convergence = False
    while not convergence:
        iteration_count += 1  # Increment the counter at the start of each iteration
        new_probabilities = probabilities.copy()
        
        for node in G.nodes():
            neighbor_sum = sum(probabilities[neighbor] * G[node][neighbor].get('weight', 1) for neighbor in G.neighbors(node))
            weighted_degree = sum(G[node][neighbor].get('weight', 1) for neighbor in G.neighbors(node))
            new_probabilities[node] = restart_component[node] + (1 - restart_prob) * neighbor_sum / weighted_degree

        convergence = np.all(np.abs(np.array(list(new_probabilities.values())) - np.array(list(probabilities.values()))) < tolerance)
        probabilities = new_probabilities
        
        # Logging for each iteration
        print(f"Iteration {iteration_count}: Convergence not reached")
   
    # After the loop completes, print the number of iterations
    print(f"Total number of iterations: {iteration_count}")
    
    return probabilities

       

# Example usage
# G = your graph
# node_importance = {node: importance_score, ...}  # Define node importance scores
# restart_prob = 0.8
# tolerance = 1e-6
# rwr_scores = run_modified_rwr(G, node_importance, restart_prob, tolerance)




#### PPR (With node importance considered) --------

def run_personalized_pagerank(G, personalization, alpha=0.85, tolerance=1e-6):
    """
    Run Personalized PageRank on a given graph.

    Parameters:
    G (networkx.Graph): The graph on which to run Personalized PageRank.
    personalization (dict): Personalization vector indicating the relative importance of each node.
    alpha (float): Damping factor representing the probability to continue the random walk.
    tolerance (float): The convergence tolerance.

    Returns:
    dict: Node importance scores from Personalized PageRank.
    """
    # Initialize probabilities
    N = len(G)
    if personalization is None:
        personalization = {n: 1 / N for n in G}
    else:
        s = sum(personalization.values())
        personalization = {k: v / s for k, v in personalization.items()}

    pagerank = dict.fromkeys(G, 1.0 / N)
    dangling_weights = personalization

    # For an undirected graph
    dangling_nodes = [n for n in G if G.degree(n, weight='weight') == 0.0]
    # For directed graph
    # dangling_nodes = [n for n in G if G.out_degree(n, weight='weight') == 0.0]
    

    # Personalized PageRank algorithm
    iteration_count = 0
    convergence = False
    while not convergence:
        iteration_count += 1
        pagerank_last = pagerank
        pagerank = dict.fromkeys(pagerank_last.keys(), 0)
        danglesum = alpha * sum(pagerank_last[n] for n in dangling_nodes)

        for n in pagerank:
            for nbr in G.neighbors(n):
                weight = G[n][nbr].get('weight', 1)  # Default weight to 1 if not specified
                pagerank[nbr] += alpha * pagerank_last[n] * weight / G.degree(n, weight='weight')
                
            pagerank[n] += danglesum * dangling_weights.get(n, 0) + (1.0 - alpha) * personalization.get(n, 0)

        err = sum(abs(pagerank[n] - pagerank_last[n]) for n in pagerank)
        convergence = err < N * tolerance

        # Logging for each iteration
        print(f"Iteration {iteration_count}: Convergence not reached")

    # After the loop completes, print the number of iterations
    print(f"Total number of iterations: {iteration_count}")

    return pagerank

# Example usage
# G = your undirected graph
# personalization = {node: score for node, score in zip(gene_ids, gamma_scores)}
# alpha = 0.85  # Commonly used damping factor in PageRank
# tolerance = 1e-6
# pagerank_scores = run_personalized_pagerank(G, personalization, alpha, tolerance)

