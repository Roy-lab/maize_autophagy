
# source activate scanpy
import networkx as nx
import matplotlib.pyplot as plt
import igraph
import numpy as np
import pandas as pd
import argparse

# input is a weighted edge list


def load_data(expr_name):
    # This is a weighted edge list
    G = nx.read_edgelist(expr_name, nodetype=str, data=(('weight',float),))
    return G


def get_cells(G):
    # This is a weighted edge list
    cell_counts=G.number_of_nodes()
    return cell_counts


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('expr', help="Network input")
    args = parser.parse_args()


    #print('Loading edge list')
    G = load_data(args.expr)


    cell_counts=get_cells(G)
    print(f'The number of nodes are {cell_counts}')


if __name__ == '__main__':
    main()

