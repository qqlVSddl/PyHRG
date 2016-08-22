#!/usr/bin/env python

import networkx as nx
import numpy as np

import optparse
import os
import sys
bindir = os.path.abspath(os.path.dirname(sys.argv[0]))
libdir = os.path.dirname(bindir) + "/lib"
sys.path.append(libdir)
sys.setrecursionlimit(100000)

from hrg import Dendrogram

def print_status(*l):
    print("\t".join([str(x) for x in l]))


def main():
    parser = optparse.OptionParser(
        description='Fits a hierarchical random graph (HRG) model to a network.  Saves the model to a file in graph markup language (GML) format.',
        prog='hrg-fit.py',
        usage='%prog [options] GRAPH_EDGELIST_FILE')

    parser.add_option('-s', '--num-steps', action='store', type=int,
        default=100000,
        help='The number of MCMC steps to take (default=100000).')

    parser.add_option('-t', '--nodetype', action='store', type='choice',
        choices=[int,str],
        default=int,
        help='The type of the nodes in the edgelist file; "int" or "str" (default="int")')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        return 1

    #G=nx.read_edgelist(args[0], nodetype=options.nodetype)
    G = nx.read_gml(args[0], label = 'id')

    D=Dendrogram.from_graph(G)

    bestL=initL=D.graph['L']
    bestI=0

    print_status("step", "L", "best L", "MC step", "deltaL")

    for i in range(1, options.num_steps):
        taken=D.monte_carlo_move()
        t = ''
        if taken:
            t = '*'
        if D.graph['L'] > bestL:
            bestL=D.graph['L']
            bestI=i
            #print_status("["+str(i)+"]", "%.3f" % bestL, "%.3f" % bestL, t, "%.3f"%D.deltaL)
        elif i % 4096 == 0:
            print_status("["+str(i)+"]", "%.3f" % D.graph['L'], "%.3f" % bestL, t, "%.3f"%D.deltaL)

        if i % 10 == 0:
            sys.stdout.flush()
    
    print(bestL)
    print('equilibrium')
    
    #Find equilibrium. You can have different convergence criteria.
    flag_eq = False
    while (not flag_eq):
        old_bestL = bestL
        for i in range(1, 10000):
            D.monte_carlo_move()
            if D.graph['L'] > bestL:
                bestL = D.graph['L']
        print_status("["+str(i)+"]", "%.3f" % D.graph['L'], "%.3f" % bestL, t, "%.3f"%D.deltaL)
        print(bestL - old_bestL)
        if (bestL - old_bestL) < 1:
            flag_eq = True
    
    print('sample')

    #sample dendrograms at regular intervals.
    n = len(G.nodes())
    matrix = np.zeros((n, n))
    for i in range(0, 10):
        tmp = D.save_prob_matrix()
        matrix = matrix + tmp
        print(i)
        for j in range(1, 1000):
            D.monte_carlo_move()
    matrix = matrix / 10   
    
    
    np.savetxt('test.txt', matrix, fmt = '%.6f', newline = '\r\n')
    return 0

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        pass
