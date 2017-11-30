'''
System Biology Module 3 Part B: Exercise B.2 B.3
Author: Weizhou Xing
please read Exercise_B23.ipynb jupyter notebook for more detailed report.
'''

import networkx as NX
import numpy as np
import matplotlib.pyplot as plt



def read_MR_txt(path='fructoseanaerobicfluxgraph.txt'):
    with open(path) as f:
        txt = f.read()
    lines = txt.split('\n')
    nl = []
    el = []

    for line in lines:
        ll = line.split(';')
        if len(ll)==1:
            nl = nl + ll
        elif len(ll)==3:
            el.append(ll)

    return nl, el


def draw_MR_Graph(nl, el):
    '''draw metabolite/reaction graph
    nl: node list
    el: edge list
    '''
    G = NX.Graph()

    rNodes = []
    mNodes = []
    for n in nl:
        G.add_node(n)
        if n[0] == 'R':
            rNodes.append(n)
        else:
            mNodes.append(n)

    nodes = {}
    nodes['R'] = rNodes
    nodes['M'] = mNodes

    nodeColor = {'R':'#0DFF0D', 'M':'#FF0D0D'}

    for e in el:
        if e[0][0] == 'M' and e[1][0] == 'R':
            direction = 'reactant' # M -> R
        else:
            direction = 'product' # R -> M

        G.add_edge(e[0], e[1],
                   weight=float(e[2]),
                   direction=direction)

    fig = plt.figure()
    pos = NX.spring_layout(G, iterations=5000)
    for nodetype in ['R', 'M']:
        for node in nodes[nodetype]:
            total_weight = 0
            for neighbour in G.neighbors(node):
                edge_data = G.get_edge_data(node, neighbour)
                total_weight = total_weight + edge_data['weight']

            NX.draw_networkx_nodes(G, pos,
                                   nodelist=[node],
                                   node_color=nodeColor[nodetype],
                                   node_size=total_weight)

    NX.draw_networkx_edges(G, pos)
    labels = {}
    for node in G.nodes():
        labels[node] = node
    NX.draw_networkx_labels(G, pos, labels,
                            font_color='black',
                            font_family='sans-serif',
                            font_size=9)
    fig.suptitle('Fructose Anaerobic')
    plt.show()
    return G


def my_dijkstra(source, G, weight=1):
    # initialize Green and Red set
    Green = [source]
    Red = list(G.nodes())
    Red.remove(source)
    dist = {}

    for n in Red:
        #  infinity
        dist[n] = float('inf')
    dist[source] = 0

    while Red:
        # min_dist = float('inf')
        for n in G.neighbors(Green[-1]):
            if n in Red:
                if isinstance(weight, int):
                    cost = weight
                elif isinstance(weight, str):
                    edge_data = G.get_edge_data(Green[-1], n)
                    cost = edge_data[weight]

                if dist[n] > dist[Green[-1]] + cost:
                    dist[n] = dist[Green[-1]] + cost

        x = min(Red, key=lambda x: dist[x])
        # choose minimal in Red set
        Green.append(x)
        Red.remove(x)

    return dist


def degree_stats(G):
    total_degree = 0
    numN = 0
    degree_distr = {}
    for t in G.adjacency():
        degree = len(t[1])
        if degree_distr.has_key(degree):
            degree_distr[degree] += 1
        else:
            degree_distr[degree] = 1
        total_degree += len(t[1])
        numN += 1

    draw_bar(degree_distr)
    mean_degree = float(total_degree) / numN
    print "Average degree", mean_degree

    # calculate variance
    var = 0.0
    for t in G.adjacency():
        degree = len(t[1])
        var += (degree - mean_degree)**2
    var = var / numN
    print "variance:", var

    adj = G.adj
    s = sorted(adj, key=lambda x: len(adj[x]))
    hubs = list()
    for n in s:
        hubs.append((n, len(adj[n])))
    print hubs
    winner = hubs[len(hubs)-1]
    print "Node", winner[0], "has the greatest degree:", winner[1]


def draw_bar(distr):
    N = max(distr.keys()) + 1
    for i in xrange(N):
        if not distr.has_key(i):
            distr[i] = 0
    # change dict into array
    y = []
    for x in xrange(N):
        y.append(distr[x])

    fig, ax = plt.subplots()
    rects = ax.bar(np.arange(N), y, 0.5, color="b")
    ax.set_ylabel("number of occurences")
    ax.set_xlabel("degree")
    ax.set_title("Degree Distribution")
    plt.show()


def compartment_analysis(G):
    freq_table = {}
    for node in G.nodes():
        # a metablite
        if node[0] == 'M':
            compartment = node[-1]

        else:
            d = {}
            for n in G.neighbors(node):
                if d.has_key(n[-1]):
                    d[n[-1]] += 1
                else:
                    d[n[-1]] = 1
            # compartment of reaction
            compartment = max(d, key=lambda x:d[x])

        if freq_table.has_key(compartment):
            freq_table[compartment].append(node)
        else:
            freq_table[compartment] = [node]

    fig, ax = plt.subplots()
    ind = np.arange(len(freq_table))
    rects = ax.bar(ind, [len(x) for x in freq_table.values()], 0.5, color="r")
    ax.set_xticks(ind)
    ax.set_xticklabels(freq_table.keys())

    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')
    ax.set_ylabel("number of R/M")
    ax.set_xlabel("compartment")
    ax.set_title("Compartment analysis")

    return freq_table


if __name__ = '__main__':
    # exercise B.2
    txt = "fructoseanaerobicfluxgraph.txt"
    nl, el = read_MR_txt(path=txt)
    G = draw_MR_Graph(nl, el)

    d1 = my_dijkstra('R_DM_atp_c_', G, weight='weight')
    d2 = my_dijkstra('R_DM_atp_c_', G, weight=1)

    print 'Source node R_DM_atp_c_ weighted distance:'
    for k,v in d1.items():
        print k,':', v

    print '\nSource node R_DM_atp_c_ unweighted distance:'
    for k,v in d2.items():
        print k, ':', v

    # use the dijkstra within networkx to compare with my result
    for node in G.nodes():
        print "R_DM_atp_c_ ->", node, "(unweighted):",\
            NX.algorithms.dijkstra_path_length(G, "R_DM_atp_c_", node, weight=1)
        print "R_DM_atp_c_ ->", node, "(weighted):",\
            NX.algorithms.dijkstra_path_length(G, "R_DM_atp_c_", node)
        print ""

    # graph degree statistics
    degree_stats(G)

    # exercise B.3 frequency table
    freq_table = compartment_analysis(G)
    node_num = 0
    for k,v in freq_table.items():
        print k,":", v, "\n"
        node_num += len(v)

    print node_num == len(G.nodes())
