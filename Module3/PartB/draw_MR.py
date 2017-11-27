#!/home/joey/Apps/SysBio2/bin/python
''' Author: Weizhou Xing
    System biology module 3
    Part B: draw M/R graph from description txt file
    Usage: python draw_MR filename.txt
'''


import networkx as NX
import matplotlib.pyplot as plt
import sys

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
    plt.show(fig)
    return G


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


def main(path='fructoseanaerobicfluxgraph.txt'):
    # # unconmment below to construct graph by computing again
    # import recon2_2PHfructoseanaerobic as phlib
    # tests_path = os.path.dirname(__file__)
    # model_path = os.path.join(tests_path, 'models')
    # model_path = os.path.normpath(model_path)
    # name = "recon_2.2"
    # print name
    # model_filename = os.path.join(model_path, name + '.xml')
    # sbml = phlib.read_sbml(model_filename)

    # # nodes list, edges list
    # nl, el = phlib.max_fluxes(sbml)

    nl, el = read_MR_txt(path)
    G = draw_MR_Graph(nl, el)


if __name__ == '__main__':
    main(sys.argv[1])
