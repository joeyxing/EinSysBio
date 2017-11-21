import networkx as NX
import pylab as P

def makegraphfromfile(nodesfilename, edgesfilename):
    atoms = {}

    atomcolors={}
    atomcolors['C']='#909090'
    atomcolors['O']='#FF0D0D'
    atomcolors['H']='#FFFFFF'
    
    atomsizes={}
    atomsizes['C']=1000
    atomsizes['O']=400
    atomsizes['H']=200

    G=NX.Graph()
    #G=nx.DiGraph()

    f=open(nodesfilename)
    lines = f.readlines()
    f.close()
    for line in lines:
        [nr, atomtype] = line.split()
        if not atomtype in atoms.keys():
            atoms[atomtype] = []
        G.add_node(nr)
        atoms[atomtype].append(nr)

    f=open(edgesfilename)
    lines = f.readlines()
    f.close()
    for line in lines:
        inds = line.split()
        G.add_edge(inds[0],inds[1])
    return G, atoms, atomcolors, atomsizes

G, atoms, atomcolors, atomsizes=makegraphfromfile('fenol2.dat', 'fenol1.dat')

fig=P.figure()
pos=NX.spring_layout(G, iterations=5000) 
for atomtype in atoms.keys():
        NX.draw_networkx_nodes(G,pos,nodelist=atoms[atomtype],
                               node_color=atomcolors[atomtype],
                               node_size=atomsizes[atomtype])
NX.draw_networkx_edges(G,pos)

labels = {}
for node in G.nodes():
    labels[node] = node
NX.draw_networkx_labels(G,pos,labels,font_color='black',
                            font_family='sans-serif')
    
P.show()
    

def find_all_paths(nw, e, path, G):
    "path is a list representing a path from s to nw"

    npaths = []

    path = path + [nw]
    print "+++", path
    if nw == e:
       print "We have succes!!!", path
       return [path]
    for node in G.neighbors(nw):
        if node not in path:
           newpaths = find_all_paths(node, e, path, G)
           for newpath in newpaths:
               npaths.append(newpath) 
    return npaths
    
print find_all_paths('1', '10', [], G)
