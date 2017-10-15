from rdkit import Chem
from rdkit.Chem import AllChem, Draw

import networkx as nx
from eden_chem.io.rdkitutils import nx_to_rdkit



def nx_to_pos(graph):
    '''
    takes nx graph
    returns pos dictionary that is usable by the eden draw functions {nodeid:(X,Y)}
    ## note that I had to scale the pos dict by .1, otherwise pyplot complains
    '''
    chem=nx_to_rdkit(graph)
    chem.UpdatePropertyCache(strict=False)
    AllChem.Compute2DCoords(chem)
    conf = chem.GetConformer(0)
    rawpos = [(conf.GetAtomPosition(i).x,conf.GetAtomPosition(i).y) for i in range(conf.GetNumAtoms())]
    #note that this only works because we dd the nodes in order when making an rdkit object
    pos={n:p for n,p in zip(graph.nodes(),rawpos)}
    return pos



def set_coordinates(chemlist):
    for m in chemlist:
        if m:
            # updateprops fixes "RuntimeError: Pre-condition Violation"
            m.UpdatePropertyCache(strict=False)
            AllChem.Compute2DCoords(m)
        else:
            raise Exception('''set coordinates failed..''')


def get_smiles_strings(graphs):
    compounds = map(nx_to_rdkit, graphs)
    return map(Chem.MolToSmiles, compounds)


def nx_to_image(graphs, n_graphs_per_line=5, size=250, title_key=None, titles=None):
    # we want a list of graphs
    if isinstance(graphs, nx.Graph):
        raise Exception("give me a list of graphs")
    # make molecule objects
    compounds = map(nx_to_rdkit, graphs)
    # print compounds

    # take care of the subtitle of each graph
    if title_key:
        legend = [g.graph.get(title_key, 'N/A') for g in graphs]
    elif titles:
        legend = titles
    else:
        legend = map(str, range(len(graphs)))
    return compounds_to_image(compounds, n_graphs_per_line=n_graphs_per_line, size=size, legend=legend)


def compounds_to_image(compounds, n_graphs_per_line=5, size=250, legend=None):
    # calculate coordinates:
    set_coordinates(compounds)
    # make the image
    return Draw.MolsToGridImage(compounds, molsPerRow=n_graphs_per_line, subImgSize=(size, size), legends=legend)
