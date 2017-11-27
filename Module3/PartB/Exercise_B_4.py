import itertools
import numpy
import os
import re

import gurobipy

import libsbml
import scipy.sparse as sparse

import recon2_2PHfructoseanaerobic as phlib


def read_block_list(path=''):
    '''
    read file return blocked reaction list
    '''
    with open(path) as f:
        block_list = f.readlines()

    for i,x in enumerate(block_list):
        block_list[i] = x.strip()
    return block_list


def block_reactions(sbml, rl):
    '''
    set kineticLaw of blocked reactions to 0 so that block
    the reaction on both side
    '''
    model = sbml.getModel()

    for rID in rl:
        reaction = model.getReaction(rID)
        kineticLaw = reaction.getKineticLaw()
        kineticLaw.getParameter('LOWER_BOUND').setValue(0)
        kineticLaw.getParameter('UPPER_BOUND').setValue(0)


def my_max_fluxes(sbml):
    '''
    Written to mimic neilswainston matlab function maxFluxes
    '''
    eps=1e-10

    objective = 'DM_atp_c_'
    print ''

    carbon_sources=[
                'EX_fru(e)',
        ]
    model = sbml.getModel()
    listOfNodes=[]
    listOfEdges=[]
    for normoxic in [False]:
        for carbon_source in carbon_sources:
            v, f_opt = my_max_flux(sbml, carbon_source, objective, normoxic, [])
            print '%s (%s):\t%g' % (carbon_source,
                                    'normoxic' if normoxic else 'anaerobic',
                                    f_opt)
            for j, reaction in enumerate(model.getListOfReactions()):
               if abs(v[j])>1e-10:
                  rID = reaction.getId()
                  txt, txt_formula, innodes, outnodes=phlib.mydisplay_reaction_and_formula_Nodes(rID, sbml, v[j]>0)
                  reactionnode=phlib.format_for_ID_SBML(rID)
                  if reactionnode not in listOfNodes:
                     listOfNodes.append(reactionnode)
                  for n in innodes+outnodes:
                      if not n in listOfNodes: listOfNodes.append(n)
                  for n in innodes:
                      if not (n, reactionnode, abs(v[j])) in listOfEdges:
                         listOfEdges.append((n, reactionnode, str(abs(v[j]))))
                  for n in outnodes:
                      if not (reactionnode, n) in listOfEdges:
                         listOfEdges.append((reactionnode, n, str(abs(v[j]))))
                  print '%s\t%s\t%g\t[%s,\t%s]' % ('reaction', reactionnode, abs(v[j]), txt, txt_formula)
    nodef=open('fructoseanaerobicfluxgraph2.txt', 'w')
    nodef.write("\n".join(listOfNodes))
    nodef.write("\n")
    nodef.write("\n".join([";".join(r) for r in listOfEdges]))
    nodef.close()

    return (listOfNodes, listOfEdges)


def my_max_flux(sbml, carbon_sourcel, objective, normoxic, media):
    '''
    Written to mimic neilswainston matlab function maxFlux
    '''
    #set_infinite_bounds(sbml)
    # block import reactions
    phlib.block_all_imports(sbml)

    bl = read_block_list(path='block')
    block_reactions(sbml, bl)

    # define carbon source
    phlib.set_import_bounds(sbml, carbon_sourcel, 50)
    # define media
    phlib.set_import_bounds(sbml, media, phlib.INF)
    if normoxic:
        phlib.set_import_bounds(sbml, 'EX_o2(e)', phlib.INF)
    # specify objective and maximise
    phlib.change_objective(sbml, objective)
    # avoid infinities
    obj_max = 1e6
    phlib.change_rxn_bounds(sbml, objective, obj_max, 'u')
    v, f_opt = phlib.optimize_cobra_model(sbml)
    if f_opt > 0.9 * obj_max:
        f_opt = phlib.INF
    return v, f_opt


def main():
    model_filename = 'models/recon_2.2.xml'
    sbml = phlib.read_sbml(model_filename)
    block_list = read_block_list(path='block')
    my_max_fluxes(sbml)

if __name__ == '__main__':
    main()
