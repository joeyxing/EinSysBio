'''This code has been adapted from recon_2.2_test.py
   (see Metabolomics (2016) 12:109
        DOI 10.1007/s11306-016-1051-4
        Recon 2.2: from reconstruction to model of human metabolism)
   Peter Hilbers
'''


import itertools
import numpy
import os
import re

import gurobipy

import libsbml
import scipy.sparse as sparse


INF = float('inf')
NAN = float('nan')


def run_model():
    '''Runs some tests on the recon_2.2 model that is
       stored in the folder d.models'''
    tests_path = os.path.dirname(__file__)
    model_path = os.path.join(tests_path, 'models')
    model_path = os.path.normpath(model_path)
    name="recon_2.2"
    print name
    model_filename = os.path.join(model_path, name + '.xml')
    sbml = read_sbml(model_filename)
    model_stats(sbml, display_errors=True)
    model_balancing(sbml, display_errors=True)
    max_fluxes(sbml)


def model_stats(sbml, display_errors=False):
    '''Statistics on the number, and type, of species and reactions
    '''
    model = sbml.getModel()

    # species statistics by SBO term
    print '\n%g\t%s' % (model.getNumSpecies(), 'species')

    # species statistics by type
    appears_in_reaction = []
    for reaction in model.getListOfReactions():
        for reactant in itertools.chain(reaction.getListOfReactants(),
                                        reaction.getListOfProducts()):
            sID = reactant.getSpecies()
            if sID not in appears_in_reaction:
                appears_in_reaction.append(sID)
    nM, nB, nE = 0, 0, 0
    for species in model.getListOfSpecies():
        if species.getId() not in appears_in_reaction:
            nE += 1
        elif species.getBoundaryCondition():
            nB += 1
            #print species.getId(), species.getName()
        else:
            nM += 1
    print '%g\t%s' % (nM, 'variable')
    print '%g\t%s' % (nB, 'fixed')
    print '%g\t%s' % (nE, 'non-reactants')

    print '\n%g\t%s' % (model.getNumReactions(), 'reactions')

    nB = 0
    rID_list = get_source_reactions(sbml)
    print '%g\t%s' % (model.getNumReactions() - len(rID_list),
                      'non-source/sink')
    print '%g\t%s' % (len(rID_list), 'source/sink')



def model_balancing(sbml, display_errors=False):
    '''
    Checks elemental balancing for all reactions in model
    '''
    model = sbml.getModel()
    rID_list = get_source_reactions(sbml)  # list of source/sink reactions
    num_balanced, num_imbalanced, num_unknown = 0, 0, 0

    for reaction in model.getListOfReactions():
        rID = reaction.getId()
        if rID not in rID_list:
            formula, unknown = elementally_balance_reaction(rID, sbml)
            if unknown:
                num_unknown += 1
                if display_errors:
                    print '\n%s\t%s\t%s\t[%s]' % ('reaction', rID, 'unknown',
                                                  formula)
            elif formula:
                num_imbalanced += 1
                if display_errors:
                    print '\n%s\t%s\t%s\t[%s]' % ('reaction', rID,
                                                  'unbalanced', formula)
            else:
                num_balanced += 1
            if display_errors and (unknown or formula):
               txt, txt_formula= display_reaction_and_formula(rID, sbml)
               print txt + '\n' + txt_formula

    print ''
    print '%g\t%s' % (num_balanced, 'reactions balanced')
    print '%g\t%s' % (num_imbalanced, 'reactions unbalanced')
    print '%g\t%s' % (num_unknown, 'reactions unknown')


def get_source_reactions(sbml):
    '''Determine source and sink reactions'''
    model = sbml.getModel()

    rID_list = []

    # strip out format used in recon 2.1
    species = model.getSpecies('M_carbon_e')
    if species:
        species.setBoundaryCondition(True)

    for reaction in model.getListOfReactions():
        nS, nP = 0, 0
        #outs=[]
        for reactant in reaction.getListOfReactants():
            sID = reactant.getSpecies()
            if not model.getSpecies(sID).getBoundaryCondition():
                nS += 1
        for product in reaction.getListOfProducts():
            sID = product.getSpecies()
            if not model.getSpecies(sID).getBoundaryCondition():
                nP += 1
        if (nS == 0) or (nP == 0):
            rID_list.append(reaction.getId())
    return rID_list


def max_fluxes(sbml):
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
            v, f_opt = max_flux(sbml, carbon_source, objective, normoxic, [])
            print '%s (%s):\t%g' % (carbon_source,
                                    'normoxic' if normoxic else 'anaerobic',
                                    f_opt)
            for j, reaction in enumerate(model.getListOfReactions()):
               if abs(v[j])>1e-10:
                  rID = reaction.getId()
                  txt, txt_formula, innodes, outnodes=mydisplay_reaction_and_formula_Nodes(rID, sbml, v[j]>0)
                  reactionnode=format_for_ID_SBML(rID)
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
    nodef=open('EX_fru(e)False.txt', 'w')
    nodef.write("\n".join(listOfNodes))
    nodef.write("\n")
    nodef.write("\n".join([";".join(r) for r in listOfEdges]))
    nodef.close()

    return (listOfNodes, listOfEdges)


def max_flux(sbml, carbon_sourcel, objective, normoxic, media):
    '''
    Written to mimic neilswainston matlab function maxFlux
    '''
    #set_infinite_bounds(sbml)
    # block import reactions
    block_all_imports(sbml)
    # define carbon source
    set_import_bounds(sbml, carbon_sourcel, 50)
    # define media
    set_import_bounds(sbml, media, INF)
    if normoxic:
        set_import_bounds(sbml, 'EX_o2(e)', INF)
    # specify objective and maximise
    change_objective(sbml, objective)
    # avoid infinities
    obj_max = 1e6
    change_rxn_bounds(sbml, objective, obj_max, 'u')
    v, f_opt = optimize_cobra_model(sbml)
    if f_opt > 0.9 * obj_max:
        f_opt = INF
    return v, f_opt


def read_sbml(filename):
    '''
    Read an SBML file from specified path.

    '''
    reader = libsbml.SBMLReader()
    sbml = reader.readSBMLFromFile(filename)
    return sbml


def block_all_imports(sbml):
    '''
    Written to mimic neilswainston matlab function blockAllImports
    '''
    model = sbml.getModel()

    for rID in get_source_reactions(sbml):
        reaction = model.getReaction(rID)
        nR, nP = 0, 0
        for reactant in reaction.getListOfReactants():
            sID = reactant.getSpecies()
            if not model.getSpecies(sID).getBoundaryCondition():
                nR += 1
        for product in reaction.getListOfProducts():
            sID = product.getSpecies()
            if not model.getSpecies(sID).getBoundaryCondition():
                nP += 1
        kineticLaw = reaction.getKineticLaw()
        if (nR == 1) and (nP == 0):
            kineticLaw.getParameter('LOWER_BOUND').setValue(0)
        if (nR == 0) and (nP == 1):
            kineticLaw.getParameter('UPPER_BOUND').setValue(0)


def change_rxn_bounds(sbml, rxn_name_list, value, bound_type='b'):
    '''
    Written to mimic the matlab function changeRxnBounds from
    http://opencobra.sf.net/
    '''
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(value, (int, float, long, complex)):
        value = [value] * len(rxn_name_list)
    if isinstance(bound_type, str):
        bound_type = [bound_type] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print 'reaction %s not found' % rID
        else:
            kineticLaw = reaction.getKineticLaw()
            if bound_type[index] in ['l', 'b']:
                kineticLaw.getParameter('LOWER_BOUND').setValue(value[index])
            if bound_type[index] in ['u', 'b']:
                kineticLaw.getParameter('UPPER_BOUND').setValue(value[index])


def change_objective(sbml, rxn_name_list, objective_coeff=1):
    '''
    Written to mimic the matlab function changeObjective from
    http://opencobra.sf.net/
    '''
    model = sbml.getModel()
    for reaction in model.getListOfReactions():
        kineticLaw = reaction.getKineticLaw()
        kineticLaw.getParameter('OBJECTIVE_COEFFICIENT').setValue(0)
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(objective_coeff, (int, float, long, complex)):
        objective_coeff = [objective_coeff] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print 'reaction %s not found' % rID
        else:
            kineticLaw = reaction.getKineticLaw()
            kineticLaw.getParameter('OBJECTIVE_COEFFICIENT').setValue(
                objective_coeff[index])


def format_for_SBML_ID(txt):
    '''
    Written to mimic the matlab function formatForSBMLID from
    http://opencobra.sf.net/
    '''
    txt = 'R_' + txt
    for symbol, replacement in [
            ('-', '_DASH_'),
            ('/', '_FSLASH_'),
            ('\\', '_BSLASH_'),
            ('(', '_LPAREN_'),
            (')', '_RPAREN_'),
            ('[', '_LSQBKT_'),
            (']', '_RSQBKT_'),
            (',', '_COMMA_'),
            ('.', '_PERIOD_'),
            ('\'', '_APOS_'),
            ('&', '&amp'),
            ('<', '&lt'),
            ('>', '&gt'),
            ('"', '&quot')]:
        txt = txt.replace(symbol, replacement)
    return txt

def format_for_ID_SBML(txt):
    '''
    Written to mimic the matlab function formatForSBMLID from
    http://opencobra.sf.net/
    '''
    for symbol, replacement in [
            ('-', '_DASH_'),
            ('/', '_FSLASH_'),
            ('\\', '_BSLASH_'),
            ('(', '_LPAREN_'),
            (')', '_RPAREN_'),
            ('[', '_LSQBKT_'),
            (']', '_RSQBKT_'),
            (',', '_COMMA_'),
            ('.', '_PERIOD_'),
            ('\'', '_APOS_'),
            ('&', '&amp'),
            ('<', '&lt'),
            ('>', '&gt'),
            ('"', '&quot')]:
        txt = txt.replace(replacement, symbol)
    return txt


def optimize_cobra_model(sbml):
    '''
    Replicate Cobra command optimizeCbModel(model,[],'one').
    '''
    bound = INF
    cobra = convert_sbml_to_cobra(sbml, bound)

    N, L, U = cobra['S'], list(cobra['lb']), list(cobra['ub'])
    f, b = list(cobra['c']), list(cobra['b'])
    v_sol, f_opt, _ = easy_lp(f, N, b, L, U, one=False)
    return v_sol, f_opt


def convert_sbml_to_cobra(sbml, bound=INF):
    '''
    Get Cobra matrices from SBML model.
    '''
    model = sbml.getModel()
    S = sparse.lil_matrix((model.getNumSpecies(), model.getNumReactions()))
    lb, ub, c, b, rev, sIDs = [], [], [], [], [], []
    for species in model.getListOfSpecies():
        sIDs.append(species.getId())
        b.append(0.)
    sIDs = [species.getId() for species in model.getListOfSpecies()]
    for j, reaction in enumerate(model.getListOfReactions()):
        for reactant in reaction.getListOfReactants():
            sID = reactant.getSpecies()
            s = reactant.getStoichiometry()
            if not model.getSpecies(sID).getBoundaryCondition():
                i = sIDs.index(sID)
                S[i, j] = S[i, j] - s
        for product in reaction.getListOfProducts():
            sID = product.getSpecies()
            s = product.getStoichiometry()
            if not model.getSpecies(sID).getBoundaryCondition():
                i = sIDs.index(sID)
                S[i, j] = S[i, j] + s

        kinetic_law = reaction.getKineticLaw()
        rxn_lb = kinetic_law.getParameter('LOWER_BOUND').getValue()
        rxn_ub = kinetic_law.getParameter('UPPER_BOUND').getValue()
        rxn_c = kinetic_law.getParameter('OBJECTIVE_COEFFICIENT').getValue()
        rxn_rev = reaction.getReversible()
        if rxn_lb < -bound:
            rxn_lb = -bound
        if rxn_ub > bound:
            rxn_ub = bound
        if rxn_lb < 0:
            rxn_rev = True
        lb.append(rxn_lb)
        ub.append(rxn_ub)
        c.append(rxn_c)
        rev.append(rxn_rev)
    lb, ub, c, b = numpy.array(lb), numpy.array(ub), numpy.array(c), \
        numpy.array(b)
    rev = numpy.array(rev)
    cobra = {'S': S, 'lb': lb, 'ub': ub, 'c': c, 'b': b, 'rev': rev}
    return cobra


def easy_lp(f, a, b, vlb, vub, one=False):
    '''
    Optimize lp using Gurobi.
    '''
    # create gurobi model
    lp = gurobipy.Model()
    lp.Params.OutputFlag = 0
    lp.Params.FeasibilityTol = 1e-9  # as per Cobra
    lp.Params.OptimalityTol = 1e-9  # as per Cobra
    rows, cols = a.shape
    # add variables to model
    for j in xrange(cols):
        LB = vlb[j]
        if LB == -INF:
            LB = -gurobipy.GRB.INFINITY
        UB = vub[j]
        if UB == INF:
            UB = gurobipy.GRB.INFINITY
        lp.addVar(lb=LB, ub=UB, obj=f[j])
    lp.update()
    lpvars = lp.getVars()
    # iterate over the rows of S adding each row into the model
    S = a.tocsr()
    for i in xrange(rows):
        start = S.indptr[i]
        end = S.indptr[i + 1]
        variables = [lpvars[j] for j in S.indices[start:end]]
        coeff = S.data[start:end]
        expr = gurobipy.LinExpr(coeff, variables)
        lp.addConstr(lhs=expr, sense=gurobipy.GRB.EQUAL, rhs=b[i])
    lp.update()
    lp.ModelSense = -1
    lp.optimize()

    u = numpy.empty(len(f))
    v = numpy.empty(len(f))
    u[:] = NAN
    v[:] = NAN
    f_opt = NAN
    conv = False
    if lp.Status == gurobipy.GRB.OPTIMAL:
        f_opt = lp.ObjVal
        conv = True
        v = [var.x for var in lp.getVars()]


    if f_opt == -0.0:
        f_opt = 0.0

    return v, f_opt, conv


def get_reaction_by_id(sbml, rID):
    '''Gets the reaction by id.'''
    model = sbml.getModel()
    reaction = model.getReaction(rID)
    if not reaction:
        # try cobra replacements
        rID = format_for_SBML_ID(rID)
        reaction = model.getReaction(rID)
    if not reaction:
        # try removing trailing underscore
        if rID[-1] == '_':
            rID = rID[:-1]
        reaction = model.getReaction(rID)
    if not reaction:
        # try adding '_in'
        reaction = model.getReaction(rID + '_in')
    if not reaction:
        # try known alternatives
        rID_map = {
            'R_DM_atp_c': 'R_HKt',  # alternative ATPase
            # alternative C10:0
            'R_EX_HC02175_LPAREN_e_RPAREN': 'R_EX_dca_LPAREN_e_RPAREN_',
            # alternative C12:0
            'R_EX_HC02176_LPAREN_e_RPAREN': 'R_EX_ddca_LPAREN_e_RPAREN_',
            # alternative C22:0
            'R_EX_docosac': 'R_EX_docosac_LPAREN_e_RPAREN_',
        }
        if rID in rID_map:
            rID = rID_map[rID]
            reaction = get_reaction_by_id(sbml, rID)
    return reaction


def set_import_bounds(sbml, rxn_name_list, value):
    '''Sets the import bounds.'''
    model = sbml.getModel()
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(value, (int, float, long, complex)):
        value = [value] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print 'reaction %s not found' % rID
        else:
            nR, nP = 0, 0
            for reactant in reaction.getListOfReactants():
                sID = reactant.getSpecies()
                if not model.getSpecies(sID).getBoundaryCondition():
                    nR += 1
            for product in reaction.getListOfProducts():
                sID = product.getSpecies()
                if not model.getSpecies(sID).getBoundaryCondition():
                    nP += 1
            kineticLaw = reaction.getKineticLaw()
            val = abs(value[index])
            if (nR == 0) and (nP == 1):
                kineticLaw.getParameter('UPPER_BOUND').setValue(val)
            elif (nR == 1) and (nP == 0):
                kineticLaw.getParameter('LOWER_BOUND').setValue(-val)
            else:
                print 'reaction %s not import' % rID


def set_infinite_bounds(sbml):
    '''
    Set default bounds to INF, rather than 1000 (say)
    '''
    model = sbml.getModel()
    for reaction in model.getListOfReactions():
        kineticLaw = reaction.getKineticLaw()
        param = kineticLaw.getParameter('LOWER_BOUND')
        if param.getValue() < -100:
            param.setValue(-INF)
        param = kineticLaw.getParameter('UPPER_BOUND')
        if param.getValue() > 100:
            param.setValue(INF)


def formula_to_map(formula):
    '''
    Tranform a formula string 'C6H12O6' into a map dictionary {'C':6, 'H':12,
    'O':6}
    '''
    formula_map = {}
    # do not parse bracketing
    if ('(' in formula) or (')' in formula):
        return formula_map
    # replace FULLR notation with R
    for R in ['FULLR3', 'FULLR2', 'FULLR']:
        formula = formula.replace(R, 'R')
    m = re.findall('([A-Z][a-z]*)([\\-]?[\\d]*[\\.]?[\\d]*)', formula)
    for duple in m:
        X = duple[0]
        n = duple[1]
        if not n:
            n = 1.0
        else:
            n = float(n)
        if X not in formula_map:
            formula_map[X] = 0.0
        formula_map[X] = formula_map[X] + n
    # do not parse n
    if 'n' in formula_map:
        formula_map = {}

    return formula_map


def map_to_formula(formula_map):
    '''
    Tranform a map dictionary {'C':6, 'H':12, 'O':6} into a formula string
    'C6H12O6'
    '''
    formula = ''
    for X in formula_map.keys():
        if formula_map[X] == 0:
            del formula_map[X]
        elif formula_map[X] == int(formula_map[X]):
            formula_map[X] = int(formula_map[X])
    # force some elements to appear first
    for X in ['C', 'H']:
        if X in formula_map:
            n = formula_map.pop(X)
            formula = formula + X
            if n != 1:
                formula = formula + str(n)
    for X in sorted(formula_map.keys()):
        n = formula_map[X]
        formula = formula + X
        if n != 1:
            formula = formula + str(n)

    return formula


def elementally_balance_reaction(rID, sbml):
    '''
    Determine the overall elemental imbalance of a reaction
    '''
    reaction = sbml.getModel().getReaction(rID)

    unknown = False
    overall_formula = ''
    for reactant in reaction.getListOfReactants():
        sID = reactant.getSpecies()
        formula = get_formula(sID, sbml)
        stoich = reactant.getStoichiometry()
        if not formula_to_map(formula) and formula not in ['.']:
            unknown = True
        overall_formula = add_formulae(overall_formula, formula, -stoich)
    for reactant in reaction.getListOfProducts():
        sID = reactant.getSpecies()
        formula = get_formula(sID, sbml)
        stoich = reactant.getStoichiometry()
        if not formula_to_map(formula) and formula not in ['.']:
            unknown = True
        overall_formula = add_formulae(overall_formula, formula, stoich)

    return overall_formula, unknown


def add_formulae(formula, formula_new, stoich=1.0):
    '''
    Calculate formula + stoich * formula_new
    e.g. formula = C6H12O6, formula_new = H2O, stoich = -6 -> C6
    '''
    formula_map = formula_to_map(formula)
    formula_map_new = formula_to_map(formula_new)

    for X in formula_map_new:
        if X not in formula_map:
            formula_map[X] = 0
        formula_map[X] = formula_map[X] + stoich * formula_map_new[X]
    return map_to_formula(formula_map)


def get_formula(sID, sbml):
    '''
    Get the formula of species with ID sID
    '''
    return get_notes_field(sID, 'FORMULA', sbml)


def get_notes_field(eID, name, sbml):
    '''
    Gets the notes field.
    '''
    element = sbml.getModel().getElementBySId(eID)
    notes = element.getNotesString()
    f = re.search(name + ':([^<]+)', notes)
    return f.group(1).strip() if f is not None else ''


def display_reaction_and_formula(rID, sbml):
    '''
    Display reaction, with formulae below each reactant.
    '''
    model = sbml.getModel()
    reaction = model.getReaction(rID)

    txt = ''
    txt_formula = ''

    for index, reactant in enumerate(reaction.getListOfReactants()):
        if index > 0:
            txt += '+ '
            txt_formula += '+ '
        stoich = reactant.getStoichiometry()
        if stoich == int(stoich):
            stoich = int(stoich)
        if stoich != 1:
            txt += str(stoich) + ' '
            txt_formula += str(stoich) + ' '
        sID = reactant.getSpecies()
        txt += sID
        formula = get_formula(sID, sbml)
        if not formula:
            formula = '?'
        txt_formula += formula
        len_txt = len(txt)
        len_txt_formula = len(txt_formula)
        if len_txt > len_txt_formula:
            txt_formula = txt_formula.ljust(len_txt)
        elif len_txt < len_txt_formula:
            txt = txt.ljust(len_txt_formula)
        txt += ' '
        txt_formula += ' '

    if reaction.getReversible():
        txt += '<-> '
        txt_formula += '<-> '
    else:
        txt += '--> '
        txt_formula += '--> '

    for index, reactant in enumerate(reaction.getListOfProducts()):
        if index > 0:
            txt += '+ '
            txt_formula += '+ '
        stoich = reactant.getStoichiometry()
        if stoich == int(stoich):
            stoich = int(stoich)
        if stoich != 1:
            txt += str(stoich) + ' '
            txt_formula += str(stoich) + ' '
        sID = reactant.getSpecies()
        txt += sID
        formula = get_formula(sID, sbml)
        if not formula:
            formula = '?'
        txt_formula += formula
        len_txt = len(txt)
        len_txt_formula = len(txt_formula)
        if len_txt > len_txt_formula:
            txt_formula = txt_formula.ljust(len_txt)
        elif len_txt < len_txt_formula:
            txt = txt.ljust(len_txt_formula)
        if index < reaction.getNumProducts() - 1:
            txt += ' '
            txt_formula += ' '

    return txt, txt_formula


def genTxtAndTxtFormulaAndNodes(sbml, l=[]):
    txt = ''
    txt_formula = ''
    nodesl=[]
    for index, reactant in enumerate(l):
        if index > 0:
            txt += '+ '
            txt_formula += '+ '
        stoich = reactant.getStoichiometry()
        if stoich == int(stoich):
            stoich = int(stoich)
        if stoich != 1:
            txt += str(stoich) + ' '
            txt_formula += str(stoich) + ' '
        sID = reactant.getSpecies()
        txt += sID
        formula = get_formula(sID, sbml)
        nodesl.append(sID)
        if not formula:
            formula = '?'
        txt_formula += formula
        len_txt = len(txt)
        len_txt_formula = len(txt_formula)
        if len_txt > len_txt_formula:
            txt_formula = txt_formula.ljust(len_txt)
        elif len_txt < len_txt_formula:
            txt = txt.ljust(len_txt_formula)
        txt += ' '
        txt_formula += ' '
    return txt, txt_formula, nodesl

def mydisplay_reaction_and_formula_Nodes(rID, sbml, forward=True):
    '''
    Display reaction, with formulae below each reactant.
    '''
    model = sbml.getModel()
    reaction = model.getReaction(rID)

    itxt = '--> '
    itxt_formula = '--> '

    rtxt, rtxt_formula, rnodes=genTxtAndTxtFormulaAndNodes(sbml, reaction.getListOfReactants())
    ptxt, ptxt_formula, pnodes=genTxtAndTxtFormulaAndNodes(sbml, reaction.getListOfProducts())
    if forward:
       txt=rtxt+itxt+ptxt
       txt_formula=rtxt_formula+itxt_formula+ptxt_formula
       innodes=rnodes
       outnodes=pnodes
    else:
       txt=ptxt+itxt+rtxt
       txt_formula=ptxt_formula+itxt_formula+rtxt_formula
       outnodes=rnodes
       innodes=pnodes
    return txt, txt_formula, innodes, outnodes


if __name__ == '__main__':
    run_model()

    print 'DONE!'
