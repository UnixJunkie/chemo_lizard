#!/usr/bin/env python

# type atoms of a molecule a la atom pairs
# (nb. pi electrons if > 0, elt. symbol, nbHA neighbors)

from __future__ import print_function

import rdkit
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AtomPairs import Pairs

PeriodicTable = Chem.GetPeriodicTable()

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = words[1]
            yield (Chem.MolFromSmiles(smile), name)

if len(sys.argv) != 2:
    print("usage: %s input.smi" % sys.argv[0])
    sys.exit(1)

def nb_heavy_atom_neighbors(a):
    neighbors = a.GetNeighbors()
    res = 0
    for n in neighbors:
        if n.GetAtomicNum() != 1:
            res = res + 1
    return res

def type_atom(a):
    nb_pi_electrons = Pairs.Utils.NumPiElectrons(a)
    symbol = PeriodicTable.GetElementSymbol(a.GetAtomicNum())
    nbHA = nb_heavy_atom_neighbors(a)
    res = ""
    if nb_pi_electrons > 0:
        res = "%d%s%d" % (nb_pi_electrons, symbol, nbHA)
    else:
        res = "%s%d" % (symbol, nbHA)
    return res

def main():
    input_smi = sys.argv[1]
    print("#name\tlogP\tMR\tMW\tHBA\HBD\RotB\tTPSA")
    for mol, name in RobustSmilesMolSupplier(input_smi):
        if mol is not None:
            logP = Descriptors.MolLogP(mol)
            molMR = Descriptors.MolMR(mol)
            molW = Descriptors.MolWt(mol)
            nbA = Descriptors.NumHAcceptors(mol)
            nbD = Descriptors.NumHDonors(mol)
            nbRotB = Descriptors.NumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            print("%s %f %f %f %d %d %d %f" %
                  (name, logP, molMR, molW, nbA, nbD, nbRotB, tpsa))

if __name__ == '__main__':
    main()
