#!/usr/bin/env python

# encode molecules from a SMILES file to a CSV files with
# a few molecular descriptors

from __future__ import print_function

import rdkit
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = words[1]
            yield (name, Chem.MolFromSmiles(smile))

if len(sys.argv) != 3:
    print("usage: %s input.smi output.csv" % sys.argv[0])
    sys.exit(1)

def main():
    input_smi = sys.argv[1]
    output_csv = sys.argv[2]
    output = open(output_csv, 'w')
    output.write("#logP\tmolMR\tmolW\tnbA\tnbD\tnbRotB\tTPSA\n")
    for name, mol in RobustSmilesMolSupplier(input_smi):
        if mol is None:
            continue
        logP = Descriptors.MolLogP(mol)
        molMR = Descriptors.MolMR(mol)
        molW = Descriptors.MolWt(mol)
        nbA = Descriptors.NumHAcceptors(mol)
        nbD = Descriptors.NumHDonors(mol)
        nbRotB = Descriptors.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)
        # FBR: TODO append counted atom pairs
        output.write("%f\t%f\t%f\t%d\t%d\t%d\t%f\n" %
                     (logP, molMR, molW, nbA, nbD, nbRotB, tpsa))
    output.close()

if __name__ == '__main__':
    main()
