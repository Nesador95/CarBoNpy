import argparse
import CarBoN_Input_Processor as carbon 
import pandas as pd

###############################################################################
# Parser setup
###############################################################################

parser = argparse.ArgumentParser()

parser.add_argument("reac", help="path to KIDA reactions file",
                    type=str)

parser.add_argument("spec", help="path to KIDA species file",
                    type=str)

args = parser.parse_args()

###############################################################################
# Parser Execution with KIDA class Network tool
###############################################################################
reac =  args.reac
spec =  args.spec

kida_file = carbon.Kida(r"" + reac, r"" + spec)

kida_file.read_species()
kida_file.read_reactions()

reactions = kida_file.reactions_dataframe()
species = kida_file.species_dataframe()
dictionary = kida_file.species_dictionary()

print(reactions)

###############################################################################
# These are the paths used
###############################################################################

# data/kida_reac_C_O_only.dat
# data/kida_spec_C_O_only.dat