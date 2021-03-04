from mso.objectives.mol_functions import qed_score, penalize_macrocycles, docking_score, penalize_molecular_weight
from mso.objectives.emb_functions import distance_score
from mso.optimizer import BasePSOptimizer
from mso.objectives.scoring import ScoringFunction
from mso.objectives.mol_functions import qed_score
from cddd.inference import InferenceModel
import mso 
import cddd
import pandas as pd 
import numpy as np 


if __name__ == "__main__":

    init_smiles = "CCSCn1cc(C(=O)O)c(=O)n1C(=O)SC" #"c1ccccc1"
    infer_model = InferenceModel()

    #smiles_embedding = infer_model.seq_to_emb([init_smiles, ])

    scoring_functions = [ScoringFunction(func=docking_score, name='docking', weight=60, 
                                        is_mol_func=True, 
                                        additional_args={'receptor': "/data1/zlzzheng/apps/mso/notebooks/5xq0_A_rec.pdb_mgltools.pdbqt", 
                                                        'pocket': [39.929, 6.848, -49.476], 
                                                        'exe': "idock",
                                                        'verbose': False,
                                                        'output_dir': "/data1/zlzzheng/apps/mso/notebooks/docking",
                                                        }),
                        ScoringFunction(func=qed_score, name="qed", weight=60, is_mol_func=True), 
                        ScoringFunction(func=penalize_macrocycles, weight=10, name="marcocycles", is_mol_func=True), 
                        ScoringFunction(func=penalize_molecular_weight, weight=30, name="pmw", is_mol_func=True)
                        ]

    opt = BasePSOptimizer.from_query(
                                    init_smiles=init_smiles,
                                    num_part=20,
                                    num_swarms=2,
                                    inference_model=infer_model,
                                    scoring_functions=scoring_functions)


    opt.run(200)

    df = opt.best_fitness_history
    df.to_csv("best_fitness_history.csv", header=True, index=True)
