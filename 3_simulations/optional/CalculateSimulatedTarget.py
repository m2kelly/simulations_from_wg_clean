import pickle
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

class CalculateSimulatedTarget:
    def __init__(self,matrix_generator,input_file,output_dir):
        self.output_path=Path(output_dir)
        self.input_file = input_file
        self.matrix = matrix_generator

#from neutral evolution simulator
    def load_pickle_records_to_df(self,path):
        rows = []
        with open(path, "rb") as f:
            while True:
                try:
                    obj = pickle.load(f)
                except EOFError:
                    break

                if isinstance(obj, list):     # came from chunked writer
                    rows.extend(obj)
                else:                          # came from per-record writer
                    rows.append(obj)

        return pd.DataFrame(rows,columns=['gene','seq','trinucs']).set_index('gene') 

    #if cannot run parallel
    def run(self):
        self.output_path.mkdir(parents=True, exist_ok=True)
        #call instance of mutation matrix generatpr 
        df = self.load_pickle_records_to_df(self.input_file)
        #just want to use innit and one method 

        for gene, row in df.iterrows():
            gene_df =self.matrix.generate_mutations(row['seq'], row['trinucs'], 1) #1 is for non syn
            gene_df.to_csv(self.output_path / f"{gene}")

    #or if want to run in parallel
    def process_gene(self,args):
        gene, seq, trinucs, output_path, matrix = args
        gene_df = matrix.generate_mutations(seq, trinucs, 1)
        gene_df.to_csv(output_path / f"{gene}")
        return gene

    def run_parallel(self):
        self.output_path.mkdir(parents=True, exist_ok=True)
        #call instance of mutation matrix generatpr 
        df = self.load_pickle_records_to_df(self.input_file)
        tasks = [
            (gene, row['seq'], row['trinucs'], self.output_path, self.matrix)
            for gene, row in df.iterrows()
        ]

        with ProcessPoolExecutor(max_workers=15) as executor:
            executor.map(self.process_gene, tasks)
               