import pandas as pd
import numpy as np
import glob
from pathlib import Path
from MultiplyTargets import MultiplyTargets


class MultiplyTargetsNoMn(MultiplyTargets):
    def __init__(
        self,
        signature: str,
        input_glob: str,
        output_dir: str,
        gene_strand_file: str,
        min_possible_muts: int = 1,
        tsb_sig_dir_pos: str = "/home/maria/signatures/TSB_signatures_scaled",
        tsb_sig_dir_neg: str = "/home/maria/signatures/TSB_signatures_scaled_rc",
        non_tsb_sig_dir: str = "/home/maria/signatures/COSMIC_sigs_uncollapsed_sum1",
    ):
        super().__init__(
            signature=signature,
            input_glob=input_glob,
            output_dir=output_dir,
            gene_strand_file=gene_strand_file,
            min_possible_muts=min_possible_muts,
            tsb_sig_dir_pos=tsb_sig_dir_pos,
            tsb_sig_dir_neg=tsb_sig_dir_neg,
            non_tsb_sig_dir=non_tsb_sig_dir,
        )

    def _find_target_sizes_of_gene(self, gene_df: pd.DataFrame, strand: str) -> tuple[float, int]:
        """
        Returns (x, n) where x = l_n(or -1 if n < threshold), and n = total possible muts.
        """
        # Multiply with alignment on both axes — safer than relying on column order
        if strand == "+":
            target_df = gene_df.mul(self.sig_df_pos, fill_value=0)
        else:
            target_df = gene_df.mul(self.sig_df_neg, fill_value=0)

        l_n = float(np.nansum(target_df.values))
        n = int(np.nansum(gene_df.values))

        if n >= self.min_possible_muts:
            return l_n, n
        return -1.0, n

    

    def calc_for_sim(self) -> None:
        """
        Process all gene files for a single species, write outputs:
        - l_{sig} : CSV without header, two columns [gene, l_n]
        - M_n_non_syn_muts : total N across genes (int)
        """
        
        gene_file_list = glob.glob(self.input_glob)
        output_directory = Path(self.output_dir)
        output_directory.mkdir(parents=True, exist_ok=True)

        genes, target_sizes = [], []
        total_N = 0

        for gene_file in gene_file_list:
            gene = gene_file[gene_file.rindex("/") + 1 :]
            strand = self._strand_for_gene(gene)
            if strand is None:
                # Skip if gene not in strand table
                continue

            gene_df = self._safe_read_gene_df(gene_file)
            if gene_df is None:
                continue

            # Compute per-gene
            x, n = self._find_target_sizes_of_gene(gene_df, strand)
            if n >= self.min_possible_muts and x != -1:
                genes.append(gene)
                target_sizes.append(x)
                total_N += n  # only count genes we accept, mirroring your logic

        # Write outputs
        if genes:
            pd.DataFrame({"gene": genes, "l_n": target_sizes}).to_csv(
                output_directory / f"l_{self.sig}",
                index=False,
                header=False)

        with open(output_directory / "M_n_non_syn_muts", "w") as f:
            f.write(str(total_N))





    