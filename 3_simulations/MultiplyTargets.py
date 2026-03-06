import pandas as pd
import numpy as np
import glob
from pathlib import Path


class MultiplyTargets:
    """
    Reusable calculator for one signature across many species.
    Mirrors your original logic, with robustness for empty/missing files.
    """

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

        self.sig = signature
        self.gene_strand_file = gene_strand_file
        self.tsb_sig_dir_pos = tsb_sig_dir_pos
        self.tsb_sig_dir_neg = tsb_sig_dir_neg
        self.non_tsb_sig_dir = non_tsb_sig_dir
        self.input_glob = input_glob
        self.output_dir = output_dir
        self.min_possible_muts = min_possible_muts
        self.SIGNATURES_NO_TSB = {
            "SBS10d","SBS22a","SBS22b","SBS23","SBS40a","SBS40b","SBS40c","SBS42",
            "SBS86","SBS87","SBS88","SBS89","SBS90","SBS91","SBS92","SBS93","SBS94",
            "SBS96","SBS97","SBS98","SBS99"
        }
        

        # Load strand map once per signature (per process)
        self.gene_strand_df = pd.read_csv(self.gene_strand_file, index_col=0)

        # Load signature vectors (pos/neg)
        if self.sig in self.SIGNATURES_NO_TSB:
            sig_path = f"{self.non_tsb_sig_dir}/{self.sig}"
            self.sig_df_pos = pd.read_csv(sig_path, index_col=0)
            self.sig_df_neg = self.sig_df_pos  # same for both strands
        else:
            self.sig_df_pos = pd.read_csv(f"{self.tsb_sig_dir_pos}/{self.sig}", index_col=0)
            self.sig_df_neg = pd.read_csv(f"{self.tsb_sig_dir_neg}/{self.sig}", index_col=0)

        # Ensure numeric, and align by index labels during mul
        #self.sig_df_pos = self.sig_df_pos.apply(pd.to_numeric, errors="coerce").fillna(0)
        #self.sig_df_neg = self.sig_df_neg.apply(pd.to_numeric, errors="coerce").fillna(0)

    @staticmethod
    def _safe_read_gene_df(path: str) -> pd.DataFrame | None:
        try:
            df = pd.read_csv(path, index_col=0)
        except Exception:
            return None
        if df.empty:
            return None
        # Ensure numeric
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
        return df

    def _find_target_sizes_of_gene(self, gene_df: pd.DataFrame, strand: str) -> tuple[float, int]:
        """
        Returns (x, n) where x = l_n/n (or -1 if n < threshold), and n = total possible muts.
        """
        # Multiply with alignment on both axes — safer than relying on column order
        if strand == "+":
            target_df = gene_df.mul(self.sig_df_pos, fill_value=0)
        else:
            target_df = gene_df.mul(self.sig_df_neg, fill_value=0)

        l_n = float(np.nansum(target_df.values))
        n = int(np.nansum(gene_df.values))

        if n >= self.min_possible_muts:
            return l_n / n, n
        return -1.0, n

    def _strand_for_gene(self, gene: str) -> str | None:
        try:
            s = self.gene_strand_df.loc[gene, "strand"]
            # If multiindex row, pick first
            if isinstance(s, pd.Series):
                s = s.iloc[0]
            return s
        except Exception:
            return None

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





    