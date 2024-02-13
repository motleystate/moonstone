"""Adapt Phylogenetic Tree to counts dataframe"""

import pandas as pd
import re


def generate_translation_dictionary(
    new_otu_id_name_ser: pd.Series
):
    """
    Args:
        - new_otu_id_name_ser: pd.Series issued from kraken 2 count dataframe with only new_otu_id_name column
    """
    level = "species"  # only works with species for now
    new_otu_id_name_ser.index = new_otu_id_name_ser.index.get_level_values(level)
    new_otu_id_name_ser = new_otu_id_name_ser[~new_otu_id_name_ser.index.str.contains("(", regex=False)].astype(str)
    new_otu_id_name_ser.index = new_otu_id_name_ser.index.str.replace("'", "''")
    dic_translate_tree = new_otu_id_name_ser.reset_index().set_index(new_otu_id_name_ser.name).to_dict()[level]
    return dic_translate_tree


def replacement(
    matchobj,
    dic_translate_tree: dict,
    quotechr: str
) -> str:
    s = matchobj.group(0).split(", ")[-1].rstrip("*"+quotechr)
    if s in dic_translate_tree.keys():
        return quotechr+dic_translate_tree[s]+quotechr
    else:
        return matchobj.group(0)


def replacing_labels(
    tree_string: str,
    dic_translate_tree: dict,
    quotechr: str = "'"
):
    return re.sub(r"'[^,]*, [0-9]*\*?'", lambda match: replacement(match, dic_translate_tree, quotechr), tree_string)


def adapt_phylogenetic_tree_to_counts_df(
    new_otu_id_name_ser: pd.Series,
    tree: str,
    output_tree_file: str = None,
    quotechr: str = "'"
):
    """
    Translate phylogenetic tree labels to names present in a counts dataframe using the txid as key
    Args:
        - new_otu_id_name_ser: pd.Series issued from count dataframe with only new_otu_id_name column
          ('NCBI_taxonomy_ID' for Kraken2, 'NCBI_tax_id' for Metaphlan3)
        - tree: path to the tree file to adapt or tree as a string. The format of the tree leaves labels should be
          '{species name}, {txid}' or '{species name}, {txid}*'
        - output_tree_file: path to the output adapted tree file.
          If None, then function return the adaptated tree as a string
        - quotechr: quote character used as delimiter of labels in tree
    """
    try:
        infile = open(tree, "r")
        T = infile.read()
        infile.close()
    except FileNotFoundError:
        T = tree

    dic_translate_tree = generate_translation_dictionary(new_otu_id_name_ser)
    T = replacing_labels(T, dic_translate_tree, quotechr)

    if output_tree_file:
        outfile = open(output_tree_file, "w")
        outfile.write(T)
        outfile.close()
    else:
        return T
