"""The goal of this script is to look through the processed metadata file
and to determine the type of variable in each column.

- String, Integer or Float
- If its a number: is the variable continuous or categorical
- Determine the number of categories (binary or more)
- In the case of continuous variable, determine normality
- For different groups, determine if they are balanced.

A second part will be to match the type of statistical analysis that could be performed.
"""

import  handleMetadata

parser = argparse.ArgumentParser(description='MetaData Description Script')
parser.add_argument("metadata", type=str, help="Clinical data input file")

read_file = handleMetadata.Inputs
