import logging
module_logger = logging.getLogger(__name__)
"""The goal of this script is to look through the processed metadata file
and to determine the type of variable in each column.

- String, Integer or Float
- If its a number: is the variable continuous or categorical
- Determine the number of categories (binary or more)
- In the case of continuous variable, determine normality
- For different groups, determine if they are balanced.

A second part will be to match the type of statistical analysis that could be performed.

@TODO
this needs to be moved to appropriate module of moonstone, to be defined
"""


class Variables(object):
    def __init__(self):
        self.logger = module_logger
        self.logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')

    def inspect(self, metadata_df):
        self.logger.info('Looking closer at the metadata structure.')
        column_dict = {c[0]: c[1] for c in enumerate(metadata_df.columns)}
        return column_dict

# import handleMetadata

# parser = argparse.ArgumentParser(description='MetaData Description Script')
# parser.add_argument("metadata", type=str, help="Clinical data input file")

# read_file = handleMetadata.Inputs
