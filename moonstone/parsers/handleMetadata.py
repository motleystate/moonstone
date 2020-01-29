"""
Here is where the clinical data for the samples will be opened.
This file is likely to be highly variable depending on the project
    and will likely need several versions of the import code.
At the least, samples names/references should match with the 'counts' file

@TODO
It is likely that redundacy exists with the MetadataParser. We would need to fuse both part of code.
"""
import pandas as pd
import logging
module_logger = logging.getLogger(__name__)


class Inputs(object):
    def __init__(self, metafile):
        self.logger = module_logger
        self.logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.metafile = metafile

    def openmeta(self, outdir_path):
        self.logger.info('Opening Metadata file: ' + self.metafile)
        df = pd.read_csv(self.metafile, sep=',', index_col=0)
        df.rename_axis("sample", inplace=True)
        df.to_csv(path_or_buf=outdir_path + '/' + 'importedMetaData.csv')
        self.logger.info('Success!')
        # variables_dict = DataAnalysisType.Variables.inspect(df, df)
        return df
