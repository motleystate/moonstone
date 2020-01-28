from moonstone.parsers.base import BaseParser


class Picrust2PathwaysParser(BaseParser):
    """
    Format is the following:

    pathways    sample_1    sample_2    sample_3
    pathway_1   14.3    123.4   12
    pathway_2   94.1    1231.1  124.2
    pathway_4   09.4    15.5    12.2
    """

    def __init__(self, *args, **kwargs):
        parsing_options = {
            'index_col': 0
        }
        super().__init__(*args, **kwargs, parsing_options=parsing_options)
