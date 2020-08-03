from moonstone.parsers.base import BaseParser


class Picrust2PathwaysParser(BaseParser):
    """
    Format is the following:

    +-----------+----------+----------+
    | pathways  | sample_1 | sample_2 |
    +===========+==========+==========+
    | pathway_1 | 14.3     | 123.4    |
    +-----------+----------+----------+
    | pathway_2 | 94.1     | 1232.1   |
    +-----------+----------+----------+
    """

    def __init__(self, *args, **kwargs):
        parsing_options = {
            'index_col': 0
        }
        super().__init__(*args, **kwargs, parsing_options=parsing_options)
