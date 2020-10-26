from moonstone.parsers.base import BaseParser


class GeneCountsParser(BaseParser):
    """
    Common way of representing gene counts per sample in a matrix.

    Format is the following:

    +-----------+----------+----------+
    | genes     | sample_1 | sample_2 |
    +===========+==========+==========+
    | gene_1    | 3        | 19       |
    +-----------+----------+----------+
    | gene_2    | 9        | 10       |
    +-----------+----------+----------+
    """

    def __init__(self, *args, **kwargs):
        parsing_options = {
            'index_col': 0
        }
        super().__init__(*args, **kwargs, parsing_options=parsing_options)
