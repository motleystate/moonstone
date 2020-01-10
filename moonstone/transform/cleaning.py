from slugify import slugify

from .base import TransformBase


class StringCleaning(TransformBase):

    def remove_trailing_spaces(self, col_name):
        self.df[col_name] = self.df[col_name].str.strip()
        self.historize(self.remove_trailing_spaces.__name__)

    def to_slug(self, col_name):
        self.df[col_name] = self.df[col_name].apply(slugify)
        self.historize(self.to_slug.__name__)
