import math
from typing import Union, List


class ColumnDef(object):
    """ Represents a data value from the Illumina InterOp output """
    def __init__(self, name: str, field: Union[str, List[str]], **kwargs):
        self.name = name

        # Allow multiple fields in one 'cell'
        if type(field) is str:
            self.fields = [field]
        else:
            self.fields = field

        # Default precision and scale to 2 and 1, respectively
        self.precision = kwargs.get('precision', 2)
        self.scale = kwargs.get('scale', 1)

    def get_value_from_row(self, row) -> str:
        """
        The attributes are a proxy of the c implementation, the values are actually a method getter
        This method works with values that have mean/std deviations, ranges, and multiple values
        """

        if len(self.fields) > 1:
            # Special case where there are multiple fields
            values = []
            for field in self.fields:
                values.append(str(self._format(getattr(row, field)().mean())))
            return ' / '.join(values)

        value_method = getattr(row, self.fields[0])
        if value_method is None:
            return 'N/A'

        value = value_method()

        if hasattr(value, 'mean'):
            # Special case when the value has a mean, add standard deviation
            return f'{self._format(value.mean())} +/- {self._format(value.stddev())}'
        elif hasattr(value, 'error_cycle_range'):
            # Special case when there is a cycle range
            error_range = value.error_cycle_range()
            first = error_range.first_cycle()
            last = error_range.last_cycle()
            if first == last:
                return self._format(first)
            return f'{self._format(first)} - {self._format(last)}'
        else:
            return self._format(value)

    def _format(self, value):
        value = value / self.scale

        if not math.isnan(value):
            value = round(value, self.precision)
        return value
