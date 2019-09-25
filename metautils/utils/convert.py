def pandas_to_python_type(pandas_dtype):
    """
    Return the name of type for Python from pandas dtype

    if the type does not exist, return the pandas dtype
    if the type exists, return the python type with pandas dtype between brackets
    """
    python_type = None
    if pandas_dtype == 'object':
        python_type = 'str'
    elif pandas_dtype == 'int64':
        python_type = 'int'
    elif pandas_dtype == 'float64':
        python_type = 'float'
    if python_type:
        return python_type
    return pandas_dtype
