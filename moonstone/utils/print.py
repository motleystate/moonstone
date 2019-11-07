def print_multi_col(str_list, line_length=80):
    """
    Allow to print several items from a list in the same line
    until it reaches the line_length limit
    """
    print()
    current_line_len = 0
    for el in str_list:
        el = f"  {el}"
        current_line_len += len(el)
        if current_line_len > line_length:
            print()
            current_line_len = len(el)
        print(el, end='')
    print()
