from string import ascii_lowercase, digits

from slugify import slugify


def generate_color_code(string: str) -> str:
    """Generate a color code from a string."""
    string = slugify(string)
    ref_char = ascii_lowercase + digits
    elements = string.split("-")
    if len(elements) > 1:
        string = elements[0][:3] + elements[1][:3]
    value_mapping = {ref_char[i]: i for i in range(0, len(ref_char))}
    value_mapping.update({" ": 0})
    color = list("#968579")
    indices = list([1, 4, 2, 3, 6, 5])
    iteration = 0
    for letter in string[:6]:
        letter_value = value_mapping[letter.lower()]
        color[indices[iteration]] = str(
            hex(letter_value + int(color[indices[iteration]]))
        )[-1]
        iteration += 1
    return "".join(color)
