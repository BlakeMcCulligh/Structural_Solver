"""
Handels validating data types of inputs.
"""

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def validate_float(new_value: str) -> bool:
    """
    Validates if the input string can be converted to a float.

    :param new_value: value to be validated.
    :type new_value: str
    :return: if the input string can be converted to a float.
    :rtype: bool
    """

    if new_value == "" or new_value == "-":
        return True
    try:
        float(new_value)
        return True
    except ValueError:
        return False

def validate_index(new_value: str) -> bool:
    """
    Validates if the input string can be converted to an integer.

    :param new_value: value to be validated.
    :type new_value: str
    :return: if the input string can be converted to an integer.
    :rtype: bool
    """

    if new_value == "":
        return True
    try:
        int(new_value)
        return True
    except ValueError:
        return False

def validate_bool(new_value: str) -> bool:
    """
    Validates if the input string can be converted to a boolean or is in the proses of typing in a boolean.

    :param new_value: value to be validated.
    :type new_value: str
    :return: if the input string can be converted to a boolean or is in the proses of typing in a boolean.
    :rtype: bool
    """

    if (new_value == "" or new_value == "F" or new_value == "Fa" or new_value == "Fal" or new_value == "Fals"
            or new_value == "T" or new_value == "Tr" or new_value == "Tru"):
        return True
    try:
        bool(new_value)
        return True
    except ValueError:
        return False