
def validate_float(new_value):
    if new_value == "" or new_value == "-":
        return True
    try:
        float(new_value)
        return True
    except ValueError:
        return False

def validate_index(new_value):
    if new_value == "":
        return True
    try:
        int(new_value)
        return True
    except ValueError:
        return False

def validate_bool(new_value):
    if (new_value == "" or new_value == "F" or new_value == "Fa" or new_value == "Fal" or new_value == "Fals"
            or new_value == "T" or new_value == "Tr" or new_value == "Tru"):
        return True
    try:
        bool(new_value)
        return True
    except ValueError:
        return False