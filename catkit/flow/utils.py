import importlib
supported_properties = [
    'energy', 'forces', 'stress',
    'charges', 'magmom', 'magmoms']


def str_to_class(classname):
    split_class = classname.split('.')
    modulename = '.'.join(split_class[:-1])
    classname = split_class[-1]
    c = getattr(importlib.import_module(modulename), classname)
    return c
