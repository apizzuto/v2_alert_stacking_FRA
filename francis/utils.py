import os
import matplotlib as mpl
mpl.use('Agg')

def get_francis_path():
    import francis
    return os.path.join(francis.__file__, '../')

def initialize_mpl_style():
    base = get_francis_path()
    print("Setting matplotlib style")
    return base + 'francis_plots.mplstyle'