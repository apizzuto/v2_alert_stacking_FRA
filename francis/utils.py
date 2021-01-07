import matplotlib as mpl
mpl.use('Agg')

def get_francis_path():
    import francis
    return os.path.join(francis.__file__, '../')

def initialize_mpl_style():
    base = get_francis_path()
    return base + 'francis_plots.mplstyle'