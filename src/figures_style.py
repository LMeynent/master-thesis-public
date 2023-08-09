import matplotlib.pyplot as plt
import seaborn as sns

def initialise_figures(context='notebook'):
    sns.set_theme(context=context, style='whitegrid', 
                  palette='blend:#8b0000,#fcdc5a',
                  rc={
                      'font.family': 'sans-serif',
                      'font.sans-serif': 'Fira Sans Condensed',
                      'axes.spines.right': False,
                      'axes.spines.top': False,
                      'xtick.bottom': True,
                      'ytick.left': True
                  })
    
def get_color(x, n):
    palette = sns.color_palette(palette='blend:#8b0000,#fcdc5a', as_cmap=True)
    return palette(x/(n-1))

def get_categorical_palette(n=5):
    if n <= 6:
        return ['#8b0000', '#fcdc5a', '#3B3561', '#6A8D73', '#78D5D7', '#c96ba5'][:n]
    
    raise NotImplementedError()