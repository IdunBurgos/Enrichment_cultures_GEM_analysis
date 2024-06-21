from matplotlib import colors

def cmap_community():
    # red mid color FF7074
    # green mid color 50C878
    cs = [ "red","#F88484","#F4CCCC","#00A36C","#9FD787","#D9EAD3","#6497b1","#b9d7f2"]
    cmap = colors.ListedColormap( cs)
    return cmap
