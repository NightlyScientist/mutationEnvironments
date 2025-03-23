# %%
from matplotlib import pyplot as plt
from matplotlib import colors

def colorbar_plot(cmap, custom_norm, imgPath, position="top"):
    fig, ax = plt.subplots(figsize=(6, 6))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=custom_norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label="", location=position)
    cbar.ax.tick_params(labelsize=11)
    #fig.savefig(f"{imgPath}/colorbar_{cmap}.svg", transparent=True)


colorbar_plot("PuOr", colors.TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=0.5), None, position="top")

# %%
if __name__ == "__main__":
    imgPath = input("Enter the path for the images: ")

    colorbar_plot("coolwarm", colors.TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=0.5), imgPath)
    colorbar_plot("turbo", colors.TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=0.5), imgPath)