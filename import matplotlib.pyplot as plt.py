import matplotlib.pyplot as plt
import mplcursors

# Create a figure and some data to plot
fig, ax = plt.subplots()
x = [1, 2, 3]
y = [4, 5, 6]
ax.plot(x, y)

# Create a HoverTool with an annotation box
tool = mplcursors.cursor(hover=True)

# Set the width of the annotation box using the width parameter of the AnnotationBbox class
tool.annotation_kwargs = {"bbox": {"facecolor": "white", "edgecolor": "black", "alpha": 0.8, "boxstyle": "round", "pad": 0.5},
                        "arrowprops": {"arrowstyle": "->"},
                        "xybox": (0, 0), "xycoords": "data", "boxcoords": "offset points",
                        "box_alignment": (0, 0), "xytext": (0, 0), "textcoords": "offset points",
                        "ha": "left", "va": "bottom", "fontsize": 8, "zorder": 10,
                        "bbox_transform": fig.transFigure,
                        "width": 200}

plt.show()
