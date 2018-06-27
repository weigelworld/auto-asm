"""
plots module

handles common plotting functionality of auto-asm
"""

from typing import List, Tuple, Union

from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


AxisLimit = Union[int, float, complex, None]
AxisLimits = Tuple[AxisLimit, AxisLimit]


def setup_single_figure(title: str = None, xlabel: str = None, ylabel: str = None,
                             xlim: AxisLimits = None, ylim: AxisLimits = None) -> (Figure, Axes):
    fig = Figure()
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlim:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim:
        ax.set_ylim(ylim[0], ylim[1])

    fig.set_tight_layout(True)

    return fig, ax


def read_length_histogram(read_lengths: List[int], assembly_id: str = None, readset_id: str = None) -> Figure:
    """
    plots a histogram of read lengths given the reads
    """

    if assembly_id and readset_id:
        read_length_hist_title = "{asm} - readset{id}".format(asm=assembly_id, id=readset_id)
    else:
        read_length_hist_title = ""

    read_length_figure, read_length_axes = setup_single_figure(
        title=read_length_hist_title,
        xlabel='Read Length (bp)',
        ylabel='Frequency'
    )

    read_length_axes.hist(read_lengths)

    return read_length_figure
