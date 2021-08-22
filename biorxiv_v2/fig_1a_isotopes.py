import matplotlib.pyplot as plt
import numpy as np
import proplot as pplt

from functions.formatting import *
from functions.base import *
from functions.logging import write_log

write_log(__file__)

output = 'save'
output_dir = current_dir / 'figures' / 'Fig_1'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_1a_isotopes'

h1 = np.array([6.9, 9.1, 5.9, 3.0, 0.7, 0.3])
h1 / h1.max()

h2 = np.array([2.1, 7.2, 9.2, 8.8, 5.5, 2.9, 0.9, 0.2])
h2 / h2.max()

sclf = 2
fig_width = 25*sclf
fig_height = 25*sclf

h1 = np.array([0.75824176, 1.        , 0.64835165, 0.32967033, 0.07692308, 0.03296703])
h2 = np.array([0.22826087, 0.7826087 , 1.        , 0.95652174, 0.59782609, 0.31521739,
               0.09782609, 0.02173913])


def lorentzian(x, x0, gamma):
    return np.divide(1, np.pi*gamma* (1 + ( (x-x0) / gamma  )**2))


x = np.linspace(0, 10, num=1000)

y1 = np.zeros_like(x)
for i, height in enumerate(h1):
    y1 += height*lorentzian(x, i+1, 0.01)

y2 = np.zeros_like(x)
for i, height in enumerate(h2):
    y2 += height*lorentzian(x, i+1, 0.01)

fig, axes = pplt.subplots(nrows=2, width=fig_width / 25.4, height=fig_height/25.4, wspace=0.01, hspace=0.05, sharey=False)
axes[0].plot(x, y1 / y1.max(), linewidth=1, color='k')
axes[0].format(ylabel='t$_{0}$')
axes[1].plot(x, y2 / y1.max(), linewidth=1, color='k')
axes[1].format(ylabel='t$_{1}$')
axes.format(xticks=[], yticks=[], xlabel='m/z')


if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')
