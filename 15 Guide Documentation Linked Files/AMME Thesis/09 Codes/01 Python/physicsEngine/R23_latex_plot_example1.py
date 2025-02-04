# import matplotlib.pyplot as plt   # side-stepping mpl's backend
# import chart_studio.plotly as py
# import plotly.tools as tls
# from plotly.graph_objs import *
# # matplotlib inline
# py.sign_in("IPython.Demo", "1fw3zw2o13")
# fig1 = plt.figure()

# import matplotlib.pyplot as plt
# import numpy as np

# x = np.linspace(-2.0, 2.0, 10000) # The x-values
# sigma = np.linspace(0.4, 1.0, 4) # Some different values of sigma

# # Here we evaluate a Gaussians for each sigma
# gaussians = [(2*np.pi*s**2)**-0.5 * np.exp(-0.5*x**2/s**2) for s in sigma]

# ax = plt.axes()

# for s,y in zip(sigma, gaussians):
#     ax.plot(x, y, lw=1.25, label=r"$sigma = %3.2f$"%s)

# formula = r"$y(x)=frac{1}{sqrt{2pisigma^2}}e^{-frac{x^2}{2sigma^2}}$"

# ax.text(0.05, 0.80, formula, transform=ax.transAxes, fontsize=20)
# ax.set_xlabel(r"$x$", fontsize=18)
# ax.set_ylabel(r"$y(x)$", fontsize=18)
# ax.legend()
# plt.show()


# import plotly.express as px

# fig = px.line(x=["a","b","c"], y=[1,3,2], title="sample figure")
# print(fig)
# fig.show()


# import plotly.graph_objects as go
# import numpy as np

# x_theo = np.linspace(-4, 4, 100)
# sincx = np.sinc(x_theo)
# x = [-3.8, -3.03, -1.91, -1.46, -0.89, -0.24, -0.0, 0.41, 0.89, 1.01, 1.91, 2.28, 2.79, 3.56]
# y = [-0.02, 0.04, -0.01, -0.27, 0.36, 0.75, 1.03, 0.65, 0.28, 0.02, -0.11, 0.16, 0.04, -0.15]

# fig = go.Figure()
# fig.add_trace(go.Scatter(
#     x=x_theo, y=sincx,
#     name='sinc(x)'
# ))
# fig.add_trace(go.Scatter(
#     x=x, y=y,
#     mode='markers',
#     name='measured',
#     error_y=dict(
#         type='constant',
#         value=0.1,
#         color='purple',
#         thickness=1.5,
#         width=3,
#     ),
#     error_x=dict(
#         type='constant',
#         value=0.2,
#         color='purple',
#         thickness=1.5,
#         width=3,
#     ),
#     marker=dict(color='purple', size=8)
# ))
# fig.show()




# import chart_studio.plotly as py
# import plotly.graph_objects as go
# import pandas as pd

# df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/gapminder2007.csv')

# fig = go.Figure(go.Scatter(x=df.gdpPercap, y=df.lifeExp, text=df.country, mode='markers', name='2007'))
# fig.update_xaxes(title_text='GDP per Capita', type='log')
# fig.update_yaxes(title_text='Life Expectancy')
# py.iplot(fig, filename='pandas-multiple-scatter')


# from matplotlib import pyplot as plt

# import numpy as np
# import matplotlib

# X = np.linspace(1, 4, 1000)
# Y = X
# Y2 = X**2
# Y_log = np.log(X)
# Y_cos = np.cos(np.pi*X)

# matplotlib.rcParams['font.size'] = 16

# plt.plot(X, Y2, X, Y, X, Y_log)
# plt.xlabel('$x \in [1, \infty)$', fontsize=14)
# plt.ylabel('$y$',fontsize=14)
# plt.legend(['$y=x^2$', '$y=x$', '$y=\ln\;x$'], fontsize=12)
# plt.title('Comparison of $\mathcal{O}(x^2)$, $\mathcal{O}(x)$, and $\mathcal{O}(\ln\;x)$', fontsize=16)
# plt.text(3.4, 8, 'Rapid growth\nwhen $y \sim \mathcal{O}(x^2)$', fontsize=14, color='gray')

# ax = plt.gca()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# plt.plot(X, Y2)
# plt.title(r'$y=x^2$');



# import matplotlib.pyplot as plt
# import numpy as np

# plt.rcParams['text.usetex'] = True


# t = np.linspace(0.0, 1.0, 100)
# s = np.cos(4 * np.pi * t) + 2

# fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
# ax.plot(t, s)

# ax.set_xlabel(r'\textbf{time (s)}')
# ax.set_ylabel('\\textit{Velocity (\N{DEGREE SIGN}/sec)}', fontsize=16)
# ax.set_title(r'\TeX\ is Number $\displaystyle\sum_{n=1}^\infty'
#              r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')





# import matplotlib.pyplot as plt
# import numpy as np

# # plt.rcParams['text.usetex'] = True

# fig, ax = plt.subplots()
# # interface tracking profiles
# N = 500
# delta = 0.6
# X = np.linspace(-1, 1, N)
# ax.plot(X, (1 - np.tanh(4 * X / delta)) / 2,    # phase field tanh profiles
#         X, (1.4 + np.tanh(4 * X / delta)) / 4, "C2",  # composition profile
#         X, X < 0, "k--")                        # sharp interface

# # legend
# ax.legend(("phase field", "level set", "sharp interface"),
#           shadow=True, loc=(0.01, 0.48), handlelength=1.5, fontsize=16)

# # the arrow
# ax.annotate("", xy=(-delta / 2., 0.1), xytext=(delta / 2., 0.1),
#             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"))
# ax.text(0, 0.1, r"$\delta$",
#         color="black", fontsize=24,
#         horizontalalignment="center", verticalalignment="center",
#         bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

# # Use tex in labels
# ax.set_xticks([-1, 0, 1])
# ax.set_xticklabels(["$-1$", r"$\pm 0$", "$+1$"], color="k", size=20)

# # Left Y-axis labels, combine math mode and text mode
# ax.set_ylabel(r"\bf{phase field} $\phi$", color="C0", fontsize=20)
# ax.set_yticks([0, 0.5, 1])
# ax.set_yticklabels([r"\bf{0}", r"\bf{.5}", r"\bf{1}"], color="k", size=20)

# # Right Y-axis labels
# ax.text(1.02, 0.5, r"\bf{level set} $\phi$",
#         color="C2", fontsize=20, rotation=90,
#         horizontalalignment="left", verticalalignment="center",
#         clip_on=False, transform=ax.transAxes)

# # Use multiline environment inside a `text`.
# # level set equations
# eq1 = (r"\begin{eqnarray*}"
#        r"|\nabla\phi| &=& 1,\\"
#        r"\frac{\partial \phi}{\partial t} + U|\nabla \phi| &=& 0 "
#        r"\end{eqnarray*}")
# ax.text(1, 0.9, eq1, color="C2", fontsize=18,
#         horizontalalignment="right", verticalalignment="top")

# # phase field equations
# eq2 = (r"\begin{eqnarray*}"
#        r"\mathcal{F} &=& \int f\left( \phi, c \right) dV, \\ "
#        r"\frac{ \partial \phi } { \partial t } &=& -M_{ \phi } "
#        r"\frac{ \delta \mathcal{F} } { \delta \phi }"
#        r"\end{eqnarray*}")
# ax.text(0.18, 0.18, eq2, color="C0", fontsize=16)

# ax.text(-1, .30, r"gamma: $\gamma$", color="r", fontsize=20)
# ax.text(-1, .18, r"Omega: $\Omega$", color="b", fontsize=20)

# plt.show()




import matplotlib.pyplot as plt
import numpy as np

# Settings
A = 6  # Want figures to be A6
plt.rc('figure', figsize=[46.82 * .5**(.5 * A), 33.11 * .5**(.5 * A)])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Create some data
x = np.arange(10)
y = x**2

# Plot
plt.plot(x, y, label=r'Dataset 1: R = t²')
plt.title(r'\LaTeX\ in Labels')
plt.xlabel(r'Time, $t$ [\textmu s]')
plt.text(7, 70, r'$R = t^2$')
plt.legend(fontsize=8)
plt.show()