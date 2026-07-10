import matplotlib.pyplot as plt
from feynman import Diagram

fig = plt.figure(figsize=(10.0, 10.0))
ax = fig.add_axes([0, 0, 1, 1], frameon=False)
gluon_style = dict(style="linear loopy", xamp=0.025, yamp=0.035, nloops=7)

diagram = Diagram(ax)
inB = diagram.vertex(xy=(0, 1), marker="")
inQ = diagram.vertex(xy=(0, 0), marker="")
vB = diagram.vertex(xy=(0.5, 0.7))
vg = diagram.vertex(xy=(0.5, 0.3))
outQ = diagram.vertex(xy=(1, 0.7), marker="")
outG = diagram.vertex(xy=(1, 0.3), marker="")

b = diagram.line(inB, vB, style="wiggly")
inQ = diagram.line(inQ, vg)
diagram.line(vg, vB)
g = diagram.line(vg, outG, **gluon_style)
outQ = diagram.line(vB, outQ)

# b.text(r"$b(q)$",t=.4,y=.05,fontsize=40)
# inQ.text(r"$q(p_1)$",t=.3,fontsize=40)
# outQ.text(r"$q'(p_2)$",t=.6,fontsize=40)
# g.text(r"$g(k)$",t=.6,fontsize=40)

diagram.plot()
plt.savefig("nlo-q-rt.pdf")
