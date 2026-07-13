from feynman import Diagram, RegularBubbleOperator, Vertex

diagram = Diagram(figsize=(6.0, 6.0))


def lepton(xy) -> Vertex:
    v = diagram.vertex(xy=xy)
    diagram.line(diagram.vertex(xy=(0.0, v.y), marker=""), v).text(
        "$l(k)$", y=0.04, t=0.3, fontsize=30
    )
    diagram.line(v, diagram.vertex(xy=(1.0, v.y), marker="")).text(
        "$l'(k-q)$", y=0.04, fontsize=30
    )
    return v


def pdf(xy: tuple[float, float]) -> Vertex:
    k = 24
    op = RegularBubbleOperator(k, xy, 0.1, 0.25)
    op.text(r"$\phi_{q}(\xi)$")
    diagram.add_operator(op)
    diagram.line(
        diagram.vertex(xy=(0.0, op.vertices[k // 4 + 1].y), marker=""),
        op.vertices[k // 4 + 1],
    )
    diagram.line(
        diagram.vertex(xy=(0.0, op.vertices[k // 4].y), marker=""), op.vertices[k // 4]
    )
    diagram.line(
        diagram.vertex(xy=(0.0, op.vertices[k // 4 - 1].y), marker=""),
        op.vertices[k // 4 - 1],
    ).text("$h(P)$", y=0.04, t=0.3, fontsize=30)
    diagram.line(
        op.vertices[k // 4 * 3 + 1],
        diagram.vertex(xy=(0.85, op.vertices[k // 4 * 3 + 1].y + 0.05), marker=""),
    )
    v = diagram.vertex(xy=(0.85, op.vertices[k // 4 * 3].y), marker="")
    diagram.line(op.vertices[k // 4 * 3], v)
    v.text("$X$", x=0.05, y=0.0, fontsize=30)
    diagram.line(
        op.vertices[k // 4 * 3 - 1],
        diagram.vertex(xy=(0.85, op.vertices[k // 4 * 3 - 1].y - 0.05), marker=""),
    )
    return op.vertices[k - 2]


# incoming splittings
pos_up = (0.4, 0.8)
inB = lepton(pos_up)
inq = pdf((0.4, 0.15))
v = diagram.vertex(xy=(0.6, 0.5))
diagram.line(inB, v, style="wiggly", phase=1.5).text(
    "$V(q)$", y=-0.2, t=0.4, fontsize=30
)
diagram.line(inq, v).text(r"$q(\xi P)$", y=0.2, t=0.15, fontsize=30)
outq = diagram.vertex(xy=(1.0, 0.5), marker="")
diagram.line(v, outq).text(r"$q'(\xi P+q)$", y=0.04, t=0.15, fontsize=30)

diagram.ax.hlines(0.67, 0, 1, linestyle="dashed", color="red")


diagram.plot()
diagram.savefig("LODIS.pdf")
