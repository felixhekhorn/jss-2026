from feynman import Diagram, RegularBubbleOperator, Vertex

diagram = Diagram(figsize=(6.0, 6.0))


def lepton(xy) -> Vertex:
    v = diagram.vertex(xy=xy)
    diagram.line(diagram.vertex(xy=(0.0, v.y), marker=""), v).text(
        "$e^-(k)$", y=0.04, t=0.3, fontsize=30
    )
    diagram.line(v, diagram.vertex(xy=(1.0, v.y * 1.1), marker="")).text(
        "$e^-(k')$", y=0.05, t=0.6, fontsize=30
    )
    return v


def pdf(xy: tuple[float, float]) -> Vertex:
    k = 24
    op = RegularBubbleOperator(k, xy, 0.1, 0.25)
    op.text("$f$")
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
    ).text("$p(P)$", y=0.04, t=0.3, fontsize=30)
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
    return op.vertices[0]


# incoming splittings
pos_up = (0.4, 0.8)
inB = lepton(pos_up)
inG = pdf((0.5, 0.15))
diagram.line(inB, inG, style="wiggly", phase=1.5).text(
    "$\\gamma(q)$", y=0.05, fontsize=30
)

diagram.plot()
diagram.savefig("DIS.pdf")
