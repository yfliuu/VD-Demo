function drawGraph(vertices, edges) {
    var canvas = d3.select("canvas").node();
    ctx = canvas.getContext("2d"),
        width = canvas.width,
        height = canvas.height;

    ctx.clearRect(0, 0, width, height);

    ctx.beginPath();
    for (var i = 0; i < vertices.length; i++) drawSite(vertices[i]);
    ctx.fillStyle = "#000";
    ctx.fill();
    ctx.strokeStyle = "#fff";
    ctx.stroke();

    ctx.beginPath();
    for (var i = 0; i < chs.length; i++)
        drawEdgeBetween(edges[i][0], edges[i][1]);
    ctx.strokeStyle = "rgba(0,0,0,0.2)";
    ctx.stroke();

}

function drawSite(site) {
    ctx.moveTo(site.x + 2.5, site.y);
    ctx.arc(site.x, site.y, 2.5, 0, 2 * Math.PI, false);
}

function drawLink(link) {
    ctx.moveTo(link.source.x, link.source.y);
    ctx.lineTo(link.target.x, link.target.y);
}

function drawEdgeBetween(p1, p2) {
    ctx.moveTo(p1.x, p1.y);
    ctx.lineTo(p2.x, p2.y);
}

function clw(p1, p2, p3) {
    return DCEL.orientation(p1, p2, p3) > 0;
}

function ccw(p1, p2, p3) {
    return DCEL.orientation(p1, p2, p3) < 0;
}

function slope(p1, p2, abs=false) {
    var slope = (p2.y - p1.y) / (p2.x - p1.x);
    if (abs) return Math.abs(slope);
    return slope;
}

function _make_edges_for_convex_hull(pts) {
    var edges = [];
    for (var i = 0; i < pts.length - 1; i++) {
        edges[i] = [pts[i], pts[i + 1]];
    }
    edges.push([pts[pts.length - 1], pts[0]]);
    return edges
}
