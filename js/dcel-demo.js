var canvas = d3.select("canvas");
var ctx = canvas.getContext("2d");
var width = canvas.width;
var height = canvas.height;

var sitesList = d3.range(100).map(function(d) { return [Math.random() * width, Math.random() * height];  });

var G = DCEL();
G.set_vertices(tuple2vert(sitesList));
CH_G = G.convex_hull();

// G is a DCEL
function draw_DCEL(G) {
    var vertices = G.vertices();
    var edges = G.edges();

    ctx.clearRect(0, 0, width, height);

    for(var i = 0; i < vertices.length; i++) {
        draw_vertex(vertices[i]);
    }

    for(var i = 0; i < edges.length; i++) {
        draw_edge(edges[i]);
    }
}

function draw_vertex(v) {
    ctx.moveTo(v.x + 2.5, v.y);
    ctx.arc(v.x, v.y, 2.5, 0, 2 * Math.PI, false);
}

function draw_edge(e) {
    // ctx.moveTo(link.source[0], link.source[1]);
    // ctx.lineTo(link.target[0], link.target[1]);
}
