class Face {
    halfedge = null;
}

class HalfEdge {
    // src and tgt are Vertices
    twin = null;
    target = null;  // Source vertex
    prev = null;    // Previous halfedge
    next = null;    // Next halfedge
    face = null;    // halfedgeident face

    constructor() {

    }

    get_src() {
        return this.twin.target;
    }
}

class Vertex {
    x = null;
    y = null;
    halfedge = null;   // This is the outgoing edge

    constructor(x, y) {
        this.x = x;
        this.y = y;
    }
}

class DCEL {
    vertices = [];
    faces = [];
    halfedges = [];
    outerface = null;
    virtualv = null;

    constructor() {
        this.outerface = Face();
        this.virtualv = Vertex(Infinity, Infinity);
    }

    sort_points_by_x(S) {
        var L = S;
        L.sort(function(a, b) {
            a.x - b.x;
        });
        return L;
    }

    // Precondition: the vertex should completely lie in the face of h
    add_vertex_at(v, h) {
        var u = h.target;
        var h1 = HalfEdge();
        var h2 = HalfEdge();
        var f = h.face;
        v.halfedge = h2;
        h1.twin = h2;
        h2.twin = h1;
        h1.target = v;
        h2.target = u;
        h1.face = f;
        h2.face = f;
        h1.next = h2;
        h2.next = h.next;
        h1.prev = h;
        h2.prev = h1;
        h.next = h1;
        h2.next.prev = h2;
    }

    // v is incident to h.face but not adjacent to h.target
    // the open line segment uv lies completely in h.face.
    split_face(h, v) {
        var f1 = Face();
        var f2 = Face();
        var h1 = HalfEdge();
        var h2 = HalfEdge();
        f1.halfedge = h1;
        f2.halfedge = h2;
        h1.twin = h2;
        h2.twin = h1;
        h1.target = v;
        h2.target = u;
        h2.next = h.next;
        h2.next.prev = h2;
        h1.prev = h;
        h.next = h1;
        var i = h2;
        while (true) {
            i.face = f2;
            if (i.target == v) { break; }
            i = i.next;
        }
        h1.next = i.next;
        h1.next.prev = h1;
        i.next = h2;
        h2.prev = i;
        i = h1;
        do {
            i.face = f1;
            i = i.next;
        } while (i.target != u);

        // delete face f
    }

    split_edge(h) {

    }

    join_face(h) {

    }

    vertices() {
        return this.vertices;
    }

    edges() {

    }

    // A function that returns the orientation of 3 points.
    // 0: colinear. >0: Clockwise. <0: Counter-Clockwise
    static orientation(p1, p2, p3) {
        return (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y);
    }


    // Return the convex hull points of this graph.
    // TODO: Check if the use of orientation is correct.
    // TODO: Modify algorithm so that the return value is actually a graph instead of points
    static convex_hull_pts(points) {
        points.sort(function(a, b) {
            return a.x == b.x ? a.y - b.y : a.x - b.x;
        });

        var lower = [];
        for (var i = 0; i < points.length; i++) {
            while (lower.length >= 2 && this.orientation(lower[lower.length - 2], lower[lower.length - 1], points[i]) <= 0) {
                lower.pop();
            }
            lower.push(points[i]);
        }

        var upper = [];
        for (var i = points.length - 1; i >= 0; i--) {
            while (upper.length >= 2 && this.orientation(upper[upper.length - 2], upper[upper.length - 1], points[i]) <= 0) {
                upper.pop();
            }
            upper.push(points[i]);
        }

        upper.pop();
        lower.pop();
        return lower.concat(upper);
    }
}


// Helper function. Convert a list of points [[0,1],[1,2],...[3,2]] to [{x:0,y:1},...{}]
function tuple2vert(plist) {
    return plist.map(function(d) { return Vertex(d[0], d[1]); })
}

function main() {
    // This block is redundant with the one in lib
    var canvas = d3.select("canvas").node();
    ctx = canvas.getContext("2d"),
        width = canvas.width,
        height = canvas.height;

    var n = 100;
    var pts = d3.range(n)
        .map(function(d) { return new Vertex(Math.random() * width, Math.random() * height); });

    // var chs = DCEL.convex_hull_pts(sitesList);
    var n1 = Math.floor(n);
    var S1 = pts.slice(0, n1);
    var S2 = pts.slice(n1);

    var p1 = DCEL.convex_hull_pts(S1);
    var p2 = DCEL.convex_hull_pts(S2);

    var edges1 = _make_edges_for_convex_hull(p1);
    var edges2 = _make_edges_for_convex_hull(p2);
    var bgs = Voronoi.bridges(p1, p2);

    drawGraph(p1.concat(p2), edges1.concat(edges2).concat(bgs));
}


main();
