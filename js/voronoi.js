// Voronoi.js, source code of CPSC-620 course project
// Author: Yifan Liu
//
// The assumption of entire algorithm: Points differ in x coordinates

// Constant indicating the type of Line
var LINE_SEG = 0;
var LINE_RAY = 1;
var LINE_LINE = 2;

// Represents a segment, ray or line.
// HalfEdge of DCEL extends this class.
// Can be built using from_pts and from_kb, where k and b are the parameters
// of y=kx+b.
class Line {
    k = null;
    b = null;
    p1 = null;
    p2 = null;
    x_domain = [-Infinity, Infinity];

    static from_pts(p1, p2, LINE_TYPE) {
        var l = Line();
        l.p1 = p1; l.p2 = p2;
        if (LINE_TYPE == LINE_SEG) {
            l.x_domain[0] = min([p1.x, p2.x]);
            l.x_domain[1] = max([p1.x, p2.x]);
        } else if (LINE_TYPE == LINE_RAY) {
            // from p1 to infinity
            l.x_domain[0] = p1.x;
        }

        kb = Line.pts2formula(p1, p2);
        l.k = kb[0]; l.b = kb[1];
        return l;
    }

    static from_kb(k, b, x_dom_min, x_dom_max) {
        var l = Line();
        l.k = k; l.b = b;
        x_domain[0] = x_dom_min;
        x_domain[1] = x_dom_max;

        l.p1 = Vertex(x_dom_min, k * x_dom_min + b);
        l.p2 = Vertex(x_dom_max, k * x_dom_max + b);

        if (l.p1.x == Infinity || l.p1.x == -Infinity) l.p1 = null;
        if (l.p2.x == Infinity || l.p2.x == -Infinity) l.p2 = null;

        return l;
    }

    static bisector_of(pt1, pt2, x_min, x_max) {
        var k = (v1.x - v2.x) / (v2.y - v1.y);
        var midp = midpoint(v1, v2);
        var b = midp.y - k * midp.x;

        return Line.from_kb(k, b, x_min, x_max);
    }

    get_pt(x) {
        return Vertex(x, this.k * x + this.b);
    }

    intersection(line) {
        if (this.k == line.k) return null;
        var x = (this.b - line.b) / (line.k - this.k);
        var y = this.k * x + this.b;

        if ((x >= this.x_domain[0] && x <= this.domain[1]) && (x >= line.x_domain[0] && x <= line.domain[1])) {
            return Vertex(x, y);
        }
        return null;
    }

    // Return a y=kx+b formula for two points p, q
    static pts2formula(p, q) {
        if (p.x == q.x) return [Infinity, p.x];
        var k = (p.y - q.y) / (p.x - q.x);
        var b = p.y - k * p.x;
        return [k, b];
    }
}

class Voronoi extends DCEL {
    // Every site uniquely maps to a face of DCEL.
    // Sites are also class Vertexes but contains id.
    // ids will be set in constructor.
    sites = new Map();

    // S is a set of Vertices
    constructor(S) {
        for (var i = 0; i < S.length; i++) {
            S[i].id = i;
            this.sites.set(i, S[i]);
        }
        rec_construct_vor(S);
    }

    vertices() {
        return this.vertices;
    }

    sites() {
        return this.sites;
    }

    rec_construct_vor(pts) {
        var n = pts.length;
        if (n <= 2) {
            // Trivial case
        } else {
            // Split & Merge
            // TODO: This sorting should be done only once
            pts.sort(function(a, b) { return a.x - b.x });

            var n1 = Math.floor(n);
            var S1 = pts.slice(0, n1);
            var S2 = pts.slice(n1);
            var vorDiag1 = this.rec_construct_vor(S1);
            var vorDiag2 = this.rec_construct_vor(S2);
            return this.merge_voronoi(vorDiag1, vorDiag2);
        }
    }

    merge_voronoi(vorL, vorR) {
        // Construct convex hull of two subgraphs
        var chsL = DCEL.convex_hull_pts(vorL.sites());
        var chsR = DCEL.convex_hull_pts(vorR.sites());
        var sigma_points = [];

        // Find lower common support
        var bdgs = this.bridges(chsL, chsR);
        var ubdg = bdgs[0]; // Upper bridge
        var lbdg = bdgs[1]; // Lower bridge

        // First find the bisector of the lower bridge and work our way up
        var sigma_low_inf = Line.bisector_of(lbdg[0], lbdg[1], -Infinity, Infinity);
        var cur_bdg = lbdg;
        var cur_bisector = sigma_low_inf;

        var sitel = cur_bdg[0]; var sitel_ind_in_ch = chsL.indexOf(sitel);
        var siter = cur_bdg[1]; var siter_ind_in_ch = chsR.indexOf(siter);
        var indl = 0;
        var indr = 0;
        var il = 0;
        var ir = 0;

        // while current working bridge is not upper bridge
        while (!arr_eq(cur_bdg, ubdg)) {
            var int_l = null;
            var int_r = null;

            // Traverse all voronoi half edges defined by sitel and siter
            face_of_sitel = this.site_map[]
            do {
                var edge = vedgesl[il];
                var int_l = cur_bisector.intersection(edge);
                if (int_l != null) {
                    break;
                }
                il = (il + 1) % vlenl;
            } while (il != indl);

            do {
                var edge = vedgesr[ir];
                var int_r = cur_bisector.intersection(edge);
                if (int_r != null) {
                    break;
                }
                ir = (ir + 1) % vlenr;
            } while (ir != indr);

            if (int_l.y < int_r.y) {
                sigma_points.push(int_r);
                siter = chsR[(chsR.length + siter_ind_in_ch - 1) % chsR.length];
                ind_r = 0;
                ir = 0;

                // Probably this is not required
                ind_l = 0;
                il = 0;
            } else {
                sigma_points.push(int_l);
                sitel = chsL[(sitel_ind_in_ch + 1) % chsL.length];
                ind_l = 0;
                il = 0;

                ind_r = 0;
                ir = 0;
            }
        }

    }

    // Return a collection of edges, in ccw order,
    // of site p.
    vedges_of_site(p) {
        // TODO: How would you implement it without using DCEL?
        // For each site, store a list of vertices that are closest to it
        // So at the merge stage, you have to figure out some ways to merge the hull
        // sites as well
        return site_map[]
    }

    static bridges(p1, p2) {
        // p1 and p2 are two convex polygons
        // p1 should be completely on the left of p2
        // returns the upper and lower bridge of p1 and p2
        var u, v, ui, vi, ind1, ind2;
        var u_tpl = max(x_comp, p1); u = u_tpl[0]; ind1 = u_tpl[1];
        var v_tpl = min(x_comp, p2); v = v_tpl[0]; ind2 = v_tpl[1];
        var ui = ind1;
        var vi = ind2;
        var n1 = p1.length;
        var n2 = p2.length;
        var done = false;

        while (!done) {
            done = true;
            while (clw(p1[ui], p2[vi], p2[(n2 + vi - 1) % n2])) {
                vi = (n2 + vi - 1) % n2;
                done = false;
            }
            while (ccw(p2[vi], p1[ui], p1[(ui + 1) % n1])) {
                ui = (ui + 1) % n1;
                done = false;
            }
        }
        var lower_bridge = [p1[ui], p2[vi]];

        ui = ind1;
        vi = ind2;
        done = false;

        while (!done) {
            done = true;
            while (clw(p2[vi], p1[ui], p1[(n1 + ui - 1) % n1])) {
                ui = (n1 + ui - 1) % n1;
                done = false;
            }
            while (ccw(p1[ui], p2[vi], p2[(vi + 1) % n2])) {
                vi = (vi + 1) % n2;
                done = false;
            }
        }
        var upper_bridge = [p1[ui], p2[vi]];
        return [upper_bridge, lower_bridge];
    }
}

function min(fn, iterable) {
    var _min = iterable[0];
    var _min_i = 0;
    for (var i = 0; i < iterable.length; i++) {
        if (fn(iterable[i], _min) < 0) {
            _min = iterable[i];
            _min_i = i;
        }
    }
    return [_min, _min_i];
}

function max(fn, iterable) {
    var _max = iterable[0];
    var _max_i = 0;
    for (var i = 0; i < iterable.length; i++) {
        if (fn(iterable[i], _max) > 0) {
            _max = iterable[i];
            _max_i = i;
        }
    }
    return [_max, _max_i];
}

function x_comp(a, b) { return a.x - b.x; }


class Face {
    halfedge = null;
}

class HalfEdge extends Line {
    // src and tgt are Vertices
    twin = null;
    target = null;  // Source vertex
    prev = null;    // Previous halfedge
    next = null;    // Next halfedge
    face = null;    // halfedgeident face

    get_src() {
        return this.twin.target;
    }
}

class Vertex {
    x = null;
    y = null;
    halfedge = null;   // This is the outgoing edge
    id = null;

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

    static sort_points_by_x(S) {
        var L = S;
        L.sort(function(a, b) {
            return a.x - b.x;
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
    var canvas = d3.select("canvas").node();
    ctx = canvas.getContext("2d"),
        width = canvas.width,
        height = canvas.height;

    var n = 100;
    var random_pts = d3.range(n)
        .map(function(d) { return new Vertex(Math.random() * width, Math.random() * height); });

    // var chs = DCEL.convex_hull_pts(sitesList);
    var pts = DCEL.sort_points_by_x(random_pts);
    var n1 = Math.floor(n / 2);
    var S1 = pts.slice(0, n1);
    var S2 = pts.slice(n1);

    var p1 = DCEL.convex_hull_pts(S1);
    var p2 = DCEL.convex_hull_pts(S2);

    var edges1 = _make_edges_for_convex_hull(p1);
    var edges2 = _make_edges_for_convex_hull(p2);
    var bgs = Voronoi.bridges(p1, p2);

    drawGraph(pts, edges1.concat(edges2).concat(bgs));
}

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
    for (var i = 0; i < edges.length; i++)
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

function clw(p1, p2, p3, report_colinear=false) {
    if (report_colinear) return DCEL.orientation(p1, p2, p3) >= 0;
    return DCEL.orientation(p1, p2, p3) > 0;
}

function ccw(p1, p2, p3, report_colinear=false) {
    if (report_colinear) return DCEL.orientation(p1, p2, p3) <= 0;
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

function arr_eq(arr1, arr2) {
    if (arr1.length != arr2.length) return false;
    for (var i = 0; i < arr1.length; i++) {
        if (arr1[i] != arr2[i])
            return false;
    }
    return true;
}

function midpoint(v1, v2) {
    return Vertex(
        (v1.x + v2.x) / 2,
        (v1.y + v2.y) / 2
    );
}

// Return the perpendicular bisector of two vertex.
// Return the k and b of formula y=kx+b
function bisector(v1, v2) {
    var k = (v1.x - v2.x) / (v2.y - v1.y);
    var midp = this.midpoint(v1, v2);
    var b = midp.y - k * midp.x;
    return [k, b];
}


// Find intersection of two lines linea and lineb.
function line_intersection(linea, lineb) {
    if (linea[2] == this.LINE_PTS) lina = pts2formula(linea[0], linea[1]).concat([this.LINE_KB, ])
    var k = linea[0];
    var b = linea[1];
    var kprime = lineb[0];
    var bprime = lineb[1];

    if (k == kprime) return null;
    var x = (b - bprime) / (kprime - k);
    var y = k * x + b;

    return Vertex(x, y);
}

// Find intersection of seg1 and seg2.
// seg1 and seg2 may be segment or ray or line.
// Ray is denoted [Vertex1, Vertex2, Infinity], from vertex1 to vertex2.
// And line, if also denoted this way, should be [Vertex1, Vertex2, Infinity, Infinity]
function ray_intersection(seg1, seg2) {
    var line_int_pt = line_intersection(pts2formula(seg1[0], seg1[1]), pts2formula(seg2[0], seg2[1]));
    if (line_int_pt == null) return null;

    if (pt_valid(seg1, p) && pt_valid(seg2, p)) return line_int_pt;
    return null;
}

// Return if a point on seg is valid
// p has to sit in the line denoted by seg
function pt_valid(seg, p) {
    // If this is a line, then valid
    if (seg[3] != undefined) return true;

    // If this is a ray, make sure the point does not lie beyond the source point
    // *p A|---->B make sure A does not sit between p and B
    if (seg[2] != undefined) {
        return !in_the_middle_of_line(p, seg[1], seg[0])
    }

    // Then seg is a segment. Make sure p is in between AB.
    return in_the_middle_of_line(seg[0], seg[1], p);
}

// Return if p is on the segment AB.
function in_the_middle_of_line(A, B, P) {
    var vec_ap = [P.x - A.x, P.y - A.y];
    var vec_bp = [P.x - B.x, P.y - A.y];

    return dot(vec_ap, vec_bp) < 0;
}

function dot(p, q) {
    return p[0] * q[0] + p[1] * q[1];
}

main();

