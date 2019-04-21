// Voronoi.js, source code of CPSC-620 course project
// Author: Yifan Liu
//
// The assumption of entire algorithm: Points differ in x coordinates

// Constant indicating the type of Line
var LINE_SEG = 0;
var LINE_RAY = 1;
var LINE_LINE = 2;

var HALF_EDGE_AUTOID = 0;
var SITE_AUTOID = 0;
var FACE_AUTOID = 0;

var INFINITY_APPROXIMATE = 10000;
var CLOSE_THRESHOLD = 0.000001;

var REACH_NONE = 0;
var REACH_BACKWARD = 1;
var REACH_FORWARD = 2;
var REACH_BOTH = 3;

var vertices = [];
var p = [];
var left = [];
var right = [];
var leftzip = [];
var rightzip = [];
var leftedge = [];
var rightedge = [];
var tos = 0;

var canvas = d3.select("canvas").node();
ctx = canvas.getContext("2d"),
    width = canvas.width,
    height = canvas.height;

// Represents a segment, ray or line.
// HalfEdge of DCEL extends this class.
// Can be built using from_pts and from_kb, where k and b are the parameters
// of y=kx+b.
class Line {
    k = null;
    b = null;
    p1 = null;
    p2 = null;
    reach = null;
    x_min = -INFINITY_APPROXIMATE;
    x_max = INFINITY_APPROXIMATE;

    constructor() {
        this.k = null;
        this.b = null;
        this.p1 = null;
        this.p2 = null;
        this.x_min = -INFINITY_APPROXIMATE;
        this.x_max = INFINITY_APPROXIMATE;
    }

    static copy(l) {
        var l = new Line();
        this.k = l.k;
        this.b = l.b;
        this.p1 = l.p1;
        this.p2 = l.p2;
        this.reach = l.reach
        this.x_min = l.x_min;
        this.x_max = l.x_max;
        return l;
    }

    is_length0() {
        return (this.p1.x == this.p2.x && this.p1.y == this.p2.y);
    }

    // Cut the line from p
    cut_from(p, from, to, extend=false) {
        if (p.x * this.k + this.b - p.y > 0.00001) {
            console.log('cut_from: p is not on the line');
            return null;
        }
        if (extend) {

        } else {
            if (p.x <= this.x_min || p.x >= this.x_max) return this;
            from = from == 1 ? this.p1 : from;
            to = to == 2 ? this.p2 : to;
            console.log('fromx and tox', from.x, to.x);
            this.p1 = from;
            this.p2 = to;
            this.x_min = from.x;
            this.x_max = to.x;
            return this;
        }
    }

    // Reverse the ray
    reverse() {
        var x_min = this.x_min; var x_max = this.x_max;
        this.x_min = x_max; this.x_max = x_min;

        if (this.x_min == INFINITY_APPROXIMATE) this.x_min = -INFINITY_APPROXIMATE;
        if (this.x_max == -INFINITY_APPROXIMATE) this.x_max = INFINITY_APPROXIMATE;

        this.p1 = new Vertex(this.x_min, this.x_min * this.k + this.b);
        this.p2 = new Vertex(this.x_max, this.x_max * this.k + this.b);
        return this;
    }

    static from_pts(p1, p2, LINE_TYPE=LINE_SEG) {
        var l = new Line();
        var kb = Line.pts2formula(p1, p2);
        l.p1 = p1; l.p2 = p2;
        l.k = kb[0]; l.b = kb[1];

        if (LINE_TYPE == LINE_SEG) {
            l.x_min = Math.min([p1.x, p2.x]);
            l.x_max = Math.max([p1.x, p2.x]);
        } else if (LINE_TYPE == LINE_RAY) {
            if (p1.x < p2.x) {
                l.x_min = p1.x;
            } else {
                l.x_max = p1.x;
            }
        } // Else: LINE_LINE. No need to modify domains.

        l.p1 = new Vertex(l.x_min, l.k * l.x_min + l.b);
        l.p2 = new Vertex(l.x_max, l.k * l.x_max + l.b);
        return l;
    }

    static perpdiv(p1, p2)
    {
        var midx = (p1.x + p2.x) / 2;
        var midy = (p1.y + p2.y) / 2;
        var dx = p2.x - p1.x;
        var dy = p2.y - p1.y;
        var div = sqrt(dx*dx + dy*dy);
        var ux = dx / div;
        var uy = dy / div;
        var lx = ux * DISTANT;
        var ly = uy * DISTANT;
        var r = new Line();
        r.p1.x = midx + ly;
        r.p1.y = midy - lx;
        r.p2.x = midx - ly;
        r.p2.y = midy + lx;
        r.reach = REACH_BOTH;
        return r;
    }

    static perpthrough(c, p1, p2) {
        var dx = p2.x - p1.x;
        var dy = p2.y - p1.y;
        var div = Math.sqrt(dx*dx + dy*dy);
        var ux = dx / div;
        var uy = dy / div;
        var lx = ux * DISTANT;
        var ly = uy * DISTANT;

        var r = new Line();
        r.p1.x = c.x;
        r.p1.y = c.y;
        r.p2.x = c.x - ly;
        r.p2.y = c.y + lx;
        r.reach = REACH_FORWARD;
        return r;
    }

    static from_kb(k, b, x_dom_min, x_dom_max) {
        var l = new Line();
        l.k = k; l.b = b;
        l.p1 = new Vertex(x_dom_min, x_dom_min * k + b);
        l.p2 = new Vertex(x_dom_max, x_dom_max * k + b);
        return l;
    }

    static bisector_of(v1, v2, x_min=-Infinity, x_max=Infinity) {
        var k = (v1.x - v2.x) / (v2.y - v1.y);
        var midp = midpoint(v1, v2);
        var b = midp.y - k * midp.x;

        return Line.from_kb(k, b, -INFINITY_APPROXIMATE, INFINITY_APPROXIMATE);
    }

    get_pt(x) {
        return Vertex(x, this.k * x + this.b);
    }

    static intersect(r, a, b) {
        var  p1 = a.p1; var p2 = a.p2;
        var  q1 = b.p1; var q2 = b.p2;

        var  dpy = p2.y - p1.y;
        var  dpx = p2.x - p1.x;
        var  dqy = q2.y - q1.y;
        var  dqx = q2.x - q1.x;

        var  dy = p1.y - q1.y;
        var  dx = p1.x - q1.x;

        var  div = dqy * dpx - dpy * dqx;
        var  pnum = dqx * dy - dqy * dx;
        var  qnum = dpx * dy - dpy * dx;

        var  up, uq, ep, eq;

        if(div == 0) {

            return false;
        }

        up = pnum / div;
        uq = qnum / div;

        ep = CLOSE_THRESHOLD / Math.sqrt(dpx * dpx + dpy * dpy);
        eq = CLOSE_THRESHOLD / Math.sqrt(dqx * dqx + dqy * dqy);

        if((a.reach & REACH_BACKWARD || up >= -ep) &&
            (a.reach & REACH_FORWARD  || up <= 1+ep) &&
            (b.reach & REACH_BACKWARD || uq >= -eq) &&
            (b.reach & REACH_FORWARD  || uq <= 1+eq)
        ) {
            r.x = p1.x + up * dpx;
            r.y = p1.y + up * dpy;
            return true;
        }

        return false;
    }

    intersection(line) {
        if (this.k == line.k) return null;
        var x = (this.b - line.b) / (line.k - this.k);
        var y = this.k * x + this.b;

        if ((x >= this.x_min && x <= this.x_max) && (x >= line.x_min && x <= line.x_max)) {
            return new Vertex(x, y);
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
    id = null;

    constructor() {
        this.id = FACE_AUTOID++;
    }
}

class HalfEdge extends Line {
    // src and tgt are Vertices
    twin = null;
    target = null;  // Target vertex
    source = null;  // Source vertex, for convenience
    prev = null;    // Previous halfedge
    next = null;    // Next halfedge
    face = null;    // halfedgeident face
    divline = null;
    id = null;

    constructor() {
        super();
        this.id = HALF_EDGE_AUTOID++;
    }

    pair_with(edge) {
        this.twin = edge;
        edge.twin = this;
    }

    static copy(source, target, divline) {
        var edge = new HalfEdge();
        edge.source = source;
        edge.target = target;
        edge.divline = divline;
        return edge;
    }

    static from_line(line) {
        var he = new HalfEdge();
        he.k = line.k;
        he.b = line.b;
        he.p1 = line.p1;
        he.p2 = line.p2;
        he.x_domain = line.x_domain;
        return he;
    }

    get_src() {
        return this.twin.target;
    }

    eq(edge) {
        return self.id == edge.id;
    }

    is_void() {
        return this.divline.is_length0();
    }

}

class Vertex {
    x = null;
    y = null;
    halfedge = null;   // This is the outgoing edge
    id = null;
    next = null;
    prev = null;
    head_edge = null;
    tail_edge = null;

    constructor(x, y) {
        this.x = x;
        this.y = y;
        this.id = SITE_AUTOID++;
    }

    eq(p) {
        return this.x == p.x && this.y == p.y;
    }

    is_virtual() {
        return this.x == INFINITY_APPROXIMATE || this.x == -INFINITY_APPROXIMATE ||
            this.y == INFINITY_APPROXIMATE || this.y == -INFINITY_APPROXIMATE;
    }

    to_drawing() {
        if (!this.is_virtual) return this;
        var x = this.x < 0 ? -INFINITY_APPROXIMATE : INFINITY_APPROXIMATE;
        return new Vertex(x, x * this.k + this.b);
    }
}

class DCEL {
    vertices = [];
    faces = [];
    halfedges = [];
    outerface = null;
    virtualv = null;  // This is the starting virtual vertex.
    virtualvs = [];    // For simplicity, we maintain a list of all virtual vertices. We don't merge them into one.

    constructor() {
        this.outerface = new Face();
        this.virtualv = new Vertex(Infinity, Infinity);

        // We connect virtual vertex to itself as a starting point
        // var h1 = new HalfEdge();
        // var h2 = new HalfEdge();
        // h1.twin = h2;
        // h2.twin = h1;
        // h1.target = this.virtualv;
        // h2.target = this.virtualv;
        // h1.next = h2;
        // h2.next = h1;
        // h1.prev = h2;
        // h2.prev = h1;
        // h1.face = this.outerface;
        // h2.face = this.outerface;
        // this.virtualv.halfedge = h1;
        // this.outerface.halfedge = h1;

        // this.halfedges.push(h1);
        // this.halfedges.push(h2);
        // this.faces.push(this.outerface);
        // this.vertices.push(this.virtualv);
        this.virtualvs.push(this.virtualv);
    }

    static sort_points_by_x(S) {
        var L = S;
        L.sort(function(a, b) {
            return a.x - b.x;
        });
        return L;
    }

    // Precondition: the vertex should completely lie in the face of h
    add_vertex_at(h, v, line=null) {
        var u = h.target;
        var h1, h2;
        if (line) {
            h1 = HalfEdge.from_line(line);
            h2 = HalfEdge.from_line(line);
        } else {
            h1 = new HalfEdge();
            h2 = new HalfEdge();
        }
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
        v.halfedge = h2;

        if (v.x == Infinity) this.virtualvs.push(v); else this.vertices.push(v);
        this.halfedges.push(h1);
        this.halfedges.push(h2);
    }

    // v is incident to h.face but not adjacent to h.target
    // the open line segment uv lies completely in h.face.
    split_face(h, v) {
        var f = h.face;
        var f1 = new Face();
        var f2 = new Face();
        var h1 = new HalfEdge();
        var h2 = new HalfEdge();
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

        // Remove f
        this.faces.filter(function(obj){ obj.id != f.id; });
        this.faces.push(f1);
        this.faces.push(f2);
        this.halfedges.push(h1);
        this.halfedges.push(h2);
    }

    // Split an edge by a vertex v.
    split_edge(h, w) {
        if (w.x * h.k + h.b != w.y) {
            Console.log('split_edge: point not on edge!');
            Console.log('point: ' + w);
            Console.log('line: ' + h);
        }
        var h1 = new HalfEdge();
        var h2 = new HalfEdge();
        var k1 = new HalfEdge();
        var k2 = new HalfEdge();
        var u = h.target;
        var v = h.twin.target;
        h1.next = h2;
        h1.prev = h.prev;
        h1.target = w;
        h1.twin = k2;
        h1.face = h.face;
        h.prev.next = h1;
        k2.target = v;
        k2.next = h.twin.next;
        k2.twin = h1;
        k2.prev = k1;
        k2.face = h.twin.face;
        h.twin.next.prev = k2;
        h2.twin = k1;
        h2.target = u;
        h2.next = h.next;
        h2.prev = h1;
        h2.face = h.face;
        h.next.prev = h2;
        k1.target = w;
        k1.twin = h2;
        k1.prev = h.twin.prev;
        k1.next = k2;
        k1.face = h.twin.face;
        h.twin.prev.next = k1;

        this.halfedges.push(h1);
        this.halfedges.push(h2);
        this.halfedges.push(k1);
        this.halfedges.push(k2);
    }

    join_face(h) {

    }

    vertices() {
        return this.vertices;
    }

    // A function that returns the orientation of 3 points.
    // 0: colinear. >0: Clockwise. <0: Counter-Clockwise
    static orientation(p1, p2, p3) {
        var r = (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y);
        if (r < 0) return -1;
        else if (r == 0) return 0;
        return 1
    }

    // Return the convex hull points of this graph.
    // TODO: Modify algorithm so that the return value is actually a graph instead of points
    static convex_hull_pts(points) {
        if (points.length <= 1) {
            return points;
        }
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

    static merge_hull(l, r) {
        var a, b, c, d, c2, d2;

        a = c = l;  b = d = r;

        c2 = vertices[c].prev;
        d2 = vertices[d].next;
        if(p[c].x == p[d].x && p[c2].x == p[d2].x) {

            c = c2;
            d = d2;
        }
        else {
            upper_line(c, d);
        }

        lower_line(a, b);

        vertices[c].next = d;
        vertices[d].prev = c;
        vertices[b].next = a;
        vertices[a].prev = b;

        // *l = a;  *r = b;
        return [a, b];
    }
}

class Voronoi extends DCEL {
    // Every site uniquely maps to a face of DCEL.
    // Sites are also class Vertexes but contains id.
    // ids will be set in constructor.
    sites = null;
    site_map = null;

    build(points) {
        var r = new Voronoi();
        for (var i = 0; i < points.length; i++) {
            this.sites.push(points[i]);
            this.vertices.push(points[i]);
        }
        for (var i = 0; i < points.length; i++) {
            this.vertices[i].next = this.vertices[i].prev = i;
            this.vertices[i].head_edge = this.vertices[i].tail_edge = null;
        }
        vertices = this.vertices;
        p = this.vertices.slice(0);
        var _size = this.sites.length * 2 + 1;
        leftzip = new Array(_size);
        rightzip = new Array(_size);
        leftedge = new Array(_size);
        rightedge = new Array(_size);
        left = new Array(_size);
        right = new Array(_size);
        if (this.sites.length >= 2) Voronoi.voronoi_rec(0, this.sites.length - 1);

        for (var i = 0; i < this.vertices.length; i++) {
            var edge, head;
            edge = head = this.vertices[i].head_edge;
            if (edge) do {
                if (edge.target > i) l[j++] = edge.divline;
                edge = edge.next;
            } while (!edge.eq(head));
        }
        return r;
    }

    static voronoi_rec(l, r) {
        var n_points = 1 + r - l;
        var m;

        if (n_points <= 3) {
            Voronoi.voronoi_base(l, r, n_points);
            return;
        }

        m = (l + r) / 2;
        Voronoi.voronoi_rec(l, m);
        Voronoi.voronoi_rec(m + 1, r);
        Voronoi.voronoi_merge(m);
    }

    // S is a set of Vertices
    constructor() {
        super();
        this.sites = [];
    }

    _build(S) {
        S.sort(function(a, b) { return a.x - b.x });
        return Voronoi.rec_construct_vor(S);
    }

    edges() {
        return this.halfedges;
    }

    static voronoi_base(l, r, n_points) {
        var q = p[l];
        var lineA, lineB, lineC;

        if (n_points == 3) {
            lineA = Line.perpdiv(q[0],q[1]);
            lineB = Line.perpdiv(q[1],q[2]);
            var c = new Vertex();
            if(Line.intersect(c,lineA,lineB)) {
                if (c.y <= q[1].y) {
                    lineA = Line.perpthrough(c, q[0], q[1]);
                    lineB = Line.perpthrough(c, q[1], q[2]);
                    lineC = Line.perpthrough(c, q[2], q[0]);
                    construct3(l, l+1, l+2, lineA, lineB, lineC);
                }
                else {
                    lineA = Line.perpthrough(c, q[1], q[0]);
                    lineB = Line.perpthrough(c, q[2], q[1]);
                    lineC = Line.perpthrough(c, q[0], q[2]);
                    construct3(l, l+2, l+1, lineC, lineB, lineA);
                }
            }
            else {
                construct2(l, l+1, l+2, lineA, lineB);
            }
        }
        else {
            lineA = Line.perpdiv(q[0], q[1]);
            construct1(l, l+1, lineA);
        }
    }

    static trivial_voronoi_3(p1, p2, p3) {
        // p1 p2 p3 should be sorted
        var v = new Voronoi();
        var mid12 = midpoint(p1, p2);
        var mid13 = midpoint(p1, p3);
        var mid23 = midpoint(p2, p3);
        var vvert = Line.bisector_of(p1, p2).intersection(Line.bisector_of(p2, p3));
        var bis12, bis13, bis23;

        bis12 = Line.from_pts(vvert, mid12, LINE_RAY);
        bis13 = Line.from_pts(vvert, mid13, LINE_RAY);
        bis23 = Line.from_pts(vvert, mid23, LINE_RAY);
        if (!in_triangle(p1, p2, p3, vvert)) {
            var lst = [distance_square(p1, p2), distance_square(p1, p3), distance_square(p2, p3)];
            var tpl = max(function(a, b) { return a - b; }, lst);
            if (tpl[1] == 0) bis12 = bis12.reverse();
            if (tpl[1] == 1) bis13 = bis13.reverse();
            if (tpl[1] == 2) bis23 = bis23.reverse();
        }
        v.sites.set(p1.id, p1);
        v.sites.set(p2.id, p2);
        v.sites.set(p3.id, p3);
        v.vertices.push(vvert);

        var h1 = HalfEdge.from_line(bis12);
        var h2 = HalfEdge.from_line(bis12);
        h1.twin = h2;
        h1.target = h1.p1.is_virtual() ? h1.p1 : h1.p2;
        h1.next = h2;
        h1.prev = h2;
        h1.face = new Face();
        h1.face.halfedge = h1;
        h1.target.halfedge = h1;

        h2.twin = h1;
        h2.target = !h1.p1.is_virtual() ? h1.p1 : h1.p2;
        h2.next = h1;
        h2.prev = h1;
        h2.face = new Face();
        h2.face.halfedge = h2;
        h2.target.halfedge = h2;

        vvert.halfedge = h1;

        v.halfedges.push(h1);
        v.halfedges.push(h2);

        v.add_vertex_at(h2, bis23.p2, bis23);
        v.add_vertex_at(h2.next.twin, bis13.p2, bis13);

        v.site_map.set(p1.id, h1.face);
        v.site_map.set(p2.id, h2.face);
        v.site_map.set(p3.id, h2.next.twin.face);

        v.faces.push(h1.face);
        v.faces.push(h2.face);

        return v;
    }

    static trivial_voronoi_2(p1, p2) {
        var v = new Voronoi();
        var l = Line.bisector_of(p1, p2);

        var h1 = HalfEdge.from_line(l);
        var h2 = HalfEdge.from_line(l);
        v.sites.set(p1.id, p1);
        v.sites.set(p2.id, p2);

        v.halfedges.push(h1);
        v.halfedges.push(h2);

        h1.target = h1.p1;
        h1.twin = h2;
        h1.next = h2;
        h1.prev = h2;
        h1.face = new Face();
        h1.face.halfedge = h1;

        h2.target = h1.p2;
        h2.twin = h1;
        h2.next = h1;
        h2.prev = h1;
        h2.face = new Face();
        h2.face.halfedge = h2;

        v.faces.push(h1.face);
        v.faces.push(h2.face);

        v.site_map.set(p1.id, h1.face);
        v.site_map.set(p2.id, h2.face);
        return v;
    }

    static trivial_voronoi_1(p) {
        var v = new Voronoi();
        v.sites.set(p.id, p);
        v.site_map.set(p.id, v.outerface);
        return v;
    }

    static rec_construct_vor(pts) {
        var n = pts.length;
        if (n <= 3) {
            // Trivial case
            if (n == 1) return Voronoi.trivial_voronoi_1(pts[0]);
            else if (n == 2) return Voronoi.trivial_voronoi_2(pts[0], pts[1]);
            return Voronoi.trivial_voronoi_3(pts[0], pts[1], pts[2]);
        } else {
            // Split & Merge
            // TODO: This sorting should be done only once
            // pts.sort(function(a, b) { return a.x - b.x });

            var n1 = Math.floor(n / 2);
            var S1 = pts.slice(0, n1);
            var S2 = pts.slice(n1);
            var vorDiag1 = Voronoi.rec_construct_vor(S1);
            var vorDiag2 = Voronoi.rec_construct_vor(S2);
            return Voronoi.merge_voronoi(vorDiag1, vorDiag2);
        }
    }

    static _merge_voronoi(vorL, vorR) {
        // Construct convex hull of two subgraphs
        var chsL = DCEL.convex_hull_pts(Array.from(vorL.sites.values()));
        var chsR = DCEL.convex_hull_pts(Array.from(vorR.sites.values()));
        var sigma_points = [];

        // Find lower common support
        var bdgs = Voronoi.bridges(chsL, chsR);
        var ubdg = bdgs[0]; // Upper bridge
        var lbdg = bdgs[1]; // Lower bridge

        drawGraph(Array.from(vorL.sites.values()), [], vorL.edges());
        drawGraph(Array.from(vorR.sites.values()), [], vorR.edges());
        drawOutStandingEdges(bdgs);

        // First find the bisector of the lower bridge and work our way up
        var sigma_upper_inf = Line.bisector_of(ubdg[0], ubdg[1], -Infinity, Infinity);
        var cur_bdg = ubdg;
        var cur_bisector = sigma_upper_inf;

        var sitel = cur_bdg[0]; var sitel_ind_in_ch = chsL.indexOf(sitel);
        var siter = cur_bdg[1]; var siter_ind_in_ch = chsR.indexOf(siter);
        var edgel = null;
        var edger = null;
        var sigma_chain = new DCEL();
        var sigma_stack = [];
        var il = sitel_ind_in_ch;
        var ir = siter_ind_in_ch;

        // We're going to construct the sigma chain by adding vertices
        // to halfedges. The starting halfedge we'll use is the circle edge
        // of virtual vertex.
        var current_sigma_chain_halfedge = sigma_chain.virtualv.halfedge;

        //----------------
        var i = 0;

        //----------------

        // while current working bridge is not upper bridge
        while (!arr_eq(cur_bdg, lbdg)) {
            var int_l = null; // Intersection of current sigma chain part with some left voronoi edge
            var int_r = null; // Intersection of current sigma chain part with some right voronoi edge
            cur_bdg = [sitel, siter];
            console.log('sitel', sitel.x, sitel.y);
            console.log('siter', siter.x, siter.y);

            i++;
            if (i > 7) break;

            // Traverse all voronoi half edges defined by sitel and eiter
            var face_of_sitel = vorL.site_map.get(sitel.id);
            var traverse_posl = face_of_sitel.halfedge;
            edgel = edgel == null ? traverse_posl : edgel;  // If edgel is not null, this means we didn't switch site at this side last round
            do {
                // int_l = cur_bisector.intersection(edgel);
                int_l = cur_bisector.intersection(traverse_posl);
                if (int_l != null) break;
                edgel = edgel.next;
            } while (edgel.eq(traverse_posl));

            var face_of_siter = vorR.site_map.get(siter.id);
            var traverse_posr = face_of_siter.halfedge;
            edger = edger == null ? traverse_posr : edger;  // If edger is not null, this means we didn't switch site at this side last round
            do {
                // int_r = cur_bisector.intersection(edger);
                int_r = cur_bisector.intersection(traverse_posr);
                if (int_r != null) break;
                edger = edger.next;
            } while (edger.eq(traverse_posr));

            if (int_l.y < int_r.y) {
                // sigma_chain.add_vertex_at(int_r, current_sigma_chain_halfedge)
                // The false below means right side
                sigma_stack.push([int_r, false, edger]);
                // current_sigma_chain_halfedge = int_r.halfedge.twin;
                siter = chsR[(chsR.length + ir - 1) % chsR.length];
                edger = null;
                ir = (chsR.length + ir - 1) % chsR.length;
            } else {
                // sigma_chain.add_vertex_at(int_l, current_sigma_chain_halfedge)
                sigma_stack.push([int_l, true, edgel]);
                // current_sigma_chain_halfedge = int_l.halfedge.twin;
                sitel = chsL[(il + 1) % chsL.length];
                edgel = null;
                il = (il + 1) % chsL.length;
            }
            cur_bisector = Line.bisector_of(sitel, siter);
        }
        // Add final infinite part

        // Merge vorL, sigma_chain and vorR
        return Voronoi.clip_edges_and_merge(vorL, vorR, sigma_stack);
    }

    static starting_edgea(a, b) {
        var t;
        var d, prevd;
        var edge = new HalfEdge();
        var end = new HalfEdge();

        end = vertices[a].tail_edge;
        prevd = clw(p[b], p[a], p[end.target]);
        edge = end.prev;
        while(!edge.eq(end)) {
            if(!edge.is_void()) {
                t = edge.target;
                d = clw(p[b], p[a], p[t]);
                if(!prevd && d) return edge;
                prevd = d;
            }
            edge = edge.prev;
        }

        return edge;
    }

    static starting_edgeb(b, a) {
        var t;
        var d, prevd;
        var edge = new HalfEdge();
        var end = new HalfEdge();

        end = vertices[b].head_edge;
        prevd = ccw(p[a], p[b], p[end.target]);
        edge = end.next;
        while(!edge.eq(end)) {
            if(!edge.is_void()) {
                t = edge.target;
                d = ccw(p[a], p[b], p[t]);
                if(!prevd && d) return edge;
                prevd = d;
            }
            edge = edge.next;
        }

        return edge;
    }

    static voronoi_merge(m) {
        // line_t zip;
        var zip = new Line();
        var  a, b, ac, bc;
        var x = new Vertex(Infinity, Infinity);
        var xa = new Vertex(Infinity, Infinity);
        var xb = new Vertex(Infinity, Infinity);
        var edgea = new HalfEdge();
        var edgeb = new HalfEdge();
        var starta = new HalfEdge();
        var startb = new HalfEdge();

        a = m;  b = m + 1;
        var tmp = merge_hull(a, b);
        a = tmp[0]; b = tmp[1];

        zip = Line.perpdiv(p[a],p[b]);

        edgea = starta = starting_edgea(a, b);
        edgeb = startb = starting_edgeb(b, a);
        for(;;) {

            ac = 0;
            do {
                if(!edgea.is_void()) {

                    if(ccw(p[b], p[a], p[edgea.target], true)) break;

                    if(Line.intersect(xa, zip, edgea.divline)) {
                        ac = 1;
                        break;
                    }
                }

                edgea = edgea.prev;
            } while(!edgea.eq(starta));

            bc = 0;
            do {
                if(!edgeb.is_void()) {
                    if(clw(p[a], p[b], p[edgeb.target], true)) break;
                    if(intersect(xb, zip, edgeb.divline)) {
                        bc = 1;
                        break;
                    }
                }
                edgeb = edgeb.next;
            } while(!edgeb.eq(startb));

            if(!(ac || bc)) break;

            if(bc == 0 || (ac && xa.y <= xb.y)) {

                if(snap_point(xa, edgea.divline)) {

                    var czip = zip;
                    fix_p2(czip, xa);
                    if(dist2(czip.p1, czip.p2) < CLOSE * CLOSE) {
                        czip = moveline(zip, xa);
                    }
                    rclip(edgea.divline, czip, xa);
                    rdelete(a, czip, xa);
                }
                else {

                    if(dist2(zip.p1, xa) < CLOSE * CLOSE) xa = zip.p1;
                    rclip(edgea.divline, zip, xa);
                    rdelete(a, zip, xa);
                }
                fix_p2(zip,xa);
                push_zipline(a, b, edgea, edgeb, zip);
                x = xa;
                a = edgea.target;
                edgea = starta = starting_edgea(a, b);
            }
            else {
                if(snap_point(xb, edgeb.divline)) {
                    var czip = zip;
                    fix_p2(czip, xb);
                    if(dist2(czip.p1, czip.p2) < CLOSE * CLOSE) {
                        czip = moveline(zip, xb);
                    }
                    lclip(edgeb.divline, czip, xb);
                    ldelete(b, czip, xb);
                }
                else {
                    if(dist2(zip.p1, xb) < CLOSE * CLOSE) xb = zip.p1;
                    lclip(edgeb.divline, zip, xb);
                    ldelete(b, zip, xb);
                }
                fix_p2(zip, xb);
                push_zipline(a, b, edgea, edgeb, zip);
                x = xb;
                b = edgeb.target;
                edgeb = startb = starting_edgeb(b, a);
            }

            zip = Line.perpthrough(x, p[a], p[b]);
        }

        push_zipline(a, b, edgea, edgeb, zip);

        pop_ziplines();
    }

    // Clip edges of vorL that's on the right of sigma
    // and vorR that's on the left side of sigma
    static clip_edges_and_merge(vorL, vorR, sigma) {
        vorL.virtualv = vorR.virtualv;
        vorL.sites = new Map([vorL.sites, vorR.sites]);
        vorL.vertices = vorL.vertices.concat(vorR.vertices);
        vorL.halfedges = vorL.halfedges.concat(vorR.halfedges);
        var prev_p = vorL.virtualv;
        sigma.push(vorL.virtualv);
        for (var i = 0; i < sigma.length; i++) {
            var tpl = sigma[i]; var intp = tpl[0];
            var is_left = tpl[1] == true; var edge = tpl[2];
            var cur_vor = is_left ? vorL : vorR;
            // var other_vor = is_left ? vorR : vorL;

            cur_vor.split_edge(edge, intp);
            cur_vor.split_face(edge, intp);
            prev_p = intp;
        }

        return vorL;
    }

    // Return a collection of edges, in ccw order,
    // of site p.
    // vedges_of_site(p) {
    // For each site, store a list of vertices that are closest to it
    // So at the merge stage, you have to figure out some ways to merge the hull
    // sites as well
    // return site_map[]
    // }

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

// Helper function. Convert a list of points [[0,1],[1,2],...[3,2]] to [{x:0,y:1},...{}]
function tuple2vert(plist) {
    return plist.map(function(d) { return Vertex(d[0], d[1]); })
}

function main() {
    var n = 4;
    var random_pts = d3.range(n)
        .map(function(d) { return new Vertex(Math.random() * width, Math.random() * height); });

    var pts = DCEL.sort_points_by_x(random_pts);
    var vor = new Voronoi().build(pts);

    drawVoronoi(vor);
}

function drawVoronoi(vor) {
    var pts = Array.from(vor.sites.values());
    drawGraph(pts, vor.vertices, vor.edges());
}

function drawGraph(sites, vertices, edges) {
    // ctx.clearRect(0, 0, width, height);

    // Draw site
    ctx.beginPath();
    for (var i = 0; i < sites.length; i++) drawSite(sites[i]);
    ctx.fillStyle = "#000";
    ctx.fill();
    ctx.strokeStyle = "#fff";
    ctx.stroke();

    ctx.beginPath();
    for (var i = 0; i < edges.length; i++)
        drawEdge(edges[i]);
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

function drawEdge(edge) {
    drawEdgeBetween(edge.p1, edge.p2);
}

function drawEdgeBetween(p1, p2) {
    ctx.moveTo(p1.x, p1.y);
    ctx.lineTo(p2.x, p2.y);
}

function drawOutStandingEdges(edges) {
    ctx.beginPath();
    for (var i = 0; i < edges.length; i++)
        drawEdgeBetween(edges[i][0], edges[i][1]);
    ctx.strokeStyle = "rgba(0,0,0,0.2)";
    ctx.stroke();
}

function clw(p1, p2, p3, report_colinear=false) {
    if (report_colinear) return DCEL.orientation(p1, p2, p3) >= 0;
    return DCEL.orientation(p1, p2, p3) > 0;
}

function ccw(p1, p2, p3, report_colinear=false) {
    if (report_colinear) return DCEL.orientation(p1, p2, p3) <= 0;
    return DCEL.orientation(p1, p2, p3) < 0;
}

function in_triangle(p1, p2, p3, v) {
    return DCEL.orientation(p1, p2, v) == DCEL.orientation(p2, p3, v) &&
        DCEL.orientation(p2, p3, v) == DCEL.orientation(p3, p1, v);
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
    return new Vertex(
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

function dist2(p1, p2) {
    return (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y);
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

function snap_point(q, l) {
    var dx1 = q.x - l.p1.x;
    var dy1 = q.y - l.p1.y;
    var dx2 = q.x - l.p2.x;
    var dy2 = q.y - l.p2.y;

    if(!(l.reach & REACH_BACKWARD != 0) &&
       dx1 * dx1 + dy1 * dy1 <= CLOSE * CLOSE) {
        q = l.p1;
        return true;
    }
    else if(!(l.reach & REACH_FORWARD) &&
            dx2 * dx2 + dy2 * dy2 <= CLOSE * CLOSE) {
        q = l.p2;
        return true;
    }
    return false;
}

function fix_p1(l, q) {
    l.reach &= ~REACH_BACKWARD;
    if(l.reach & REACH_FORWARD) {
        var dx = q.x - l.p1.x;
        var dy = q.y - l.p1.y;
        l.p2.x += dx;
        l.p2.y += dy;
    }
    l.p1 = q;
}

function fix_p2(l, q) {
    l.reach &= ~REACH_FORWARD;
    if(l.reach & REACH_BACKWARD) {
        var dx = q.x - l.p2.x;
        var dy = q.y - l.p2.y;
        l.p1.x += dx;
        l.p1.y += dy;
    }
    l.p2 = q;
}

function moveline(l, c) {
    var dx = c.x - l.p1.x;
    var dy = c.y - l.p1.y;
    l.p1 = c;
    l.p2.x += dx;
    l.p2.y += dy;
    return l;
}

function rclip(line, zip, x) {
    if(cross(line, zip) < 0) {
        fix_p1(line, x);
    }
    else {
        fix_p2(line, x);
    }
}

function lclip(line, zip, x) {
    if(cross(line, zip) > 0) {
        fix_p1(line, x);
    }
    else {
        fix_p2(line, x);
    }
}

function rdelete(v, zip, x) {
    var edge = new HalfEdge();
    var head = new HalfEdge();
    var p1 = zip.p1;
    var p2 = zip.p2;
    var l0 = new Line();
    l0.p1 = x;
    l0.p2 = x;
    l0.reach = REACH_NONE;
    edge = head = vertices[v].head_edge;
    do {
        if(!edge.is_void()) {

            if(clw(p1, p2, edge.divline.p1, true)
                && clw(p1, p2, edge.divline.p2, true)) {
                edge.divline = l0;
            }
        }
        edge = edge.next;
    } while(!edge.eq(head));
}

function ldelete(v, zip, x) {
    var edge = new HalfEdge();
    var head = new HalfEdge();
    var p1 = zip.p1;
    var p2 = zip.p2;
    var l0 = new Line();
    l0.p1 = x;
    l0.p2 = x;
    l0.reach = REACH_NONE;
    edge = head = vertices[v].head_edge;
    do {
        if(!is_void_edge(edge)) {

        if(ccw(p1, p2, edge.divline.p1, true)
           && ccw(p1, p2, edge.divline.p2, true)) {
            edge.divline = l0;
        }
        }
        edge = edge.next;
    } while(!edge.eq(head));
}

function cross(l1, l2) {
    return (l1.p2.x - l1.p1.x) * (l2.p2.y - l2.p1.y)
         - (l1.p2.y - l1.p1.y) * (l2.p2.x - l2.p1.x);
}

function push_zipline(a, b, edgea, edgeb, zip) {
    var l = Line.copy(zip);
    var e1 = HalfEdge.copy(a, b, l);
    var e2 = HalfEdge.copy(a, b, l);
    e1.pair_with(e2);
    left.push(a);
    right.push(b);
    leftedge.push(edgea);
    rightedge.push(edgeb);
    leftzip.push(e1);
    rightzip.push(e1);
    tos++;
}

function pop_ziplines() {
    var i, last;

    last = tos - 1;
    if(last != 0) while(left[last-1] == left[last]) last--;

    for(i = 0; i < last; i++) {
        add_edge_after(leftedge[i], leftzip[i]);
    }

    for(i = tos-1; i >= last; i--) add_edge_before(vertices[left[i]].head_edge, leftzip[i]);

    last = tos-1;
    if(last != 0) while(right[last-1] == right[last]) last--;
    for(i = 0; i < last; i++) {
        add_edge_before(rightedge[i], rightzip[i]);
    }
    for(i = tos-1; i >= last; i--) add_edge_after(vertices[right[i]].tail_edge, rightzip[i]);

    last = tos-1;
    vertices[left[0]].tail_edge = leftzip[0];
    vertices[right[0]].head_edge = rightzip[0];
    vertices[left[last]].head_edge = leftzip[last];
    vertices[right[last]].tail_edge = rightzip[last];

    clean_edgelist(left[0]);
    i = 1;
    while(i < tos) {
        if(left[i] != left[i-1]) clean_edgelist(left[i]);
        i++;
    }

    clean_edgelist(right[0]);
    i = 1;
    while(i < tos) {
        if(right[i] != right[i-1]) clean_edgelist(right[i]);
        i++;
    }

    tos = 0;
}

function delete_line(l) {
    n_lines--;
}

function remove_edge(edge) {
    var v = edge.source;
    var prev = edge.prev;
    var next = edge.next;
    prev.next = next;
    next.prev = prev;
    if(vertices[v].head_edge == edge) vertices[v].head_edge = next;
    if(vertices[v].tail_edge == edge) vertices[v].tail_edge = prev;
}

function delete_edge_pair(edge) {
    remove_edge(edge.partner);
    remove_edge(edge);
    // delete_line(edge.divline);
}

function assign_edge(source, edge) {
    edge.next = edge.prev = edge;
    vertices[source].head_edge = vertices[source].tail_edge = edge;
}

function assign_edge2(source, edge1, edge2) {
    edge1.next = edge1.prev = edge2;
    edge2.next = edge2.prev = edge1;
    vertices[source].head_edge = edge1;
    vertices[source].tail_edge = edge2;
}

function add_edge(edge_before, edge_after, edge) {
    edge.next = edge_after;
    edge.prev = edge_before;
    edge_before.next = edge;
    edge_after.prev = edge;
}

function add_edge_after(edge_before, edge) {
    add_edge(edge_before, edge_before.next, edge);
}

function add_edge_before(edge_after, edge) {
    add_edge(edge_after.prev, edge_after, edge);
}

function clean_edgelist(v) {
    var edge = new HalfEdge();
    var next = new HalfEdge();
    var head = vertices[v].head_edge;
    if(head.is_void()) {
        do {
            head = head.next;
        } while(head.is_void());
        vertices[v].head_edge = head;
    }
    edge = head.next;
    while(!edge.eq(head)) {
        next = edge.next;
        if(edge.is_void()) delete_edge_pair(edge);
        edge = next;
    }
}

function linkhull3(v1, v2, v3) {
    vertices[v1].prev = v3;
    vertices[v1].next = v2;
    vertices[v2].prev = v1;
    vertices[v2].next = v3;
    vertices[v3].prev = v2;
    vertices[v3].next = v1;
}

function linkhull2(v1, v2) {
    vertices[v1].prev = v2;
    vertices[v1].next = v2;
    vertices[v2].prev = v1;
    vertices[v2].next = v1;
}

function construct3(v1, v2, v3, lA, lB, lC) {
    var sA = Line.copy(lA);
    var sB = Line.copy(lB);
    var sC = Line.copy(lC);
    var e12 = HalfEdge.copy(v1, v2, sA);
    var e21 = HalfEdge.copy(v2, v1, sA);
    var e23 = HalfEdge.copy(v2, v3, sB);
    var e32 = HalfEdge.copy(v3, v2, sB);
    var e31 = HalfEdge.copy(v3, v1, sC);
    var e13 = HalfEdge.copy(v1, v3, sC);
    e12.pair_with(e21);
    e23.pair_with(e32);
    e31.pair_with(e13);
    assign_edge2(v1, e12, e13);
    assign_edge2(v2, e23, e21);
    assign_edge2(v3, e31, e32);
    linkhull3(v1, v2, v3);
}

function construct2(v1, v2, v3, lA, lB) {
    var sA = Line.copy(lA);
    var sB = Line.copy(lB);
    var e12 = HalfEdge.copy(v1, v2, sA);
    var e21 = HalfEdge.copy(v2, v1, sA);
    var e23 = HalfEdge.copy(v2, v3, sB);
    var e32 = HalfEdge.copy(v3, v2, sB);
    e12.pair_with(e21);
    e23.pair_with(e32);
    assign_edge(v1, e12);
    assign_edge2(v2, e21, e23);
    assign_edge(v3, e32);
    linkhull3(v1, v3, v2);
}

function construct1(v1, v2, l) {
    var s = Line.copy(l);
    var e12 = HalfEdge.copy(v1, v2, s);
    var e21 = HalfEdge.copy(v2, v1, s);
    e12.pair_with(e21);
    assign_edge(v1, e12);
    assign_edge(v2, e21);
    linkhull2(v1, v2);
}

function store_lines(r) {
    var i, j;
    var l;

    // l = r.l = malloc(n_lines * sizeof(line_t));
    j = 0;
    for (i = 0; i < vertices.length; i++) {
        var edge, head;
        edge = head = vertices[i].head_edge;
        if (edge) do {
            if (edge.target > i) l[j++] = edge.divline;
            edge = edge.next;
        } while (!edge.eq(head));
    }
    // r.n_lines = j;
}

main();

