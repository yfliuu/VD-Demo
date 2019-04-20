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

// Represents a segment, ray or line.
// HalfEdge of DCEL extends this class.
// Can be built using from_pts and from_kb, where k and b are the parameters
// of y=kx+b.
class Line {
    k = null;
    b = null;
    p1 = null;
    p2 = null;
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
    target = null;  // Source vertex
    prev = null;    // Previous halfedge
    next = null;    // Next halfedge
    face = null;    // halfedgeident face
    id = null;

    constructor() {
        super();
        this.id = HALF_EDGE_AUTOID++;
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
}

class Vertex {
    x = null;
    y = null;
    halfedge = null;   // This is the outgoing edge
    id = null;

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
}

class Voronoi extends DCEL {
    // Every site uniquely maps to a face of DCEL.
    // Sites are also class Vertexes but contains id.
    // ids will be set in constructor.
    sites = null;
    site_map = null;

    // S is a set of Vertices
    constructor() {
        super();
        this.sites = new Map([]);
        this.site_map = new Map([]);
    }

    build(S) {
        return Voronoi.rec_construct_vor(S);
    }

    edges() {
        return this.halfedges;
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
            pts.sort(function(a, b) { return a.x - b.x });

            var n1 = Math.floor(n / 2);
            var S1 = pts.slice(0, n1);
            var S2 = pts.slice(n1);
            var vorDiag1 = Voronoi.rec_construct_vor(S1);
            var vorDiag2 = Voronoi.rec_construct_vor(S2);
            return Voronoi.merge_voronoi(vorDiag1, vorDiag2);
        }
    }

    static merge_voronoi(vorL, vorR) {
        console.log(vorL);
        console.log(vorR);
        // Construct convex hull of two subgraphs
        var chsL = DCEL.convex_hull_pts(Array.from(vorL.sites.values()));
        var chsR = DCEL.convex_hull_pts(Array.from(vorR.sites.values()));
        var sigma_points = [];

        // Find lower common support
        var bdgs = Voronoi.bridges(chsL, chsR);
        var ubdg = bdgs[0]; // Upper bridge
        var lbdg = bdgs[1]; // Lower bridge

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

        // We're going to construct the sigma chain by adding vertices
        // to halfedges. The starting halfedge we'll use is the circle edge
        // of virtual vertex.
        var current_sigma_chain_halfedge = sigma_chain.virtualv.halfedge;

        // while current working bridge is not upper bridge
        while (!arr_eq(cur_bdg, lbdg)) {
            var int_l = null; // Intersection of current sigma chain part with some left voronoi edge
            var int_r = null; // Intersection of current sigma chain part with some right voronoi edge
            cur_bdg = [sitel, siter];

            console.log(cur_bdg);

            // Traverse all voronoi half edges defined by sitel and siter
            var face_of_sitel = vorL.site_map.get(sitel.id);
            var traverse_posl = face_of_sitel.halfedge;
            edgel = edgel == null ? traverse_posl : edgel;  // If edgel is not null, this means we didn't switch site at this side last round
            do {
                int_l = cur_bisector.intersection(edgel);
                if (int_l != null) break;
                edgel = edgel.next;
            } while (edgel.eq(traverse_posl));

            var face_of_siter = vorR.site_map.get(siter.id);
            var traverse_posr = face_of_siter.halfedge;
            edger = edger == null ? traverse_posr : edger;  // If edger is not null, this means we didn't switch site at this side last round
            do {
                int_r = cur_bisector.intersection(edger);
                if (int_r != null) break;
                edger = edger.next;
            } while (edger.eq(traverse_posr));

            if (int_l.y < int_r.y) {
                // sigma_chain.add_vertex_at(int_r, current_sigma_chain_halfedge)
                // The false below means right side
                sigma_stack.push([int_r, false, edger]);
                current_sigma_chain_halfedge = int_r.halfedge.twin;
                siter = chsR[(chsR.length + siter_ind_in_ch - 1) % chsR.length];
                edgel = null;
            } else {
                // sigma_chain.add_vertex_at(int_l, current_sigma_chain_halfedge)
                sigma_stack.push([int_l, true, edgel]);
                current_sigma_chain_halfedge = int_l.halfedge.twin;
                sitel = chsL[(sitel_ind_in_ch + 1) % chsL.length];
                edger = null;
            }
        }
        // Add final infinite part

        // Merge vorL, sigma_chain and vorR
        console.log(sigma_stack);
        return Voronoi.clip_edges_and_merge(vorL, vorR, sigma_stack);
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
    var canvas = d3.select("canvas").node();
    ctx = canvas.getContext("2d"),
        width = canvas.width,
        height = canvas.height;

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
    var canvas = d3.select("canvas").node();
    ctx = canvas.getContext("2d"),
        width = canvas.width,
        height = canvas.height;

    ctx.clearRect(0, 0, width, height);

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

function distance_square(p1, p2) {
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

main();

