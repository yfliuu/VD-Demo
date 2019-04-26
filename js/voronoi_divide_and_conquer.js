// CPSC-620 Course Project
// Author: Yifan Liu
// Voronoi Diagram, divide and conquer
//
// Part of the code was borrowed from http://www.cosc.canterbury.ac.nz/tad.takaoka/alg/geometry/voronoi.c

var TEST_NUM_PTS = 100;
var EXTEND = 100000;
var DISTANT = 100000;
var CLOSE = 0.000001;

var REACH_NONE = 0;
var REACH_BACKWARD = 1;
var REACH_FORWARD = 2;
var REACH_BOTH = 3;

var MAX_POINTS = 10000;

var canvas = d3.select("canvas").node();
ctx = canvas.getContext("2d"),
    width = canvas.width,
    height = canvas.height;
var WIN_WIDTH  = width;
var WIN_HEIGHT = height;
var CLIP_WIDTH = (5 * WIN_WIDTH);
var CLIP_HEIGHT = (5 * WIN_HEIGHT);

class Point {
    constructor() {
        this.x = null;
        this.y = null;
    }

    static val_close(v1, v2) {
        return Math.abs(v1 - v2) <= CLOSE;
    }

    static from_coord(x, y) {
        var pt = new Point();
        pt.x = x;
        pt.y = y;
        return pt;
    }

    copy() {
        var newp = new Point();
        newp.x = this.x;
        newp.y = this.y;
        return newp;
    }

    static point_cmp(p1, p2) {
        var d = p1.x - p2.x;
        if (d == 0) d = p1.y - p2.y;
        if (d == 0) return 0;
        return d > 0 ? 1 : -1;
    }

    static comp_points(p1, p2) {
        var d = p1.x - p2.x;
        if (d == 0) return p1.y - p2.y;
        return d;
    }

    inside() {
        return this.x > 0 && this.x < CLIP_WIDTH && this.y > 0 && this.y < CLIP_HEIGHT ? 1 : 0;
    }

    dist2(p2) {
        var dx = p2.x - this.x;
        var dy = p2.y - this.y;
        return dx * dx + dy * dy;
    }
}

class Line {
    constructor() {
        this.p1 = new Point();
        this.p2 = new Point();
        this.reach = null;
    }

    static from_pts(p1x, p1y, p2x, p2y) {
        var l = new Line();
        l.p1 = new Point();
        l.p2 = new Point();
        l.p1.x = p1x;
        l.p1.y = p1y;
        l.p2.x = p2x;
        l.p2.y = p2y;
        return l;
    }

    copy() {
        var newl = new Line();
        newl.p1 = this.p1.copy();
        newl.p2 = this.p2.copy();
        newl.reach = this.reach;
        return newl;
    }


    other(p1) {
        if (Point.val_close(p1.x, this.p1.x) && Point.val_close(p1.y == this.p1.y)) {
            return this.p2;
        } else if (Point.val_close(p1.x, this.p2.x) && Point.val_close(p1.y == this.p2.y)) {
            return this.p1;
        } else {
            return null;
        }
    }

    outside_point() {
        if (this.p1.inside()) return this.p2;
        return this.p1;
    }

    fix_p1(q) {
        this.reach &= ~REACH_BACKWARD;
        if (this.reach & REACH_FORWARD) {
            var dx = q.x - this.p1.x;
            var dy = q.y - this.p1.y;
            this.p2.x += dx;
            this.p2.y += dy;
        }
        this.p1 = q;
    }

    fix_p2(q) {
        this.reach &= ~REACH_FORWARD;
        if (this.reach & REACH_BACKWARD) {
            var dx = q.x - this.p2.x;
            var dy = q.y - this.p2.y;
            this.p1.x += dx;
            this.p1.y += dy;
        }
        this.p2 = q;
    }

    rclip(zip, x) {
        if (this.cross(zip) < 0) this.fix_p1(x);
        else this.fix_p2(x);
    }

    lclip(zip, x) {
        if (this.cross(zip) > 0) this.fix_p1(x);
        else this.fix_p2(x);
    }

    is_length0() {
        return (this.p1.x == this.p2.x && this.p1.y == this.p2.y);
    }

    static new_line(l) {
        var r = new Line();
        r = Object.assign(r, l);
        n_lines++;
        return r;
    }

    move(c) {
        var l = this.copy();
        l.p1 = c;
        l.p2.x += c.x - l.p1.x;
        l.p2.y += c.y - l.p1.y;
        return l;
    }

    static bisector(p1, p2) {
        var midx = (p1.x + p2.x) / 2;
        var midy = (p1.y + p2.y) / 2;
        var dx = p2.x - p1.x;
        var dy = p2.y - p1.y;
        var div = Math.sqrt(dx * dx + dy * dy);

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

    static bisector_through(c, p1, p2) {
        var dx = p2.x - p1.x;
        var dy = p2.y - p1.y;
        var div = Math.sqrt(dx * dx + dy * dy);
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

    extend() {
        var dx = this.p2.x - this.p1.x;
        var dy = this.p2.y - this.p1.y;
        var midx = (this.p1.x + this.p2.x) / 2;
        var midy = (this.p1.y + this.p2.y) / 2;
        var r = Math.sqrt(dx * dx + dy * dy);
        var ux = dx / r;
        var uy = dy / r;
        var lx = ux * EXTEND;
        var ly = uy * EXTEND;
        if (this.reach & REACH_BACKWARD) {
            this.p1.x = midx - lx;
            this.p1.y = midy - ly;
        }
        if (this.reach & REACH_FORWARD) {
            this.p2.x = midx + lx;
            this.p2.y = midy + ly;
        }
        return this.copy();
    }

    clip() {
        var border = [];
        border.push(Line.from_pts(0, 0, CLIP_WIDTH, 0));
        border.push(Line.from_pts(CLIP_WIDTH, 0, CLIP_WIDTH, CLIP_HEIGHT));
        border.push(Line.from_pts(CLIP_WIDTH, CLIP_HEIGHT, 0, CLIP_HEIGHT));
        border.push(Line.from_pts(0, CLIP_HEIGHT, 0, 0));
        var x = new Point();
        var s1, s2;
        var i, result;

        s1 = (this.reach & REACH_BACKWARD) ? 0 : this.p1.inside();
        s2 = (this.reach & REACH_FORWARD)  ? 0 : this.p2.inside();
        if (s1 && s2) return true;

        if (s1 || s2) {
            for (i = 0; i < 4; i++) {
                if (intersect(x, border[i],this)) {
                    this.rclip(border[i], x.copy());
                    return true;
                }
            }
        }
        result = 0;
        for (i = 0; i < 4; i++) {
            if (intersect(x, border[i], this)) {
                this.rclip(border[i], x.copy());
                result |= 1;
            }
        }
        return result;
    }

    cross(l2) {
        return (this.p2.x - this.p1.x) * (l2.p2.y - l2.p1.y)
            - (this.p2.y - this.p1.y) * (l2.p2.x - l2.p1.x);
    }
}

class HalfEdge {
    constructor() {
        this.source = null;
        this.target = null;
        this.divline = new Line();
        this.twin = null;
        this.next = null;
        this.prev = null;
    }

    static new_edge(source, target, divline) {
        var edge = new HalfEdge();
        edge.source = source;
        edge.target = target;
        edge.divline = divline;
        return edge;
    }

    is_void() {
        return this.divline.is_length0();
    }

    twin_with(e2) {
        this.twin = e2;
        e2.twin = this;
    }
}

class Vertex {
    next = null;
    prev = null;
    head_edge = null;
    tail_edge = null;
    p = null;

    constructor() {
        this.next = null;
        this.prev = null;
        this.head_edge =  null;
        this.tail_edge = null;
        this.p = null;
    }
}

class Cell {
    id = null;
    point = null;
    color = null;
    edges = null;
    vertices = null;
    vindexes = null;

    constructor() {
        this.id = null;
        this.point = null;
        this.edges = [];
        this.color = 0;
        this.vertices = [];
        this.vindexes = [];
    }

    intersect(cell) {
        for (var i = 0; i < cell.vindexes.length; i++) {
            if (cell.vindexes.has(cell.vindexes[i])) return true;
        }
        return false;
    }
}

class Voronoi {
    constructor(points, n) {
        var rpoints;
        var rvertices;
        var i;

        this.n_points = n;
        this.n_lines = 0;
        this.l = 0;

        rpoints = this.p = new Array(n);
        for (i = 0; i < n; i++) rpoints[i] = points[i];
        rpoints.sort(Point.point_cmp);

        this.n_vertices = remove_duplicates(rpoints, n);
        rvertices = this.vertices = new Array(this.n_vertices);
        for (i = 0; i < this.n_vertices; i++) {
            rvertices[i] = new Vertex();
            rvertices[i].next = rvertices[i].prev = i;
            rvertices[i].head_edge = rvertices[i].tail_edge = 0;
        }
    }

    cells() {
        var remaining = new Set();
        for (var i = 0; i < this.n_vertices; i++) remaining.add(i);
        var cells = []

        while (remaining.size > 0) {
            var ind = remaining.values().next().value;
            var vert = this.vertices[ind];
            var edge = vert.head_edge;
            var cell = new Cell();

            do {
                if (!edge.is_void()) {
                    var target_point = new Point();
                    cell.id = cells.length;
                    cell.edges.push(edge);
                    if (intersect(target_point, edge.divline, edge.next.divline)) {
                        cell.vertices.push(target_point);
                        if (edge == edge.next.next) {
                            cell.vertices.push(edge.divline.outside_point());
                            cell.vertices.push(edge.next.divline.outside_point());
                            cell.vertices.push(target_point);
                            edge = edge.next;
                        }
                    } else {
                        var endp = edge.divline.p1.inside() ? edge.divline.p2 : edge.divline.p1;
                        var next = edge.next.divline.p1.inside() ? edge.next.divline.p2 : edge.next.divline.p1;
                        cell.vertices.push(endp);
                        cell.vertices.push(next);
                    }
                    remaining.delete(edge.source);
                }
                edge = edge.next;
            } while (edge != vert.head_edge);
            cells.push(cell);
        }
        return cells;
    }

    static voronoi(points, n) {
        var r;
        r = new Voronoi(points, n);
        p = r.p;
        vertices = r.vertices;

        var k = 2 * r.n_vertices + 1;
        leftzip = new Array(k);
        rightzip = new Array(k);
        leftedge = new Array(k);
        rightedge = new Array(k);
        left = new Array(k);
        right = new Array(k);
        tos = 0;

        n_lines = 0;
        if (r.n_vertices >= 2) Voronoi.rec_voronoi(0, r.n_vertices - 1);
        r.store_lines();

        return r;
    }

    static voronoi_merge(m) {
        var zip = new Line();
        var a, b;
        var ac, bc;
        var x = new Point();
        var xa = new Point();
        var xb = new Point();
        var edgea = new HalfEdge();
        var edgeb = new HalfEdge();
        var starta = new HalfEdge();
        var startb = new HalfEdge();

        a = m;  b = m + 1;
        var tpl = merge_hull(a, b);
        a = tpl[0]; b = tpl[1];

        zip = Line.bisector(p[a], p[b]);

        edgea = starta = starting_edgea(a, b);
        edgeb = startb = starting_edgeb(b, a);
        while (true) {
            ac = 0;
            do {
                if (!edgea.is_void()) {
                    if (ccweq(p[b], p[a], p[edgea.target])) break;
                    if (intersect(xa, zip.copy(), edgea.divline.copy())) {
                        ac = 1;
                        break;
                    }
                }
                edgea = edgea.prev;
            } while (edgea != starta);

            bc = 0;
            do {
                if (!edgeb.is_void()) {
                    if (clweq(p[a], p[b], p[edgeb.target])) break;
                    if (intersect(xb, zip.copy(), edgeb.divline.copy())) {
                        bc = 1;
                        break;
                    }
                }
                edgeb = edgeb.next;
            } while (edgeb != startb);

            if (!(ac || bc)) break;

            if (bc == 0 || (ac && xa.y <= xb.y)) {
                if (snap_point(xa, edgea.divline.copy())) {
                    var czip = zip.copy();
                    czip.fix_p2(xa.copy());
                    if (czip.p1.dist2(czip.p2) < CLOSE * CLOSE) {
                        czip = zip.move(xa);
                    }
                    edgea.divline.rclip(czip.copy(), xa.copy());
                    rdelete(a, czip.copy(), xa.copy());
                }
                else {
                    if (zip.p1.dist2(xa) < CLOSE * CLOSE) xa = zip.p1.copy();
                    edgea.divline.rclip(zip.copy(), xa.copy());
                    rdelete(a, zip.copy(), xa.copy());
                }
                zip.fix_p2(xa.copy());
                push_zipline(a, b, edgea, edgeb, zip.copy());
                x = xa.copy();
                a = edgea.target;
                edgea = starta = starting_edgea(a, b);
            }
            else {
                if (snap_point(xb, (edgeb.divline))) {
                    var czip = zip.copy();
                    czip.fix_p2(xb.copy());
                    if (czip.p1.dist2(czip.p2) < CLOSE * CLOSE) {
                        czip = zip.move(xb);
                    }
                    lclip(edgeb.divline, czip.copy(), xb.copy());
                    ldelete(b, czip.copy(), xb.copy());
                }
                else {
                    if (zip.p1.dist2(xb) < CLOSE * CLOSE) xb = zip.p1.copy();
                    edgeb.divline.lclip(zip.copy(), xb.copy());
                    ldelete(b, zip.copy(), xb.copy());
                }
                zip.fix_p2(xb.copy());
                push_zipline(a, b, edgea, edgeb, zip.copy());
                x = xb.copy();
                b = edgeb.target;
                edgeb = startb = starting_edgeb(b, a);
            }

            zip = Line.bisector_through(x, p[a], p[b]);
        }

        push_zipline(a, b, edgea, edgeb, zip.copy());
        pop_ziplines();
    }

    static rec_voronoi(l, r) {
        var n_points = 1 + r - l;
        var m;

        if (n_points <= 3) {
            Voronoi.voronoi_base(l, r, n_points);
            return;
        }

        m = Math.floor((l + r) / 2);
        Voronoi.rec_voronoi(l, m);
        Voronoi.rec_voronoi(m + 1, r);
        Voronoi.voronoi_merge(m);
    }

    static voronoi_base(l, r, n_points) {
        var q = p[l];
        var lineA, lineB, lineC;

        if (n_points == 3) {
            lineA = Line.bisector(p[l],p[l + 1]);
            lineB = Line.bisector(p[l + 1],p[l + 2]);
            var c = new Point();
            if (intersect(c, lineA, lineB)) {
                if (c.y <= p[l + 1].y) {
                    lineA = Line.bisector_through(c, p[l], p[l + 1]);
                    lineB = Line.bisector_through(c, p[l + 1], p[l + 2]);
                    lineC = Line.bisector_through(c, p[l + 2], p[l]);
                    Voronoi.construct3(l, l + 1, l + 2, lineA, lineB, lineC);
                }
                else {
                    lineA = Line.bisector_through(c, p[l + 1], p[l]);
                    lineB = Line.bisector_through(c, p[l + 2], p[l + 1]);
                    lineC = Line.bisector_through(c, p[l], p[l + 2]);
                    Voronoi.construct3(l, l + 2, l + 1, lineC, lineB, lineA);
                }
            }
            else {
                Voronoi.construct2(l, l + 1, l + 2, lineA, lineB);
            }
        }
        else {
            lineA = Line.bisector(p[l], p[l + 1]);
            Voronoi.construct1(l, l + 1, lineA);
        }
    }

    store_lines() {
        var i, j;
        var l;

        l = this.l = new Array(n_lines);
        j = 0;
        for (i = 0; i < this.n_vertices; i++) {
            var edge, head;
            edge = head = vertices[i].head_edge;
            if (edge) do {
                if (edge.target > i) l[j++] = (edge.divline);
                edge = edge.next;
            } while (edge != head);
        }
        this.n_lines = j;
    }

    static construct3(v1, v2, v3, lA, lB, lC) {
        var sA = Line.new_line(lA);
        var sB = Line.new_line(lB);
        var sC = Line.new_line(lC);
        var e12 = HalfEdge.new_edge(v1, v2, sA);
        var e21 = HalfEdge.new_edge(v2, v1, sA);
        var e23 = HalfEdge.new_edge(v2, v3, sB);
        var e32 = HalfEdge.new_edge(v3, v2, sB);
        var e31 = HalfEdge.new_edge(v3, v1, sC);
        var e13 = HalfEdge.new_edge(v1, v3, sC);
        e12.twin_with(e21);
        e23.twin_with(e32);
        e31.twin_with(e13);
        assign_edge2(v1, e12, e13);
        assign_edge2(v2, e23, e21);
        assign_edge2(v3, e31, e32);
        Voronoi.linkhull3(v1, v2, v3);
    }

    static construct2(v1, v2, v3, lA, lB) {
        var sA = Line.new_line(lA);
        var sB = Line.new_line(lB);
        var e12 = HalfEdge.new_edge(v1, v2, sA);
        var e21 = HalfEdge.new_edge(v2, v1, sA);
        var e23 = HalfEdge.new_edge(v2, v3, sB);
        var e32 = HalfEdge.new_edge(v3, v2, sB);
        e12.twin_with(e21);
        e23.twin_with(e32);
        assign_edge(v1, e12);
        assign_edge2(v2, e21, e23);
        assign_edge(v3, e32);
        Voronoi.linkhull3(v1, v3, v2);
    }

    static construct1(v1, v2, l) {
        var s = Line.new_line(l);
        var e12 = HalfEdge.new_edge(v1, v2, s);
        var e21 = HalfEdge.new_edge(v2, v1, s);
        e12.twin_with(e21);
        assign_edge(v1, e12);
        assign_edge(v2, e21);
        Voronoi.linkhull2(v1, v2);
    }

    static linkhull3(v1, v2, v3) {
        vertices[v1].prev = v3;
        vertices[v1].next = v2;
        vertices[v2].prev = v1;
        vertices[v2].next = v3;
        vertices[v3].prev = v2;
        vertices[v3].next = v1;
    }

    static linkhull2(v1, v2) {
        vertices[v1].prev = v2;
        vertices[v1].next = v2;
        vertices[v2].prev = v1;
        vertices[v2].next = v1;
    }
}

var input_points = new Array(MAX_POINTS);
var n_input_points = 0;

var vertices;
var p;

var leftzip, rightzip;
var left, right;
var leftedge, rightedge;

var tos;

var n_lines;

function turn(c, p1, p2) {
    return (p2.y - c.y) * (p1.x - c.x) - (p1.y - c.y) * (p2.x - c.x);
}

function clw(c, p1, p2) {
    return (turn(c, p1, p2) < 0);
}

function clweq(c, p1, p2) {
    return (turn(c, p1, p2) <= 0);
}

function ccw(c, p1, p2) {
    return (turn(c, p1, p2) > 0);
}

function ccweq(c, p1, p2) {
    return (turn(c, p1, p2) >= 0);
}

// If a and b intersects, the result will be stored in r
function intersect(r, a, b) {
    var p1 = a.p1, p2 = a.p2;
    var q1 = b.p1, q2 = b.p2;

    var dpy = p2.y - p1.y;
    var dpx = p2.x - p1.x;
    var dqy = q2.y - q1.y;
    var dqx = q2.x - q1.x;

    var dy = p1.y - q1.y;
    var dx = p1.x - q1.x;

    var div = dqy * dpx - dpy * dqx;
    var pnum = dqx * dy - dqy * dx;
    var qnum = dpx * dy - dpy * dx;

    var up, uq, ep, eq;

    if (div == 0) return false;

    up = pnum / div;
    uq = qnum / div;

    ep = CLOSE / Math.sqrt(dpx * dpx + dpy * dpy);
    eq = CLOSE / Math.sqrt(dqx * dqx + dqy * dqy);

    if ((a.reach & REACH_BACKWARD || up >= -ep) &&
        (a.reach & REACH_FORWARD  || up <= 1 + ep) &&
        (b.reach & REACH_BACKWARD || uq >= -eq) &&
        (b.reach & REACH_FORWARD  || uq <= 1 + eq)
    ) {
        r.x = p1.x + up * dpx;
        r.y = p1.y + up * dpy;
        return true;
    }

    return false;
}

function snap_point(q, l) {
    var dx1 = q.x - l.p1.x;
    var dy1 = q.y - l.p1.y;
    var dx2 = q.x - l.p2.x;
    var dy2 = q.y - l.p2.y;

    if (!(l.reach & REACH_BACKWARD) &&
        dx1 * dx1 + dy1 * dy1 <= CLOSE * CLOSE) {
        q = l.p1;
        return true;
    }
    else if (!(l.reach & REACH_FORWARD) &&
        dx2 * dx2 + dy2 * dy2 <= CLOSE * CLOSE) {
        q = l.p2;
        return true;
    }
    return false;
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
    if (vertices[v].head_edge == edge) vertices[v].head_edge = next;
    if (vertices[v].tail_edge == edge) vertices[v].tail_edge = prev;
}

function delete_edge_pair(edge) {
    remove_edge(edge.twin);
    remove_edge(edge);
    delete_line(edge.divline);
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
    var edge, next;
    var head = vertices[v].head_edge;
    if (head.is_void()) {
        do {
            head = head.next;
        } while (head.is_void());
        vertices[v].head_edge = head;
    }
    edge = head.next;
    while (edge != head) {
        next = edge.next;
        if (edge.is_void()) delete_edge_pair(edge);
        edge = next;
    }
}

function push_zipline(a, b, edgea, edgeb, zip) {
    var l = Line.new_line(zip);
    var e1 = HalfEdge.new_edge(a, b, l);
    var e2 = HalfEdge.new_edge(b, a, l);
    e1.twin_with(e2);
    left[tos] = a;
    right[tos] = b;
    leftedge[tos] = edgea;
    rightedge[tos] = edgeb;
    leftzip[tos] = e1;
    rightzip[tos] = e2;
    tos++;
}

function pop_ziplines() {
    var i, last;

    last = tos - 1;
    if (last != 0) while (left[last - 1] == left[last]) last--;

    for (i = 0; i < last; i++) {
        add_edge_after(leftedge[i], leftzip[i]);
    }

    for (i = tos - 1; i >= last; i--) add_edge_before(vertices[left[i]].head_edge, leftzip[i]);

    last = tos - 1;
    if (last != 0) while (right[last - 1] == right[last]) last--;
    for (i = 0; i < last; i++) {
        add_edge_before(rightedge[i], rightzip[i]);
    }
    for (i = tos - 1; i >= last; i--) add_edge_after(vertices[right[i]].tail_edge, rightzip[i]);

    last = tos - 1;
    vertices[left[0]].tail_edge = leftzip[0];
    vertices[right[0]].head_edge = rightzip[0];
    vertices[left[last]].head_edge = leftzip[last];
    vertices[right[last]].tail_edge = rightzip[last];

    clean_edgelist(left[0]);
    i = 1;
    while (i < tos) {
        if (left[i] != left[i - 1]) clean_edgelist(left[i]);
        i++;
    }

    clean_edgelist(right[0]);
    i = 1;
    while (i < tos) {
        if (right[i] != right[i - 1]) clean_edgelist(right[i]);
        i++;
    }

    tos = 0;
}

function rdelete(v, zip, x) {
    var edge = new HalfEdge();
    var head = new HalfEdge();
    var p1 = zip.p1.copy();
    var p2 = zip.p2.copy();
    var l0 = new Line();
    l0.p1 = x;
    l0.p2 = x;
    l0.reach = REACH_NONE;
    edge = head = vertices[v].head_edge;
    do {
        if (!edge.is_void()) {
            if (clweq(p1, p2, edge.divline.p1)
                && clweq(p1, p2, edge.divline.p2)) {
                edge.divline = l0;
            }
        }
        edge = edge.next;
    } while (edge != head);
}

function ldelete(v, zip, x) {
    var edge = new HalfEdge();
    var head = new HalfEdge();
    var p1 = zip.p1.copy();
    var p2 = zip.p2.copy();
    var l0 = new Line();
    l0.p1 = x;
    l0.p2 = x;
    l0.reach = REACH_NONE;
    edge = head = vertices[v].head_edge;
    do {
        if (!edge.is_void()) {
            if (ccweq(p1, p2, edge.divline.p1)
                && ccweq(p1, p2, edge.divline.p2)) {
                edge.divline = l0;
            }
        }
        edge = edge.next;
    } while (edge != head);
}

function starting_edgea(a, b) {
    var t;
    var d, prevd;
    var edge, end;

    end = vertices[a].tail_edge;
    prevd = clw(p[b], p[a], p[end.target]);
    edge = end.prev;
    while (edge != end) {
        if (!edge.is_void()) {
            t = edge.target;
            d = clw(p[b], p[a], p[t]);
            if (!prevd && d) return edge;
            prevd = d;
        }
        edge = edge.prev;
    }

    return edge;
}

function starting_edgeb(b, a) {
    var t;
    var d, prevd;
    var edge, end;

    end = vertices[b].head_edge;
    prevd = ccw(p[a], p[b], p[end.target]);
    edge = end.next;
    while (edge != end) {
        if (!edge.is_void()) {
            t = edge.target;
            d = ccw(p[a], p[b], p[t]);
            if (!prevd && d) return edge;
            prevd = d;
        }
        edge = edge.next;
    }

    return edge;
}

function lower_line(l, r) {
    var old;
    var a = l;
    var b = r;
    var a2 = vertices[a].next;
    var b2 = vertices[b].prev;

    while (clw(p[a], p[b], p[b2])) {
        b = b2;
        b2 = vertices[b].prev;
    }
    while (true) {
        old = a;
        while (ccw(p[b], p[a], p[a2])) {
            a = a2;
            a2 = vertices[a].next;
        }
        if (a == old) break;
        old = b;
        while (clw(p[a], p[b], p[b2])) {
            b = b2;
            b2 = vertices[b].prev;
        }
        if (b == old) break;
    }
    return [a, b];
}

function upper_line(l, r) {
    var old;
    var a = l;
    var b = r;
    var a2 = vertices[a].prev;
    var b2 = vertices[b].next;

    while (ccw(p[a], p[b], p[b2])) {
        b = b2;
        b2 = vertices[b].next;
    }
    while (true) {
        old = a;
        while (clw(p[b], p[a], p[a2])) {
            a = a2;
            a2 = vertices[a].prev;
        }
        if (a == old) break;
        old = b;
        while (ccw(p[a], p[b], p[b2])) {
            b = b2;
            b2 = vertices[b].next;
        }
        if (b == old) break;
    }

    while (ccweq(p[a], p[b], p[b2])
        && p[a].dist2(p[b2]) >= p[a].dist2(p[b])) {
        b = b2;
        b2 = vertices[b].next;
    }
    while (clweq(p[b], p[a], p[a2])
        && p[b].dist2(p[a2]) >= p[b].dist2(p[a])) {
        a = a2;
        a2 = vertices[a].prev;
    }
    return [a, b];
}

function merge_hull(l, r) {
    var a, b, c, d, c2, d2;
    a = c = l;  b = d = r;
    c2 = vertices[c].prev;
    d2 = vertices[d].next;
    if (p[c].x == p[d].x && p[c2].x == p[d2].x) {
        c = c2;
        d = d2;
    }
    else {
        var tpl = upper_line(c, d);
        c = tpl[0]; d = tpl[1];
    }

    var tpl = lower_line(a, b);
    a = tpl[0]; b = tpl[1];

    vertices[c].next = d;
    vertices[d].prev = c;
    vertices[b].next = a;
    vertices[a].prev = b;

    return [a, b];
}

function remove_duplicates(p, n) {
    var i, j;
    var s, t;

    t = p[0];
    i = 1; j = 1;
    while (j < n) {
        s = p[j];
        if (Point.comp_points(s, t) != 0) {
            p[j] = p[i];
            t = p[i++] = s;
        }
        j++;
    }
    return i;
}

function showpoint(p) {
    draw_site(p.x, p.y);
}

function showline(line) {
    var l = line.extend();
    if (l.clip()) {
        draw_edge(l.p1, l.p2);
    }
}

function draw_edge(p1, p2) {
    ctx.moveTo(p1.x, p1.y);
    ctx.lineTo(p2.x, p2.y);
}

function draw_diagram() {
    var r = Voronoi.voronoi(input_points, n_input_points);
    var i;

    ctx.clearRect(0, 0, width, height);

    draw_cell_of_voronoi(r);

    ctx.beginPath();
    for (i = 0; i < r.n_points; i++) showpoint(r.p[i]);
    ctx.fillStyle = "#000";
    ctx.fill();
    ctx.strokeStyle = "#fff";
    ctx.stroke();

    ctx.beginPath();
    for (i = 0; i < r.n_lines; i++) showline(r.l[i]);
    ctx.strokeStyle = "#000";
    ctx.stroke();
}

function gen_points(n) {
    n_input_points = 0;
    var random_pts = d3.range(n)
        .map(function(d) { return [Math.random() * width, Math.random() * height]; });
    for (var i = 0; i < random_pts.length; i++) {
        input_points[n_input_points++] = Point.from_coord(random_pts[i][0], random_pts[i][1]);
    }
}

function gen_random_color_mapping(n) {
    var result = new Object();

    for (var i = 0; i < n; i++) {
        var v1 = Math.floor(Math.random() * 255);
        var v2 = Math.floor(Math.random() * 255);
        var v3 = Math.floor(Math.random() * 255);
        result[i] = 'rgb(' + v1.toString() + ',' + v2.toString() + ',' + v3.toString() + ')';
    }
    return result;
}

function draw_cell_of_voronoi(r) {
    var cells = r.cells();
    var color_mapping = gen_random_color_mapping(r.n_vertices);
    for (var i = 0; i < cells.length; i++) {
        var cell = cells[i];
        var start = cell.vertices[0];
        ctx.beginPath();
        ctx.moveTo(start.x, start.y);
        for (var j = 0; j < cell.vertices.length; j++) {
            var vert = cell.vertices[j];
            ctx.lineTo(vert.x, vert.y);
        }
        ctx.closePath();
        ctx.fillStyle = color_mapping[i];
        ctx.fill();
    }
}

function redraw(n) {
    gen_points(n);
    draw_diagram();
}

function draw_site(x, y) {
    ctx.moveTo(x + 2.5, y);
    ctx.arc(x, y, 2.5, 0, 2 * Math.PI, false);
}

function configConfirmCallback() {
    var n = parseInt(document.getElementById("nPts").value);
    if(n <= 1 || n > MAX_POINTS) {
        errDiv.innerHTML = '<font color="red">The number of points should be between 2~' + MAX_POINTS.toString() +  '.</font>';
        errDiv.style.display = "block";
    } else if (Number.isNaN(n)){
        errDiv.innerHTML = '<font color="red">Please input something.</font>';
        errDiv.style.display = "block";
    } else {
        redraw(n);
        errDiv.style.display = "none"
    }
}

var nPtsInput = document.getElementById("nPts");
nPtsInput.addEventListener("keyup", function(event) {
    if (event.keyCode === 13) {
        event.preventDefault();
        configConfirmCallback();
    }
});

function main() {
    gen_points(TEST_NUM_PTS);
    draw_diagram();
}

main();

