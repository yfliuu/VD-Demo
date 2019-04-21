var TEST_NUM_PTS = 10;
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

    copy() {
        var newp = new Point();
        newp.x = this.x;
        newp.y = this.y;
        return newp;
    }
};

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
};

class Edge {
    constructor() {
        this.source = null;
        this.target = null;
        this.divline = new Line();
        this.twin = null;
        this.next = null;
        this.prev = null;
    }
};

class Vertex {
    next = null;
    prev = null;
    head_edge = null;
    tail_edge = null;

    constructor() {
        this.next = null;
        this.prev = null;
        this.head_edge =  null;
        this.tail_edge = null;
    }
};

class Voronoi_result {
    constructor() {
        this.vertices = null;
        this.p = null;
        this.l = null;
        this.n_vertices = null;
        this.n_points = null;
        this.n_lines = null;
    }
};

var input_points = new Array(MAX_POINTS);
var n_input_points = 0;

var vertices;
var p;

var leftzip, rightzip;
var left, right;
var leftedge, rightedge;

var tos;

var n_lines;

function dist2(p1, p2) {
    var dx = p2.x - p1.x;
    var dy = p2.y - p1.y;
    return dx * dx + dy * dy;
}

function turn(c, p1, p2) {
    return (p2.y - c.y) * (p1.x - c.x) - (p1.y - c.y) * (p2.x - c.x);
}

function clockwise(c, p1, p2) {
    return (turn(c, p1, p2) < 0);
}
function clockwiseeq(c, p1, p2) {
    return (turn(c, p1, p2) <= 0);
}
function anticlockwise(c, p1, p2) {
    return (turn(c, p1, p2) > 0);
}
function anticlockwiseeq(c, p1, p2) {
    return (turn(c, p1, p2) >= 0);
}

function cross(l1, l2) {
    return (l1.p2.x - l1.p1.x) * (l2.p2.y - l2.p1.y)
        - (l1.p2.y - l1.p1.y) * (l2.p2.x - l2.p1.x);
}

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

    if(div == 0) {
        return 0;
    }

    up = pnum / div;
    uq = qnum / div;

    ep = CLOSE / Math.sqrt(dpx * dpx + dpy * dpy);
    eq = CLOSE / Math.sqrt(dqx * dqx + dqy * dqy);

    if((a.reach & REACH_BACKWARD || up >= -ep) &&
        (a.reach & REACH_FORWARD  || up <= 1+ep) &&
        (b.reach & REACH_BACKWARD || uq >= -eq) &&
        (b.reach & REACH_FORWARD  || uq <= 1+eq)
    )
    {
        r.x = p1.x + up * dpx;
        r.y = p1.y + up * dpy;
        return 1;
    }

    return 0;
}

function perpdiv(p1, p2) {
    var midx = (p1.x + p2.x) / 2;
    var midy = (p1.y + p2.y) / 2;
    var dx = p2.x - p1.x;
    var dy = p2.y - p1.y;
    var div = Math.sqrt(dx*dx + dy*dy);
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

function perpthrough(c, p1, p2) {
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

function point_cmp(p1, p2) {
    var d = p1.x - p2.x;
    if(d == 0) d = p1.y - p2.y;
    if(d == 0) return 0;
    return d > 0 ? 1 : -1;
}

function comp_points(p1, p2) {
    var d = p1.x - p2.x;
    if(d == 0) return p1.y - p2.y;
    return d;
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

function snap_point(q, l) {
    var dx1 = q.x - l.p1.x;
    var dy1 = q.y - l.p1.y;
    var dx2 = q.x - l.p2.x;
    var dy2 = q.y - l.p2.y;

    if(!(l.reach & REACH_BACKWARD) &&
        dx1 * dx1 + dy1 * dy1 <= CLOSE * CLOSE) {
        q = l.p1;
        return 1;
    }
    else if(!(l.reach & REACH_FORWARD) &&
        dx2 * dx2 + dy2 * dy2 <= CLOSE * CLOSE) {
        q = l.p2;
        return 1;
    }
    return 0;
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

function is_length0(line) {
    return (line.p1.x == line.p2.x && line.p1.y == line.p2.y);
}

function new_line(l) {
    var r = new Line();
    r = Object.assign(r, l);
    n_lines++;
    return r;
}

function delete_line(l) {
    n_lines--;
}

function new_edge(source, target, divline) {
    var edge = new Edge();
    edge.source = source;
    edge.target = target;
    edge.divline = divline;
    return edge;
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

function pair(e1, e2) {
    e1.twin = e2;
    e2.twin = e1;
}

function is_void_edge(edge) {
    return is_length0(edge.divline);
}

function clean_edgelist(v) {
    var edge, next;
    var head = vertices[v].head_edge;
    if(is_void_edge(head)) {
        do {
            head = head.next;
        } while(is_void_edge(head));
        vertices[v].head_edge = head;
    }
    edge = head.next;
    while(edge != head) {
        next = edge.next;
        if(is_void_edge(edge)) delete_edge_pair(edge);
        edge = next;
    }
}

function alloc_zipstack(k) {
    leftzip = new Array(k);
    rightzip = new Array(k);
    leftedge = new Array(k);
    rightedge = new Array(k);
    left = new Array(k);
    right = new Array(k);
    tos = 0;
}

function push_zipline(a, b, edgea, edgeb, zip) {
    var l = new_line(zip);
    var e1 = new_edge(a, b, l);
    var e2 = new_edge(b, a, l);
    pair(e1, e2);
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

    last = tos-1;
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

function rdelete(v, zip, x) {
    var edge = new Edge();
    var head = new Edge();
    var p1 = zip.p1.copy();
    var p2 = zip.p2.copy();
    var l0 = new Line();
    l0.p1 = x;
    l0.p2 = x;
    l0.reach = REACH_NONE;
    edge = head = vertices[v].head_edge;
    do {
        if(!is_void_edge(edge)) {
            if(clockwiseeq(p1, p2, edge.divline.p1)
                && clockwiseeq(p1, p2, edge.divline.p2)) {
                edge.divline = l0;
            }
        }
        edge = edge.next;
    } while(edge != head);
}

function ldelete(v, zip, x) {
    var edge = new Edge();
    var head = new Edge();
    var p1 = zip.p1.copy();
    var p2 = zip.p2.copy();
    var l0 = new Line();
    l0.p1 = x;
    l0.p2 = x;
    l0.reach = REACH_NONE;
    edge = head = vertices[v].head_edge;
    do {
        if(!is_void_edge(edge)) {
            if(anticlockwiseeq(p1, p2, edge.divline.p1)
                && anticlockwiseeq(p1, p2, edge.divline.p2)) {
                edge.divline = l0;
            }
        }
        edge = edge.next;
    } while(edge != head);
}

function starting_edgea(a, b) {
    var t;
    var d, prevd;
    var edge, end;

    end = vertices[a].tail_edge;
    prevd = clockwise(p[b], p[a], p[end.target]);
    edge = end.prev;
    while(edge != end) {
        if(!is_void_edge(edge)) {
            t = edge.target;
            d = clockwise(p[b], p[a], p[t]);
            if(!prevd && d) return edge;
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
    prevd = anticlockwise(p[a], p[b], p[end.target]);
    edge = end.next;
    while(edge != end) {
        if(!is_void_edge(edge)) {
            t = edge.target;
            d = anticlockwise(p[a], p[b], p[t]);
            if(!prevd && d) return edge;
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

    while(clockwise(p[a], p[b], p[b2])) {
        b = b2;
        b2 = vertices[b].prev;
    }
    for(;;) {
        old = a;
        while(anticlockwise(p[b], p[a], p[a2])) {
            a = a2;
            a2 = vertices[a].next;
        }
        if(a == old) break;
        old = b;
        while(clockwise(p[a], p[b], p[b2])) {
            b = b2;
            b2 = vertices[b].prev;
        }
        if(b == old) break;
    }
    l = a;
    r = b;
    return [a, b];
}

function upper_line(l, r) {
    var old;
    var a = l;
    var b = r;
    var a2 = vertices[a].prev;
    var b2 = vertices[b].next;

    while(anticlockwise(p[a], p[b], p[b2])) {
        b = b2;
        b2 = vertices[b].next;
    }
    for(;;) {
        old = a;
        while(clockwise(p[b], p[a], p[a2])) {
            a = a2;
            a2 = vertices[a].prev;
        }
        if(a == old) break;
        old = b;
        while(anticlockwise(p[a], p[b], p[b2])) {
            b = b2;
            b2 = vertices[b].next;
        }
        if(b == old) break;
    }

    while(anticlockwiseeq(p[a], p[b], p[b2])
        && dist2(p[a], p[b2]) >= dist2(p[a],p[b])) {
        b = b2;
        b2 = vertices[b].next;
    }
    while(clockwiseeq(p[b], p[a], p[a2])
        && dist2(p[b],p[a2]) >= dist2(p[b],p[a])) {
        a = a2;
        a2 = vertices[a].prev;
    }

    l = a;
    r = b;
    return [a, b];
}

function merge_hull(l, r) {
    var a, b, c, d, c2, d2;

    a = c = l;  b = d = r;

    c2 = vertices[c].prev;
    d2 = vertices[d].next;
    if(p[c].x == p[d].x && p[c2].x == p[d2].x) {
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

// function showhull(p);
// function showstate(l, r);
// function showline(l);
// function cycle_edges(i);
function voronoi_merge(m) {
    var zip = new Line();
    var a, b;
    var ac, bc;
    var x = new Point();
    var xa = new Point();
    var xb = new Point();
    var edgea = new Edge();
    var edgeb = new Edge();
    var starta = new Edge();
    var startb = new Edge();

    a = m;  b = m + 1;
    var tpl = merge_hull(a, b);
    a = tpl[0]; b = tpl[1];

    zip = perpdiv(p[a], p[b]);

    edgea = starta = starting_edgea(a, b);
    edgeb = startb = starting_edgeb(b, a);
    for(;;) {
        ac = 0;
        do {
            if(!is_void_edge(edgea)) {
                if(anticlockwiseeq(p[b], p[a], p[edgea.target])) break;
                if(intersect(xa, zip.copy(), edgea.divline.copy())) {
                    ac = 1;
                    break;
                }
            }
            edgea = edgea.prev;
        } while(edgea != starta);

        bc = 0;
        do {
            if(!is_void_edge(edgeb)) {
                if(clockwiseeq(p[a], p[b], p[edgeb.target])) break;
                if(intersect(xb, zip.copy(), edgeb.divline.copy())) {
                    bc = 1;
                    break;
                }
            }
            edgeb = edgeb.next;
        } while(edgeb != startb);

        if(!(ac || bc)) break;

        if(bc == 0 || (ac && xa.y <= xb.y)) {
            if(snap_point(xa, edgea.divline.copy())) {
                var czip = zip.copy();
                fix_p2(czip, xa.copy());
                if(dist2(czip.p1, czip.p2) < CLOSE * CLOSE) {
                    czip = moveline(zip.copy(), xa.copy());
                }
                rclip(edgea.divline, czip.copy(), xa.copy());
                rdelete(a, czip.copy(), xa.copy());
            }
            else {
                if(dist2(zip.p1, xa) < CLOSE * CLOSE) xa = zip.p1.copy();
                rclip(edgea.divline, zip.copy(), xa.copy());
                rdelete(a, zip.copy(), xa.copy());
            }
            fix_p2(zip, xa.copy());
            push_zipline(a, b, edgea, edgeb, zip.copy());
            x = xa.copy();
            a = edgea.target;
            edgea = starta = starting_edgea(a, b);
        }
        else {
            if(snap_point(xb, (edgeb.divline))) {
                var czip = zip.copy();
                fix_p2(czip, xb.copy());
                if(dist2(czip.p1, czip.p2) < CLOSE * CLOSE) {
                    czip = moveline(zip.copy(), xb.copy());
                }
                lclip(edgeb.divline, czip.copy(), xb.copy());
                ldelete(b, czip.copy(), xb.copy());
            }
            else {
                if(dist2(zip.p1, xb) < CLOSE * CLOSE) xb = zip.p1.copy();
                lclip(edgeb.divline, zip.copy(), xb.copy());
                ldelete(b, zip.copy(), xb.copy());
            }
            fix_p2(zip, xb.copy());
            push_zipline(a, b, edgea, edgeb, zip.copy());
            x = xb.copy();
            b = edgeb.target;
            edgeb = startb = starting_edgeb(b, a);
        }

        zip = perpthrough(x, p[a], p[b]);
    }

    push_zipline(a, b, edgea, edgeb, zip.copy());

    pop_ziplines();
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
    var sA = new_line(lA);
    var sB = new_line(lB);
    var sC = new_line(lC);
    var e12 = new_edge(v1, v2, sA);
    var e21 = new_edge(v2, v1, sA);
    var e23 = new_edge(v2, v3, sB);
    var e32 = new_edge(v3, v2, sB);
    var e31 = new_edge(v3, v1, sC);
    var e13 = new_edge(v1, v3, sC);
    pair(e12, e21);
    pair(e23, e32);
    pair(e31, e13);
    assign_edge2(v1, e12, e13);
    assign_edge2(v2, e23, e21);
    assign_edge2(v3, e31, e32);
    linkhull3(v1, v2, v3);
}

function construct2(v1, v2, v3, lA, lB) {
    var sA = new_line(lA);
    var sB = new_line(lB);
    var e12 = new_edge(v1, v2, sA);
    var e21 = new_edge(v2, v1, sA);
    var e23 = new_edge(v2, v3, sB);
    var e32 = new_edge(v3, v2, sB);
    pair(e12,e21);
    pair(e23,e32);
    assign_edge(v1, e12);
    assign_edge2(v2, e21, e23);
    assign_edge(v3, e32);
    linkhull3(v1, v3, v2);
}

function construct1(v1, v2, l) {
    var s = new_line(l);
    var e12 = new_edge(v1, v2, s);
    var e21 = new_edge(v2, v1, s);
    pair(e12, e21);
    assign_edge(v1, e12);
    assign_edge(v2, e21);
    linkhull2(v1, v2);
}

function voronoi_base(l, r, n_points) {
    var q = p[l];
    // var q = p.slice(l);
    var lineA, lineB, lineC;

    if(n_points == 3) {
        lineA = perpdiv(p[l],p[l + 1]);
        lineB = perpdiv(p[l + 1],p[l + 2]);
        var c = new Point();
        if(intersect(c,lineA,lineB)) {
            if(c.y <= p[l + 1].y) {
                lineA = perpthrough(c, p[l], p[l + 1]);
                lineB = perpthrough(c, p[l + 1], p[l + 2]);
                lineC = perpthrough(c, p[l + 2], p[l]);
                construct3(l, l+1, l+2, lineA, lineB, lineC);
            }
            else {
                lineA = perpthrough(c, p[l + 1], p[l]);
                lineB = perpthrough(c, p[l + 2], p[l + 1]);
                lineC = perpthrough(c, p[l], p[l + 2]);
                construct3(l, l+2, l+1, lineC, lineB, lineA);
            }
        }
        else {
            construct2(l, l+1, l+2, lineA, lineB);
        }
    }
    else {
        lineA = perpdiv(p[l], p[l + 1]);
        construct1(l, l+1, lineA);
    }
}

function voronoi_rr(l, r) {
    var n_points = 1 + r - l;
    var m;

    if (n_points <= 3) {
        voronoi_base(l,r,n_points);
        return;
    }

    m = Math.floor((l + r) / 2);
    voronoi_rr(l, m);
    voronoi_rr(m + 1, r);
    voronoi_merge(m);
}

function remove_duplicates(p, n) {
    var i, j;
    var s, t;

    t = p[0];
    i = 1; j = 1;
    while(j < n) {
        s = p[j];
        if(comp_points(s,t) != 0) {
            p[j] = p[i];
            t = p[i++] = s;
        }
        j++;
    }
    return i;
}

function voronoi_result_alloc(points, n) {
    var r;
    var rpoints;
    var rvertices;
    var i;

    // r = new Array(1);
    r = new Voronoi_result();

    r.n_points = n;
    r.n_lines = 0;

    rpoints = r.p = new Array(n);
    for(i = 0; i < n; i++) rpoints[i] = points[i];
    rpoints.sort(point_cmp);

    r.n_vertices = remove_duplicates(rpoints, n);
    rvertices = r.vertices = new Array(r.n_vertices);
    for(i = 0; i < r.n_vertices; i++) {
        rvertices[i] = new Vertex();
        rvertices[i].next = rvertices[i].prev = i;
        rvertices[i].head_edge = rvertices[i].tail_edge = 0;
    }

    return r;
}

function store_lines(r) {
    var i, j;
    var l;

    l = r.l = new Array(n_lines);
    j = 0;
    for(i = 0; i < r.n_vertices; i++) {
        var edge, head;
        edge = head = vertices[i].head_edge;
        if(edge) do {
            if(edge.target > i) l[j++] = (edge.divline);
            edge = edge.next;
        } while(edge != head);
    }
    r.n_lines = j;
}

function set_globals(r) {
    p = r.p;
    vertices = r.vertices;
}

function voronoi(points, n) {
    var r;
    r = voronoi_result_alloc(points, n);
    set_globals(r);

    alloc_zipstack(2*r.n_vertices + 1);

    n_lines = 0;
    if(r.n_vertices >= 2) voronoi_rr(0, r.n_vertices-1);
    store_lines(r);

    return r;
}

function inside(p) {
    return p.x > 0 && p.x < CLIP_WIDTH && p.y > 0 && p.y < CLIP_HEIGHT ? 1 : 0;
}

function clip_line(l) {
    var border = [];
    border.push(Line.from_pts(0, 0, CLIP_WIDTH, 0));
    border.push(Line.from_pts(CLIP_WIDTH, 0, CLIP_WIDTH, CLIP_HEIGHT));
    border.push(Line.from_pts(CLIP_WIDTH, CLIP_HEIGHT, 0, CLIP_HEIGHT));
    border.push(Line.from_pts(0, CLIP_HEIGHT, 0, 0));
    var x = new Point();
    var s1, s2;
    var i, result;

    s1 = (l.reach & REACH_BACKWARD) ? 0 : inside(l.p1);
    s2 = (l.reach & REACH_FORWARD)  ? 0 : inside(l.p2);
    if(s1 && s2) return 1;

    if(s1 || s2) {
        for(i = 0; i < 4; i++) {
            if(intersect(x,border[i],l)) {
                // rclip(l, border[i], x);
                rclip(l, border[i], x.copy());
                return 1;
            }
        }
    }
    result = 0;
    for(i = 0; i < 4; i++) {
        if(intersect(x,border[i],l)) {
            rclip(l, border[i], x.copy());
            result |= 1;
        }
    }
    return result;
}

function extend(l) {
    var dx = l.p2.x - l.p1.x;
    var dy = l.p2.y - l.p1.y;
    var midx = (l.p1.x + l.p2.x) / 2;
    var midy = (l.p1.y + l.p2.y) / 2;
    var r = Math.sqrt(dx * dx + dy * dy);
    var ux = dx / r;
    var uy = dy / r;
    var lx = ux * EXTEND;
    var ly = uy * EXTEND;
    if(l.reach & REACH_BACKWARD) {
        l.p1.x = midx - lx;
        l.p1.y = midy - ly;
    }
    if(l.reach & REACH_FORWARD) {
        l.p2.x = midx + lx;
        l.p2.y = midy + ly;
    }
    return l;
}

function showpoint(p) {
    drawSite(p.x, p.y);
}

function showline(line) {
    var l = extend(line);
    if(clip_line(l)) {
        drawEdgeBetween(l.p1, l.p2);
    }
}

function drawEdgeBetween(p1, p2) {
    ctx.moveTo(p1.x, p1.y);
    ctx.lineTo(p2.x, p2.y);
}

function draw_diagram() {
    var i;
    var r;
    r = voronoi(input_points, n_input_points);

    ctx.clearRect(0, 0, width, height);
    ctx.beginPath();
    for(i = 0; i < r.n_points; i++) { showpoint(r.p[i]); }
    ctx.fillStyle = "#000";
    ctx.fill();
    ctx.strokeStyle = "#fff";
    ctx.stroke();

    ctx.beginPath();
    for(i = 0; i < r.n_lines; i++) { showline(r.l[i]); }
    ctx.strokeStyle = "rgba(0,0,0,0.2)";
    ctx.stroke();
}

function add_point(x, y) {
    var p = new Vertex();
    p.x = x;
    p.y = y;
    input_points[n_input_points++] = p;
}

function clear_points() {
    n_input_points = 0;
}

function pointset1() {
    var n = TEST_NUM_PTS;
    var random_pts = d3.range(n)
        .map(function(d) { return [Math.random() * width, Math.random() * height]; });
    for(var i = 0; i < random_pts.length; i++) {
        add_point(random_pts[i][0], random_pts[i][1]);
    }
}

function drawSite(x, y) {
    ctx.moveTo(x + 2.5, y);
    ctx.arc(x, y, 2.5, 0, 2 * Math.PI, false);
}

function main() {
    pointset1();
    draw_diagram();
    return 0;
}

main();

