var canvas = d3.select("canvas").node();
ctx = canvas.getContext("2d"),
    width = canvas.width,
    height = canvas.height;

var n = 10;

var X = new Array(n), Y = new Array(n), C = new Array(n);

function draw(mt=2) {
    ctx.fillStyle = "white"; ctx.fillRect(0, 0, width, height);

    for (var i = 0; i < n; i++) {
        X[i] = randgp(width); Y[i] = randgp(height); C[i] = randhclr();
    }

    // Brute force
    // For every pixel on the drawing plane, compare the distance
    // with every site and find the closest one
    for (var x = 0; x < width; x++) {
        for (var y = 0; y < height; y++) {
            var _min = Infinity;
            var _min_ind = -1;
            for (var i = 0; i < n; i++) {
                var d = metric(X[i] - x, Y[i] - y, mt);
                if (d < _min) {
                    _min = d;
                    _min_ind = i;
               }
            }
            ctx.fillStyle = C[_min_ind];
            ctx.fillRect(x, y, 1, 1);
        }
    }

    // sites are represented as little black square
    ctx.fillStyle = "black";
    for (var i = 0; i < n; i++) {
        ctx.fillRect(X[i], Y[i], 3, 3);
    }
}

function metric(x, y, mt=2) {
    if (mt == 1) { return Math.sqrt(x * x + y * y); }
    if (mt == 2) { return Math.abs(x) + Math.abs(y); }
}

function randgp(max) {
    return Math.floor(Math.random() * max);
}

function randhclr() {
    return "#"+
        ("00"+randgp(256).toString(16)).slice(-2) +
        ("00"+randgp(256).toString(16)).slice(-2) +
        ("00"+randgp(256).toString(16)).slice(-2);
}

draw();
