var canvas = d3.select("canvas").on("touchmove mousemove", moved).node(),
    context = canvas.getContext("2d"),
    width = canvas.width,
    height = canvas.height;

var sitesList = d3.range(100)
    .map(function(d) { return [Math.random() * width, Math.random() * height]; });

var voronoi = d3.voronoi()
    .extent([[-1, -1], [width + 1, height + 1]]);

var nPtsInput = document.getElementById("nPts");
nPtsInput.addEventListener("keyup", function(event) {
    if (event.keyCode === 13) {
        event.preventDefault();
        configConfirmCallback();
    }
});

var showDelaunay = true;
redraw(sitesList);

function moved() {
    sitesList[0] = d3.mouse(this);
    redraw(sitesList);
}

$("#distMetricButton :input").change(function() {
    console.log($this);
});

function configConfirmCallback() {
    var n = parseInt(document.getElementById("nPts").value);
    if(n <= 1 || n > 300) {
        errDiv.text = "The number of points should between 2~300.";
        errDiv.style.display = "block";
    } else if (Number.isNaN(n)){
        errDiv.text("Please input something");
        errDiv.style.display = "block";
    } else {
        sitesList = d3.range(n)
            .map(function(d) { return [Math.random() * width, Math.random() * height]; });
        redraw(sitesList);
    }
}

function triToggleBtnClick() {
    showDelaunay = !showDelaunay;
    redraw(siteList);
}

function redraw(sites) {
    var diagram = voronoi(sites),
        links = diagram.links(),
        polygons = diagram.polygons();

    context.clearRect(0, 0, width, height);
    context.beginPath();
    drawCell(polygons[0]);
    context.fillStyle = "#f00";
    context.fill();

    context.beginPath();
    for (var i = 0, n = polygons.length; i < n; ++i) drawCell(polygons[i]);
    context.strokeStyle = "#000";
    context.stroke();

    if (showDelaunay) {
        context.beginPath();
        for (var i = 0, n = links.length; i < n; ++i) drawLink(links[i]);
        context.strokeStyle = "rgba(0,0,0,0.2)";
        context.stroke();
    }

    context.beginPath();
    drawSite(sites[0]);
    context.fillStyle = "#fff";
    context.fill();

    context.beginPath();
    for (var i = 1, n = sites.length; i < n; ++i) drawSite(sites[i]);
    context.fillStyle = "#000";
    context.fill();
    context.strokeStyle = "#fff";
    context.stroke();
}

function drawSite(site) {
    context.moveTo(site[0] + 2.5, site[1]);
    context.arc(site[0], site[1], 2.5, 0, 2 * Math.PI, false);
}

function drawLink(link) {
    context.moveTo(link.source[0], link.source[1]);
    context.lineTo(link.target[0], link.target[1]);
}

function drawCell(cell) {
    if (!cell) return false;
    context.moveTo(cell[0][0], cell[0][1]);
    for (var j = 1, m = cell.length; j < m; ++j) {
        context.lineTo(cell[j][0], cell[j][1]);
    }
    context.closePath();
    return true;
}
