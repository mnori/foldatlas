// IDEAS

// Max window size:
// 	http://bl.ocks.org/raffazizzi/3691274

// Disable clicking outside the brush
// 	http://stackoverflow.com/questions/18036836/disable-clearing-of-d3-js-brush

(D3nome = function(config) { this.init(config); }).prototype = {

	// config should include chromosome data?
	init: function(config) {

		this.maxBrushRange = 2500000;
		this.minBrushRange = 25000;

		// Set this to config.initialBPRange
		// Storing the extent locally is a hack that kills 2 birds with 1 stone
		//  - lets us constrain the box size
		//  - lets us disable the clear functionality
		this.brushExtent = [0, this.minBrushRange];

		this.config = config;
		this.chromosomes = this.config.chromosomes;

		// select the first chromosome
		this.selectedChromosome = config.selectedChromosome; 

		// function for fetching and parsing data into the right format.
		this.fetchData = this.config.fetchData 

		// total dimensions of the browser 
		// TODO use config values
		this.totDims = {x: 898, y: 300}

		// height of nav bar area 
		this.navHeight = 40;

		this.brush = null;
		this.zoom = null;

		this.viewXScale = null;
		this.viewYScale = null;
		this.viewXAxis = null;

		this.navXScale = null;
		this.navYScale = null;
		this.navXAxis = null;
		
		this.viewElement = null;
		this.navBoxNode = null;

		this.data = null;

		this.xFormat = function (d) {
	        var prefix = d3.formatPrefix(d);
	        return prefix.scale(d) + prefix.symbol + "b";
	    }

		this.draw();
	},

	draw: function() {
		var buf = 
			"Select your chromosome:"+
	        "<select id=\"chromosome-selector\">";

	    var chromosomes = this.chromosomes;

	    for (var i = 0; i < chromosomes.length; i++) {
	    	var selected = i == this.selectedChromosome ? "selected=\"selected\"" : "";
    	buf += 
	    		"<option value=\""+chromosomes[i].id+"\" "+selected+">"+
    		 		chromosomes[i].id+" ("+chromosomes[i].length+" bp)"+
    			"</option>";
	    }

		buf += "</select><br />";
		buf += "<svg id=\"d3nome-canvas\"></svg>"

		$(this.config.container).html(buf);
		
		this.initViewer()
	},

	// Set up the chromosome scrollbar.
	initViewer: function() {

		var totDims = this.totDims;

		var navDims = {x: totDims.x, y: this.navHeight};
		// dimensions of view genome browser area
		var viewDims = {x: 898, y: (totDims.y - navDims.y)};

		var bp = this.chromosomes[this.selectedChromosome].length

		var svg = d3.select("#d3nome-canvas")
			.attr("width", totDims.x)
		    .attr("height", totDims.y)

		// view scales
		this.viewXScale = d3.scale.linear()
	        .domain([0, bp]) // this should correspond to bp
	        .range([0, navDims.x])

	    this.viewYScale = d3.scale.linear()
	        .range([navDims.y, 0]);

	    // nav scales
		this.navXScale = d3.scale.linear()
	        .domain([0, bp]) // this should correspond to bp
	        .range([0, navDims.x]);

	    this.navYScale = d3.scale.linear()
	        .range([navDims.y, 0]);

	    // nav axis
       	this.navXAxis = d3.svg.axis()
		    .scale(this.navXScale)
		    .orient("bottom")
		    .ticks(8)
			.tickFormat(this.xFormat);

		// view axis
		this.viewXAxis = d3.svg.axis()
		    .scale(this.viewXScale)
		    .orient("bottom")
		    .ticks(8)
		    .tickFormat(this.xFormat);

		// create view element

		this.zoom = d3.behavior.zoom()
			.x(this.viewXScale)
			.size([viewDims.x, viewDims.y])
			.on('zoom', $.proxy(this.onZoom, this))
			.on('zoomend', $.proxy(this.onZoomEnd, this))

		// console.log(this.zoom.size()); // size

		this.viewElement = svg.append("g")
			.attr("class", "d3nome-view")
			.attr("width", viewDims.x)
			.attr("height", viewDims.y)
			.attr("transform", "translate("+0+","+0+")")
			// .call(this.zoom)

		// view chart area - add x axis
		this.viewElement.append("g")
		    .attr("class", "x axis")
		    .attr("width", navDims.x)
		    .attr("height", navDims.y)
		    .attr("transform", "translate("+0+","+(viewDims.y - navDims.y)+")")
		    .call(this.viewXAxis)

		// Append an invisible overlay rectangle to recieve zoom commands
		svg.append("rect")
			.attr("class", "d3nome-overlay")
			.attr("width", totDims.x)
			.attr("height", totDims.y)
			.call(this.zoom);

		// create the viewport, i.e. the brush	
		this.brush = d3.svg.brush()
		    .x(this.navXScale)

		    // attach brush size change event handler
		    .on("brush", $.proxy(this.onBrush, this))

		    // on mouseup, load data.
		    .on("brushend", $.proxy(function() {

				// get the bp coords, use them to fetch the data
		    	// var extent = this.brush.extent(); 
		    	// var chrID = this.chromosomes[this.selectedChromosome].id
		    	// var start = Math.round(extent[0]);
		    	// var end = Math.round(extent[1]);

		    	this.loadData();
			}, this));

		// navbar - add brush
		this.navBoxNode = svg.append("g")
		    .attr("class", "d3nome-viewport")
		    .call(this.brush);
		this.navBoxNode
		    .selectAll("rect")
		    .attr("height", navDims.y)
		    .attr("transform", "translate("+0+","+viewDims.y+")")

	 	// navbar - add x axis
		svg.append("g")
		    .attr("class", "x axis")
		    .attr("id", "navbar-x-axis")
		    .attr("width", navDims.x)
		    .attr("height", navDims.y)
		    .attr("transform", "translate(0," + viewDims.y + ")")
		    .call(this.navXAxis)

		this.updateBrush(this.brushExtent);

		this.loadData(this.chromosomes[this.selectedChromosome].id, this.brushExtent[0], this.brushExtent[1]);
	},

	onBrush: function() {

    	// if brush is empty, use the this.brushExtent - i.e. the previous brush value
    	// otherwise grab the new extent from the brush.

    	// using this.brushExtent when empty disables the clearing behaviour
    	var domain = this.brush.empty() ? this.brushExtent : this.brush.extent()

		domain[0] = Math.round(domain[0])
    	domain[1] = Math.round(domain[1]);

    	// enforce minimum and maximum constraints on the extent

    	domain = this.restrictDomainSize(this.brushExtent, domain);

    	// store the extent data for comparison later.
    	this.brushExtent = domain;
    	this.updateBrush(domain);

    	// redraw the data
    	this.drawData();
	},

	restrictDomainSize: function(oldDomain, newDomain) {
		var newRange = newDomain[1] - newDomain[0];
    	if (newRange > this.maxBrushRange) { 
    		if (newDomain[0] < oldDomain[0]) { // dragged backwards
    			newDomain[0] = oldDomain[1] - this.maxBrushRange;
    			newDomain[1] = oldDomain[1];
    		} else { // dragged forward
    			newDomain[0] = oldDomain[0]
    			newDomain[1] = oldDomain[0] + this.maxBrushRange;
    		}
    	} else if (newRange < this.minBrushRange) {
    		if (newDomain[0] > oldDomain[0]) { // dragged backwards
    			newDomain[0] = oldDomain[1] - this.minBrushRange;
    			newDomain[1] = oldDomain[1];
    		} else { // dragged forward
    			newDomain[0] = oldDomain[0]
    			newDomain[1] = oldDomain[0] + this.minBrushRange;
    		}
    	}
    	return newDomain;
	},

	onZoom: function() {
		// this is all about prevent dragging past the boundary.
		var bp = this.chromosomes[this.selectedChromosome].length

		var domain = this.viewXScale.domain();

		// enforce min/max zoom size
		// if domain is too big, make it smaller

		domain = restrictDomainSize();

		console.log("this.brushExtent", this.brushExtent);
		console.log("domain: "+domain);

		// enforce chromosome boundaries
		if (domain[0] < 0 || domain[1] > bp) {
			if (domain[0] < 0) {
				var offset = -domain[0];
				domain = [domain[0] + offset, domain[1] + offset];

			} else if (domain[1] > bp) {
				var offset = domain[1] - bp; // positive number
				domain = [domain[0] - offset, domain[1] - offset];
			}
			// Translate the zoom by x scaled according to tot number of nucleotides
			var scaledTransX = Math.round(domain[0] / bp);
			this.zoom.translate([scaledTransX, 0]);
		}

	    // update the brush's version of the extent
	    this.brush.extent(domain);

	    // update the view X scale domain with the new data
		this.viewXScale.domain(domain);

		// keep track of our local extent data
		// Rename to navBoundaries
		this.brushExtent = domain;

		// update the view X axis as well
		this.viewElement.select(".x.axis").call(this.viewXAxis);

		// tell d3 to redraw the brush - this is important!
	    this.navBoxNode.call(this.brush);

		this.drawData();
	},

	onZoomEnd: function() {
		this.loadData();
	},

	updateBrush: function(domain) {
	    // update the brush's version of the extent
	    this.brush.extent(domain);

	    // update the view X scale domain with the new data
		this.viewXScale.domain(domain);

		// must update the zoom as well
		this.zoom.x(this.viewXScale);

		// update the view X axis as well
		this.viewElement.select(".x.axis").call(this.viewXAxis);

		// tell d3 to redraw the brush - this is important!
	    this.navBoxNode.call(this.brush);
	},

	loadData: function() {
		// fetch int value of chr

		var chrID = this.chromosomes[this.selectedChromosome].id
    	var start = Math.round(this.brushExtent[0]);
    	var end = Math.round(this.brushExtent[1]);

		var chrNum = chrID[3];
		var url = this.config.dataUrl+"?chr="+chrNum+"&start="+start+"&end="+end;

		$.ajax({
			url: url,
			context: this
		}).done($.proxy(function(results) {
			this.parseData($.parseJSON(results));
		}, this));
	},

	parseData: function(data) {
		this.data = data;

		// do cool stuff with the data here. 
		// maybs arrange into lanes for overlapping transcripts

		this.drawData();
	},

	drawData: function() {
		// Remove old elements
		this.viewElement.selectAll("rect").remove()

		// Add new elements

		var element = this.viewElement.selectAll(".d3nome-view")
		element
			.data(this.data).enter() // select missing nodes
			.append("rect")
			.attr("class", "d3nome-feature")

			.attr("x", $.proxy(function(d, i) { return this.viewXScale(d.start); }, this))
			.attr("y", function(d, i) { return 50; })

		    .attr("width", $.proxy(function(d, i) { 
		    	return (this.viewXScale(d.end) - this.viewXScale(d.start)); 
		    }, this))
		    .attr("height", 30)
	}
}

