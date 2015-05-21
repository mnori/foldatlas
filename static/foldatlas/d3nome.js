// IDEAS

// Max window size:
// 	http://bl.ocks.org/raffazizzi/3691274

// Disable clicking outside the brush
// 	http://stackoverflow.com/questions/18036836/disable-clearing-of-d3-js-brush

(D3nome = function(config) { this.init(config); }).prototype = {

	// config should include chromosome data?
	init: function(config) {

		this.maxBrushRange = 2500000;
		this.minBrushRange = 1000000;

		// Set this to config.initialBPRange
		// Storing the extent locally is a hack that kills 2 birds with 1 stone
		//  - lets us constrain the box size
		//  - lets us disable the clear functionality
		this.brushExtent = [0, 1000000];

		this.config = config;
		this.chromosomes = this.config.chromosomes;

		// select the first chromosome
		this.selectedChromosome = config.selectedChromosome; 

		// function for fetching and parsing data into the right format.
		this.fetchData = this.config.fetchData 
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
		buf = "<svg id=\"d3nome-canvas\"></svg>"

		$(this.config.container).html(buf);
		
		this.initViewer()
	},

	// Set up the chromosome scrollbar.
	initViewer: function() {

		// total dimensions of the browser
		var totDims = {x: 898, y: 300}

		// dimensions of nav bar area 
		var navDims = {x: 898, y: 50};

		// dimensions of view genome browser area
		var viewDims = {x: 898, y: (totDims.y - navDims.y)}

		var bp = this.chromosomes[this.selectedChromosome].length

		var svg = d3.select("#d3nome-canvas")
			.attr("width", totDims.x)
		    .attr("height", totDims.y)

		// view scales
		var viewXScale = d3.scale.linear()
	        .domain([0, bp]) // this should correspond to bp
	        .range([0, navDims.x]);

	    var viewYScale = d3.scale.linear()
	    	.domain([0, 1]) // this doesn"t matter - there is no y data
	        .range([navDims.y, 0]);

	    // nav scales
		var navXScale = d3.scale.linear()
	        .domain([0, bp]) // this should correspond to bp
	        .range([0, navDims.x]);

	    var navYScale = d3.scale.linear()
	        .range([navDims.y, 0]);

	    // nav axis
       	var navXAxis = d3.svg.axis()
		    .scale(navXScale)
		    .orient("bottom");

		// view axis
		var mainXAxis = d3.svg.axis()
		    .scale(viewXScale)
		    .orient("bottom");

		// view area - add element
		var viewElement = svg.append("g")
			.attr("class", "d3nome-view")
			.attr("width", viewDims.x)
			.attr("height", viewDims.y)
			.attr("transform", "translate("+0+","+0+")")

		// view chart area - add x axis
		viewElement.append("g")
		    .attr("class", "x axis")
		    .attr("width", navDims.x)
		    .attr("height", navDims.y)
		    .attr("transform", "translate("+0+","+(viewDims.y - navDims.y)+")")
		    .call(mainXAxis)

		// create the viewport, i.e. the brush	
		var brush = d3.svg.brush()
		    .x(navXScale)

		    // attach brush size change event handler
		    .on("brush", $.proxy(function() {

		    	// problem is that brush.extent[1] seems to get stuck on the starting point.
		    	var domain = brush.empty() ? this.brushExtent : brush.extent()

				domain[0] = Math.round(domain[0])
		    	domain[1] = Math.round(domain[1]);

		    	var newRange = domain[1] - domain[0];
		    	if (newRange > this.maxBrushRange) { 
		    		if (domain[0] < this.brushExtent[0]) { // dragged backwards
		    			domain[0] = this.brushExtent[1] - this.maxBrushRange;
		    			domain[1] = this.brushExtent[1];
		    		} else { // dragged forward
		    			domain[0] = this.brushExtent[0]
		    			domain[1] = this.brushExtent[0] + this.maxBrushRange;
		    		}
		    	} else if (newRange < this.minBrushRange) {
		    		if (domain[0] > this.brushExtent[0]) { // dragged backwards
		    			domain[0] = this.brushExtent[1] - this.minBrushRange;
		    			domain[1] = this.brushExtent[1];
		    		} else { // dragged forward
		    			domain[0] = this.brushExtent[0]
		    			domain[1] = this.brushExtent[0] + this.minBrushRange;
		    		}
		    	}

		    	this.brushExtent = domain;

			    // update the extent
			    brush.extent(domain);

			    // update the view scale domain with the new data
				viewXScale.domain(domain);

				// update the view axis
				viewElement.select(".x.axis").call(mainXAxis);

				// tell d3 to redraw the brush - this is important!
			    navBoxNode.call(brush);

		    }, this))

		    // attach mouseup event handler
		    .on("brushend", $.proxy(function() {

				// get the bp coords, use them to fetch the data
		    	var extent = brush.extent(); 
		    	var chrID = this.chromosomes[this.selectedChromosome].id
		    	var start = Math.round(extent[0]);
		    	var end = Math.round(extent[1]);

		    	this.loadData(chrID, start, end);
			}, this));

		// navbar - add brush
		var navBoxNode = svg.append("g")
		    .attr("class", "d3nome-viewport")
		    .call(brush);

		navBoxNode
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
		    .call(navXAxis)

		// update the extent
	    brush.extent(this.brushExtent);

	    // update the view scale domain with the new data
		viewXScale.domain(this.brushExtent);

		// update the view axis
		viewElement.select(".x.axis").call(mainXAxis);

		// tell d3 to redraw the brush - this is important!
		// note - this breaks the clear handler
	    navBoxNode.call(brush);
	},

// 	drawBar: function(domain) {
		

// 		var x_range = d3.scale.linear()
// 			.domain([0, range.length])
// 			.range([0, width]); 

// 		svg.selectAll("rect.items").remove();

// 		svg.selectAll("rect.items")
// 		.data(range)
// 		.enter().append("svg:rect")
// 		.attr("class", "items")
// 		.attr("x", function(d, i) {return x_range(i);})
// 		.attr("y", 0)
// 		.attr("width",  width/range.length-2)
// 		.attr("height", 100)
// 		.attr("fill", function(d) {return d})
// 		.attr("title", function(d) {return d});
// }
	// }

	loadData: function(chrID, start, end) {
		// fetch int value of chr

		var chrNum = chrID[3];
		var url = this.config.dataUrl+"?chr="+chrNum+"&start="+start+"&end="+end;

		$.ajax({
			url: url,
			context: this
		}).done($.proxy(function(results) {
			this.parseAndDisplayData($.parseJSON(results));
		}, this));
	},

	parseAndDisplayData: function(data) {
		console.log("Data grabbed");
		// console.log("DATA", data);
	}
}

