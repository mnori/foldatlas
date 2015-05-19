// IDEAS

// Max window size:
// 	http://bl.ocks.org/raffazizzi/3691274

// Disable clicking outside the brush
// 	http://stackoverflow.com/questions/18036836/disable-clearing-of-d3-js-brush

(D3nome = function(config) { this.init(config); }).prototype = {

	// config should include chromosome data?
	init: function(config) {
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
		
		this.initScrollbar()
	},

	// Set up the chromosome scrollbar.
	initScrollbar: function() {
		var navDims = {x: 898, y: 50};
		var bp = this.chromosomes[this.selectedChromosome].length

		var svg = d3.select("#d3nome-canvas")
			.attr("width", navDims.x)
		    .attr("height", navDims.y)

		var navXScale = d3.scale.linear()
	        .domain([0, bp]) // this should correspond to bp
	        .range([0, navDims.x]);

	    var navYScale = d3.scale.linear()
	        // .domain([0, 1]) // this doesn"t matter - there is no y data
	        .range([navDims.y, 0]);

       	var navXAxis = d3.svg.axis()
		    .scale(navXScale)
		    .orient("bottom");

		// create the viewport, i.e. the brush	
		var brush = d3.svg.brush()
		    .x(navXScale)
		    .on("brushend", $.proxy(function() {

				// get the bp coords, use them to fetch the data
		    	var extent = brush.extent(); 
		    	var chrID = this.chromosomes[this.selectedChromosome].id
		    	var start = Math.round(extent[0]);
		    	var end = Math.round(extent[1]);

		    	this.loadData(chrID, start, end);
			}, this));

		// add x axis to navbar
		svg.append("g")
		    .attr("class", "x axis")
		    .attr("width", navDims.x)
		    .attr("height", navDims.y)
		    .attr("transform", "translate(0," + 0 + ")")
		    .call(navXAxis)

		// add the viewport brush element
		svg.append("g")
		    .attr("class", "d3nome-viewport")
		    .call(brush)
		    .selectAll("rect")
		    .attr("height", navDims.y);
	},

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
		console.log(data);
	}
}

