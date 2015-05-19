
(D3nome = function(config) { this.init(config); }).prototype = {

	// config should include chromosome data?
	init: function(config) {
		this.config = config;
		this.selectedChromosome = config.selectedChromosome; // select the first chromosome
		this.draw();
	},

	draw: function() {
		var buf = 
			"Select your chromosome:"+
	        "<select id=\"chromosome-selector\">";

	    var chromosomes = this.config.chromosomes;

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

		var bp = this.config.chromosomes[this.selectedChromosome].length

		console.log("bp: "+bp);

		var svg = d3.select("#d3nome-canvas")
			.attr("width", navDims.x)
		    .attr("height", navDims.y)

		var navXScale = d3.scale.linear()
	        .domain([0, bp]) // this should correspond to bp
	        .range([0, navDims.x]);

	    var navYScale = d3.scale.linear()
	        .domain([0, 1]) // this doesn"t matter - there is no y data
	        .range([navDims.y, 0]);

       	var navXAxis = d3.svg.axis()
		    .scale(navXScale)
		    .orient("bottom");

		// create the viewport, i.e. the brush	
		var brush = d3.svg.brush()
		    .x(navXScale)

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
	redrawChart: function() {

	}
}

