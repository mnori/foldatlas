
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

		buf += "</select>";
		buf += "<div id=\"d3nome-scrollbar\"></div>"

		this.initScrollbar()

		$(this.config.container).html(buf);
	}

	// set up the D3 scrollbar
	initScrollbar: function() {

		var navDims = {x: 500, y: 100};
		var maxD

		var navXScale = d3.time.scale()
	        .domain([0, 1000])
	        .range([0, navDims.x]);

	    var navYScale = d3.scale.linear()
	        .domain([0, 1]) // this doesn't matter - there is no y data
	        .range([navDims.y, 0]);

		var viewport = d3.svg.brush()
		    .x(navXScale)
		    .on("brush", function () {

		    	// this updates the main chart's x scale
		        // xScale.domain(viewport.empty() ? navXScale.domain() : viewport.extent());

		        // redraws the main chart's x scale
		        // redrawChart();
		    });
	},

	redrawChart: function() {

	}
}

