/**
 * D3nome: D3-based genome browser
 * Matthew Norris 2015
 */
(D3nome = function(config) { this.init(config); }).prototype = {

	// config should include chromosome data?
	init: function(config) {

		this.labelHeight = 20;
		
		this.intronBulge = 10;
		this.intronBulgeOffset = 20;

		// total height of lane excluding the margin
		this.transcriptHeight = 15;
		this.transcriptLaneMargin = 10;
		this.transcriptLaneHeight = this.intronBulge + this.transcriptHeight + this.labelHeight

		this.geneHeight = 15;
		this.geneLaneMargin = 10;

		this.maxBrushRange = 2500000;
		this.minBrushRange = 25000;
		this.initRange = 250000;

		// If zoomed out more than this, it will just fetch gene boundaries rather than features
		this.simpleThreshold = 100000;

		// If zoomed out further than this, no data will be shown.
		this.blankThreshold = 2500000;

		// Set this to config.initialBPRange
		// Storing the extent locally is a hack that kills 2 birds with 1 stone
		//  - lets us constrain the box size
		//  - lets us disable the clear functionality
		this.navBoundaries = [0, this.minBrushRange];

		this.config = config;
		this.chromosomes = this.config.chromosomes;

		// select the first chromosome
		this.selectedChromosome = config.selectedChromosome; 

		// function for fetching and parsing data into the right format.
		this.fetchData = this.config.fetchData;

		// total dimensions of the browser 
		// TODO use config values
		this.horizMargin = 20;

		this.innerSvgDims = {x: 902, y: 300}
		this.totSvgDims = {x: 902 + (this.horizMargin * 2), y: 300};

		// Left / right padding
		// Prevents labels being chopped off at the sides

		this.initialSvgDims = {x: this.totSvgDims.x, y: this.totSvgDims.y}

		this.minHeight = this.totSvgDims.y;

		this.navDims = null;
		this.viewDims = null;

		// height of nav bar area (axis is the same as this)
		this.navHeight = 30;

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

		this.initialContainerHeight = null;
		this.gridElement = null;

		// Detailed data about gene features
		this.featureData = null;

		// Generalised data about genes
		this.geneData = null;

		this.xFormat = function (d) {
	        var prefix = d3.formatPrefix(d);
	        return prefix.scale(d) + prefix.symbol + "b";
	    }
		this.draw();
	},

	draw: function() {
		// Draw canvas and its resize bar
		var buf = 	
			"<div id=\"d3nome-canvas-container\" style=\"width: "+this.totSvgDims.x+"px;\">"+
				"<svg id=\"d3nome-canvas\"></svg>"+
			"</div>"+
			"<div id=\"d3nome-resize-bar\" class=\"ui-resizable-handle ui-resizable-s\" style=\""+
					"width: "+this.innerSvgDims.x+"px; "+
					"left: "+this.horizMargin+"px;\">"+
				"..."+
			"</div>"+
			"<div class=\"d3nome-flank-overlay l\"></div>"+
			"<div class=\"d3nome-flank-overlay r\"></div>";

		// Draw chromosome select menu
		buf += 
	        "<select id=\"d3nome-chromosome-selector\">";

	    var chromosomes = this.chromosomes;

	    for (var i = 0; i < chromosomes.length; i++) {
	    	var selected = i == this.selectedChromosome ? "selected=\"selected\"" : "";
	    	var len = Math.round(chromosomes[i].length / 1000000)+" Mb";
    		buf += 
	    		"<option value=\""+i+"\" "+selected+">"+
    		 		chromosomes[i].id+" ("+len+")"+
    			"</option>";
	    }
		buf += 	"</select>";

		// Add the HTML to the container
		$(this.config.container).html(buf);
		$(this.config.container).css({
			"width": this.totSvgDims.x+"px",
			"margin-left": (-this.horizMargin)+"px",
			"margin-right": (-this.horizMargin)+"px"
		});

		// Initialise the chromosome selector menu
		this.initChromosomeSelector();

		// Initialise the viewer SVG
		this.initViewer();

		// Add resizer
		this.initialContainerHeight = $(this.config.container).height()

		$(this.config.container).resizable({
			handles: {
				s: "#d3nome-resize-bar"
			}, 
			minHeight: this.minHeight
		}).bind({resize: $.proxy(function(event, ui) {
			var newContainerHeight = ui.size.height;
			var heightDiff = newContainerHeight - this.initialContainerHeight;
			var newSvgHeight = this.initialSvgDims.y + heightDiff;
			this.totSvgDims.y = newSvgHeight;

			

			this.setSvgDims();

			// must also set the overlay height
			this.setOverlayDims();
			
			// reset the x grid - could also do initViewer(), but that's rather slow
			this.calcDims();
			this.setXGrid();

			// draw data so that lines appear in the background
			this.drawData();

		}, this)});
	},

	initChromosomeSelector: function() {

		// Transform native element to fancy jquery version
		$("#d3nome-chromosome-selector").selectmenu({width: 150});

		$("#d3nome-chromosome-selector").on(
			"selectmenuselect",
			$.proxy(function(event, ui) {
				var chrInd = $(event.target).val();
				if (chrInd != this.selectedChromosome) { // only jump if the chromosome has changed.
					this.jumpToPosition(chrInd, [0, this.initRange], false);
				}
			}, this)
		);

		$("#d3nome-chromosome-selector-button").css({"right": (10 + this.horizMargin)+"px"});
	},

	chrIDToInd: function(chrID) {
		for (var i = 0; i < this.chromosomes.length; i++) {
			if (chrID == this.chromosomes[i].id) {
				return i;
			}
		}
	},

	// coords - 2 element array
	jumpToPosition: function(chrInd, coords, updateMenu) {

		var oldChrInd = this.selectedChromosome;
		this.selectedChromosome = chrInd;

		if (oldChrInd != this.selectedChromosome) { 
			// New chromosome selected. Must whole thing
			this.draw();
		}

		// make sure coords are within the chromosome's range
		var chrLen = this.chromosomes[this.selectedChromosome].length;
		var coords = [
			coords[0] < 0 ? 0 : coords[0],
			coords[1] >= chrLen ? chrLen - 1 : coords[1]
		];

		if (updateMenu) {
			// if selected chromosome has changed, must update the chromosome selector menu
			// $("#d3nome-chromosome-selector").selectmenu("value", chrInd);
			$("#d3nome-chromosome-selector").val(chrInd);
			$("#d3nome-chromosome-selector").selectmenu("refresh", true);
		}

		// update the view
		this.navBoundaries = coords;
		this.updateBrush(this.navBoundaries);
		this.loadData();
	},

	// Change all the overlay dimensions to reflect the new height of the plotting area
	setOverlayDims: function() {
		var styleStr = 
			"width: "+this.viewDims.x+"px; "+
			"height: "+(this.viewDims.y - this.navHeight)+"px; "+
			"top: "+(this.navHeight * 2)+"px; "+
			"left: "+(this.horizMargin)+"px;";
		d3.select("#d3nome-overlay").attr("style", styleStr)
		d3.select("#d3nome-underlay").attr("style", styleStr)

		// Set size of left flank cover - this prevents genes poking out the sides
		$(".d3nome-flank-overlay.l").css({
			"top": (this.navHeight * 2)+"px",
			"left": "0px",
			"width": this.horizMargin+"px",
			"height": (this.viewDims.y - this.navHeight)+"px"
		});

		// Same but for the right hand side
		$(".d3nome-flank-overlay.r").css({
			"top": (this.navHeight * 2)+"px",
			"left": (this.viewDims.x + this.horizMargin + 1)+"px",
			"width": this.horizMargin+"px",
			"height": (this.viewDims.y - this.navHeight)+"px"
		});
	},

	calcDims: function() {
		this.navDims = {x: this.totSvgDims.x - (this.horizMargin * 2), y: this.navHeight};
		this.viewDims = {x: this.totSvgDims.x - (this.horizMargin * 2), y: (this.totSvgDims.y - this.navDims.y)};
	},

	setSvgDims: function() {
		// grab the elements
		var canvasContainer = d3.select("#d3nome-canvas-container")
		var svg = d3.select("#d3nome-canvas")

		// set HTML width and height (not actually necessary)
		// svg.attr("width", this.totSvgDims.x).attr("height", this.totSvgDims.y);
		// canvasContainer.attr("width", this.totSvgDims.x).attr("height", this.totSvgDims.y);

		// set style width and height
		var styleStr = 
			"width: "+this.totSvgDims.x+"px; "+
			"height: "+this.totSvgDims.y+"px; ";
		svg.attr("style", styleStr);
		canvasContainer.attr("style", styleStr);
	},

	// Set up the chromosome scrollbar.
	initViewer: function() {

		$("#d3nome-canvas").empty();

		this.calcDims();

		var bp = this.chromosomes[this.selectedChromosome].length

		this.setSvgDims();

		// view scales
		this.viewXScale = d3.scale.linear()
	        .domain([0, bp]) // this should correspond to bp
	        .range([0, this.navDims.x])

	    this.viewYScale = d3.scale.linear()
	        .range([this.navDims.y, 0]);

	    // nav scales
		this.navXScale = d3.scale.linear()
	        .domain([0, bp]) // this should correspond to bp
	        .range([0, this.navDims.x]);

	    this.navYScale = d3.scale.linear()
	        .range([this.navDims.y, 0]);

	    // nav axis
       	this.navXAxis = d3.svg.axis()
		    .scale(this.navXScale)
		    .orient("top")
		    .ticks(8)
			.tickFormat(this.xFormat);

		// chart axis lines
		this.navXAxis = d3.svg.axis()
		    .scale(this.navXScale)
		    .orient("top")
		    .ticks(8)
			.tickFormat(this.xFormat);

		// view x axis
		this.viewXAxis = d3.svg.axis()
		    .scale(this.viewXScale)
		    .orient("top")
		    .ticks(8)
		    .tickFormat(this.xFormat);

		this.zoom = d3.behavior.zoom()
			.x(this.viewXScale)
			.size([this.viewDims.x, this.viewDims.y])
			.on('zoom', $.proxy(this.onZoom, this))
			.on('zoomend', $.proxy(this.onZoomEnd, this))

		this.zoomTranslate = this.zoom.translate();
		this.zoomScale = this.zoom.scale();

		var svg = d3.select("#d3nome-canvas");
		this.viewElement = svg.append("g")
			.attr("class", "d3nome-view")
			.attr("width", this.viewDims.x)
			.attr("height", this.viewDims.y)
			.attr("transform", "translate("+0+","+this.navDims.y+")")

		// view chart area - add x axis
		this.viewElement.append("g")
		    .attr("class", "d3nome-x d3nome-axis")
		    .attr("width", this.navDims.x)
		    .attr("height", this.navDims.y)
		    .attr("transform", "translate("+this.horizMargin+","+(this.navDims.y)+")")
		    .call(this.viewXAxis)

		// add grid lines too
		this.setXGrid();

		// This is a hack for sending mouse events behind the overlay.
		// see http://www.vinylfox.com/forwarding-mouse-events-through-layers/
		var getElementBehind = function(element) {

			// Mouse position within the offset div
			var point = d3.mouse(element);

			// Absolute position of the overlay
			var overlayOffset = $("#d3nome-overlay").offset()

			// Get total offset, must take window scroll offsets into account
			var totOffset = {
				x: point[0] + overlayOffset.left - $(window).scrollLeft(),
				y: point[1] + overlayOffset.top - $(window).scrollTop()
			}

			// Put the event behind the overlay
			$("#d3nome-overlay").hide()
			var elementOut = $(document.elementFromPoint(totOffset.x, totOffset.y))
			$("#d3nome-overlay").show() 
			return elementOut; // .trigger(eventName);
		};		

		d3.select("#d3nome-canvas-container")
			.append("div")
			.attr("id", "d3nome-underlay")

		d3.select("#d3nome-canvas-container").append("div")
			.attr("id", "d3nome-overlay")
			.call(this.zoom)

			// this gets us a nice pointer when the user hovers over a gene label
			// also decorates the gene labels
			.on("mousemove", function(ev) {
				$(".d3nome-transcript-label-hover").removeClass("d3nome-transcript-label-hover");
				var element = getElementBehind(this); 
				if (element.hasClass("d3nome-transcript-label")) {
					$(this).addClass("d3nome-overlay-hover");
					element.addClass("d3nome-transcript-label-hover");
				} else {
					$(this).removeClass("d3nome-overlay-hover")
				}
			})

			// sends the click event to the correct gene label element
			.on("click", function() { 
				var element = getElementBehind(this); 
				element.trigger("click");
			});

		// Append the underlay - this is where the labels will go

		this.setOverlayDims();

		// create the viewport, i.e. the brush	
		this.brush = d3.svg.brush()
		    .x(this.navXScale)

		    // attach brush size change event handler
		    .on("brush", $.proxy(this.onBrush, this))

		    // on mouseup, load data.
		    .on("brushend", $.proxy(function() {
		    	this.loadData();
			}, this));

		// navbar - add brush
		this.navBoxNode = svg.append("g")
		    .attr("class", "d3nome-viewport")
		    .call(this.brush);
		this.navBoxNode
		    .selectAll("rect")
		    .attr("height", this.navDims.y)
		    .attr("transform", "translate("+this.horizMargin+","+0+")")

	 	// navbar - add x axis
		svg.append("g")
		    .attr("class", "d3nome-x d3nome-axis")
		    .attr("id", "navbar-x-axis")
		    .attr("width", this.navDims.x)
		    .attr("height", this.navDims.y)
		    .attr("transform", "translate("+this.horizMargin+","+this.navDims.y+")")
		    .call(this.navXAxis)

		this.updateBrush(this.navBoundaries);

		// this.loadData(this.chromosomes[this.selectedChromosome].id, this.navBoundaries[0], this.navBoundaries[1]);
		this.loadData();
	},

	setXGrid: function() {
		this.viewElement.selectAll(".d3nome-x.d3nome-grid").remove()
		this.viewElement.selectAll(".d3nome-grid-bg").remove()

		// Light grey background for the main gene display area
		this.viewElement.append("rect")
			.attr("class", "d3nome-grid-bg")
			.attr("width", this.viewDims.x)
			.attr("height",	this.viewDims.y)
			.attr("x", this.horizMargin)
			.attr("y", this.navDims.y + 1);

		// Add the gridlines on top of the grey background
		this.viewXGrid = d3.svg.axis()
		    .scale(this.viewXScale)
		    .orient("top")
		    .ticks(8)
		    .tickFormat("")
		    .tickSize(this.viewDims.y - this.navDims.y)
		    .outerTickSize(0);

		this.gridElement = this.viewElement.append("g")
		    .attr("class", "d3nome-x d3nome-grid")
		    .attr("transform", "translate("+this.horizMargin+","+(this.viewDims.y + 1)+")")
			.call(this.viewXGrid);
	},

	onBrush: function() {

    	// if brush is empty, use the this.navBoundaries - i.e. the previous brush value
    	// otherwise grab the new extent from the brush.

    	// using this.navBoundaries when empty disables the clearing behaviour
    	var domain = this.brush.empty() ? this.navBoundaries : this.brush.extent()

		domain[0] = Math.round(domain[0])
    	domain[1] = Math.round(domain[1]);

    	// store the extent data for comparison later.
    	this.navBoundaries = domain;
    	this.updateBrush(domain);

    	// redraw the data
    	this.drawData();
	},

	onZoom: function() {
		// this is all about prevent dragging past the boundary.
		var bp = this.chromosomes[this.selectedChromosome].length

		var domain = this.viewXScale.domain();
		var range = domain[1] - domain[0];

		if (domain[0] <= 0 && domain[1] >= bp) {
			range = domain[1] - domain[0];
			domain = [0, bp];
			this.zoom.translate([0, 0]);
			this.zoom.scale(this.zoom.scale() * (range / bp));

		} else if (domain[0] <= 0) {
			var offset = -domain[0];
			domain = [domain[0] + offset, domain[1] + offset];
			this.zoom.translate([0, 0]);

		} else if (domain[1] >= bp) {
			var offset = domain[1] - bp; // positive number
			domain = [domain[0] - offset, domain[1] - offset];
			var scaledTransX = domain[0] / bp;
			this.zoom.translate([scaledTransX, 0]);
		}

	    // update the brush's version of the extent
	    this.brush.extent(domain);

	    // update the view X scale domain with the new data
		this.viewXScale.domain(domain);

		// keep track of our local extent data
		// Rename to navBoundaries
		this.navBoundaries = domain;

		// update the view X axis as well
		this.viewElement.select(".d3nome-x.d3nome-axis").call(this.viewXAxis);
		this.viewElement.select(".d3nome-x.d3nome-grid").call(this.viewXGrid); // .. and the grid

		// tell d3 to redraw the brush - this is important!
	    this.navBoxNode.call(this.brush);

		this.drawData();
	},

	onZoomEnd: function() {
		this.loadData();
		this.zoom.x(this.viewXScale); // this fixes a shit load of issues!
	},

	updateBrush: function(domain) {
	    // update the brush's version of the extent
	    this.brush.extent(domain);

	    // update the view X scale domain with the new data
		this.viewXScale.domain(domain);

		// must update the zoom as well
		this.zoom.x(this.viewXScale);

		// update the view X axis/grid as well
		this.viewElement.select(".d3nome-x.d3nome-axis").call(this.viewXAxis);
		this.viewElement.select(".d3nome-x.d3nome-grid").call(this.viewXGrid);

		// tell d3 to redraw the brush - this is important!
	    this.navBoxNode.call(this.brush);
	},

	// Get boundaries specified by the user's zoom
	getBounds: function() {
		var start = Math.round(this.navBoundaries[0])
		var end = Math.round(this.navBoundaries[1])
		return {
			start: start,
			end: end,
			diff: end - start
		}
	},

	// Get expanded and clipped data boundaries, for API requests
	getDataBounds: function() {
		var bounds = this.getBounds();
		var diff = bounds.diff;

		// Add some extra padding.
    	bounds.start -= diff;
    	bounds.end += diff;

    	// Make sure it falls within the range
    	if (bounds.start < 0) {
    		bounds.start = 0;
    	}
    	var chrLen = this.chromosomes[this.selectedChromosome].length;
    	if (bounds.end >= chrLen) {
    		bounds.end = chrLen;
    	}
    	bounds.diff = bounds.end - bounds.start;
    	return bounds;
	},

	loadData: function() {
		var bounds = this.getBounds();
		if (bounds.diff < this.simpleThreshold) {
			this.loadFeatureData();
		} else if (bounds.diff >= this.simpleThreshold && bounds.diff < this.blankThreshold) {
			this.loadGeneData();
		} else {
			this.drawData(); // this will draw the blank
		}
	},

	drawData: function() {
		this.viewElement.selectAll(".d3nome-feature-cds rect").remove()
		this.viewElement.selectAll(".d3nome-feature-utr rect").remove()
		this.viewElement.selectAll("path.d3nome-feature-intron").remove()
		this.viewElement.selectAll("g.d3nome-transcript").remove()
		this.viewElement.selectAll("text.d3nome-zoomed-out-message").remove()
		d3.select("#d3nome-underlay").selectAll(".d3nome-transcript-label").remove()

		var bounds = this.getBounds();

		if (	this.featureData != null && 
				bounds.diff < this.simpleThreshold) {

			this.drawFeatureData();
		} else if (
				this.geneData != null && 
				bounds.diff >= this.simpleThreshold && 
				bounds.diff < this.blankThreshold) {

			this.drawGeneData();

		} else if (bounds.diff >= this.blankThreshold) {
			this.drawBlankMessage();

		} else {
			// .. maybs show loading spinner?
		}

		// // "Prefer" data sources - use the other if one is missing
		// if (this.geneData != null && (
		// 			(bounds.diff >= this.simpleThreshold) ||  // TODO less than blankThreshold
		// 			(bounds.diff < this.simpleThreshold && this.featureData == null)
		// 		)) {
		// 	this.drawGeneData();
		// } else if (
		// 		this.featureData != null && (
		// 			(bounds.diff < this.simpleThreshold) || 
		// 			(bounds.diff >= this.simpleThreshold && this.geneData == null) // TODO less than blankThreshold
		// 		)) { 
		// 	this.drawFeatureData();
		// } 
	},

	drawBlankMessage: function() {
		var element = this.viewElement // selectAll(".d3nome-view")
			.append("text")
			.attr("class", "d3nome-zoomed-out-message")
			.attr("x", this.viewDims.x / 2)
			.attr("y", this.viewDims.y / 2)
			.attr("text-anchor", "middle")
			.attr("dominant-baseline", "central") // vertical centering
			.text("Zoom in to see genes.")


	},

	// Loads simplified gene data via ajax
	loadGeneData: function() {
		var chrID = this.chromosomes[this.selectedChromosome].id
		var chrNum = chrID[3];
		var bounds = this.getDataBounds();
		var url = this.config.genesUrl+"?chr="+chrNum+"&start="+bounds.start+"&end="+bounds.end;

		$.ajax({
			url: url,
			context: this
		}).done($.proxy(function(results) {
			this.parseGeneData($.parseJSON(results));
		}, this));
	},

	// Parses simplified gene data
	parseGeneData: function(data) {
		// ... must arrange into lanes
		this.addLaneData(data);
		this.geneData = data;
		this.drawData();
	},

	drawGeneData: function() {
		var element = this.viewElement.selectAll(".d3nome-view")
		var getYPos = $.proxy(function(d) {

			// this.geneHeight = 15;
			// this.geneLaneMargin: 10;

			// d.lane = 0; // temp hack
			return this.geneLaneMargin + (d.lane * (this.geneHeight + this.geneLaneMargin));
		}, this);

		var geneGroups = element
			.data(this.geneData).enter()
			.append('g')
			.attr('class', "d3nome-transcript")
			.append("rect")
			.attr("class", function(d) { return "d3nome-gene "+d.direction; })
			.attr("x", $.proxy(function(d, i) { return this.viewXScale(d.start) + this.horizMargin; }, this))
			.attr("y", $.proxy(function(d, i) { 
				return this.navDims.y + this.geneLaneMargin + getYPos(d); 
			}, this))
			.attr("width", $.proxy(function(d, i) { 
				return (this.viewXScale(d.end + 1) - this.viewXScale(d.start)); 
			}, this))
			.attr("height", this.geneHeight)

	},

	loadFeatureData: function() {
		var chrID = this.chromosomes[this.selectedChromosome].id
		var chrNum = chrID[3];
		var bounds = this.getDataBounds();
		var url = this.config.featuresUrl+"?chr="+chrNum+"&start="+bounds.start+"&end="+bounds.end;

		$.ajax({
			url: url,
			context: this
		}).done($.proxy(function(results) {
			this.parseFeatureData($.parseJSON(results));
		}, this));
	},

	// Parses the data from the API into a format that can be easily 
	// understood by the D3 library.
	parseFeatureData: function(data) {

		// create an empty object to use as a key-value store
		var transcripts = {};

		// first pass - collect UTR and CDS sequences
		for (var i = 0; i < data.length; i++) {
			var feature = data[i];
			if (feature.feature_type == "transcript") {
				var transcriptID = feature.id;
				var direction = feature.direction;

				if (transcripts[transcriptID] === undefined) {
					transcripts[transcriptID] = {
						id: transcriptID,
						start: feature.start,
						end: feature.end,
						direction: direction,
						features: []
					};
				}
			} else if (	feature.feature_type == "CDS" || 
				feature.feature_type.indexOf("UTR") > -1) {

				var transcriptID = feature["Parent"];
				var transcript = transcripts[transcriptID];
				var direction = transcript.direction;

				// add the feature data
				transcript.features.push({
					transcriptID: feature.Parent,
					type: (feature.feature_type == "CDS") ? "cds" : "utr",
					start: feature.start,
					end: feature.end,
					direction: direction
				});
			}
		}

		// Find the non-CDS transcripts
		var nonCds = {};
		$.each(transcripts, function(transcriptID, transcript) {
			if (transcript.features.length == 0) {
				nonCds[transcriptID] = transcript;
			}
		});

		// For non-CDS transcripts, add exons as UTRs
		for (var i = 0; i < data.length; i++) {
			var feature = data[i];
			var transcriptID = feature["Parent"];
			if (	nonCds[transcriptID] !== undefined && 
					feature.feature_type == "exon") {

				var transcript = nonCds[transcriptID];
				transcript.features.push({
					transcriptID: transcriptID,
					type: "utr",
					start: feature.start,
					end: feature.end,
					direction: transcript.direction
				});
			}
		}

		// Now infer the intronic sequences from the features already stored
		var sortFeatures = function(a, b) {
			if (a.start > b.start) {
				return 1;
			} else {
				return -1;
			}
		}
		$.each(transcripts, function(transcriptID, transcript) {
			var features = transcript.features;

			// sort features by their positions.
			features.sort(sortFeatures);
			
			var prevFeature = null;
			var introns = [];

			// figure out locations of introns
			for (var i = 0; i < features.length; i++) {
				var currFeature = features[i];

				if (	prevFeature != null && 
						currFeature.start - prevFeature.end > 1) { // there is a gap

					introns.push({
						transcriptID: transcriptID,
						type: "intron",
						start: prevFeature.end + 1, // offset by 1 since coords are inclusive.
						end: currFeature.start - 1,
						direction: currFeature.direction
					});
				}
				prevFeature = currFeature;
			}

			// add the introns to the features array
			features = features.concat(introns);

			// sort features, including introns, by their positions.
			features.sort(sortFeatures);

			// set back into the data structure.
			transcript.features = features;
		});

		// TODO make this more generalised


		var featuresArray = [];
		$.each(transcripts, function(transcriptID, transcript) {
			featuresArray.push(transcript);
		});

		this.addLaneData(featuresArray);

		var dataOut = [];
		for (var i = 0; i < featuresArray.length; i++) {
			var transcript = featuresArray[i];
			// add lane metadata to each feature
			for (var j = 0; j < transcript.features.length; j++) {
				transcript.features[j].lane = transcript.lane;
			}

			// add transcript to final output
			dataOut.push(transcript);
		}

		this.featureData = dataOut;
		this.drawData();
	},

	// Add lane data to features objects. The features must contain start and end properties for it to work.
	addLaneData: function(features) {

		// Figure out whether feature overlaps another in the given lane
		var overlaps = function(feature, lane) {
			for (var i = 0; i < lane.length; i++) {
				var otherFeature = lane[i];

				// test whether feature overlaps the other feature.
				if (	feature.end > otherFeature.start && 
						feature.start < otherFeature.end) {
					return true;
				}
			}
			return false;
		}

		// Organise features into "lanes", according to overlaps.
		var featuresByLane = []
		for (var i = 0; i < features.length; i++) {
			var feature = features[i];
			var laneNum = 0;
			while(true) {
				if (featuresByLane[laneNum] === undefined) {
					featuresByLane[laneNum] = [];
				}
				if (!overlaps(feature, featuresByLane[laneNum])) {
					feature.lane = laneNum;
					featuresByLane[laneNum].push(feature);
					break;
				}
				laneNum++;
			}
		}
	},

	drawFeatureData: function() {
		var element = this.viewElement.selectAll(".d3nome-view")
		var getYPos = $.proxy(function(d) {
			return this.transcriptLaneMargin + (d.lane * (this.transcriptLaneHeight + this.transcriptLaneMargin));
		}, this);

		// Append 1 group element per transcript
		var transcriptGroups = element
			.data(this.featureData).enter()
			.append('g')
			.attr('class', "d3nome-transcript")

		// Append text as div to get HTML formatting and nice background
		d3.select("#d3nome-underlay").selectAll("#d3nome-underlay")

			.data(this.featureData).enter().append("div")
			.attr("class", "d3nome-transcript-label")
			.attr("data-transcript_id", function(d) { return d.id; })
			.attr("style", $.proxy(function(d) {
				var leftVal = this.viewXScale(d.start);
				var topVal = this.transcriptLaneMargin + getYPos(d) + this.transcriptHeight;
				var out = 	"left: "+Math.round(leftVal)+"px; "+
							"top: "+topVal+"px";
				return out;
			}, this))
			.text(function(d) { return d.id; })

		// When gene label is clicked, use the callback		
		$(".d3nome-transcript-label").bind("click", {self:this}, function(event) {
			var self = event.data.self;
			var transcriptID = $(this).data("transcript_id");
			self.config.geneClick(transcriptID); // callback
		})

		// Find features of specific type.
		function getFeatures(transcript, featureType) {
			var out = []
			var features = transcript.features;
			for (var i = 0; i < features.length; i++) {
				if (features[i].type == featureType) {
					out.push(features[i]);
				}
			}
			return out;
		}

		// UTRs
		transcriptGroups.selectAll('g.d3nome-transcript')
			.data(function(transcript) { return getFeatures(transcript, "utr"); })
			.enter()
			.append("rect")
			.attr("class", function(d) { return "d3nome-feature-utr "+d.direction; })
			.attr("x", $.proxy(function(d, i) { return this.viewXScale(d.start) + this.horizMargin; }, this))
			.attr("y", $.proxy(function(d, i) { 
				return this.navDims.y + this.transcriptLaneMargin + getYPos(d); 
			}, this))
			.attr("width", $.proxy(function(d, i) { 
				// Must add 1 here since boundaries are inclusive.
				return (this.viewXScale(d.end + 1) - this.viewXScale(d.start)); 
			}, this))
			.attr("height", this.transcriptHeight)

		// CDSs
		transcriptGroups.selectAll('g.d3nome-transcript')
			.data(function(transcript) { return getFeatures(transcript, "cds"); })
			.enter()
			.append("rect")
			.attr("class", function(d) { return "d3nome-feature-cds "+d.direction; })
			.attr("x", $.proxy(function(d, i) { return this.viewXScale(d.start) + this.horizMargin; }, this))
			.attr("y", $.proxy(function(d, i) { 
				return this.navDims.y + this.transcriptLaneMargin + getYPos(d); 
			}, this))
			.attr("width", $.proxy(function(d, i) { 
				return (this.viewXScale(d.end + 1) - this.viewXScale(d.start)); 
			}, this))
			.attr("height", this.transcriptHeight)

		// Introns - represented using bezier curves
		// see https://developer.mozilla.org/en-US/docs/Web/SVG/Tutorial/Paths
		// for explanation of drawing method
		transcriptGroups.selectAll('g.d3nome-transcript')
			.data(function(transcript) { return getFeatures(transcript, "intron"); })
			.enter()
			.append("path")
			.attr("d", $.proxy(function(d) {
				var bulge = this.intronBulge;
				var bulgeOffset = this.intronBulgeOffset;

				var yOffset = this.navDims.y + this.transcriptLaneMargin + getYPos(d) + 1;
				var startStr = (this.viewXScale(d.start) + this.horizMargin)+" "+yOffset;
				var endStr = (this.viewXScale(d.end + 1) + this.horizMargin)+" "+yOffset;

				var control1 = (this.viewXScale(d.start+bulgeOffset) + this.horizMargin)+" "+(yOffset-bulge);
				var control2 = (this.viewXScale((d.end + 1)-bulgeOffset) + this.horizMargin)+" "+(yOffset-bulge);

				return "M"+startStr+" C "+control1+", "+control2+", "+endStr;
			}, this))
			.attr("class", function(d) { return "d3nome-feature-intron "+d.direction; })
	}
}

