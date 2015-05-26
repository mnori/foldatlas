/**
 * D3nome: D3-based genome browser
 * Matthew Norris 2015
 */
(D3nome = function(config) { this.init(config); }).prototype = {

	// config should include chromosome data?
	init: function(config) {

		this.transcriptHeight = 15;
		this.labelHeight = 20;
		this.laneMargin = 10;
		this.intronBulge = 10;
		this.intronBulgeOffset = 20;

		// total height of lane excluding the margin
		this.laneHeight = this.intronBulge + this.transcriptHeight + this.labelHeight

		this.maxBrushRange = 2500000;
		this.minBrushRange = 25000;

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
		this.fetchData = this.config.fetchData 

		// total dimensions of the browser 
		// TODO use config values
		this.totSvgDims = {x: 898, y: 300}
		this.initialSvgDims = {x: this.totSvgDims.x, y: this.totSvgDims.y}

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

		buf += 	"</select><br />"+
				"<svg id=\"d3nome-canvas\"></svg>"+
				"<div id=\"d3nome-resize-bar\" class=\"ui-resizable-handle ui-resizable-s\" style=\"width: "+this.totSvgDims.x+"px;\">"+
					"..."+
				"</div>"

		// Add the HTML
		$(this.config.container).html(buf);

		// Initialise the viewer SVG
		this.initViewer()

		// Add resizer
		this.initialContainerHeight = $(this.config.container).height()

		$(this.config.container).resizable({
			handles: {"s": "#d3nome-resize-bar"}, 
		}).bind({resize: $.proxy(function(event, ui) {
			var newContainerHeight = ui.size.height;
			var heightDiff = newContainerHeight - this.initialContainerHeight;
			var newSvgHeight = this.initialSvgDims.y + heightDiff;
			this.totSvgDims.y = newSvgHeight;

			// set the new canvas height
			$("#d3nome-canvas").height(this.totSvgDims.y);

			// must also set the overlay height
			$("#d3nome-overlay").height(this.totSvgDims.y);

			// reset the x grid - could also do initViewer(), but that's rather slow
			this.calcDims();
			this.setXGrid();

		}, this)});

		
	},

	calcDims: function() {
		this.navDims = {x: this.totSvgDims.x, y: this.navHeight};
		this.viewDims = {x: 898, y: (this.totSvgDims.y - this.navDims.y)};
	},

	// Set up the chromosome scrollbar.
	initViewer: function() {

		$("#d3nome-canvas").empty();

		this.calcDims();

		var bp = this.chromosomes[this.selectedChromosome].length

		var svg = d3.select("#d3nome-canvas")
			.attr("width", this.totSvgDims.x)
		    .attr("height", this.totSvgDims.y)

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
		    .attr("transform", "translate("+0+","+(this.navDims.y)+")")
		    .call(this.viewXAxis)

		// add grid lines too
		this.setXGrid();

		// SVG version
		// svg.append("rect")
		// 	.attr("class", "d3nome-overlay")
		// 	.attr("width", totSvgDims.x)
		// 	.attr("height", totSvgDims.y)
		// 	.call(this.zoom);

		// OVERLAY
		// Uses a div overlay
		var foreignObject = this.viewElement.append("foreignObject")
			.attr("x", 0).attr("y", this.navDims.y)
		
		var htmlDoms = foreignObject.append("xhtml:body")
		    .style("margin",0)
		    .style("padding",0);

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

		// Append the overlay div
		htmlDoms.append("div")
			.attr("id", "d3nome-overlay")
			.attr("style", "width: "+this.viewDims.x+"px; height: "+this.viewDims.y+"px;")
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
		    .attr("transform", "translate("+0+","+0+")")

	 	// navbar - add x axis
		svg.append("g")
		    .attr("class", "d3nome-x d3nome-axis")
		    .attr("id", "navbar-x-axis")
		    .attr("width", this.navDims.x)
		    .attr("height", this.navDims.y)
		    .attr("transform", "translate(0,"+this.navDims.y+")")
		    .call(this.navXAxis)

		this.updateBrush(this.navBoundaries);

		this.loadData(this.chromosomes[this.selectedChromosome].id, this.navBoundaries[0], this.navBoundaries[1]);
	},

	setXGrid: function() {
		this.viewElement.selectAll(".d3nome-x.d3nome-grid").remove()
		this.viewXGrid = d3.svg.axis()
		    .scale(this.viewXScale)
		    .orient("top")
		    .ticks(8)
		    .tickFormat("")
		    .tickSize(this.viewDims.y - this.navDims.y)
		    .outerTickSize(0);

		this.gridElement = this.viewElement.append("g")
		    .attr("class", "d3nome-x d3nome-grid")
		    .attr("transform", "translate("+0+","+(this.viewDims.y + 1)+")")
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

	loadData: function() {
		// fetch int value of chr

		var chrID = this.chromosomes[this.selectedChromosome].id
    	var start = Math.round(this.navBoundaries[0]);
    	var end = Math.round(this.navBoundaries[1]);

		var chrNum = chrID[3];
		var url = this.config.dataUrl+"?chr="+chrNum+"&start="+start+"&end="+end;

		$.ajax({
			url: url,
			context: this
		}).done($.proxy(function(results) {
			this.parseData($.parseJSON(results));
		}, this));
	},

	// Parses the data from the API into a format that can be easily 
	// understood by the D3 library.
	parseData: function(data) {

		// create an empty object to use as a key-value store
		var transcripts = {};

		// first pass - collect UTR and CDS sequences
		for (var i = 0; i < data.length; i++) {
			var feature = data[i];

			if (feature.feature_type == "transcript") {
				var transcriptID = feature.id;

				if (transcripts[transcriptID] === undefined) {
					transcripts[transcriptID] = {
						id: transcriptID,
						start: feature.start,
						end: feature.end,
						features: []
					};
				}
			}

			if (	feature.feature_type == "CDS" || 
					feature.feature_type.indexOf("UTR") > -1) {

				var transcriptID = feature["Parent"];
				var transcript = transcripts[transcriptID];

				// add the feature data
				transcript.features.push({
					transcriptID: feature.Parent,
					type: (feature.feature_type == "CDS") ? "cds" : "utr",
					start: feature.start,
					end: feature.end
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
					end: feature.end
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
						end: currFeature.start - 1
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

		// Figure out whether transcript overlaps another in the given lane
		var overlaps = function(transcript, lane) {
			for (var i = 0; i < lane.length; i++) {
				var otherTranscript = lane[i];

				// test whether transcript overlaps the other transcript.
				if (	(	transcript.start 		>= otherTranscript.start && 
						 	transcript.start 		<= otherTranscript.end) ||
						(	transcript.end 			>= otherTranscript.start && 
							transcript.end 			<= otherTranscript.end) ||
						(	otherTranscript.start 	>= transcript.start  && 
						 	otherTranscript.start 	<= transcript.end) ||
						(	otherTranscript.end 	>= transcript.start && 
						 	otherTranscript.end 	<= transcript.end)) {
					return true;
				}
			}
			return false;
		}

		// Organise transcripts into "lanes", according to overlaps.
		var transcriptsByLane = []
		$.each(transcripts, function(transcriptID, transcript) {
			var laneNum = 0;
			while(true) {
				if (transcriptsByLane[laneNum] === undefined) {
					transcriptsByLane[laneNum] = [];
				}
				if (!overlaps(transcript, transcriptsByLane[laneNum])) {
					transcript.lane = laneNum;
					transcriptsByLane[laneNum].push(transcript);
					break;
				}
				laneNum++;
			}
		});

		var dataOut = [];
		$.each(transcripts, function(transcriptID, transcript) {
			// add lane metadata to the features
			for (var i = 0; i < transcript.features.length; i++) {
				transcript.features[i].lane = transcript.lane;
			}

			// add transcript to final output
			dataOut.push(transcript);
		});

		this.data = dataOut;
		this.drawData();
	},

	drawData: function() {
		this.viewElement.selectAll("rect").remove()
		this.viewElement.selectAll("path.d3nome-feature-intron").remove()
		this.viewElement.selectAll("g.d3nome-transcript").remove()
		this.viewElement.selectAll("text.d3nome-transcript-label").remove()

		var element = this.viewElement.selectAll(".d3nome-view")

		var getYPos = $.proxy(function(d) {
			return this.laneMargin + (d.lane * (this.laneHeight + this.laneMargin));
		}, this);

		// Append 1 group element per transcript
		var transcriptGroups = element
			.data(this.data).enter()
			.append('g')
			.attr('class', "d3nome-transcript")

		// Append text as ForeignObject to get HTML formatting and nice background
		var foreignObjects = transcriptGroups.append("foreignObject")
		    .attr("x", $.proxy(function(d, i) { return this.viewXScale(d.start); }, this))
		    .attr("y", $.proxy(function(d, i) { 
		    	return (this.navDims.y * 2) + this.laneMargin + getYPos(d) + this.transcriptHeight; 
		    }, this))

		var htmlDoms = foreignObjects.append("xhtml:body")
		    .style("margin",0)
		    .style("padding",0);

		htmlDoms.append("div")
			.attr("class", "d3nome-transcript-label")
			.attr("data-transcript_id", function(d) { return d.id; })
			.text(function(d) { return d.id; })
	
		// When gene label is clicked, use the callback		
		$(".d3nome-transcript-label").bind("click", {self:this}, function(event) {
			var self = event.data.self;
			var transcriptID = $(this).data("transcript_id");
			self.config.geneClick(transcriptID); // callback
		})

		// .on("mouseover", function(ev) { sendEventBehind(this, "mouseover"); })
		// .on("mouseout", function(ev) { sendEventBehind(this, "mouseout"); })

		// old method with SVG labels
		// transcriptGroups
		// 	.append("text")
		// 	.attr("class", "d3nome-transcript-label")
		//     .attr("x", $.proxy(function(d, i) { return this.viewXScale(d.start); }, this))
		//     .attr("y", $.proxy(function(d, i) { 
		//     	return (this.navDims.y * 2) + getYPos(d) + this.transcriptHeight; 
		//     }, this))
		//     .attr("dy", ".35em")
		//     .text(function(d) { return d.id; })
		//     .on("click", function() { alert("test"); })

		////////////////////////////////////////////////////

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
			.attr("class", "d3nome-feature-utr")
			.attr("x", $.proxy(function(d, i) { return this.viewXScale(d.start); }, this))
			.attr("y", $.proxy(function(d, i) { 
				return this.navDims.y + this.laneMargin + getYPos(d); 
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
			.attr("class", "d3nome-feature-cds")
			.attr("x", $.proxy(function(d, i) { return this.viewXScale(d.start); }, this))
			.attr("y", $.proxy(function(d, i) { 
				return this.navDims.y + this.laneMargin + getYPos(d); 
			}, this))
			.attr("width", $.proxy(function(d, i) { 
				return (this.viewXScale(d.end + 1) - this.viewXScale(d.start)); 
			}, this))
			.attr("height", this.transcriptHeight)

		// Introns - represented using bezier curves
		transcriptGroups.selectAll('g.d3nome-transcript')
			.data(function(transcript) { return getFeatures(transcript, "intron"); })
			.enter()
			.append("path")
			.attr("d", $.proxy(function(d) {
				// we are drawing a Bezier curve.
				// see https://developer.mozilla.org/en-US/docs/Web/SVG/Tutorial/Paths
				// for explanation.

				var bulge = this.intronBulge;
				var bulgeOffset = this.intronBulgeOffset;

				var yOffset = this.navDims.y + this.laneMargin + getYPos(d) + 1;
				var startStr = this.viewXScale(d.start)+" "+yOffset;
				var endStr = this.viewXScale(d.end + 1)+" "+yOffset;

				var control1 = this.viewXScale(d.start+bulgeOffset)+" "+(yOffset-bulge);
				var control2 = this.viewXScale((d.end + 1)-bulgeOffset)+" "+(yOffset-bulge);

				return "M"+startStr+" C "+control1+", "+control2+", "+endStr;
			}, this))
			.attr("class", "d3nome-feature-intron")
	}
}

