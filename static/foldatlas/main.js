
// Define the website's controller. All page changes must go through this badboy
function BrowserController(config) {

	var instance = this;
	
	this.nucsPerRow = 80;
	this.staticBase = config.staticBaseUrl;
	this.reactivities = {}

	// initialise: function(config) {
	// 	this.staticBase = window.browserControllerConfig;
	// },

	this.init = function() {
		instance.drawReactivities();
	}

	this.selectTranscript = function(transcriptID) {
		$("#loading-indicator").show()
		this.changeUrl(transcriptID, "/transcript/"+transcriptID)

		$.ajax({
			url: "/ajax/transcript/"+transcriptID,
			context: this
		}).done(function(results) {
			$("#loading-indicator").hide();
			$("#transcript-data").empty();
			$("#transcript-data").html(results);
			instance.drawReactivities();
		});
	}

	this.changeUrl = function(title, url) {
	    if (typeof (history.pushState) != "undefined") {
	        var obj = { Page: title, Url: url };
	        history.pushState(obj, obj.Page, obj.Url);
	    } else {
	        alert("Your browser does not support HTML5. Please upgrade it.");
	    }
	}

	this.drawReactivities = function(reactivities) {

		var data = this.getReactivitiesJson()
		if (data == null) {
			return;
		}

		var nDataRows = data.length;
		var nChartRows = Math.ceil(nDataRows / this.nucsPerRow);

		console.log("nChartRows: "+nChartRows);

		// Define chart dimensions including axis panelMargins
		var panelMargin = {top: 15, right: 60, bottom: 30, left: 70}
		var panelTotDims = {x: 898, y: 100}

		// dims without margins
		var panelDims = { 
			x: panelTotDims.x - panelMargin.left - panelMargin.right,
			y: panelTotDims.y - panelMargin.bottom - panelMargin.top
		}

		// Total dimensions of chart across all panels and margins.
		var chartDims = {
			x: panelTotDims.x,
			y: panelTotDims.y * nChartRows
		}

		// var	width = panelTotDims.x - panelMargin.left - panelMargin.right,
		// 	height = panelTotDims.y - panelMargin.bottom - panelMargin.top;

		// Init the chart's container element
		var chart = d3.select("#reactivities-chart")
			.attr("width", chartDims.x)
			.attr("height", chartDims.y)

		// Define the scales
		var yScale = d3.scale.linear()

			// range maps to the pixel dimensions.
		    .range([panelDims.y, 0])

		    // domain describes the range of data to show.
		    .domain([0, d3.max(data, function(d) { return d.reactivity; })]);

		// when there is no reactivity data, degrade gracefully
		if (isNaN(yScale.domain()[1])) { 
			
			chart.append("text")
		      .attr("transform", "translate("+(panelTotDims.x / 2)+", "+(panelTotDims.y / 2)+")")
		      .style("text-anchor", "middle")
		      .text("No reactivity data");
			return;
		}
		
		for (var rowN = 0; rowN < nChartRows; rowN++) {
			var start = rowN * this.nucsPerRow;
			var end = start + this.nucsPerRow;

			if (end > nDataRows) {
				end = nDataRows;
			}
			var nucsThisRow = end - start;

			// this is for panel positioning.
			panelYOffset = rowN * panelTotDims.y;

			// Shows nucleotide numbers
			var rangeX = parseInt(panelDims.x * ((nucsThisRow - 1) / this.nucsPerRow));

			var xScale = d3.scale.linear()
			    .range([0, rangeX], .1) 
			    // .domain([0, d3.max(data, function(d) { return d.position; })]);
			    // .domain([start, end - 1])
			    .domain([start, end - 1])

			    // custom tick format, to show nucleotides.

			// Shows nucleotide letters (BORKED)
			// var xScale = d3.scale.ordinal()
			//     .domain(data.map(function(d) { 
   //  				return d.nuc; 
			// 	}))
			//     .rangeRoundBands([0, panelDims.x]);

		   	// Create axis objects
			var xAxis = d3.svg.axis()
			    .scale(xScale)
			    .orient("bottom")
			    .ticks(nucsThisRow)

			    // this replaces the position number with the nucleotide at that
				.tickFormat(function(d, i) {
					return data[d].nuc;
				});

			var yAxis = d3.svg.axis()
			    .scale(yScale)
			    .orient("left")
	    		.ticks(3); // how many ticks to show.

			// Add x-axis objectsto the chart
			chart.append("g")
				.attr("class", "x axis")
				.attr("transform", "translate("+
					panelMargin.left+","+
					(panelYOffset + panelDims.y + panelMargin.top)
				+")")
				.call(xAxis);

			// Add x axis label
		    chart.append("text")
		        .attr("transform", "translate("+
		        	(panelMargin.left + (panelDims.x / 2))+","+
		        	(panelYOffset + panelDims.y + panelMargin.top)+
	        	")")
		        .style("text-anchor", "middle")
		        .attr("dy", "2.7em");
		        // .text("Nucleotide");

			// Add y-axis objects to the chart
			chart.append("g")
				.attr("class", "y axis")
				.attr("transform", "translate("+
					panelMargin.left+","+
					(panelYOffset + panelMargin.top)+
				")")
				.call(yAxis)

			// Add y-axis label
			chart.append("text")
		        .attr("transform", "rotate(-90)")
		        .attr("y", panelMargin.left) // this is actually X direction, because we rotated.
		        .attr("x", (-panelYOffset - (panelMargin.top + (panelDims.y / 2))))
		        .attr("dy", "-2.7em")
		        .style("text-anchor", "middle")
		        .text("Reactivity");

	        // Add length label
	        var panelDimsX = panelDims.x * (nucsThisRow / this.nucsPerRow);
			chart.append("text")
		        .attr("x", panelMargin.left + panelDimsX)
		        .attr("y", panelYOffset + panelMargin.top + panelDims.y)
		        .style("text-anchor", "left")
		        .attr("dy", "1.3em")
		        .text(end);

		    // add the actual line chart
			var lineGen = d3.svg.line()
			    .x(function(d) { return panelMargin.left + xScale(d.position); })
			    .y(function(d) { return panelYOffset + panelMargin.top + yScale(d.reactivity); });

			// must slice the data to get out the bit we're interested in
	 		chart.append("path")
				.datum(data.slice(start, end))
				.attr("class", "line")
				.attr("d", lineGen);

		} // end looping through chart rows

		

		// old code for generating bar chart
		// chart.append('svg:path')
		// 	.attr('d', lineGen(data))
		// 	.attr('stroke', 'green')
		// 	.attr('stroke-width', 2)
		// 	.attr('fill', 'none');

		// // // Define bar width, depends on n values and also width of canvas
		// var barWidth = width / data.length;

		// // Add bar elements to the chart
		// var bar = chart.selectAll("g")
		// 	.data(data)
		// 	.enter().append("g")
		// 	.attr("transform", function(d, i) { return "translate(" + (panelMargin.left + (i * barWidth)) + ", "+panelMargin.top+")"; })

		// // Add rectangles to the bars
		// bar.append("rect")
		// 	.attr("y", function(d) { return yScale(d.reactivity); })
		// 	.attr("height", function(d) { return height - yScale(d.reactivity); })
		// 	.attr("width", barWidth);	

	}

	this.getReactivitiesJson = function() {
		var html = $("#reactivities-json").html();

		if (html == undefined) { // this means no reactivity data to show
			return null;
		}
		var json = $.parseJSON(html);
		return json;
	}

	this.init()
}

$(document).ready(function () { 
    window.genoverse = new Genoverse(window.genoverseConfig); 
    window.browserController = new BrowserController(window.foldatlasConfig);
});

$("#chromosome-selector").change(function() {
	// TODO reload the page .. difficult / impossible to change the chromosome 
	// ... 
});