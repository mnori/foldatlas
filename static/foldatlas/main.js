// Create the browser controller class
var BrowserController = Class.extend({

	// Constructor
	init: function(config) {
		this.nucsPerRow = 80;
		this.staticBase = config.staticBaseUrl;
		this.reactivities = {}
		this.drawNucleotideMeasurements();
	},

	// Jump to a specific transcript page
	selectTranscript: function(transcriptID) {
		$("#loading-indicator").show()
		this.changeUrl(transcriptID, "/transcript/"+transcriptID)
		$.ajax({
			url: "/ajax/transcript/"+transcriptID,
			context: this
		}).done(function(results) {
			$("#loading-indicator").hide();
			$("#transcript-data").empty();
			$("#transcript-data").html(results);
			this.drawNucleotideMeasurements();
		});
	},

	// HTML5 change URL method
	changeUrl: function(title, url) {
	    if (typeof (history.pushState) != "undefined") {
	        var obj = { Page: title, Url: url };
	        history.pushState(obj, obj.Page, obj.Url);
	    } else {
	    	// TODO show a prettier warning
	        alert("Your browser does not support HTML5. Please upgrade it.");
	    }
	},

	getMeasurementJson: function() {
		var html = $("#nucleotide-measurement-json").html();
		if (html == undefined) { // this means no measurement data to show
			return null;
		}
		var json = $.parseJSON(html);
		return json;
	},

	drawNucleotideMeasurements: function(reactivities) {
		var data = this.getMeasurementJson()
		if (data) {
			this.drawExperiment(data[1])
			this.drawExperiment(data[2])
		}
	},

	// Visualises the measurement data.
	drawExperiment: function(experiment_data) {

		var yLabelText = (experiment_data["type"] == "dms_reactivity") 
			? "Reactivity"
			: "Occupancy";

		// Add header and graph container svg element
		// Also Add the SVG itself
		var chart_id = "nucleotide-measurement-chart_"+experiment_data["id"];
		var buf = 
			"<h2>"+experiment_data["description"]+"</h2>"+
			"<svg id=\""+chart_id+"\"></svg>"

		$("#nucleotide-measurement-charts").append(buf)

		var data = experiment_data["data"]
		if (data == null) { // can happen
			return;
		}

		var nDataRows = data.length;
		var nChartRows = Math.ceil(nDataRows / this.nucsPerRow);

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
		var chart = d3.select("#"+chart_id)
			.attr("width", chartDims.x)
			.attr("height", chartDims.y)

		// Define the scales
		var yScale = d3.scale.linear()

			// range maps to the pixel dimensions.

		    // domain describes the range of data to show.
		    .range([panelDims.y, 0])
		    .domain([0, d3.max(data, function(d) { return d.measurement; })]);

		// when there is no measurement data, degrade gracefully
		if (isNaN(yScale.domain()[1])) { 
			
			chart.append("text")
		      .attr("transform", "translate("+(panelTotDims.x / 2)+", "+(panelTotDims.y / 2)+")")
		      .style("text-anchor", "middle")
		      .text("No measurement data");
			return;
		}
		
		for (var rowN = 0; rowN < nChartRows; rowN++) { // each iteration = 1 chart row
			var start = rowN * this.nucsPerRow;
			var end = start + this.nucsPerRow;

			if (end > nDataRows) {
				end = nDataRows;
			}

			var dataSlice = data.slice(start, end);

			var nucsThisRow = end - start;

			// this is for panel positioning.
			panelYOffset = rowN * panelTotDims.y;

			// Shows nucleotide numbers
			var rangeX = parseInt(panelDims.x * ((nucsThisRow - 1) / this.nucsPerRow));

			var xScale = d3.scale.linear()
			    .range([0, rangeX], .1) 
			    .domain([start - 0.5, (end - 1) + 0.5])

		   	// Create axis objects
			var xAxis = d3.svg.axis()
			    .scale(xScale)
			    .orient("bottom")
			    .ticks(nucsThisRow)
				.tickFormat(function(d, i) { return data[d].nuc; })

			// xAxis.select("tick")
			// 	.append("rect")
			// 	.setAttr("width", 10)
			// 	.setAttr("height", 10);

			var yAxis = d3.svg.axis()
			    .scale(yScale)
			    .orient("left")
	    		.ticks(3); // how many ticks to show.

			// Add x-axis objects to the chart.
			var bgWidth = parseInt(panelDims.x / this.nucsPerRow) + 1;

			chart.append("g")
				.attr("class", "x axis")
				.attr("transform", "translate("+
					panelMargin.left+","+
					(panelYOffset + panelDims.y + panelMargin.top)
				+")")
				.call(xAxis)

				// select the X axis tick element
				.selectAll(".tick")

				// aadd a rect to it, first child means it draws in the background
				.insert("rect", ":first-child")

				// set position and dimensions of rectangle element.
				.attr("transform", "translate("+(-bgWidth / 2)+", "+10+")")
				.attr("width", bgWidth)
				.attr("height", 10)

				// highlight nucleotides with missing reactivities
				.attr("class", function(n, i) {
					return (dataSlice[i].measurement == null) ? 
						"missing-bg" : "not-missing-bg";
				})

			// Add y-axis objects to the chart
			chart.append("g")
				.attr("class", "y axis")
				.attr("transform", "translate("+
					panelMargin.left+","+
					(panelYOffset + panelMargin.top)+
				")")
				.call(yAxis);

			// Add x axis label
		    chart.append("text")
		        .attr("transform", "translate("+
		        	(panelMargin.left + (panelDims.x / 2))+","+
		        	(panelYOffset + panelDims.y + panelMargin.top)+
	        	")")
		        .style("text-anchor", "middle")
		        .attr("dy", "2.7em");
		        // .text("Nucleotide");

			// Add y-axis label
			chart.append("text")
		        .attr("transform", "rotate(-90)")
		        .attr("y", panelMargin.left) // this is actually X direction, because we rotated.
		        .attr("x", (-panelYOffset - (panelMargin.top + (panelDims.y / 2))))
		        .attr("dy", "-2.7em")
		        .style("text-anchor", "middle")
		        .text(yLabelText);

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
			    .y(function(d) { return panelYOffset + panelMargin.top + yScale(d.measurement); });

	 		chart.append("path")
				.datum(dataSlice) // get data specific to this row
				.attr("class", "line")
				.attr("d", lineGen);
		} // End looping through chart rows
	}
})



// $("#chromosome-selector").change(function() {
// 	// TODO reload the page .. difficult / impossible to change the chromosome 
// 	// ... 
// });